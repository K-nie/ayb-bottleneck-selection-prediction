#!/usr/bin/env python3
"""
fineSTRUCTURE / ChromoPainter haplotype-based population structure inference
on the AYB panel (Lawson et al. 2012 PLoS Genet).

Method
------
  1. Use the Beagle-phased VCF from script 36 (results/36_ihs/data/).
  2. Convert per-chromosome to ChromoPainter `.phase` format
     (haplotype rows, SNP positions in a P line).
  3. Build the ID file from the 95 sample names.
  4. Run `fs <proj>.cp -idfile ... -phasefiles ... -go` in "linked"
     mode using a uniform-rate recombfile (1.06 cM/Mb cowpea proxy).
  5. Parse the *_linked_tree.xml MCMC output to extract the cluster
     assignment.
  6. Plot the coancestry heatmap and report the cluster vector.

Output
------
results/37_finestructure/
    data/
        ayb.ids, Ss*.phase, Ss*.recombfile, ayb.cp ...
    tables/
        cluster_assignments.csv      sample -> fineSTRUCTURE cluster
        coancestry.csv               n x n haplotype copy matrix
    figures/
        fig81_coancestry_heatmap.png/.pdf
        fig82_finestructure_tree.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PHASED = ROOT / "results" / "36_ihs" / "data" / "ayb.phased.vcf.gz"
OUT = ROOT / "results" / "37_finestructure"
DATA = OUT / "data"
TAB = OUT / "tables"
FIG = OUT / "figures"
for d in (DATA, TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

# Java tool quirks again: keep the project working dir off the spaced path
TMP = Path("/tmp/ayb_fs")
TMP.mkdir(parents=True, exist_ok=True)

FS_BIN = str(Path.home() / "miniconda3" / "envs" / "popgen_env" / "bin" / "fs")
CM_PER_MB = 1.06


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Parse the phased VCF
# ---------------------------------------------------------------------------
print("[load] phased VCF...")
samples = None
records = {}        # chrom -> list of (pos, [allele_a, allele_b, ...])
with gzip.open(PHASED, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            cols = line.rstrip().split("\t")
            samples = cols[9:]
            continue
        f = line.rstrip().split("\t")
        chrom = f[0]; pos = int(f[1])
        gts = f[9:]
        # Phased GT looks like "0|1"; we want the two haplotypes per sample.
        hap_a, hap_b = [], []
        for g in gts:
            a, b = g.replace(".", "0").replace("/", "|").split("|")
            hap_a.append(int(a)); hap_b.append(int(b))
        records.setdefault(chrom, []).append((pos, hap_a, hap_b))

n_samp = len(samples)
n_hap = 2 * n_samp
print(f"[load] {n_samp} samples; "
      f"chromosomes = {len(records)}; total markers = "
      f"{sum(len(v) for v in records.values())}")


# ---------------------------------------------------------------------------
# Write ChromoPainter .phase files (1 per chromosome) and .recombfile
# ChromoPainter format (haploid version):
#   line 1: N (number of haplotypes)
#   line 2: S (number of SNPs)
#   line 3: P pos1 pos2 ... posS
#   lines 4..(N+3): haplotypes as concatenated 0/1 strings of length S
# ---------------------------------------------------------------------------
phase_files = []
recomb_files = []
chrs = sorted(records.keys(),
              key=lambda s: int(re.search(r"\d+", s).group()))
for c in chrs:
    rows = sorted(records[c], key=lambda r: r[0])
    positions = [r[0] for r in rows]
    S = len(positions)
    hap_strings = ["" for _ in range(n_hap)]
    for (_pos, hap_a, hap_b) in rows:
        for i in range(n_samp):
            hap_strings[2 * i]     += str(hap_a[i])
            hap_strings[2 * i + 1] += str(hap_b[i])
    pfile = TMP / f"{c}.phase"
    with open(pfile, "w") as fh:
        fh.write(f"{n_hap}\n{S}\nP " + " ".join(str(p) for p in positions) + "\n")
        for h in hap_strings:
            fh.write(h + "\n")
    phase_files.append(pfile)
    # ChromoPainter recombfile: 2 cols + header
    # `start.pos`  `recom.rate.perbp`
    # Last marker has rate -9.
    rfile = TMP / f"{c}.recombfile"
    morgan_per_bp = CM_PER_MB / 100.0 / 1e6   # cM/Mb -> Morgans/bp
    with open(rfile, "w") as fh:
        fh.write("start.pos\trecom.rate.perbp\n")
        for i, p in enumerate(positions):
            if i == len(positions) - 1:
                fh.write(f"{p}\t-9\n")
            else:
                fh.write(f"{p}\t{morgan_per_bp:.4e}\n")
    recomb_files.append(rfile)
print(f"[phase] wrote {len(phase_files)} phase + recomb files to {TMP}")

# ID file: <sample> <pop=NULL> <include=1>
id_file = TMP / "ayb.ids"
with open(id_file, "w") as fh:
    for s in samples:
        # use a safe label (alphanumeric + underscore)
        safe = re.sub(r"[^A-Za-z0-9_]", "_", s)
        fh.write(f"{safe}\t1\t1\n")


# ---------------------------------------------------------------------------
# Build + run the fineSTRUCTURE project
# Linked mode with one phase file + one recombfile per chromosome.
# `-go` makes fs run all stages: cp painting, chromocombine, fs MCMC, fs tree.
# ---------------------------------------------------------------------------
proj = TMP / "ayb_fs.cp"
# If the previous pipeline already produced the combined coancestry +
# tree, skip the long re-run and go straight to parsing.
combined = TMP / "ayb_fs_linked.chunkcounts.out"
tree_xml = TMP / "ayb_fs_linked_tree.xml"
SKIP_FS = combined.exists() and tree_xml.exists()
if SKIP_FS:
    print(f"[fs] skipping re-run — found existing {combined.name} + "
          f"{tree_xml.name}")
else:
    if proj.exists():
        proj.unlink()
    for sub in TMP.iterdir():
        if sub.is_dir() and sub.name.startswith("ayb_fs"):
            shutil.rmtree(sub)

# fs concatenates absolute log paths with stage1/ and writes to an invalid
# nested path. Workaround: pass *relative* filenames after cd'ing into TMP.
phase_arg  = ",".join(p.name for p in phase_files)
recomb_arg = ",".join(p.name for p in recomb_files)

if not SKIP_FS:
    print("[fs] stage 1: project init...")
    r = subprocess.run([FS_BIN, "ayb_fs.cp", "-n",
                        "-idfile", id_file.name,
                        "-phasefiles", phase_arg,
                        "-recombfiles", recomb_arg,
                        "-hpc", "0",
                        "-ploidy", "2",
                        "-go"],
                       capture_output=True, text=True, timeout=14400,
                       cwd=str(TMP))
    print("--- fs stdout tail ---"); print(r.stdout[-1500:])
    if r.returncode != 0:
        print("--- fs stderr tail ---"); print(r.stderr[-1500:])

    # fs `-go` reports planned commands but may need repeated invocation to
    # advance through stages on small projects. Loop up to 12 times.
    for stage_pass in range(12):
        r = subprocess.run([FS_BIN, "ayb_fs.cp", "-go"],
                           capture_output=True, text=True, timeout=14400,
                           cwd=str(TMP))
        tail = r.stdout[-600:]
        print(f"--- fs pass {stage_pass+1} stdout tail ---"); print(tail)
        if r.returncode != 0:
            print("--- fs pass stderr tail ---"); print(r.stderr[-600:])
        done_markers = ("All actions completed", "All stages complete",
                        "Project is complete", "Nothing to be done")
        if any(m in r.stdout for m in done_markers):
            print("[fs] pipeline reports complete.")
            break


# ---------------------------------------------------------------------------
# Parse outputs:
#   <proj>.chunkcounts.out   coancestry / chunk-count matrix (full)
#   <proj>_linked_tree.xml   MCMC tree with cluster assignments
# ---------------------------------------------------------------------------
# Combined coancestry produced by chromocombine sits directly under TMP.
chunk_path = TMP / "ayb_fs_linked.chunkcounts.out"
xml_path   = TMP / "ayb_fs_linked_tree.xml"
mcmc_path  = TMP / "ayb_fs_linked_mcmc.xml"
print(f"[parse] chunk_path = {chunk_path} (exists={chunk_path.exists()})")
print(f"[parse] tree xml  = {xml_path} (exists={xml_path.exists()})")

cluster_df = pd.DataFrame()

if chunk_path.exists() and chunk_path.stat().st_size > 0:
    # chunkcounts.out:
    #   line 1: "#Cfactor 0.0276..."
    #   line 2: "Recipient TSs156A TSs361 ..."     (space-delimited)
    #   line 3..: "<recipient> count count count..."
    with open(chunk_path) as fh:
        _cfactor = fh.readline()
        header = fh.readline().split()
        # the first column header is "Recipient"; the rest are donor IDs
        donor_ids = header[1:]
        rows = []
        recipients = []
        for line in fh:
            f = line.split()
            if not f:
                continue
            recipients.append(f[0])
            rows.append([float(x) for x in f[1:]])
    cc = pd.DataFrame(rows, index=recipients, columns=donor_ids)
    cc.to_csv(TAB / "coancestry.csv")
    print(f"[parse] coancestry matrix: {cc.shape}")

    # Reorder rows + columns by fineSTRUCTURE tree leaf order so the
    # cluster block-diagonal structure becomes visible. Without this the
    # heatmap looks uniformly orange and conveys no clustering.
    tree_order = None
    if xml_path.exists() and xml_path.stat().st_size > 0:
        m = re.search(r"<Tree>(.*?)</Tree>", xml_path.read_text(),
                       flags=re.DOTALL)
        if m:
            newick = m.group(1).strip()
            try:
                from Bio import Phylo
                import io
                tr = Phylo.read(io.StringIO(newick), "newick")
                tree_order = [tip.name for tip in tr.get_terminals()
                              if tip.name]
            except Exception:
                tree_order = None
    if tree_order is not None:
        # only keep samples present in cc; drop any tree-order names that
        # aren't in the matrix (the underscore-safe IDs may differ slightly)
        keep = [s for s in tree_order if s in cc.index and s in cc.columns]
        if len(keep) >= cc.shape[0] // 2:
            cc_ord = cc.loc[keep, keep]
        else:
            cc_ord = cc
    else:
        cc_ord = cc

    # ----- Fig 81. coancestry heatmap (log-transformed, tree-ordered) -----
    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    M = cc_ord.values.astype(float)
    # log1p compresses the diagonal vs. off-diagonal contrast (self-copy is
    # always large in CP output) so we can see the off-diagonal structure
    im = ax.imshow(np.log1p(M), cmap="magma", aspect="equal")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("log(1 + chunkcount)")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(f"Donor (n={cc_ord.shape[1]}, fs-tree order)")
    ax.set_ylabel(f"Recipient (n={cc_ord.shape[0]}, fs-tree order)")
    ax.set_title("ChromoPainter coancestry — AYB 95-line panel")
    save(fig, "fig81_coancestry_heatmap")
    print("[fig] fig81_coancestry_heatmap (tree-ordered)")
else:
    print("[warn] no chunkcounts.out found — fs pipeline did not finish")
    cc = pd.DataFrame()

# ----- parse the MAP cluster assignments from the tree XML -----
# fineSTRUCTURE encodes Pop assignments as
#   <Pop>(s1,s2)(s3)(s4,s5,s6)...</Pop>
# in the <Iteration> block of the tree XML. Find the LAST iteration (final
# MAP partition) and parse parenthetical groups.
if xml_path.exists() and xml_path.stat().st_size > 0:
    txt = xml_path.read_text()
    pop_blocks = re.findall(r"<Pop>(.*?)</Pop>", txt, flags=re.DOTALL)
    K_match = re.search(r"<K>(\d+)</K>", txt)
    K_value = int(K_match.group(1)) if K_match else None
    if pop_blocks:
        last_pop = pop_blocks[-1].strip()
        # split on `)(`, then strip the leading `(` and trailing `)`
        groups = re.findall(r"\(([^)]+)\)", last_pop)
        rows = []
        for idx, g in enumerate(groups, start=1):
            for s in re.split(r"[,\s]+", g):
                if s:
                    rows.append({"sample": s, "cluster": f"K{idx:02d}"})
        cluster_df = pd.DataFrame(rows)
        cluster_df.to_csv(TAB / "cluster_assignments.csv", index=False)
        print(f"[parse] K = {K_value}; "
              f"{cluster_df['cluster'].nunique()} clusters parsed; "
              f"{cluster_df.shape[0]} samples assigned")

        # ----- Fig 82. dendrogram of fineSTRUCTURE tree -----
        # Newick is in the <Tree>...</Tree> block.
        m = re.search(r"<Tree>(.*?)</Tree>", txt, flags=re.DOTALL)
        if m:
            newick = m.group(1).strip()
            newick_path = DATA / "ayb_fs_linked.newick"
            newick_path.write_text(newick + "\n")
            # let ete3 or biopython draw it; biopython is already in the env
            try:
                from Bio import Phylo
                import io
                tr = Phylo.read(io.StringIO(newick), "newick")
                fig, ax = plt.subplots(figsize=(8, 14),
                                       constrained_layout=True)
                Phylo.draw(tr, axes=ax, do_show=False,
                           branch_labels=lambda c: "",
                           label_func=lambda c: c.name if c.is_terminal()
                                                       else "")
                ax.set_title(f"fineSTRUCTURE coancestry tree "
                             f"(MAP K = {K_value})")
                save(fig, "fig82_finestructure_tree")
                print("[fig] fig82_finestructure_tree")
            except Exception as e:
                print(f"[warn] tree draw failed: {e}")

# Mirror the whole /tmp tree into OUT/data for the record
print(f"[copy] mirroring /tmp/ayb_fs into {DATA}...")
for f in TMP.glob("*"):
    if f.is_file():
        try:
            shutil.copy(f, DATA / f.name)
        except Exception:
            pass

print(f"\nOutputs in: {OUT}")
