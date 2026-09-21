#!/usr/bin/env python3
"""
Beagle 5 phasing + selscan iHS extended-haplotype-homozygosity scan for
recent positive selection (Voight et al. 2006 PLoS Biol).

Method:
  1. Build a VCF for the AYB-anchored markers (95 diploid samples).
  2. Phase with Beagle 5.5 (Browning et al. 2018).
  3. Per-chromosome selscan iHS run with a uniform 1.06 cM/Mb genetic map
     (Lonardi 2019 cowpea proxy).
  4. Combine + standardise iHS scores across chromosomes (selscan norm).
  5. Render a Manhattan-style |iHS| > 2 selection-signal plot.

Output:
  results/36_ihs/
    data/
      ayb.vcf.gz            unphased input
      ayb.phased.vcf.gz     phased
      *.ihs.out             per-chromosome selscan output
    tables/
      ihs_combined.csv      standardised iHS per SNP
    figures/
      fig80_ihs_manhattan.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
OUT = ROOT / "results" / "36_ihs"
DATA = OUT / "data"
TAB = OUT / "tables"
FIG = OUT / "figures"
for d in (DATA, TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
CM_PER_MB = 1.06   # Lonardi 2019 cowpea proxy
POPGEN_BIN = Path.home() / "miniconda3" / "envs" / "popgen_env" / "bin"
BEAGLE = str(POPGEN_BIN / "beagle")
SELSCAN = str(POPGEN_BIN / "selscan")
NORM = str(POPGEN_BIN / "norm")

# Java tools (Beagle) cannot tolerate spaces in CLI args. Operate from a
# space-free /tmp working directory and mirror results back into OUT/data.
TMP = Path("/tmp/ayb_ihs")
TMP.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def write_vcf(samples, snps, dosage_matrix, vcf_path):
    """Write a basic VCF (sample x SNP); diploid GT field; dosage rounded."""
    with gzip.open(vcf_path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        fh.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        for c in sorted(set(snps["chr"])):
            fh.write(f"##contig=<ID={c}>\n")
        cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                "INFO", "FORMAT"] + samples
        fh.write("\t".join(cols) + "\n")
        for j in range(len(snps)):
            chrom = snps["chr"].iloc[j]
            pos = int(snps["pos"].iloc[j])
            rs = snps["rs"].iloc[j]
            ref = snps["ref"].iloc[j]; alt = snps["alt"].iloc[j]
            geno_strs = []
            for i in range(len(samples)):
                d = dosage_matrix[i, j]
                if np.isnan(d):
                    geno_strs.append("./.")
                else:
                    d_int = int(round(d))
                    if d_int == 0:
                        geno_strs.append("0/0")
                    elif d_int == 1:
                        geno_strs.append("0/1")
                    else:
                        geno_strs.append("1/1")
            row = [chrom, str(pos), rs, ref, alt, ".", "PASS", ".",
                   "GT"] + geno_strs
            fh.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Build VCF from filtered + AYB-anchored dosage
# ---------------------------------------------------------------------------
print("[load] HapMap + filtered dosage...")
hm = pd.read_csv(HAPMAP, low_memory=False)
META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
        "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META]
calls = hm[sample_cols].astype(str)
ref_alt = hm["alleles"].str.split("/", expand=True)
ref_alt.columns = ["ref", "alt"]
dosage = np.full((len(hm), len(sample_cols)), np.nan)
for i in range(len(hm)):
    ref, alt = ref_alt.iloc[i]
    arr = calls.iloc[i].values
    dosage[i, arr == ref + ref] = 0.0
    dosage[i, (arr == ref + alt) | (arr == alt + ref)] = 1.0
    dosage[i, arr == alt + alt] = 2.0
dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)
call_rate_m = dosage.notna().sum(axis=1) / dosage.shape[1]
maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                 1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]
keep_m = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values
dosage = dosage.loc[keep_m, keep_s]
ref_alt_f = ref_alt.loc[keep_m]
ref_alt_f.index = dosage.index

# AYB-anchored only
anc = pd.read_csv(ANCHOR).dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
anc["snp_pos_ayb"] = anc["snp_pos_ayb"].astype(int)
anc = anc[anc["rs"].isin(dosage.index)].copy()
anc = anc.sort_values(["chr_ayb", "snp_pos_ayb"]).reset_index(drop=True)

snps = pd.DataFrame({
    "rs": anc["rs"].values,
    "chr": anc["chr_ayb"].values,
    "pos": anc["snp_pos_ayb"].values,
    "ref": ref_alt_f.loc[anc["rs"], "ref"].values,
    "alt": ref_alt_f.loc[anc["rs"], "alt"].values,
})
# Beagle 5 rejects duplicate (chrom, pos). Keep the first record per locus
# and report how many we drop so the count is auditable in the log.
n_pre = len(snps)
snps = snps.drop_duplicates(subset=["chr", "pos"], keep="first").reset_index(drop=True)
n_post = len(snps)
if n_pre != n_post:
    print(f"[dedup] dropped {n_pre - n_post} duplicate-position markers "
          f"({n_pre} -> {n_post})")
samples = dosage.columns.tolist()
# transpose to samples x snps, aligned to the deduplicated rs order
D = dosage.loc[snps["rs"]].T.values
print(f"[load] {D.shape[0]} samples x {D.shape[1]} markers")

vcf_path = TMP / "ayb.vcf.gz"
print(f"[vcf] writing {vcf_path}...")
write_vcf(samples, snps, D, vcf_path)
print(f"[vcf] {vcf_path.stat().st_size/1024:.0f} KB written")
# keep a copy in OUT/data for the record
shutil.copy(vcf_path, DATA / "ayb.vcf.gz")


# ---------------------------------------------------------------------------
# Phase with Beagle 5 (operates entirely under /tmp/ayb_ihs to dodge spaces)
# ---------------------------------------------------------------------------
phased_tmp = TMP / "ayb.phased.vcf.gz"
print(f"[beagle] phasing → {phased_tmp}...")
r = subprocess.run([BEAGLE,
                    f"gt={vcf_path}",
                    f"out={TMP/'ayb.phased'}",
                    "nthreads=4"],
                   capture_output=True, text=True, timeout=3600,
                   cwd=str(TMP))
print("--- beagle stdout tail ---"); print(r.stdout[-500:])
if r.returncode != 0:
    print("--- beagle stderr tail ---"); print(r.stderr[-500:])
    raise SystemExit("Beagle failed")
phased_path = DATA / "ayb.phased.vcf.gz"
shutil.copy(phased_tmp, phased_path)


# ---------------------------------------------------------------------------
# Per-chromosome selscan iHS run
# ---------------------------------------------------------------------------
print("[selscan] per-chromosome iHS...")
# selscan accepts plain VCF (gunzipped); we decompress per-chrom
import gzip as _gz

chrom_files = {}
with _gz.open(phased_path, "rt") as fh:
    header = []
    for line in fh:
        if line.startswith("#"):
            header.append(line); continue
        c = line.split("\t", 1)[0]
        chrom_files.setdefault(c, []).append(line)

n_total = 0
for c, lines in chrom_files.items():
    out_vcf = TMP / f"{c}.phased.vcf"
    with open(out_vcf, "w") as fh:
        fh.writelines(header); fh.writelines(lines)
    # genetic map: chrom, pos_bp, cM, recomb_rate (cM/Mb)
    map_path = TMP / f"{c}.map"
    map_rows = []
    for ln in lines:
        f = ln.split("\t")
        pos = int(f[1])
        cm = pos * CM_PER_MB / 1e6
        map_rows.append(f"{c}\t.\t{cm:.6f}\t{pos}")
    map_path.write_text("\n".join(map_rows) + "\n")
    # selscan iHS — run from /tmp working dir.
    # DArTseq density on AYB ~500 kb/marker, much sparser than the WGS
    # 1-5 kb density selscan defaults assume. Relax --max-gap to 5 Mb and
    # --max-extend to 20 Mb so EHH can decay before selscan gives up. Limit
    # of statistical interpretation is reported in the methods.
    print(f"  {c}: {len(lines)} variants...")
    r = subprocess.run(
        [SELSCAN, "--ihs", "--vcf", str(out_vcf), "--map", str(map_path),
         "--out", str(TMP / f"{c}.ihs"), "--threads", "4",
         "--max-gap", "5000000",
         "--max-extend", "20000000"],
        capture_output=True, text=True, timeout=1800,
        cwd=str(TMP))
    if r.returncode != 0:
        print(f"    selscan failed: {r.stderr[-300:]}")
        continue
    n_total += len(lines)

print(f"[selscan] total variants processed: {n_total}")


# ---------------------------------------------------------------------------
# Normalise iHS across chromosomes
# ---------------------------------------------------------------------------
ihs_files = sorted(TMP.glob("*.ihs.out"))
print(f"[norm] normalising {len(ihs_files)} per-chromosome files...")
if ihs_files:
    # selscan's `norm` tool standardises iHS within MAF bins, pooled across
    # input files. Output is files named <orig>.100bins.norm
    r = subprocess.run([NORM, "--ihs", "--files",
                        *(str(f) for f in ihs_files)],
                       capture_output=True, text=True, timeout=1800,
                       cwd=str(TMP))
    print("--- norm stdout tail ---"); print(r.stdout[-400:])
    # mirror per-chrom and normalised files back into OUT/data
    for f in TMP.glob("*.ihs.out*"):
        shutil.copy(f, DATA / f.name)
    norm_files = sorted(TMP.glob("*.ihs.out.100bins.norm"))
    print(f"[norm] {len(norm_files)} normalised files")
    rows = []
    for f in norm_files:
        if f.stat().st_size == 0: continue
        d = pd.read_csv(f, sep="\t", header=None,
                        names=["rs", "pos", "freq", "ihh1", "ihh0",
                               "uniH_or_iHS", "norm_iHS", "crit"])
        c = re.match(r"(Ss\d+)\.", f.name).group(1)
        d["chr"] = c
        rows.append(d)
    if rows:
        ihs = pd.concat(rows, ignore_index=True)
        ihs.to_csv(TAB / "ihs_combined.csv", index=False)
        sig = (ihs["norm_iHS"].abs() > 2).sum()
        print(f"[norm] {len(ihs):,} markers with normalised iHS; "
              f"{sig} with |iHS| > 2")
    else:
        ihs = pd.DataFrame()
else:
    ihs = pd.DataFrame()


# ---------------------------------------------------------------------------
# Manhattan plot of |iHS|
# ---------------------------------------------------------------------------
if not ihs.empty:
    chrs = sorted(ihs["chr"].unique(),
                  key=lambda s: int(re.search(r"\d+", s).group()))
    offsets, mids, x_cursor = {}, [], 0
    xs = []
    ihs = ihs.sort_values(["chr", "pos"]).reset_index(drop=True)
    for c in chrs:
        sub = ihs[ihs["chr"] == c]
        offsets[c] = x_cursor
        xs.extend((sub["pos"].values + x_cursor).tolist())
        mids.append(x_cursor + (sub["pos"].max() - sub["pos"].min()) / 2)
        x_cursor += sub["pos"].max() + 5_000_000
    ihs["x"] = xs

    fig, ax = plt.subplots(figsize=(11, 3.8), constrained_layout=True)
    palette = [WONG["blue"], "#7c7c7c"]
    for i, c in enumerate(chrs):
        sub = ihs[ihs["chr"] == c]
        ax.scatter(sub["x"], sub["norm_iHS"].abs(), s=10,
                   color=palette[i % 2], alpha=0.7, edgecolor="none")
    ax.axhline(2, color=WONG["vermillion"], linewidth=0.7, linestyle="--",
               label="|iHS| = 2")
    ax.set_xticks(mids); ax.set_xticklabels(chrs, fontsize=8)
    ax.set_ylabel("|standardised iHS|")
    ax.set_title(f"selscan iHS scan ({len(ihs):,} markers; uniform "
                 f"{CM_PER_MB} cM/Mb map)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, axis="y", alpha=0.4)
    save(fig, "fig80_ihs_manhattan")
    print("[fig] fig80_ihs_manhattan")

print(f"\nOutputs in: {OUT}")
