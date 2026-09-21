#!/usr/bin/env python3
"""
Folded site-frequency spectrum + Stairway Plot 2 (Liu & Fu 2020 Nat Commun)
demographic inference for the AYB panel.

Stairway Plot 2 reconstructs Ne(t) from the SFS without requiring an
outgroup (folded SFS is used). Pair with the LD-derived current Ne = 659 in
script 11 for a two-time-scale view.

This driver:
  1. Build the folded SFS from the AYB-anchored dosage matrix.
  2. Write a Stairway Plot blueprint file with sensible defaults.
  3. Build the StairwayPlot Java pipeline scripts.
  4. Run the Stairway Plot Java executor.
  5. Plot the resulting Ne(t) curve from the .summary output.

Outputs (results/33_stairway/)
------------------------------
data/
    sfs.txt                 folded SFS used as input
    ayb.blueprint           Stairway Plot config
    stairway_run/           Stairway Plot working directory
tables/
    ne_through_time.csv     Stairway Plot summary curve
figures/
    fig79_stairway_ne_time.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import os
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
STAIRWAY = ROOT / "refs" / "tools" / "stairway_plot_v2.1.2" / "stairway_plot_es"
OUT = ROOT / "results" / "33_stairway"
DATA = OUT / "data"
TAB = OUT / "tables"
FIG = OUT / "figures"
for d in (DATA, TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
GEN_TIME_YEARS = 1.0           # one-year cycle for AYB landraces (a defensible
                                # default; can be parameterised later)
MUTATION_RATE = 7e-9            # per-site, per-generation (Ossowski et al. 2010
                                # for Arabidopsis — common legume proxy)
TOTAL_GENOME_BP = 649_800_000   # AYB chromosome-scale assembly span


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Build filtered dosage + folded SFS
# ---------------------------------------------------------------------------
print("[load] HapMap...")
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

# AYB-anchored only — they're the unbiased subset, no cowpea proxy bias
anc = pd.read_csv(ANCHOR).dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
dosage = dosage.loc[dosage.index.intersection(anc["rs"])]
n_samples = dosage.shape[1]
n_alleles = 2 * n_samples
print(f"[load] {dosage.shape[0]} markers x {n_samples} samples "
      f"({n_alleles} alleles)")


# ---------------------------------------------------------------------------
# Folded SFS
# ---------------------------------------------------------------------------
print("[sfs] computing folded SFS...")
# alt-allele count per marker (treat each sample as 2 alleles; round any
# imputed values back to integer)
ct = dosage.fillna(0).round().astype(int).sum(axis=1).values
ct = np.clip(ct, 0, n_alleles)
# fold to minor-allele count
fold = np.minimum(ct, n_alleles - ct)
# folded SFS has bins 1..(n_alleles/2)
max_minor = n_alleles // 2
sfs = np.zeros(max_minor + 1, dtype=int)
for m in fold:
    if 1 <= m <= max_minor:
        sfs[m] += 1
print(f"[sfs] folded SFS sum = {sfs.sum()} segregating sites")

# write SFS as space-separated counts (Stairway Plot format: bins 1..max
# without the monomorphic 0 bin)
sfs_path = DATA / "sfs.txt"
with open(sfs_path, "w") as fh:
    fh.write(" ".join(str(int(x)) for x in sfs[1:]) + "\n")
print(f"[sfs] wrote {sfs_path}")


# ---------------------------------------------------------------------------
# Stairway Plot blueprint
# Format reference: stairway_plot_v2.1.2/two-epoch_fold.blueprint
# ---------------------------------------------------------------------------
blueprint_path = DATA / "ayb.blueprint"
n_seq = n_alleles
L = TOTAL_GENOME_BP
# nrand: breakpoint values used by Stairway. Recommended:
# floor((n_seq-2)/4), floor((n_seq-2)/2), floor((n_seq-2)*3/4), n_seq-2.
nrand = [(n_seq - 2) // 4, (n_seq - 2) // 2, 3 * (n_seq - 2) // 4, n_seq - 2]
project_dir = str(DATA / "stairway_run").rstrip("/")

# ensure project_dir is clean
import shutil
if Path(project_dir).exists():
    shutil.rmtree(project_dir)

with open(blueprint_path, "w") as fh:
    fh.write(f"""#input setting
popid: AYB_95
nseq: {n_seq}
L: {L}
whether_folded: true
SFS: {" ".join(str(int(x)) for x in sfs[1:])}
#smallest_size_of_SFS_bin_used_for_estimation: 1
#largest_size_of_SFS_bin_used_for_estimation: {max_minor}
pct_training: 0.67
nrand: {" ".join(str(int(x)) for x in nrand)}
project_dir: {project_dir}
stairway_plot_dir: {STAIRWAY}
ninput: 200
#random_seed: 6
#output setting
mu: {MUTATION_RATE}
year_per_generation: {GEN_TIME_YEARS}
#plot setting
plot_title: AYB_95
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12
""")
print(f"[bp] wrote {blueprint_path}")


# ---------------------------------------------------------------------------
# Run Stairway Builder (generates per-job shell script) + Stairway plot
# Both steps are Java; the JAR ships in stairway_plot_es/
#
# Java's classpath resolution cannot tolerate spaces in the path. Skip the
# in-project Java pipeline if the off-tree /tmp run has already produced
# the final summary; in that case the result is reused as-is.
# ---------------------------------------------------------------------------
TMP_PROJECT = Path("/tmp/ayb_stairway/stairway_data/stairway_run")
TMP_SUMMARY = TMP_PROJECT / "AYB_95.final.summary"

if TMP_SUMMARY.exists():
    print(f"[run] skipping Java pipeline — found {TMP_SUMMARY}")
else:
    java_classpath = str(STAIRWAY)
    print(f"[run] Stairway Plot 2 (Java)...")
    env = os.environ.copy()
    env["PATH"] = "/usr/bin:" + env.get("PATH", "")

    # Step 1: builder turns blueprint into a per-experiment shell script
    r = subprocess.run(
        ["java", "-cp", java_classpath, "Stairbuilder", str(blueprint_path)],
        capture_output=True, text=True, env=env)
    print("--- builder stdout tail ---")
    print(r.stdout[-600:])
    if r.returncode != 0:
        print("--- builder stderr ---"); print(r.stderr[-600:])
        raise SystemExit("Stairbuilder failed")

    # Step 2: run the per-experiment shell script
    shell_script = Path(f"{blueprint_path}.sh")
    if not shell_script.exists():
        raise SystemExit(f"Stairbuilder did not create {shell_script}")
    print(f"[run] executing {shell_script}...")
    r2 = subprocess.run(["bash", str(shell_script)],
                        capture_output=True, text=True,
                        cwd=DATA, env=env, timeout=3600)
    print("--- run stdout tail ---")
    print(r2.stdout[-600:])
    if r2.returncode != 0:
        print("--- run stderr ---"); print(r2.stderr[-600:])

# ---------------------------------------------------------------------------
# Parse summary + plot Ne(t) ourselves (matplotlib version of stairway's PDF)
# ---------------------------------------------------------------------------
# Stairway Plot's Java pipeline cannot tolerate spaces in classpath paths,
# so the actual pipeline was driven from /tmp/ayb_stairway/. Look there
# first for the .final.summary; fall back to the in-project location.
TMP_PROJECT = Path("/tmp/ayb_stairway/stairway_data/stairway_run")
summary_files = list(TMP_PROJECT.glob("*final.summary"))
if not summary_files:
    summary_files = list(Path(project_dir).glob("*final.summary"))
print(f"[parse] found summary files: {[str(s) for s in summary_files]}")
if summary_files:
    summary = pd.read_csv(summary_files[0], sep="\t")
    summary.to_csv(TAB / "ne_through_time.csv", index=False)
    print(f"[parse] {len(summary)} time points")

    # Stairway's summary file leads with thousands of rows where the
    # cumulative mutation rate (and therefore year) is essentially zero
    # (down to ~10^-310). These rows compress the plot into illegible
    # range. Filter to year >= 1 generation before plotting.
    plot_df = summary[summary["year"] >= 1.0].copy()
    print(f"[plot] {len(plot_df)} rows with year >= 1 yr")
    print(plot_df.head().to_string())

    fig, ax = plt.subplots(figsize=(8, 5.0), constrained_layout=True)
    ax.plot(plot_df["year"], plot_df["Ne_median"],
            color=WONG["blue"], linewidth=1.8, label="median Ne")
    if "Ne_2.5%" in plot_df.columns and "Ne_97.5%" in plot_df.columns:
        ax.fill_between(plot_df["year"],
                         plot_df["Ne_2.5%"], plot_df["Ne_97.5%"],
                         color=WONG["blue"], alpha=0.18,
                         label="95 % CI")
    ax.axhline(659, color=WONG["vermillion"], linewidth=0.8,
                linestyle="--", alpha=0.7,
                label="LD-derived current Ne = 659")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Years before present "
                  "(μ = 7e-9 / bp / yr, generation = 1 yr, L = 112 kb DArT)")
    ax.set_ylabel(r"Effective population size $N_e$")
    ax.set_title("AYB demographic history — Stairway Plot 2")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, which="both", alpha=0.4)
    save(fig, "fig79_stairway_ne_time")
    print("[fig] fig79_stairway_ne_time")
else:
    print("[warn] no summary file found")

print(f"\nOutputs in: {OUT}")
