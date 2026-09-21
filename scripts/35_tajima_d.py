#!/usr/bin/env python3
"""
Windowed Tajima's D, nucleotide diversity (π), and Watterson's θ across the
AYB-anchored genome.

Method
------
For each chromosome (Ss01-Ss11), sliding 250 kb windows with 100 kb step.
Per window with ≥ 4 segregating SNPs:
  - π = mean pairwise nucleotide differences per site
  - θ_W = S / a1 where a1 = Σ_{i=1}^{n-1} 1/i (Watterson 1975)
  - Tajima's D = (π − θ_W) / sqrt(Var)  (Tajima 1989)

Per-chromosome Manhattan-style plot of D values; significant D < -2 or > +2
flagged (rough selection / balancing signal heuristic).

Outputs
-------
results/35_tajima_d/
  tables/
    tajima_d_windowed.csv    chrom, mid_bp, n_snps, pi, theta_W, tajima_D
  figures/
    fig77_tajima_d_manhattan.png/.pdf
    fig78_pi_vs_theta.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
OUT = ROOT / "results" / "35_tajima_d"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
WINDOW_BP = 250_000
STEP_BP = 100_000
MIN_SNPS = 4


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def tajima_D(seg_S: int, pi: float, n: int) -> float:
    """Tajima 1989 D statistic for n samples, S segregating sites, π per-site."""
    if seg_S < 1 or n < 2:
        return np.nan
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / i**2 for i in range(1, n))
    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = 2.0 * (n**2 + n + 3.0) / (9.0 * n * (n - 1.0))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2.0) / (a1 * n) + a2 / a1**2
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    theta_W = seg_S / a1
    var = e1 * seg_S + e2 * seg_S * (seg_S - 1.0)
    if var <= 0:
        return np.nan
    return (pi - theta_W) / np.sqrt(var)


# ---------------------------------------------------------------------------
# Build filtered dosage + AYB-anchored positions
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

anc = pd.read_csv(ANCHOR).dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
anc["snp_pos_ayb"] = anc["snp_pos_ayb"].astype(int)
anc = anc[anc["rs"].isin(dosage.index)]
print(f"[load] AYB-anchored markers in panel: {len(anc)}")


# ---------------------------------------------------------------------------
# Per-window Tajima's D, π, θ_W
# Treat each diploid sample as 2 alleles (count_alt = dosage); MAF-filtered.
# π per site = (2 * count_alt * (n_alleles − count_alt)) /
#              (n_alleles * (n_alleles − 1))
# Sum over sites within window -> π_window; divide by window-bp for per-site.
# ---------------------------------------------------------------------------
print(f"[scan] sliding {WINDOW_BP/1000:.0f} kb / {STEP_BP/1000:.0f} kb...")
rows = []
chrs = sorted(anc["chr_ayb"].unique(),
              key=lambda s: int(re.search(r"\d+", s).group()))

for c in chrs:
    sub = anc[anc["chr_ayb"] == c].sort_values("snp_pos_ayb").reset_index(drop=True)
    if len(sub) < MIN_SNPS:
        continue
    rs_arr = sub["rs"].values
    pos_arr = sub["snp_pos_ayb"].values
    # extract dosage subset for this chromosome
    D = dosage.loc[rs_arr]
    # mean-impute NA for π computation
    col_mean = D.mean(axis=1, skipna=True)
    D = D.apply(lambda row: row.fillna(col_mean[row.name]), axis=1)
    n_samples = D.shape[1]
    n_alleles = 2 * n_samples  # treat each sample as 2 alleles
    chrom_start = pos_arr.min()
    chrom_end = pos_arr.max()
    for win_start in range(chrom_start, chrom_end + 1, STEP_BP):
        win_end = win_start + WINDOW_BP
        in_win = (pos_arr >= win_start) & (pos_arr < win_end)
        if in_win.sum() < MIN_SNPS:
            continue
        win_D = D.iloc[in_win.nonzero()[0]]
        # alt-allele count per site
        count_alt = win_D.sum(axis=1).values
        # per-site pi (Nei & Li 1979)
        site_pi = np.where(
            n_alleles >= 2,
            (2.0 * count_alt * (n_alleles - count_alt)) /
            (n_alleles * (n_alleles - 1)),
            0.0,
        )
        S = int(in_win.sum())   # segregating sites in window (all are by definition post-MAF filter)
        pi_total = float(site_pi.sum())
        pi_persite = pi_total / WINDOW_BP
        a1 = sum(1.0 / i for i in range(1, n_alleles))
        theta_W = S / a1
        D_stat = tajima_D(S, pi_total, n_alleles)
        rows.append({
            "chr": c,
            "win_start": int(win_start),
            "win_end": int(win_end),
            "mid_bp": int((win_start + win_end) / 2),
            "n_snps": S,
            "pi_total": pi_total,
            "pi_per_site": float(pi_persite),
            "theta_W": float(theta_W),
            "tajima_D": float(D_stat),
        })

df = pd.DataFrame(rows)
df.to_csv(TAB / "tajima_d_windowed.csv", index=False)
print(f"[scan] {len(df)} windows; finite Tajima D in {df['tajima_D'].notna().sum()}")
print(f"  median D: {df['tajima_D'].median():.3f}")
print(f"  D < -2 (selective sweep / expansion): {(df['tajima_D'] < -2).sum()}")
print(f"  D > +2 (balancing / bottleneck):       {(df['tajima_D'] >  2).sum()}")


# ---------------------------------------------------------------------------
# Fig 77. Manhattan-style Tajima's D plot
# ---------------------------------------------------------------------------
df_anc = df.dropna(subset=["tajima_D"]).copy()
chrs_present = sorted(df_anc["chr"].unique(),
                      key=lambda s: int(re.search(r"\d+", s).group()))
offsets, mids, x_cursor = {}, [], 0
xs = []
df_anc = df_anc.sort_values(["chr", "mid_bp"]).reset_index(drop=True)
for c in chrs_present:
    sub = df_anc[df_anc["chr"] == c]
    offsets[c] = x_cursor
    xs.extend((sub["mid_bp"].values + x_cursor).tolist())
    mids.append(x_cursor + (sub["mid_bp"].max() - sub["mid_bp"].min()) / 2)
    x_cursor += sub["mid_bp"].max() + 5_000_000
df_anc["x"] = xs

fig, ax = plt.subplots(figsize=(11, 4.0), constrained_layout=True)
palette = [WONG["blue"], "#7c7c7c"]
for i, c in enumerate(chrs_present):
    sub = df_anc[df_anc["chr"] == c]
    ax.scatter(sub["x"], sub["tajima_D"], s=18,
               color=palette[i % 2], alpha=0.85, edgecolor="none")
ax.axhline(0, color="black", linewidth=0.6)
ax.axhline(-2, color=WONG["vermillion"], linewidth=0.6,
           linestyle="--", alpha=0.7, label="D = -2 sweep threshold")
ax.axhline(+2, color=WONG["green"], linewidth=0.6,
           linestyle="--", alpha=0.7, label="D = +2 balancing threshold")
ax.set_xticks(mids); ax.set_xticklabels(chrs_present, fontsize=8)
ax.set_ylabel("Tajima's D")
ax.set_title(f"Windowed Tajima's D across the AYB-anchored genome "
             f"({len(df_anc)} windows, {int(WINDOW_BP/1000)} kb / "
             f"{int(STEP_BP/1000)} kb)")
ax.legend(loc="upper right", fontsize=8)
ax.grid(True, axis="y", alpha=0.4)
save(fig, "fig77_tajima_d_manhattan")
print("[fig] fig77_tajima_d_manhattan")


# ---------------------------------------------------------------------------
# Fig 78. π vs θ_W per window (the Tajima D's geometry view)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
ax.scatter(df_anc["theta_W"], df_anc["pi_total"], s=24,
           color=WONG["blue"], alpha=0.7, edgecolor="none")
lim = max(df_anc["theta_W"].max(), df_anc["pi_total"].max()) * 1.05
ax.plot([0, lim], [0, lim], "k--", linewidth=0.6, alpha=0.5)
ax.set_xlabel(r"Watterson's $\theta_W$")
ax.set_ylabel(r"Nucleotide diversity $\pi$")
ax.set_title(r"$\pi$ vs $\theta_W$ per window")
ax.grid(True, alpha=0.4)
save(fig, "fig78_pi_vs_theta")
print("[fig] fig78_pi_vs_theta")


print(f"\nOutputs in: {OUT}")
