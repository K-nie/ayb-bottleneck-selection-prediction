#!/usr/bin/env python3
"""
AYB DArTseq diversity panel: marker/sample QC, PCA, clustering, GWAS power analysis.

Inputs
------
- AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv
- Wet_Chemistry_Data.xlsx (Means sheet: Tannin, Phenol, Flavinoid, Antioxidant)

Outputs (under results/01_qc_pca_power/)
----------------------------------------
tables/
  marker_qc.csv              per-SNP call rate, MAF, het, PIC, alleles, chrom, pos
  sample_qc.csv              per-sample call rate, observed het, joinable to phenotype
  filtered_marker_count.txt  filter cascade counts
  pca_coords.csv             PC scores (PC1..PC10) + cluster + phenotype join
  pca_variance.csv           variance explained per PC
  power_curve.csv            GWAS power per (MAF, h2_qtl) grid
phenotype_summary.csv        means/sd/skew/normality per trait + correlation matrix
figures/
  marker_qc_panels.png       MAF, call-rate, het, PIC histograms
  sample_qc_panels.png       sample call-rate, het, scatter
  pca_scree.png              variance-explained scree
  pca_pc1_pc2.png            PC1 vs PC2 colored by k-means cluster
  pca_pc1_pc2_pheno.png      PC1 vs PC2 colored by each phenotype (4-panel)
  power_heatmap.png          power across MAF x h2_qtl grid

Notes
-----
- AYB (Sphenostylis stenocarpa) is highly autogamous; HWE is not enforced.
- HapMap missing call code is NN; observed het = AB / (AA+AB+BB).
- Markers are filtered to: call-rate >= 0.90, MAF >= 0.05. Samples to call-rate >= 0.90.
- Power analysis uses the standard single-marker additive model with genomic-control-style
  Bonferroni threshold; treats sample size, MAF, and per-QTL variance fraction as the axes
  most relevant to a 105-line panel.

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PHENO_XLSX = ROOT / "Wet_Chemistry_Data.xlsx"

OUT = ROOT / "results" / "01_qc_pca_power"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

sns.set_context("talk")
sns.set_style("whitegrid")


# ----------------------------------------------------------------------------
# Load HapMap and convert to dosage / metadata frames
# ----------------------------------------------------------------------------
print("[load] reading HapMap...")
hm = pd.read_csv(HAPMAP, low_memory=False)
META_COLS = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
             "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META_COLS]
print(f"[load] markers={len(hm)}, samples={len(sample_cols)}")

calls = hm[sample_cols].astype(str)


def dosage_row(row: pd.Series, alleles: str) -> np.ndarray:
    """Convert one HapMap marker row to 0/1/2 dosage of the second (alternate) allele.

    HapMap calls are two-letter strings (AA / AG / GG / NN). The `alleles` field
    is `ref/alt` (e.g. C/A). Dosage counts the alt allele.
    """
    ref, alt = alleles.split("/")
    het1, het2 = ref + alt, alt + ref
    hom_ref = ref + ref
    hom_alt = alt + alt
    out = np.full(len(row), np.nan)
    arr = row.values
    out[arr == hom_ref] = 0.0
    out[(arr == het1) | (arr == het2)] = 1.0
    out[arr == hom_alt] = 2.0
    return out


print("[load] converting calls to dosage...")
dosage = np.vstack([dosage_row(calls.iloc[i], hm["alleles"].iloc[i])
                    for i in range(len(hm))])
dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)


# ----------------------------------------------------------------------------
# Per-marker QC
# ----------------------------------------------------------------------------
print("[qc] per-marker QC...")
n_samp = dosage.shape[1]
called = dosage.notna().sum(axis=1)
call_rate_m = called / n_samp
mean_dose = dosage.mean(axis=1, skipna=True)
maf = np.minimum(mean_dose / 2.0, 1.0 - mean_dose / 2.0)
het = (dosage == 1).sum(axis=1) / called.replace(0, np.nan)
# PIC for biallelic locus: 1 - sum(p_i^2) - sum_{i<j} 2 p_i^2 p_j^2
p = mean_dose / 2.0
q = 1.0 - p
pic = 1.0 - (p ** 2 + q ** 2) - (2.0 * (p ** 2) * (q ** 2))

marker_qc = pd.DataFrame({
    "rs": hm["rs#"].values,
    "alleles": hm["alleles"].values,
    "chrom": hm["chrom"].values,
    "pos": hm["pos"].values,
    "call_rate": call_rate_m.values,
    "maf": maf.values,
    "het": het.values,
    "pic": pic.values,
})
marker_qc.to_csv(TAB / "marker_qc.csv", index=False)


# ----------------------------------------------------------------------------
# Per-sample QC
# ----------------------------------------------------------------------------
print("[qc] per-sample QC...")
n_mark = dosage.shape[0]
call_rate_s = dosage.notna().sum(axis=0) / n_mark
het_s = (dosage == 1).sum(axis=0) / dosage.notna().sum(axis=0).replace(0, np.nan)
sample_qc = pd.DataFrame({
    "sample": dosage.columns,
    "call_rate": call_rate_s.values,
    "observed_het": het_s.values,
})
sample_qc.to_csv(TAB / "sample_qc.csv", index=False)


# ----------------------------------------------------------------------------
# Filter and impute
# ----------------------------------------------------------------------------
SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05

keep_samp = sample_qc.loc[sample_qc["call_rate"] >= SAMP_CR, "sample"].tolist()
keep_mark = marker_qc.loc[(marker_qc["call_rate"] >= MARK_CR) &
                          (marker_qc["maf"] >= MIN_MAF), "rs"].tolist()

X = dosage.loc[keep_mark, keep_samp]
print(f"[filter] samples kept {len(keep_samp)} / {n_samp}")
print(f"[filter] markers kept {len(keep_mark)} / {n_mark}")

# mean-impute remaining NA (per marker)
X = X.T  # samples as rows
X = X.fillna(X.mean(axis=0))

with (TAB / "filtered_marker_count.txt").open("w") as f:
    f.write(f"start_markers\t{n_mark}\n")
    f.write(f"start_samples\t{n_samp}\n")
    f.write(f"sample_callrate_>=_{SAMP_CR}\t{len(keep_samp)}\n")
    f.write(f"marker_callrate_>=_{MARK_CR}_and_maf_>=_{MIN_MAF}\t{len(keep_mark)}\n")


# ----------------------------------------------------------------------------
# Phenotype load
# ----------------------------------------------------------------------------
print("[pheno] loading wet chemistry...")
pheno = pd.read_excel(PHENO_XLSX, sheet_name="Means")
pheno = pheno.rename(columns={"Genotypes": "sample", "Flavinoid": "Flavonoid"})
TRAITS = ["Tannin", "Phenol", "Flavonoid", "Antioxidant"]

# normality + summary stats
rows = []
for t in TRAITS:
    v = pheno[t].dropna().values
    sh = stats.shapiro(v)
    rows.append({
        "trait": t, "n": len(v),
        "mean": np.mean(v), "sd": np.std(v, ddof=1),
        "median": np.median(v),
        "skew": stats.skew(v), "kurtosis": stats.kurtosis(v),
        "shapiro_W": sh.statistic, "shapiro_p": sh.pvalue,
    })
pheno_sum = pd.DataFrame(rows)
corr = pheno[TRAITS].corr(method="pearson")
pheno_sum.to_csv(OUT / "phenotype_summary.csv", index=False)
corr.to_csv(OUT / "phenotype_correlation.csv")


# ----------------------------------------------------------------------------
# PCA
# ----------------------------------------------------------------------------
print("[pca] running PCA...")
Xs = StandardScaler().fit_transform(X.values)
pca = PCA(n_components=min(10, Xs.shape[0] - 1, Xs.shape[1]))
pcs = pca.fit_transform(Xs)
var_exp = pd.DataFrame({
    "PC": [f"PC{i+1}" for i in range(pcs.shape[1])],
    "variance_explained": pca.explained_variance_ratio_,
    "cumulative": np.cumsum(pca.explained_variance_ratio_),
})
var_exp.to_csv(TAB / "pca_variance.csv", index=False)

# K-means cluster on PC1..PC5 — use small K range, pick by silhouette
from sklearn.metrics import silhouette_score
sil = {}
for k in range(2, 8):
    km = KMeans(n_clusters=k, n_init=25, random_state=0)
    lab = km.fit_predict(pcs[:, :5])
    sil[k] = silhouette_score(pcs[:, :5], lab)
best_k = max(sil, key=sil.get)
print(f"[pca] best_k by silhouette = {best_k} (scores: {sil})")
km = KMeans(n_clusters=best_k, n_init=50, random_state=0)
clusters = km.fit_predict(pcs[:, :5])

pca_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(pcs.shape[1])])
pca_df["sample"] = X.index.values
pca_df["cluster"] = clusters
pca_df = pca_df.merge(pheno, on="sample", how="left")
pca_df.to_csv(TAB / "pca_coords.csv", index=False)


# ----------------------------------------------------------------------------
# GWAS power analysis
# ----------------------------------------------------------------------------
# Single-marker linear regression. Under additive coding, the per-marker QTL variance
# fraction (h2_qtl) determines the non-centrality parameter:
#   NCP = n * h2_qtl / (1 - h2_qtl)
# Test threshold is Bonferroni-corrected by the effective marker count (we use the
# filtered count). Power is the probability that the chi-square(df=1, NCP) exceeds
# the threshold.

print("[power] computing power curves...")
n = len(keep_samp)
m_eff = len(keep_mark)
alpha = 0.05
thresh_chi2 = stats.chi2.ppf(1 - alpha / m_eff, df=1)
print(f"[power] n={n}, m_eff={m_eff}, alpha={alpha}, Bonferroni chi2 threshold={thresh_chi2:.2f}")

maf_grid = np.array([0.05, 0.10, 0.20, 0.30, 0.40, 0.50])
h2_grid = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50])

power_rows = []
for h2 in h2_grid:
    for mf in maf_grid:
        # NCP = n * h2 / (1 - h2). MAF enters indirectly: low MAF -> higher SE of the
        # estimated allele effect for the same h2. Account for this via the marker-
        # variance factor 2*MAF*(1-MAF) being smaller, which scales NCP as:
        #   NCP_eff = NCP * (2*MAF*(1-MAF)) / (2*0.5*0.5)
        # That is, normalised so MAF=0.5 -> full NCP.
        ncp_base = n * h2 / max(1e-9, 1 - h2)
        scale = (2 * mf * (1 - mf)) / 0.5
        ncp = ncp_base * scale
        power = 1 - stats.ncx2.cdf(thresh_chi2, df=1, nc=ncp)
        power_rows.append({"maf": mf, "h2_qtl": h2, "ncp": ncp, "power": power})
power_df = pd.DataFrame(power_rows)
power_df.to_csv(TAB / "power_curve.csv", index=False)

# minimum-detectable h2 at power 0.80 for each MAF
mdq_rows = []
for mf in maf_grid:
    h_search = np.linspace(0.01, 0.95, 500)
    powers = []
    for h2 in h_search:
        ncp = (n * h2 / max(1e-9, 1 - h2)) * (2 * mf * (1 - mf)) / 0.5
        p = 1 - stats.ncx2.cdf(thresh_chi2, df=1, nc=ncp)
        powers.append(p)
    powers = np.array(powers)
    idx = np.where(powers >= 0.80)[0]
    mdq_rows.append({
        "maf": mf,
        "h2_for_power_0.80": float(h_search[idx[0]]) if len(idx) else float("nan"),
    })
mdq = pd.DataFrame(mdq_rows)
mdq.to_csv(TAB / "min_detectable_h2.csv", index=False)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
print("[fig] marker QC panels...")
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
sns.histplot(marker_qc["maf"].dropna(), bins=40, ax=axes[0, 0], color="#3b6fb6")
axes[0, 0].axvline(MIN_MAF, ls="--", c="k")
axes[0, 0].set(title=f"Marker MAF (cutoff {MIN_MAF})", xlabel="MAF")

sns.histplot(marker_qc["call_rate"].dropna(), bins=40, ax=axes[0, 1], color="#3b6fb6")
axes[0, 1].axvline(MARK_CR, ls="--", c="k")
axes[0, 1].set(title=f"Marker call rate (cutoff {MARK_CR})", xlabel="Call rate")

sns.histplot(marker_qc["het"].dropna(), bins=40, ax=axes[1, 0], color="#3b6fb6")
axes[1, 0].set(title="Observed heterozygosity per marker", xlabel="H_obs")

sns.histplot(marker_qc["pic"].dropna(), bins=40, ax=axes[1, 1], color="#3b6fb6")
axes[1, 1].set(title="PIC per marker", xlabel="PIC")
fig.suptitle(f"AYB DArTseq marker QC — {len(marker_qc)} SNPs, {n_samp} samples", y=1.02)
fig.tight_layout()
fig.savefig(FIG / "marker_qc_panels.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "marker_qc_panels.pdf", bbox_inches="tight")
plt.close(fig)

print("[fig] sample QC panels...")
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
sns.histplot(sample_qc["call_rate"], bins=30, ax=axes[0], color="#c14a4a")
axes[0].axvline(SAMP_CR, ls="--", c="k")
axes[0].set(title=f"Sample call rate (cutoff {SAMP_CR})", xlabel="Call rate")
sns.scatterplot(data=sample_qc, x="call_rate", y="observed_het", ax=axes[1],
                s=50, color="#c14a4a", edgecolor="black", alpha=0.75)
axes[1].set(title="Sample call rate vs observed het",
            xlabel="Call rate", ylabel="Observed heterozygosity")
fig.tight_layout()
fig.savefig(FIG / "sample_qc_panels.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "sample_qc_panels.pdf", bbox_inches="tight")
plt.close(fig)

print("[fig] PCA scree...")
fig, ax = plt.subplots(figsize=(7, 5))
ax.bar(var_exp["PC"], var_exp["variance_explained"] * 100, color="#3b6fb6")
ax2 = ax.twinx()
ax2.plot(var_exp["PC"], var_exp["cumulative"] * 100, color="#c14a4a", marker="o")
ax2.set_ylabel("Cumulative % variance", color="#c14a4a")
ax.set_ylabel("% variance per PC")
ax.set_title(f"PCA scree — {Xs.shape[1]} SNPs after QC")
fig.tight_layout()
fig.savefig(FIG / "pca_scree.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "pca_scree.pdf", bbox_inches="tight")
plt.close(fig)

print("[fig] PCA scatter by cluster...")
fig, ax = plt.subplots(figsize=(8, 7))
palette = sns.color_palette("Set2", n_colors=best_k)
for k in range(best_k):
    m = pca_df["cluster"] == k
    ax.scatter(pca_df.loc[m, "PC1"], pca_df.loc[m, "PC2"],
               s=60, color=palette[k], edgecolor="black",
               label=f"Cluster {k+1} (n={m.sum()})")
ax.set_xlabel(f"PC1 ({var_exp.loc[0, 'variance_explained']*100:.1f} %)")
ax.set_ylabel(f"PC2 ({var_exp.loc[1, 'variance_explained']*100:.1f} %)")
ax.set_title(f"PCA — AYB 105 lines, K-means k={best_k}")
ax.legend(loc="best", frameon=True)
fig.tight_layout()
fig.savefig(FIG / "pca_pc1_pc2.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "pca_pc1_pc2.pdf", bbox_inches="tight")
plt.close(fig)

print("[fig] PCA colored by phenotype...")
fig, axes = plt.subplots(2, 2, figsize=(13, 11))
for ax, trait in zip(axes.ravel(), TRAITS):
    sc = ax.scatter(pca_df["PC1"], pca_df["PC2"], c=pca_df[trait],
                    s=55, cmap="viridis", edgecolor="black")
    plt.colorbar(sc, ax=ax, label=trait)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(trait)
fig.suptitle("PCA scores colored by trait value", y=1.01)
fig.tight_layout()
fig.savefig(FIG / "pca_pc1_pc2_pheno.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "pca_pc1_pc2_pheno.pdf", bbox_inches="tight")
plt.close(fig)

print("[fig] power heatmap...")
piv = power_df.pivot(index="h2_qtl", columns="maf", values="power")
fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(piv, annot=True, fmt=".2f", cmap="rocket_r", vmin=0, vmax=1, ax=ax,
            cbar_kws={"label": "Power"})
ax.invert_yaxis()
ax.set_xlabel("Minor allele frequency (MAF)")
ax.set_ylabel("QTL variance fraction (h$^2_{qtl}$)")
ax.set_title(f"GWAS power — n={n}, Bonferroni $\\alpha$ = {0.05/m_eff:.2e}")
fig.tight_layout()
fig.savefig(FIG / "power_heatmap.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "power_heatmap.pdf", bbox_inches="tight")
plt.close(fig)

# minimum-detectable h2 plot
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(mdq["maf"], mdq["h2_for_power_0.80"], marker="o", color="#3b6fb6", lw=2)
ax.set_xlabel("MAF")
ax.set_ylabel("QTL h$^2$ needed for 80% power")
ax.set_title(f"Minimum-detectable QTL effect — n={n}, m={m_eff}")
ax.grid(True, alpha=0.4)
fig.tight_layout()
fig.savefig(FIG / "min_detectable_h2.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "min_detectable_h2.pdf", bbox_inches="tight")
plt.close(fig)


# ----------------------------------------------------------------------------
# Console summary
# ----------------------------------------------------------------------------
print("\n=========== SUMMARY ===========")
print(f"Markers (input):  {n_mark}")
print(f"Samples (input):  {n_samp}")
print(f"Markers post-QC:  {len(keep_mark)}  (call >= {MARK_CR}, MAF >= {MIN_MAF})")
print(f"Samples post-QC:  {len(keep_samp)}  (call >= {SAMP_CR})")
print(f"Median MAF (pre-filter): {marker_qc['maf'].median():.3f}")
print(f"Median PIC (pre-filter): {marker_qc['pic'].median():.3f}")
print(f"Mean sample call rate:   {sample_qc['call_rate'].mean():.3f}")
print(f"PCA: PC1 + PC2 explain  {(var_exp.loc[0,'variance_explained'] + var_exp.loc[1,'variance_explained'])*100:.1f}%")
print(f"K-means best k (silhouette): {best_k}")
print("\nGWAS power highlights (n={}, m_eff={}, alpha=0.05 Bonferroni):".format(n, m_eff))
for mf in maf_grid:
    h2_need = mdq.loc[mdq['maf']==mf, 'h2_for_power_0.80'].values[0]
    print(f"  MAF={mf:.2f}: minimum h2_qtl for 80% power = {h2_need:.2f}")
print(f"\nOutputs written to: {OUT}")
