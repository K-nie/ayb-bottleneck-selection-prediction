"""
81 -- Per-chromosome triangular LD r^2 heatmaps for Ss05, Ss10, Ss04.

Direct visualisation of the LD architecture on the three chromosomes that
carry the paper-load-bearing signals:

  Ss05 -- carries the Cluster-2 attribution-driven sites from Paper A1.
  Ss10 -- carries ALMT4_2 at 15.4 Mb (Paper A2 four-method consensus locus).
  Ss04 -- carries the suggestive MATE-cluster oxalate locus and is the
          chromosome with the strongest Tajima D distribution at the
          chromosome scale (median D = 2.13, 55 % of windows D > 2).

For each chromosome, we compute pairwise r^2 between AYB-anchored markers
on that chromosome using the working dosage matrix (mean imputation) and
draw a rotated triangular heatmap with the marker position on the x axis.
A horizontal track below the heatmap shows marker density and flags the
A2 candidate-gene callout coordinates.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = PROJ / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = PROJ / "refs/ayb_genome/ayb_marker_anchoring.csv"
OUT = PROJ / "results/81_per_chr_ld_heatmap"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05

CHRS = ["Ss05", "Ss10", "Ss04"]
CALLOUTS = {
    "Ss10": [("ALMT4_2", 15394673)],
    "Ss04": [("MATE cluster", 67792286)],
    "Ss05": [],
}


def load_dosage():
    hm = pd.read_csv(HAPMAP, low_memory=False)
    META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sample_cols = [c for c in hm.columns if c not in META]
    ref_alt = hm["alleles"].str.split("/", expand=True)
    ref_alt.columns = ["ref", "alt"]
    dosage = np.full((len(hm), len(sample_cols)), np.nan)
    calls = hm[sample_cols].astype(str)
    for i in range(len(hm)):
        r_, a_ = ref_alt.iloc[i]
        arr = calls.iloc[i].values
        dosage[i, arr == r_ + r_] = 0.0
        dosage[i, (arr == r_ + a_) | (arr == a_ + r_)] = 1.0
        dosage[i, arr == a_ + a_] = 2.0
    dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)
    cr_m = dosage.notna().sum(axis=1) / dosage.shape[1]
    maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2,
                     1 - dosage.mean(axis=1, skipna=True) / 2)
    cr_s = dosage.notna().sum(axis=0) / dosage.shape[0]
    keep_m = ((cr_m >= MARK_CR) & (maf >= MIN_MAF)).values
    keep_s = (cr_s >= SAMP_CR).values
    dosage = dosage.loc[keep_m, keep_s]
    return dosage


def pairwise_r2(X: np.ndarray) -> np.ndarray:
    """X: samples x markers, NaN-imputed by column mean."""
    Xs = X - np.nanmean(X, axis=0)
    Xs = np.where(np.isnan(Xs), 0.0, Xs)
    norms = np.linalg.norm(Xs, axis=0)
    norms[norms == 0] = 1.0
    Xn = Xs / norms
    r = Xn.T @ Xn
    return r ** 2


LD_CMAP = LinearSegmentedColormap.from_list(
    # Haploview-style white-to-dark-red, but with a faint cream base so
    # r^2 values near the panel median (~0.012) still render visibly
    # instead of disappearing into the white background.
    "ld",
    ["#FFF3E0", "#FFD8A8", "#F4A261", "#E76F51", "#C0392B", "#7B241C"],
    N=256,
)


def draw_triangle(r2: np.ndarray, positions: np.ndarray, ax, vmax=0.5,
                    gamma: float = 0.40):
    """Render the upper triangle of `r2` as rotated diamonds. Uses a
    PowerNorm with gamma < 1 to amplify the low-r^2 end of the colour
    scale -- in a sparse-LD panel like AYB (median r^2 ~ 0.012), linear
    scaling collapses everything below 0.05 into invisible white, while
    gamma = 0.4 maps 0.012 to roughly 22 %% of the colour ramp."""
    import matplotlib.colors as mcolors
    n = r2.shape[0]
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    patches = []
    vals = []
    pos_mb = positions / 1e6
    for i in range(n):
        for j in range(i + 1, n):
            x = (pos_mb[i] + pos_mb[j]) * 0.5
            y = (pos_mb[j] - pos_mb[i]) * 0.5
            hw = max((pos_mb[j] - pos_mb[i]) * 0.25, 0.02)
            patches.append(Polygon([
                (x - hw, y), (x, y - hw), (x + hw, y), (x, y + hw),
            ]))
            vals.append(r2[i, j])
    norm = mcolors.PowerNorm(gamma=gamma, vmin=0.0, vmax=vmax)
    pc = PatchCollection(patches, cmap=LD_CMAP, norm=norm,
                          alpha=0.95, edgecolor="none")
    pc.set_array(np.array(vals))
    ax.add_collection(pc)
    ax.autoscale_view()
    return pc


def main() -> None:
    dosage = load_dosage()
    anc = pd.read_csv(ANCHOR)
    anc = anc[anc["rs"].isin(dosage.index)].copy()
    anc = anc[["rs", "chr_ayb", "snp_pos_ayb"]].dropna()
    print(f"AYB-anchored markers in QC-passed dosage: {len(anc)}")

    fig, axes = plt.subplots(len(CHRS), 1, figsize=(13.5, 3.8 * len(CHRS)),
                              gridspec_kw={"hspace": 0.55})
    if len(CHRS) == 1:
        axes = [axes]

    rows_summary = []
    for ax, chrom in zip(axes, CHRS):
        sub = anc[anc["chr_ayb"] == chrom].sort_values("snp_pos_ayb").copy()
        if sub.empty:
            ax.text(0.5, 0.5, f"No anchored markers on {chrom}",
                    ha="center", transform=ax.transAxes)
            continue
        markers = sub["rs"].tolist()
        positions = sub["snp_pos_ayb"].values.astype(float)
        X = dosage.loc[markers].T.values.astype(float)
        # column-mean impute
        col_mu = np.nanmean(X, axis=0)
        X = np.where(np.isnan(X), col_mu, X)
        r2 = pairwise_r2(X)
        # save tables
        triu_i, triu_j = np.triu_indices_from(r2, k=1)
        long = pd.DataFrame({
            "rs_i": np.array(markers)[triu_i],
            "rs_j": np.array(markers)[triu_j],
            "pos_i": positions[triu_i],
            "pos_j": positions[triu_j],
            "r2": r2[triu_i, triu_j],
        })
        long.to_csv(OUT / f"tables/ld_pairs_{chrom}.csv", index=False)
        med_r2 = float(np.median(long["r2"]))
        p95_r2 = float(np.percentile(long["r2"], 95))
        rows_summary.append({
            "chr": chrom, "n_markers": len(markers),
            "n_pairs": len(long), "r2_median": med_r2,
            "r2_p95": p95_r2,
            "span_mb": float((positions[-1] - positions[0]) / 1e6),
        })

        pc = draw_triangle(r2, positions, ax, vmax=0.5)
        ax.set_xlim(positions.min() / 1e6 - 0.5,
                    positions.max() / 1e6 + 0.5)
        ax.set_ylim(0, (positions.max() - positions.min()) / 2e6 * 1.05)
        ax.set_aspect("equal")
        ax.set_ylabel("Pair half-distance (Mb)")
        ax.set_xlabel(f"{chrom} position (Mb)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_title(f"{chrom}: n = {len(markers)} anchored markers\n"
                     f"median pair $r^2$ = {med_r2:.3f}; "
                     f"95th pct = {p95_r2:.3f}",
                     fontsize=10, pad=8)
        # callouts
        for nm, pos in CALLOUTS.get(chrom, []):
            ax.axvline(pos / 1e6, color="#222", lw=1.0, ls=":", alpha=0.85)
            ax.text(pos / 1e6, ax.get_ylim()[1] * 0.95, nm,
                     ha="center", fontsize=8, color="#222",
                     bbox=dict(facecolor="white", edgecolor="none",
                               alpha=0.8, pad=1.5))

    for k, ax in enumerate(axes):
        ax.text(-0.045, 1.10, chr(ord("a") + k), transform=ax.transAxes,
                fontsize=12, fontweight="bold", ha="right", va="bottom",
                zorder=1000)

    pd.DataFrame(rows_summary).to_csv(OUT / "tables/per_chr_summary.csv",
                                      index=False)
    print(pd.DataFrame(rows_summary))

    # one colour bar across the figure -- must use the same PowerNorm
    # the per-axis triangles use, so the colourbar ticks correspond to
    # the actual r^2 values rather than the gamma-transformed scale
    import matplotlib.colors as mcolors
    cb_norm = mcolors.PowerNorm(gamma=0.40, vmin=0.0, vmax=0.5)
    sm = plt.cm.ScalarMappable(cmap=LD_CMAP, norm=cb_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.022, pad=0.015,
                          ticks=[0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5])
    cbar.set_label(r"pairwise $r^2$")

    fig.suptitle("Per-chromosome triangular LD heatmaps: Ss05, Ss10, Ss04",
                 fontsize=12, y=1.01)
    fig.savefig(OUT / "figures/fig_per_chr_ld.png", dpi=300,
                bbox_inches="tight")
    fig.savefig(OUT / "figures/fig_per_chr_ld.pdf", bbox_inches="tight")
    plt.close(fig)
    print("\nWrote", OUT / "figures/fig_per_chr_ld.png")


if __name__ == "__main__":
    main()
