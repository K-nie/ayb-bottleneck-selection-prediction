"""
79 -- fineSTRUCTURE coancestry heatmap with phenotype side-annotation.

Cross-paper bridge B <-> A2: do the K = 24 haplotype-resolution clusters
recovered by ChromoPainter + fineSTRUCTURE track any of the 13 phenotype
BLUPs?

  - Heatmap: 95 x 95 coancestry chunk-count matrix, samples reordered to put
    cluster members adjacent.
  - Right-hand annotations: 13 columns, one per trait, coloured by per-sample
    BLUP z-score.
  - Top bar: cluster identity (K01 .. K24).
  - Statistics table: per-trait Kruskal-Wallis H of BLUP ~ K = 24 cluster;
    BH-FDR across the 13 traits.

Honest framing: significant trait ~ cluster associations at K = 24 do not
imply causal genetic structure. They imply that the breeding pedigrees that
form the panel's haplotype subdivisions also stratify trait variance. We
report cluster-level group means and explicit "n_cluster < 3 dropped" notes
so the test is not over-stated.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import kruskal

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, save_figure, panel_label
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
COANC = PROJ / "results/37_finestructure/tables/coancestry.csv"
CLUSTS = PROJ / "results/37_finestructure/tables/cluster_assignments.csv"
BLUP = PROJ / "results/00_blup_pipeline/tables/phenotype_blups.csv"
OUT = PROJ / "results/79_finestructure_x_phenotype"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

TRAITS = ["Tannin", "Phenol", "Flavonoid", "Antioxidant", "Seed_Length",
          "Seed_Width", "Seed_Thickness", "Mass_of_Seeds", "Seed_Coat_Tannin",
          "Crude_Protein", "Total_Oxalate", "Soluble_Oxalate",
          "Insoluble_Oxalate"]


def main() -> None:
    coanc_df = pd.read_csv(COANC, index_col=0)
    clusts = pd.read_csv(CLUSTS).set_index("sample")
    blup = pd.read_csv(BLUP).set_index("sample")

    # align samples
    common = sorted(set(coanc_df.index) & set(clusts.index))
    print(f"Common samples (coancestry x cluster): {len(common)}")
    coanc_df = coanc_df.loc[common, common]
    clusts = clusts.loc[common]

    # reorder samples by cluster id then by sample name within cluster
    ordered = clusts.sort_values(["cluster"]).index.tolist()
    M = coanc_df.loc[ordered, ordered].values
    clusters_ordered = clusts.loc[ordered, "cluster"].values

    # Kruskal-Wallis per trait against K = 24 cluster
    rows = []
    cluster_means = []
    for t in TRAITS:
        if t not in blup.columns:
            continue
        s = blup[t].reindex(ordered)
        df = pd.DataFrame({"cluster": clusters_ordered, "y": s.values}).dropna()
        groups = [
            g["y"].values for _, g in df.groupby("cluster")
            if len(g) >= 2
        ]
        if len(groups) < 2:
            continue
        h, p = kruskal(*groups)
        n_used = sum(len(g) for g in groups)
        rows.append({
            "trait": t, "n_used": int(n_used),
            "n_groups": int(len(groups)),
            "kruskal_H": float(h), "p_raw": float(p),
        })
        # per-cluster mean BLUP for this trait
        means = df.groupby("cluster")["y"].agg(["mean", "size"]).reset_index()
        means["trait"] = t
        cluster_means.append(means)
    stat = pd.DataFrame(rows)
    # BH-FDR
    from statsmodels.stats.multitest import multipletests
    if len(stat):
        ok = ~stat["p_raw"].isna()
        if ok.any():
            stat.loc[ok, "p_bh"] = multipletests(
                stat.loc[ok, "p_raw"], method="fdr_bh"
            )[1]
    stat = stat.sort_values("p_raw")
    stat.to_csv(OUT / "tables/kruskal_trait_vs_K24.csv", index=False)
    print("\nKruskal-Wallis trait ~ K=24 cluster:")
    print(stat)
    if cluster_means:
        pd.concat(cluster_means, ignore_index=True).to_csv(
            OUT / "tables/per_cluster_trait_means.csv", index=False
        )

    # ----- Figure -----
    n = len(ordered)
    fig = plt.figure(figsize=(14, 8.6))
    gs = fig.add_gridspec(
        2, 3, height_ratios=[0.45, 14], width_ratios=[14, 0.45, 7.5],
        hspace=0.04, wspace=0.06,
    )

    ax_top = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[1, 0], sharex=ax_top)
    ax_left_pad = fig.add_subplot(gs[1, 1])  # spacing
    ax_left_pad.axis("off")
    ax_traits = fig.add_subplot(gs[1, 2], sharey=ax_main)

    # main heatmap (log scale, off-diagonal)
    M_plot = M.copy()
    np.fill_diagonal(M_plot, np.nan)
    vmax = float(np.nanpercentile(M_plot, 99))
    vmin = float(np.nanpercentile(M_plot, 5))
    im = ax_main.imshow(
        M_plot, cmap="viridis", aspect="auto", vmin=vmin, vmax=vmax,
        interpolation="nearest",
    )
    ax_main.set_xlim(-0.5, n - 0.5); ax_main.set_ylim(n - 0.5, -0.5)
    ax_main.set_xticks([]); ax_main.set_yticks([])
    ax_main.set_xlabel("Clusters")
    ax_main.set_ylabel("Samples")
    cbar = fig.colorbar(im, ax=ax_main, fraction=0.04, pad=0.01,
                        orientation="vertical")
    cbar.set_label("Chunk-count (off-diagonal)")

    # top: cluster label bar
    cluster_levels = sorted(set(clusters_ordered))
    cmap_clust = plt.get_cmap("tab20b", len(cluster_levels))
    clust_to_idx = {c: i for i, c in enumerate(cluster_levels)}
    clust_codes = np.array([clust_to_idx[c] for c in clusters_ordered])[
        np.newaxis, :
    ]
    ax_top.imshow(clust_codes, cmap=cmap_clust, aspect="auto",
                  interpolation="nearest")
    ax_top.set_yticks([])
    ax_top.set_xticks([])
    ax_top.set_title(f"K = 24 fineSTRUCTURE clusters along the sample axis "
                     f"(n = {n} samples)", fontsize=10)
    # mark cluster boundaries
    prev = clusters_ordered[0]; cluster_starts = [0]; cluster_ids = [prev]
    for i, c in enumerate(clusters_ordered[1:], start=1):
        if c != prev:
            ax_main.axvline(i - 0.5, color="white", lw=0.45, alpha=0.7)
            ax_main.axhline(i - 0.5, color="white", lw=0.45, alpha=0.7)
            cluster_starts.append(i); cluster_ids.append(c)
            prev = c
    # label cluster ids on the top axis. White glyphs with a dark stroke halo
    # stay legible over any tab20b block colour (some blocks are pale, where
    # plain white text vanished before).
    for i, (s_, c_) in enumerate(zip(cluster_starts, cluster_ids)):
        nxt = cluster_starts[i + 1] if i + 1 < len(cluster_starts) else n
        mid = (s_ + nxt - 1) / 2.0
        if (nxt - s_) >= 1:
            ax_top.text(mid, 0.5, str(c_), ha="center", va="center",
                         fontsize=5.0, fontweight="bold", color="white",
                         rotation=90, zorder=5,
                         path_effects=[pe.withStroke(linewidth=1.0,
                                                     foreground="black")])

    # right: trait BLUP z-score matrix (samples x traits)
    blup_mat = np.full((n, len(TRAITS)), np.nan)
    for j, t in enumerate(TRAITS):
        if t not in blup.columns:
            continue
        s = blup[t].reindex(ordered).values.astype(float)
        sigma = np.nanstd(s)
        if sigma > 0:
            blup_mat[:, j] = (s - np.nanmean(s)) / sigma
    cmap_z = LinearSegmentedColormap.from_list(
        "z", ["#3050a0", "#dddddd", "#a04040"], N=256
    )
    im2 = ax_traits.imshow(
        blup_mat, aspect="auto", cmap=cmap_z, vmin=-3, vmax=3,
        interpolation="nearest",
    )
    ax_traits.set_xticks(np.arange(len(TRAITS)))
    ax_traits.set_xticklabels([t.replace("_", " ") for t in TRAITS],
                              rotation=90, ha="center", va="top", fontsize=8)
    ax_traits.tick_params(axis="x", length=2, pad=2)
    ax_traits.set_yticks([])
    ax_traits.set_title("Trait BLUP z-score", fontsize=9)
    fig.colorbar(im2, ax=ax_traits, fraction=0.04, pad=0.02, label="z-score")

    panel_label(ax_top, "a", x=-0.035, y=1.15)
    panel_label(ax_traits, "b", x=-0.14, y=1.02)

    fig.suptitle("fineSTRUCTURE coancestry (K = 24) with per-trait BLUP "
                 "z-scores across the 95-line panel", y=0.995)
    written = save_figure(fig, OUT / "figures/fig_coancestry_x_phenotype")
    plt.close(fig)
    print("\nWrote", *written)


if __name__ == "__main__":
    main()
