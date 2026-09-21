"""
88 -- Orientation figure (paper B, Figure B1): SNP-PCA + ADMIXTURE K=3.

Panel-level population structure of the 95-line S. stenocarpa DArTseq panel,
built as the reader's orientation figure. Two panels:

  A) SNP-PCA scatter (PC1 vs PC2) of the 95 accessions, each point coloured by
     its majority ADMIXTURE K=3 ancestry component. This ties the K=2 PCA split
     the panel reads on PC1 to the K=3 ADMIXTURE partition in a single view.
  B) ADMIXTURE K=3 stacked-ancestry bar plot, one bar per accession, samples
     ordered by majority component then by dominant ancestry fraction.

K=3 (not K=4) is used because the cross-validation error is minimised at K=3
(results/09_admixture/tables/admixture_cv_error.csv: CV = 0.5143 at K=3 vs
0.52122 at K=4) and the manuscript reports the panel structure as K=3 ADMIXTURE.
There are no per-accession geographic coordinates for this genebank panel, so a
sample map is not drawn; the PCA panel is the spatial orientation instead.

Inputs:
  results/01_qc_pca_power/tables/pca_coords.csv     (95 x PC1..PC10 + sample)
  results/01_qc_pca_power/tables/pca_variance.csv   (per-PC variance explained)
  results/09_admixture/tables/admixture_Q_K3.csv    (95 x Q1..Q3, sample index)

Outputs:
  results/88_orientation_pca_admixture/figures/fig_orientation_pca_admixture.{png,pdf}
  results/88_orientation_pca_admixture/tables/sample_ancestry_assignment.csv
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, CLUSTER_PAL, publishable_axes, save_figure, panel_label
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PCA = PROJ / "results/01_qc_pca_power/tables/pca_coords.csv"
PVE = PROJ / "results/01_qc_pca_power/tables/pca_variance.csv"
QK3 = PROJ / "results/09_admixture/tables/admixture_Q_K3.csv"
OUT = PROJ / "results/88_orientation_pca_admixture"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

K = 3
QCOLS = [f"Q{i}" for i in range(1, K + 1)]
# distinct hues for the K=3 ancestry components
ANC_COL = [CLUSTER_PAL[0], CLUSTER_PAL[1], CLUSTER_PAL[2]]


def main() -> None:
    pca = pd.read_csv(PCA)
    pve = pd.read_csv(PVE).set_index("PC")["variance_explained"]
    q = pd.read_csv(QK3, index_col=0)
    print(f"Loaded PCA ({len(pca)} samples), Q-K{K} ({len(q)} samples).")

    # align on sample ID
    q = q.rename_axis("sample").reset_index()
    df = pca.merge(q[["sample"] + QCOLS], on="sample", how="inner")
    if len(df) != len(pca):
        raise SystemExit(f"sample mismatch: {len(df)} merged vs {len(pca)} PCA rows")

    # majority ancestry component per sample
    df["anc"] = df[QCOLS].values.argmax(axis=1)
    df["anc_frac"] = df[QCOLS].values.max(axis=1)

    df[["sample", "PC1", "PC2"] + QCOLS + ["anc", "anc_frac"]].to_csv(
        OUT / "tables/sample_ancestry_assignment.csv", index=False)
    counts = df["anc"].value_counts().sort_index()
    print("Majority-ancestry group sizes:",
          {f"Q{i+1}": int(counts.get(i, 0)) for i in range(K)})

    pc1_pve = 100 * float(pve.get("PC1", np.nan))
    pc2_pve = 100 * float(pve.get("PC2", np.nan))

    fig = plt.figure(figsize=(12.5, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.22)

    # ---- Panel A: PCA scatter coloured by majority ancestry ----
    axA = fig.add_subplot(gs[0, 0])
    for i in range(K):
        sub = df[df["anc"] == i]
        axA.scatter(sub["PC1"], sub["PC2"], s=46, color=ANC_COL[i],
                    alpha=0.85, edgecolor="white", linewidth=0.6,
                    label=f"ADMIXTURE cluster {i+1}", zorder=3)
    axA.axhline(0, color=WONG["grey"], linewidth=0.6, alpha=0.5, zorder=1)
    axA.axvline(0, color=WONG["grey"], linewidth=0.6, alpha=0.5, zorder=1)
    axA.set_xlabel(f"PC1 ({pc1_pve:.1f}% variance)")
    axA.set_ylabel(f"PC2 ({pc2_pve:.1f}% variance)")
    axA.set_title("SNP-PCA, coloured by ADMIXTURE K = 3", loc="left")
    panel_label(axA, "a")
    axA.legend(loc="best", fontsize=8.5, frameon=False)
    publishable_axes(axA, grid="both", grid_alpha=0.18)

    # ---- Panel B: ADMIXTURE K=3 stacked bar ----
    axB = fig.add_subplot(gs[0, 1])
    order = df.sort_values(["anc", "anc_frac"],
                           ascending=[True, False]).reset_index(drop=True)
    x = np.arange(len(order))
    bottom = np.zeros(len(order))
    for i in range(K):
        vals = order[f"Q{i+1}"].values.astype(float)
        axB.bar(x, vals, bottom=bottom, width=1.0, color=ANC_COL[i],
                edgecolor="none", label=f"cluster {i+1}")
        bottom += vals
    # thin separators between majority-ancestry blocks
    block_edges = np.where(np.diff(order["anc"].values) != 0)[0]
    for e in block_edges:
        axB.axvline(e + 0.5, color="white", linewidth=1.1, zorder=5)
    axB.set_xlim(-0.5, len(order) - 0.5)
    axB.set_ylim(0, 1)
    axB.set_xticks([])
    axB.set_ylabel("Ancestry proportion")
    axB.set_xlabel(f"Accession (n = {len(order)}, ordered by majority cluster)")
    axB.set_title("ADMIXTURE K = 3 ancestry (CV-optimal)", loc="left")
    panel_label(axB, "b")
    for side in ("top", "right"):
        axB.spines[side].set_visible(False)

    fig.suptitle(
        "Panel-level population structure of the 95-line "
        "S. stenocarpa DArTseq panel", y=1.02, fontsize=12.5)

    written = save_figure(fig, OUT / "figures/fig_orientation_pca_admixture")
    plt.close(fig)
    print("Wrote:", *[str(p) for p in written], sep="\n  ")


if __name__ == "__main__":
    main()
