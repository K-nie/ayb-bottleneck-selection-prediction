"""
37fig -- house-style replot of the two fineSTRUCTURE supp figures from saved
artifacts:

  B19 (fig81_coancestry_heatmap): ChromoPainter chunk-count coancestry matrix,
      log1p-compressed, re-ordered by the fineSTRUCTURE MAP tree so cluster
      blocks fall on the diagonal.
  B13 (fig82_finestructure_tree): the MAP coancestry tree (K = 24), tips
      coloured by their MAP cluster assignment.

The parent script (37_finestructure.py) runs the full ChromoPainter +
fineSTRUCTURE MCMC pipeline on every invocation, so it is not re-run to restyle.
This reads the frozen coancestry matrix, the saved Newick, and the parsed cluster
assignments.

Inputs  (results/37_finestructure/):
  tables/coancestry.csv            95 x 95 chunk-count matrix (raw order)
  tables/cluster_assignments.csv   sample -> MAP cluster (K01..K24)
  data/ayb_fs_linked.newick        MAP coancestry tree

Outputs (results/37_finestructure/figures/):
  fig81_coancestry_heatmap.{png,pdf}
  fig82_finestructure_tree.{png,pdf}
"""
from __future__ import annotations
import io
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, save_figure
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
BASE = PROJ / "results/37_finestructure"
CC = BASE / "tables/coancestry.csv"
CLUST = BASE / "tables/cluster_assignments.csv"
NEWICK = BASE / "data/ayb_fs_linked.newick"
OUT = BASE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

K_MAP = 24


def tree_tip_order(tr) -> list[str]:
    return [t.name for t in tr.get_terminals()]


def cluster_palette(n: int) -> list:
    # 24 reasonably distinct colours from two qualitative maps.
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    return [base[i % len(base)] for i in range(n)]


def main() -> None:
    cc = pd.read_csv(CC, index_col=0)
    clusters = pd.read_csv(CLUST)
    samp2clust = dict(zip(clusters["sample"], clusters["cluster"]))
    tr = Phylo.read(io.StringIO(NEWICK.read_text().strip()), "newick")

    order = [s for s in tree_tip_order(tr) if s in cc.index and s in cc.columns]
    if len(order) < cc.shape[0] // 2:
        order = list(cc.index)
    cc_ord = cc.loc[order, order]
    print(f"coancestry {cc.shape}; ordered {cc_ord.shape} by tree tips")

    # cluster -> colour, in the order clusters first appear along the tree
    clist = sorted(clusters["cluster"].unique())
    pal = cluster_palette(len(clist))
    clust2col = dict(zip(clist, pal))

    # ---------- B19: coancestry heatmap ----------
    fig, ax = plt.subplots(figsize=(8.6, 7.6))
    M = cc_ord.values.astype(float)
    im = ax.imshow(np.log1p(M), cmap="magma", aspect="equal")
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("log(1 + chunk count)")
    # thin cluster colour strip along the top + left edges so blocks are legible
    n = len(order)
    strip = np.array([clust2col[samp2clust.get(s, clist[0])] for s in order])
    ax_top = ax.inset_axes([0, 1.005, 1, 0.018], transform=ax.transAxes)
    ax_left = ax.inset_axes([-0.028, 0, 0.018, 1], transform=ax.transAxes)
    ax_top.imshow(strip[np.newaxis, :, :], aspect="auto")
    ax_left.imshow(strip[:, np.newaxis, :], aspect="auto")
    for a in (ax_top, ax_left):
        a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values():
            s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(f"Donor (n = {n}, fineSTRUCTURE-tree order)")
    ax.set_ylabel(f"Recipient (n = {n}, fineSTRUCTURE-tree order)")
    ax.set_title(f"ChromoPainter coancestry, MAP K = {K_MAP} "
                 "(African yam bean 95-line panel)", loc="left", pad=14)
    fig.tight_layout()
    written = save_figure(fig, OUT / "fig81_coancestry_heatmap")
    plt.close(fig)
    print("Wrote", *written)

    # ---------- B13: fineSTRUCTURE tree ----------
    fig, ax = plt.subplots(figsize=(8, 14))
    Phylo.draw(tr, axes=ax, do_show=False,
               branch_labels=lambda c: "",
               label_func=lambda c: c.name if c.is_terminal() else "")
    # colour each tip label by its MAP cluster
    for txt in ax.texts:
        name = txt.get_text().strip()
        if name in samp2clust:
            txt.set_color(clust2col[samp2clust[name]])
            txt.set_fontsize(6.5)
    ax.set_title(f"fineSTRUCTURE coancestry tree (MAP K = {K_MAP})", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel("")
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    written = save_figure(fig, OUT / "fig82_finestructure_tree")
    plt.close(fig)
    print("Wrote", *written)


if __name__ == "__main__":
    main()
