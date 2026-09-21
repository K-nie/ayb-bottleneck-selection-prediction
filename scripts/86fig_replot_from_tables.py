"""
86fig -- house-style replot of B30, the per-pathway selection-statistic
enrichment forest, from the saved table written by
86_pathway_selection_enrichment.py.

The parent script parses the ~30 s Funannotate GFF and runs a permutation test
(n_perm = 2999) per pathway on every invocation, so it is not re-run to restyle;
this reads the frozen enrichment table.

Two panels of horizontal deltas (pathway mean minus genome-wide background):
  A  panel ROH fraction delta
  B  number of iHS outliers within +/- 100 kb, delta
Bars are coloured by permutation significance and direction (green = positive
enriched at p_perm < 0.05, vermillion = negative, grey = ns).

Input : results/86_pathway_selection_enrichment/tables/pathway_selection_enrichment.csv
Output: results/86_pathway_selection_enrichment/figures/fig_pathway_selection_enrichment.{png,pdf}
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, save_figure
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAB = PROJ / "results/86_pathway_selection_enrichment/tables/pathway_selection_enrichment.csv"
OUT = PROJ / "results/86_pathway_selection_enrichment/figures"
OUT.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05


def short(name: str) -> str:
    # keep the plain pathway label, drop the parenthetical gene-list gloss
    return name.split(" (")[0].strip()


def bar_colors(deltas, pvals):
    out = []
    for d, p in zip(deltas, pvals):
        if p < ALPHA and d > 0:
            out.append(WONG["green"])
        elif p < ALPHA and d < 0:
            out.append(WONG["vermillion"])
        else:
            out.append(WONG["grey"])
    return out


def panel(ax, pdf, delta_col, p_col, xlabel, title):
    pl = pdf.sort_values(delta_col).reset_index(drop=True)
    y = np.arange(len(pl))
    ax.barh(y, pl[delta_col], color=bar_colors(pl[delta_col], pl[p_col]),
            edgecolor="white")
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{short(n)}  (n={ng})"
                        for n, ng in zip(pl["pathway"], pl["n_genes"])],
                       fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left")
    ax.grid(True, alpha=0.25, axis="x")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main() -> None:
    pdf = pd.read_csv(TAB)
    print(f"{len(pdf)} pathways; ROH p<{ALPHA}: {(pdf['p_perm_roh'] < ALPHA).sum()}; "
          f"iHS p<{ALPHA}: {(pdf['p_perm_ihs'] < ALPHA).sum()}")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6))
    panel(axes[0], pdf, "roh_delta", "p_perm_roh",
          "Delta panel ROH fraction (pathway - background)",
          "A  Per-pathway ROH-density enrichment")
    panel(axes[1], pdf, "ihs_delta", "p_perm_ihs",
          "Delta iHS outliers within +/- 100 kb (pathway - background)",
          "B  Per-pathway iHS-proximity enrichment")

    handles = [
        plt.Line2D([], [], color=WONG["green"], lw=6,
                   label=f"enriched (p_perm < {ALPHA})"),
        plt.Line2D([], [], color=WONG["vermillion"], lw=6,
                   label=f"depleted (p_perm < {ALPHA})"),
        plt.Line2D([], [], color=WONG["grey"], lw=6, label="not significant"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Per-pathway selection-statistic enrichment "
                 "(curated gene sets vs genome-wide Funannotate background)",
                 fontsize=12, y=1.01)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    written = save_figure(fig, OUT / "fig_pathway_selection_enrichment")
    plt.close(fig)
    print("Wrote", *written)


if __name__ == "__main__":
    main()
