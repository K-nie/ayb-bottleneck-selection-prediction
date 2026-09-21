"""
36fig -- house-style replot of B18, the standalone selscan iHS Manhattan, from
results/36_ihs/tables/ihs_combined.csv.

B18 is the companion to B6 (which carries the Paper A2 candidate-gene callouts);
this standalone view plots |standardised iHS| across Ss01-Ss11 with the single
|iHS| = 2 outlier threshold and no callouts. The parent script (36_ihs.py) runs
Beagle phasing + per-chromosome selscan on every invocation, so it is not re-run
to restyle; this reads the frozen combined table.

Manhattan convention matched to B6/B17: alternating Wong blue / sky-blue per
chromosome, dashed threshold line, house axes.

Input : results/36_ihs/tables/ihs_combined.csv
Output: results/36_ihs/figures/fig80_ihs_manhattan.{png,pdf}
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, publishable_axes, save_figure
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAB = PROJ / "results/36_ihs/tables/ihs_combined.csv"
OUT = PROJ / "results/36_ihs/figures"
OUT.mkdir(parents=True, exist_ok=True)

CM_PER_MB = 1.06   # Lonardi 2019 cowpea proxy, uniform map
GAP = 5_000_000


def main() -> None:
    ihs = pd.read_csv(TAB).dropna(subset=["norm_iHS"]).copy()
    ihs["absiHS"] = ihs["norm_iHS"].abs()
    n_sig = int((ihs["absiHS"] > 2).sum())
    print(f"{len(ihs)} markers; |iHS|>2: {n_sig} ({100 * n_sig / len(ihs):.1f} %)")

    chrs = sorted(ihs["chr"].unique(),
                  key=lambda s: int(re.search(r"\d+", s).group()))
    ihs = ihs.sort_values(["chr", "pos"]).reset_index(drop=True)
    mids, xs, x_cursor = [], [], 0
    for c in chrs:
        sub = ihs[ihs["chr"] == c]
        xs.extend((sub["pos"].values + x_cursor).tolist())
        mids.append(x_cursor + (sub["pos"].max() - sub["pos"].min()) / 2)
        x_cursor += sub["pos"].max() + GAP
    ihs["x"] = xs

    fig, ax = plt.subplots(figsize=(11, 3.9))
    for i, c in enumerate(chrs):
        sub = ihs[ihs["chr"] == c]
        col = WONG["blue"] if (i % 2 == 0) else WONG["skyblue"]
        ax.scatter(sub["x"], sub["absiHS"], s=11, color=col, alpha=0.75,
                   edgecolors="none", rasterized=True)
    ax.axhline(2, color=WONG["vermillion"], lw=0.7, ls="--", alpha=0.8,
               label="|iHS| = 2 (outlier threshold)")
    ax.set_xticks(mids)
    ax.set_xticklabels(chrs, fontsize=8)
    ax.set_xlabel(r"$\it{S.\ stenocarpa}$ pseudo-chromosome (AYB-anchored markers)")
    ax.set_ylabel("|standardised iHS|")
    ax.set_title(f"selscan iHS scan ({len(ihs):,} markers, "
                 f"{n_sig} with |iHS| > 2; uniform {CM_PER_MB} cM/Mb map)",
                 loc="left")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9, edgecolor="none")
    publishable_axes(ax, grid="y")

    fig.tight_layout()
    written = save_figure(fig, OUT / "fig80_ihs_manhattan")
    plt.close(fig)
    print("Wrote", *written)


if __name__ == "__main__":
    main()
