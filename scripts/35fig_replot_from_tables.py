"""
35fig -- house-style replot of B17, the windowed Tajima's D Manhattan across
Ss01-Ss11, from results/35_tajima_d/tables/tajima_d_windowed.csv.

The parent script (35_tajima_d.py) recomputes the windowed SFS scan on every
run, so it is not re-run to restyle; this reads the frozen table. Window and step
sizes are read back from the table rather than hard-coded.

Manhattan convention matched to the other pop-gen scans (B6/B18): alternating
Wong blue / sky-blue per chromosome, dashed sweep (D = -2) and balancing/
bottleneck (D = +2) threshold lines, zero line, house axes.

Input : results/35_tajima_d/tables/tajima_d_windowed.csv
Output: results/35_tajima_d/figures/fig77_tajima_d_manhattan.{png,pdf}
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, publishable_axes, save_figure
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAB = PROJ / "results/35_tajima_d/tables/tajima_d_windowed.csv"
OUT = PROJ / "results/35_tajima_d/figures"
OUT.mkdir(parents=True, exist_ok=True)

GAP = 5_000_000


def main() -> None:
    df = pd.read_csv(TAB)
    df = df.dropna(subset=["tajima_D"]).copy()
    win_bp = int((df["win_end"] - df["win_start"]).median())
    step_bp = int(df.sort_values(["chr", "mid_bp"])
                    .groupby("chr")["mid_bp"].diff().median())
    print(f"{len(df)} windows; window {win_bp/1000:.0f} kb / step {step_bp/1000:.0f} kb")
    print(f"median D = {df['tajima_D'].median():.3f}; "
          f"D<-2: {(df['tajima_D'] < -2).sum()}; D>+2: {(df['tajima_D'] > 2).sum()}")

    chrs = sorted(df["chr"].unique(),
                  key=lambda s: int(re.search(r"\d+", s).group()))
    df = df.sort_values(["chr", "mid_bp"]).reset_index(drop=True)
    mids, xs, x_cursor = [], [], 0
    for c in chrs:
        sub = df[df["chr"] == c]
        xs.extend((sub["mid_bp"].values + x_cursor).tolist())
        mids.append(x_cursor + (sub["mid_bp"].max() - sub["mid_bp"].min()) / 2)
        x_cursor += sub["mid_bp"].max() + GAP
    df["x"] = xs

    fig, ax = plt.subplots(figsize=(11, 4.2))
    for i, c in enumerate(chrs):
        sub = df[df["chr"] == c]
        col = WONG["blue"] if (i % 2 == 0) else WONG["skyblue"]
        ax.scatter(sub["x"], sub["tajima_D"], s=14, color=col, alpha=0.85,
                   edgecolors="none", rasterized=True)
    ax.axhline(0, color="black", lw=0.6)
    ax.axhline(-2, color=WONG["vermillion"], lw=0.7, ls="--", alpha=0.8,
               label="$D = -2$ (recent sweep)")
    ax.axhline(+2, color=WONG["green"], lw=0.7, ls="--", alpha=0.8,
               label="$D = +2$ (bottleneck / balancing)")
    ax.set_xticks(mids)
    ax.set_xticklabels(chrs, fontsize=8)
    ax.set_xlabel(r"$\it{S.\ stenocarpa}$ pseudo-chromosome (AYB-anchored markers)")
    ax.set_ylabel("Tajima's $D$")
    ax.set_title(f"Windowed Tajima's $D$ across the genome "
                 f"({len(df)} windows, {win_bp/1000:.0f} kb / {step_bp/1000:.0f} kb)",
                 loc="left")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9, edgecolor="none")
    publishable_axes(ax, grid="y")

    fig.tight_layout()
    written = save_figure(fig, OUT / "fig77_tajima_d_manhattan")
    plt.close(fig)
    print("Wrote", *written)


if __name__ == "__main__":
    main()
