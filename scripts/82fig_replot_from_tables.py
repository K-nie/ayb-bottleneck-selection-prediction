"""
82fig -- house-style replot of B14, the leak-corrected top-quintile binary
classifier benchmark, from the saved tables written by
82_classifier_top_quintile.py.

The parent script refits 13 traits x 5 models under stratified 5 x 5 CV plus a
15 % holdout on every run, so it is not re-run to restyle; this reads the frozen
CSVs.

Three panels:
  A  CV ROC-AUC (mean across 25 folds), trait x model heatmap.
  B  CV top-5 breeder accuracy, trait x model heatmap.
  C  Final-holdout ROC-AUC for the per-trait CV winner (horizontal bars).

Inputs  (results/82_classifier_top_quintile/tables/):
  cv_aggregated.csv        per trait x model CV means
  cv_winner_per_trait.csv  best model per trait by CV AUC mean
  final_holdout.csv        holdout metrics per trait x model
Output: results/82_classifier_top_quintile/figures/fig_classifier_top_quintile.{png,pdf}
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, save_figure, panel_label
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAB = PROJ / "results/82_classifier_top_quintile/tables"
OUT = PROJ / "results/82_classifier_top_quintile/figures"
OUT.mkdir(parents=True, exist_ok=True)

TOP_K_BREEDER = 5


def nice(t: str) -> str:
    return t.replace("_", " ")


def main() -> None:
    agg = pd.read_csv(TAB / "cv_aggregated.csv")
    best = pd.read_csv(TAB / "cv_winner_per_trait.csv")
    hold = pd.read_csv(TAB / "final_holdout.csv")

    traits = sorted(agg["trait"].unique())
    models = sorted(agg["model"].unique())
    auc_mat = np.full((len(traits), len(models)), np.nan)
    topk_mat = np.full((len(traits), len(models)), np.nan)
    for i, t in enumerate(traits):
        for j, m in enumerate(models):
            sub = agg[(agg["trait"] == t) & (agg["model"] == m)]
            if not sub.empty:
                auc_mat[i, j] = sub["auc_mean"].iloc[0]
                topk_mat[i, j] = sub["topk_mean"].iloc[0]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, max(5.0, 0.45 * len(traits))))

    cmap_auc = LinearSegmentedColormap.from_list(
        "auc", ["#cccccc", "#dddddd", WONG["skyblue"], WONG["green"]], N=256)
    im0 = axes[0].imshow(auc_mat, cmap=cmap_auc, aspect="auto", vmin=0.3, vmax=0.85)
    axes[0].set_xticks(range(len(models)))
    axes[0].set_xticklabels(models, rotation=35, ha="right")
    axes[0].set_yticks(range(len(traits)))
    axes[0].set_yticklabels([nice(t) for t in traits], fontsize=9)
    axes[0].set_title("CV ROC-AUC (mean over 25 folds)", loc="left")
    panel_label(axes[0], "a")
    for i in range(auc_mat.shape[0]):
        for j in range(auc_mat.shape[1]):
            v = auc_mat[i, j]
            if not np.isnan(v):
                axes[0].text(j, i, f"{v:.2f}", ha="center", va="center",
                             fontsize=8, color="black" if v < 0.7 else "white")
    fig.colorbar(im0, ax=axes[0], fraction=0.04, pad=0.02, label="ROC-AUC")

    cmap_topk = LinearSegmentedColormap.from_list(
        "topk", ["#cccccc", "#dddddd", WONG["yellow"], WONG["purple"]], N=256)
    im1 = axes[1].imshow(topk_mat, cmap=cmap_topk, aspect="auto", vmin=0, vmax=1)
    axes[1].set_xticks(range(len(models)))
    axes[1].set_xticklabels(models, rotation=35, ha="right")
    axes[1].set_yticks(range(len(traits)))
    axes[1].set_yticklabels([""] * len(traits))
    axes[1].set_title(f"CV top-{TOP_K_BREEDER} breeder accuracy", loc="left")
    panel_label(axes[1], "b")
    for i in range(topk_mat.shape[0]):
        for j in range(topk_mat.shape[1]):
            v = topk_mat[i, j]
            if not np.isnan(v):
                axes[1].text(j, i, f"{v:.2f}", ha="center", va="center",
                             fontsize=8, color="black" if v < 0.6 else "white")
    fig.colorbar(im1, ax=axes[1], fraction=0.04, pad=0.02,
                 label=f"Top-{TOP_K_BREEDER} accuracy")

    ax = axes[2]
    if not hold.empty:
        hold_win = hold.merge(best[["trait", "model"]], on=["trait", "model"])
        hold_win = hold_win.sort_values("auc", ascending=True)
        y = np.arange(len(hold_win))
        colors = [WONG["green"] if v >= 0.7 else WONG["skyblue"] if v >= 0.5
                  else WONG["purple"] for v in hold_win["auc"]]
        ax.barh(y, hold_win["auc"], color=colors, edgecolor="white")
        ax.axvline(0.5, color="black", ls="--", lw=1, label="AUC = 0.5 (chance)")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{nice(t)} ({m})" for t, m
                            in zip(hold_win["trait"], hold_win["model"])],
                           fontsize=8)
        ax.set_xlabel("Holdout ROC-AUC")
        ax.set_title("Final-holdout AUC, per-trait CV winner", loc="left")
        panel_label(ax, "c")
        ax.set_xlim(0, 1)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        ax.grid(True, alpha=0.25, axis="x")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    fig.suptitle("Top-quintile binary classifier benchmark "
                 "(13 traits x 5 models; stratified 5 x 5 CV + 15 % holdout)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    written = save_figure(fig, OUT / "fig_classifier_top_quintile")
    plt.close(fig)
    print("Wrote", *written)


if __name__ == "__main__":
    main()
