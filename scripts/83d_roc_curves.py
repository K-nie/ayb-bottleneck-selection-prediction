"""
83d -- ROC curves per (trait x algorithm) for the Shiu Lab top-quintile
classifier benchmark.

Each Shiu Lab run stored per-sample mean prediction scores across 20
balanced reps in `runs/<trait>_<algo>/<trait>_<algo>_scores.txt`. This
script pools those mean scores into a single ROC curve per (trait, algo),
then renders a multi-panel grid (one panel per trait) with all four algos
overlaid plus the diagonal chance line and per-algo AUC printed in the
legend.

Output:
  fig_shiu_classifier_roc_curves.png/.pdf  -- 3-wide grid, 10 trait panels
                                                (Insoluble/Total Oxalate dropped
                                                -- no balanced runs).

Replaces the earlier fig_shiu_classifier_auc heatmap that lived in the
manuscript as a less-readable presentation of the same data.

Author: Benjamin Narh-Madey
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, publishable_axes, save_figure
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
RUNS = PROJ / "results/83_shiu_lab_classifier/runs"
OUT = PROJ / "results/83_shiu_lab_classifier/figures"
OUT.mkdir(parents=True, exist_ok=True)

ALGO_COLOUR = {
    "LogReg": WONG["blue"],
    "RF":     WONG["vermillion"],
    "GB":     WONG["green"],
    "SVMrbf": WONG["orange"],
}
ALGO_ORDER = ["LogReg", "RF", "GB", "SVMrbf"]

# Insoluble_Oxalate and Total_Oxalate are dropped: neither trait produced any
# balanced Shiu-Lab run with both classes present, so their ROC panels are
# empty. Soluble_Oxalate is retained -- it has full runs for all four algos.
TRAIT_ORDER = [
    "Crude_Protein", "Soluble_Oxalate",
    "Tannin", "Phenol", "Flavonoid", "Antioxidant",
    "Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds",
]


def load_scores(trait: str, algo: str) -> pd.DataFrame | None:
    """Load the per-sample mean score file. Returns DataFrame with columns
    ['sample', 'class', 'mean_score'] or None if missing.

    The Shiu Lab `_scores.txt` files have a duplicate `ID` column in the
    header, which shifts the header-to-data mapping by one for every
    subsequent column. We therefore index by POSITION, not by header name:
      col 0 -> sample ID (TSs...)
      col 1 -> true class label (0 / 1)
      col 2 -> mean prediction score across 20 balanced reps.
    """
    p = RUNS / f"{trait}_{algo}" / f"{trait}_{algo}_scores.txt"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t", dtype=str)
    if df.shape[1] < 3:
        return None
    out = pd.DataFrame({
        "sample": df.iloc[:, 0].astype(str),
        "class": pd.to_numeric(df.iloc[:, 1], errors="coerce"),
        "mean_score": pd.to_numeric(df.iloc[:, 2], errors="coerce"),
    }).dropna()
    out["class"] = out["class"].astype(int)
    return out


def main() -> None:
    n_traits = len(TRAIT_ORDER)
    n_cols = 3
    n_rows = int(np.ceil(n_traits / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(11.5, 3.2 * n_rows),
                              constrained_layout=True)
    axes_flat = axes.flatten()

    summary_rows = []
    for ax_idx, trait in enumerate(TRAIT_ORDER):
        ax = axes_flat[ax_idx]
        row, col = divmod(ax_idx, n_cols)
        # Chance diagonal: drawn once per panel, kept OUT of the legend so
        # each legend carries only the four model curves.
        ax.plot([0, 1], [0, 1], color=WONG["grey"], linewidth=0.8,
                  linestyle="--", zorder=1, label="_nolegend_")
        any_drawn = False
        for algo in ALGO_ORDER:
            scores = load_scores(trait, algo)
            if scores is None or scores["class"].nunique() < 2:
                continue
            try:
                auc = roc_auc_score(scores["class"], scores["mean_score"])
                fpr, tpr, _ = roc_curve(scores["class"], scores["mean_score"])
            except Exception:
                continue
            ax.plot(fpr, tpr,
                    color=ALGO_COLOUR[algo], linewidth=1.6, alpha=0.9,
                    label=f"{algo}  {auc:.3f}", zorder=3)
            summary_rows.append({"trait": trait, "algo": algo, "auc": auc,
                                   "n": int(len(scores))})
            any_drawn = True
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_yticks([0, 0.5, 1.0])
        # Outer-edge axis labels only, so inner panels stay uncluttered.
        # A panel is "bottom" if no drawn panel sits directly below it -- this
        # keeps x-axis labels on the last populated panel of every column even
        # when the final row is ragged (10 traits in a 3-wide grid).
        is_bottom = (ax_idx + n_cols) >= n_traits
        if col == 0:
            ax.set_ylabel("True positive rate", fontsize=9)
        else:
            ax.tick_params(labelleft=False)
        if is_bottom:
            ax.set_xlabel("False positive rate", fontsize=9)
        else:
            ax.tick_params(labelbottom=False)
        ax.set_title(trait.replace("_", " "), fontsize=10.5,
                      fontweight="bold")
        if any_drawn:
            leg = ax.legend(loc="lower right", fontsize=7.5, title="AUC",
                            title_fontsize=7.5, framealpha=0.92,
                            edgecolor="none", handlelength=1.4,
                            labelspacing=0.25, borderpad=0.3)
            leg._legend_box.align = "left"
        else:
            ax.text(0.5, 0.5, "no balanced runs",
                    transform=ax.transAxes, ha="center", fontsize=9,
                    color=WONG["grey"])
        publishable_axes(ax, grid="both", grid_alpha=0.20)

    for ax in axes_flat[len(TRAIT_ORDER):]:
        ax.axis("off")

    fig.suptitle("Top-quintile classifier ROC curves "
                  "(Shiu Lab ML-Pipeline, 5-fold CV x 20 balanced reps)",
                  fontsize=12, y=1.01)
    written = save_figure(fig, OUT / "fig_shiu_classifier_roc_curves")
    plt.close(fig)
    print("Wrote", *written)

    if summary_rows:
        sdf = pd.DataFrame(summary_rows).pivot(index="trait", columns="algo",
                                                  values="auc")
        sdf = sdf.reindex(TRAIT_ORDER).reindex(columns=ALGO_ORDER)
        sdf.to_csv(PROJ / "results/83_shiu_lab_classifier/tables/auc_from_roc.csv")
        print("\nPer-trait per-algo AUC computed directly from pooled mean scores:")
        print(sdf.round(3).to_string())


if __name__ == "__main__":
    main()
