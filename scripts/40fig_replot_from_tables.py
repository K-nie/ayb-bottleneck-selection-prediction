"""
40fig -- house-style replot of the four genomic-prediction supp figures
(B20-B23) directly from the saved result tables written by
40_ml_genomic_prediction.py.

The parent script recomputes the full nested-CV benchmark (KMeans + per-fold
model fits + SHAP) on every run, so it is not re-run just to restyle figures.
This script reads the frozen CSV outputs and re-renders them in the shared
_figstyle house style with the labelling fixes:
  - trait labels rendered with spaces, not underscores
  - the forest legend moved outside the plotting area (was overlapping the
    bottom trait rows)
  - _figstyle apply()/publishable_axes()/save_figure() for consistent 300-dpi
    PNG+PDF output

Inputs  (results/40_ml_genomic_prediction/tables/):
  cv_predictive_ability.csv   -> B20 model-comparison forest
  ml_vs_daetwyler.csv         -> B21 best-ML vs Daetwyler ceiling
  top_k_accuracy.csv          -> B22 top-10 ranking-accuracy heatmap
  cohens_d_delta_r.csv        -> B23 Cohen's d on delta-r vs RR-BLUP

Outputs (results/40_ml_genomic_prediction/figures/): overwrites the four
fig85..fig88 basenames the manifest already points at.

Author: Benjamin Narh-Madey
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, publishable_axes, save_figure
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAB = PROJ / "results/40_ml_genomic_prediction/tables"
FIG = PROJ / "results/40_ml_genomic_prediction/figures"
FIG.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "Dummy":        WONG["grey"],
    "RR-BLUP":      WONG["blue"],
    "ElasticNet":   WONG["skyblue"],
    "RandomForest": WONG["green"],
    "XGBoost":      WONG["vermillion"],
    "LightGBM":     WONG["orange"],
    "KernelRidge":  WONG["purple"],
    "MLP":          WONG["yellow"],
}
MODEL_ORDER = ["Dummy", "RR-BLUP", "ElasticNet", "RandomForest", "XGBoost",
               "LightGBM", "KernelRidge", "MLP"]


def nice(t: str) -> str:
    return str(t).replace("_", " ")


def fig_forest(cv: pd.DataFrame) -> None:
    traits = list(cv["trait"].unique())
    n_t, n_m, height = len(traits), len(MODEL_ORDER), 0.10
    y_pos = np.arange(n_t)[::-1]
    n_reps = int(cv["n_reps"].dropna().iloc[0]) if "n_reps" in cv else 0

    fig, ax = plt.subplots(figsize=(11, 8), constrained_layout=True)
    for i, m in enumerate(MODEL_ORDER):
        sub = cv[cv["model"] == m].set_index("trait").reindex(traits)
        y = y_pos + (i - (n_m - 1) / 2) * height
        ax.errorbar(sub["r_mean"], y,
                    xerr=[sub["r_mean"] - sub["r_ci_lo"],
                          sub["r_ci_hi"] - sub["r_mean"]],
                    fmt="o", color=PALETTE[m], markersize=4,
                    elinewidth=0.8, capsize=2, label=m)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([nice(t) for t in traits], fontsize=9)
    ax.set_ylim(-0.6, n_t - 0.4)
    ax.set_xlabel(f"Pearson r (nested 5$\\times$5 CV $\\times$ {n_reps} "
                  f"outer repeats)")
    ax.set_title("Genomic-prediction model comparison (AYB DArTseq panel)",
                 loc="left")
    publishable_axes(ax, grid="x")
    # Legend parked outside the axes so it never sits on the bottom traits.
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8,
              title="Model", title_fontsize=8.5, handletextpad=0.4,
              borderpad=0.4, labelspacing=0.35, frameon=False)
    print("[fig]", *save_figure(fig, FIG / "fig85_model_comparison_forest"))
    plt.close(fig)


def fig_daetwyler(ceil_df: pd.DataFrame) -> None:
    ceil_df = ceil_df.sort_values("daetwyler_ceiling").reset_index(drop=True)
    xpos = np.arange(len(ceil_df))
    fig, ax = plt.subplots(figsize=(11, 6.0), constrained_layout=True)
    ax.bar(xpos, ceil_df["daetwyler_ceiling"], width=0.6,
           color="#e0e0e0", edgecolor=WONG["grey"], linewidth=0.6,
           label="Daetwyler ceiling", zorder=1)
    ax.errorbar(xpos, ceil_df["best_r"],
                yerr=[ceil_df["best_r"] - ceil_df["best_r_lo"],
                      ceil_df["best_r_hi"] - ceil_df["best_r"]],
                fmt="o", color=WONG["vermillion"], markersize=5,
                elinewidth=0.8, capsize=2, zorder=3,
                label="Best ML model r (95 % CI)")
    for i, m in enumerate(ceil_df["best_model"].values):
        y_text = max(float(ceil_df["daetwyler_ceiling"].iloc[i]),
                     float(ceil_df["best_r_hi"].iloc[i])) + 0.03
        ax.text(xpos[i], y_text, m, rotation=90, ha="center", va="bottom",
                fontsize=7.5, color=WONG["grey"])
    ax.set_xticks(xpos)
    ax.set_xticklabels([nice(t) for t in ceil_df["trait"]], rotation=35,
                       ha="right", fontsize=9)
    ax.set_ylabel("Pearson r")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylim(top=float(ceil_df["best_r_hi"].max()) + 0.16)
    ax.set_title("Best ML model vs Daetwyler ceiling per trait", loc="left")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    publishable_axes(ax, grid="y")
    print("[fig]", *save_figure(fig, FIG / "fig86_ml_vs_daetwyler_ceiling"))
    plt.close(fig)


def fig_topk(tk: pd.DataFrame) -> None:
    traits = list(tk["trait"].unique())
    top_k = int(tk["top_k"].dropna().iloc[0]) if "top_k" in tk else 10
    pivot = tk.pivot(index="trait", columns="model",
                     values="topk_acc_mean").reindex(
        index=traits, columns=MODEL_ORDER)
    fig, ax = plt.subplots(figsize=(8.5, 7), constrained_layout=True)
    im = ax.imshow(pivot.values, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(MODEL_ORDER)))
    ax.set_xticklabels(MODEL_ORDER, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(traits)))
    ax.set_yticklabels([nice(t) for t in traits], fontsize=8)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=7, color="white" if v < 0.5 else "black")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02,
                 label=f"top-{top_k} overlap")
    ax.set_title(f"Top-{top_k} ranking accuracy (predicted vs observed)",
                 loc="left")
    print("[fig]", *save_figure(fig, FIG / "fig87_top_k_accuracy_heatmap"))
    plt.close(fig)


def fig_cohens(coh: pd.DataFrame) -> None:
    traits = list(dict.fromkeys(coh["trait"]))
    coh_p = coh.pivot(index="trait", columns="comparison",
                      values="cohens_d").reindex(index=traits)
    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    im = ax.imshow(coh_p.values, cmap="RdBu_r", vmin=-1.5, vmax=1.5,
                   aspect="auto")
    ax.set_xticks(np.arange(coh_p.shape[1]))
    ax.set_xticklabels([c.replace(" vs RR-BLUP", "") for c in coh_p.columns],
                       rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(coh_p.shape[0]))
    ax.set_yticklabels([nice(t) for t in coh_p.index], fontsize=8)
    for i in range(coh_p.shape[0]):
        for j in range(coh_p.shape[1]):
            v = coh_p.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                        fontsize=7, color="white" if abs(v) > 0.8 else "black")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02,
                 label="Cohen's d on $\\Delta$r")
    ax.set_title("Effect size on $\\Delta$r vs RR-BLUP (positive = ML better)",
                 loc="left")
    print("[fig]", *save_figure(fig, FIG / "fig88_wilcoxon_effect_size"))
    plt.close(fig)


def main() -> None:
    fig_forest(pd.read_csv(TAB / "cv_predictive_ability.csv"))
    fig_daetwyler(pd.read_csv(TAB / "ml_vs_daetwyler.csv"))
    fig_topk(pd.read_csv(TAB / "top_k_accuracy.csv"))
    coh = pd.read_csv(TAB / "cohens_d_delta_r.csv")
    if not coh.empty:
        fig_cohens(coh)


if __name__ == "__main__":
    main()
