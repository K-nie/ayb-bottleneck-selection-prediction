"""
42fig -- house-style replot of the two multi-task supp figures (B15, B28)
from the saved result tables written by 42_ml_multi_task.py.

The parent script recomputes the multi-task CV benchmark on every run, so it
is not re-run just to restyle figures. This reads the frozen CSVs and
re-renders in the shared _figstyle house style with readable group names and
spaced trait labels.

Inputs  (results/42_ml_multi_task/tables/):
  multi_task_per_trait_r.csv  -> B15 per-group forest
  delta_r_vs_single.csv       -> B28 delta-r heatmap vs single-trait RR-BLUP
Outputs (results/42_ml_multi_task/figures/): overwrites fig93 / fig94.

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
TAB = PROJ / "results/42_ml_multi_task/tables"
FIG = PROJ / "results/42_ml_multi_task/figures"
FIG.mkdir(parents=True, exist_ok=True)

GROUPS = {
    "A_oxalates":  (["Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate"],
                    "Oxalate fractions"),
    "B_prot_ox":   (["Crude_Protein", "Total_Oxalate", "Soluble_Oxalate",
                     "Insoluble_Oxalate"], "Protein + oxalates"),
    "C_seed_size": (["Seed_Length", "Seed_Width", "Seed_Thickness",
                     "Mass_of_Seeds"], "Seed size & mass"),
    "D_antiox":    (["Antioxidant", "Flavonoid", "Phenol"], "Antioxidant panel"),
}
MODEL_ORDER = ["Single RR-BLUP", "Single XGBoost", "Multi-trait GBLUP",
               "Multi-output XGB", "Shared-trunk MLP"]
PALETTE = {
    "Single RR-BLUP":    WONG["grey"],
    "Single XGBoost":    "#c9c9c9",
    "Multi-trait GBLUP": WONG["blue"],
    "Multi-output XGB":  WONG["vermillion"],
    "Shared-trunk MLP":  WONG["green"],
}


def nice(t: str) -> str:
    return str(t).replace("_", " ")


def fig_forest(df: pd.DataFrame) -> None:
    heights = [len(v[0]) for v in GROUPS.values()]
    fig, axes = plt.subplots(
        len(GROUPS), 1, figsize=(10, 1.6 + 1.35 * sum(heights)),
        constrained_layout=True,
        gridspec_kw={"height_ratios": heights})
    n_reps = int(df["n_reps"].dropna().iloc[0]) if "n_reps" in df else 0
    for ax, (g, (traits, label)) in zip(axes, GROUPS.items()):
        sub = df[df["group"] == g]
        n_t, height = len(traits), 0.13
        y_pos = np.arange(n_t)[::-1]
        for i, m in enumerate(MODEL_ORDER):
            s = sub[sub["model"] == m].set_index("trait").reindex(traits)
            y = y_pos + (i - (len(MODEL_ORDER) - 1) / 2) * height
            ax.errorbar(s["r_mean"], y,
                        xerr=[s["r_mean"] - s["r_ci_lo"],
                              s["r_ci_hi"] - s["r_mean"]],
                        fmt="o", color=PALETTE[m], markersize=5,
                        elinewidth=0.8, capsize=2, label=m)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([nice(t) for t in traits], fontsize=9)
        ax.set_ylim(-0.6, n_t - 0.4)
        n_g = int(sub["n"].iloc[0]) if len(sub) else 0
        ax.set_title(f"{label}  (n = {n_g})", fontsize=9.5, loc="left")
        publishable_axes(ax, grid="x")
        if ax is axes[-1]:
            ax.set_xlabel(f"Pearson r (5-fold CV $\\times$ {n_reps} repeats)")
        if ax is axes[0]:
            ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
                      fontsize=7.5, title="Model", title_fontsize=8,
                      labelspacing=0.35, frameon=False)
    fig.suptitle("Multi-task vs single-trait genomic prediction "
                 "(AYB DArTseq panel)", fontsize=12)
    print("[fig]", *save_figure(fig, FIG / "fig93_multi_task_forest"))
    plt.close(fig)


def fig_delta(delta: pd.DataFrame) -> None:
    if delta.empty:
        return
    grp_label = {k: v[1] for k, v in GROUPS.items()}
    delta = delta.copy()
    delta["trait_grp"] = (delta["group"].map(grp_label).fillna(delta["group"])
                          + " :: " + delta["trait"].map(nice))
    pivot = delta.pivot(index="trait_grp", columns="model",
                        values="delta_r_mean")
    fig, ax = plt.subplots(figsize=(7.5, max(4, 0.38 * len(pivot))),
                           constrained_layout=True)
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-0.10, vmax=0.10,
                   aspect="auto")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=8)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.3f}", ha="center", va="center",
                        fontsize=7, color="white" if abs(v) > 0.06 else "black")
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02,
                 label="$\\Delta$r vs single-trait RR-BLUP")
    ax.set_title("Multi-task gain over single-trait RR-BLUP "
                 "(positive = better)", loc="left", fontsize=10.5)
    print("[fig]", *save_figure(fig, FIG / "fig94_delta_r_heatmap"))
    plt.close(fig)


def main() -> None:
    fig_forest(pd.read_csv(TAB / "multi_task_per_trait_r.csv"))
    fig_delta(pd.read_csv(TAB / "delta_r_vs_single.csv"))


if __name__ == "__main__":
    main()
