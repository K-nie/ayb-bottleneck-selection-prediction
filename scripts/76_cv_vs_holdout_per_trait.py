"""
76 -- CV r vs holdout r per-trait, per-model -- the Bellot 2018 pattern.

For each (trait, model) pair, plot mean CV r on x against final holdout r on
y. Points on the y = x diagonal indicate honest CV that tracks held-out
performance; points strongly off the diagonal indicate either CV over-optimism
(Bellot 2018 small-n deep-net failure) or holdout luck (n_hold = 4 -- 14
fragility).

Colour by model (Wong palette), marker by trait grouping (oxalate / seed-size
/ other). Trait labels placed via adjustText so no overlap. Legend pulled
beneath the figure into a single horizontal block.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, adjust_labels, publishable_axes, panel_label
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
CV = PROJ / "results/40_ml_genomic_prediction/tables/cv_predictive_ability.csv"
HO = PROJ / "results/40_ml_genomic_prediction/tables/final_holdout_scores.csv"
OUT = PROJ / "results/76_cv_vs_holdout"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

OXALATE = {"Soluble_Oxalate", "Insoluble_Oxalate", "Total_Oxalate"}
SEED = {"Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds"}

MODEL_COLOUR = {
    "Dummy":        "#999999",
    "RR-BLUP":      WONG["blue"],
    "ElasticNet":   WONG["purple"],
    "RandomForest": WONG["green"],
    "XGBoost":      WONG["vermillion"],
    "LightGBM":     WONG["orange"],
    "KernelRidge":  WONG["skyblue"],
    "MLP":          WONG["yellow"],
}


def trait_marker(t: str) -> str:
    if t in OXALATE:
        return "o"
    if t in SEED:
        return "s"
    return "^"


def main() -> None:
    cv = pd.read_csv(CV)
    ho = pd.read_csv(HO)

    merged = ho.merge(
        cv[["trait", "model", "r_mean", "r_ci_lo", "r_ci_hi", "n_work"]],
        left_on=["trait", "best_model"], right_on=["trait", "model"],
        how="left",
    )
    merged.to_csv(OUT / "tables/cv_vs_holdout_winners.csv", index=False)
    print(f"Loaded {len(merged)} (trait, best-model) winner pairs.")
    print(merged[["trait", "best_model", "n_work", "cv_r_mean", "holdout_r",
                  "n_hold"]])

    all_pairs = cv.merge(
        ho[["trait", "best_model"]].rename(columns={"best_model": "winner"}),
        on="trait", how="left",
    )
    all_pairs.to_csv(OUT / "tables/cv_vs_holdout_all_models.csv", index=False)

    fig = plt.figure(figsize=(14, 7.5))
    gs = fig.add_gridspec(
        2, 2, height_ratios=[6.0, 1.2], width_ratios=[1.0, 1.0],
        hspace=0.20, wspace=0.18,
    )

    # ---- Panel A: CV vs holdout for each trait's winner ----
    ax = fig.add_subplot(gs[0, 0])
    point_handles = {}
    trait_label_texts = []
    for _, r in merged.iterrows():
        col = MODEL_COLOUR.get(r["best_model"], "#666")
        mk = trait_marker(r["trait"])
        xerr_lo = max(0.0, float(r["cv_r_mean"] - r["r_ci_lo"]))
        xerr_hi = max(0.0, float(r["r_ci_hi"] - r["cv_r_mean"]))
        h = ax.errorbar(
            r["cv_r_mean"], r["holdout_r"],
            xerr=[[xerr_lo], [xerr_hi]],
            fmt=mk, color=col, ecolor=col, ms=11,
            markeredgecolor="black", markeredgewidth=0.7,
            elinewidth=0.9, alpha=0.92, capsize=2.5, zorder=4,
        )
        point_handles[r["best_model"]] = h
        trait_label_texts.append(
            ax.text(r["cv_r_mean"], r["holdout_r"],
                    f"{r['trait'].replace('_', ' ')}",
                    fontsize=8.5, color="#222", zorder=10))

    lim = float(max(abs(merged["cv_r_mean"]).max(),
                    abs(merged["holdout_r"]).max()) * 1.20)
    ax.plot([-lim, lim], [-lim, lim], color="black",
            linestyle="--", linewidth=0.9,
            label="y = x (honest CV)", zorder=2)
    ax.axhline(0, color="#666", linewidth=0.5, zorder=1)
    ax.axvline(0, color="#666", linewidth=0.5, zorder=1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"CV r (mean, 5 $\times$ 5 $\times$ 10 nested CV)",
                  fontsize=10.5)
    ax.set_ylabel(r"Final holdout r ($n_{\mathrm{hold}}$ per trait)",
                  fontsize=10.5)
    ax.set_title("Per-trait CV-winner: CV r vs holdout r",
                  fontsize=11)
    panel_label(ax, "a")
    publishable_axes(ax, grid="both", grid_alpha=0.20)
    adjust_labels(trait_label_texts, ax=ax,
                  expand=(1.35, 1.55), force_text=(1.0, 1.2),
                  only_move={"text": "xy", "static": "xy"})

    # ---- Panel B: per-trait CV-holdout gap, sorted ----
    ax = fig.add_subplot(gs[0, 1])
    merged["gap"] = merged["cv_r_mean"] - merged["holdout_r"]
    merged_sorted = merged.sort_values("gap")
    colours = [MODEL_COLOUR.get(m, "#666")
               for m in merged_sorted["best_model"]]
    y_pos = np.arange(len(merged_sorted))
    ax.barh(y_pos, merged_sorted["gap"], color=colours, alpha=0.85,
            edgecolor="black", linewidth=0.5, zorder=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [f"{r['trait'].replace('_', ' ')} | {r['best_model']}"
         for _, r in merged_sorted.iterrows()],
        fontsize=8.5,
    )
    ax.axvline(0, color="black", linewidth=0.8, zorder=2)
    ax.set_xlabel(r"CV r minus holdout r (positive = CV over-optimism)",
                  fontsize=10.5)
    ax.set_title("CV-vs-holdout gap, sorted", fontsize=11)
    panel_label(ax, "b")
    publishable_axes(ax, grid="x", grid_alpha=0.20)

    # ---- Combined bottom legend ----
    ax_leg = fig.add_subplot(gs[1, :])
    ax_leg.axis("off")
    model_handles = [
        mpatches.Patch(facecolor=col, edgecolor="black",
                       label=m, alpha=0.85, linewidth=0.6)
        for m, col in MODEL_COLOUR.items() if m in merged["best_model"].values
    ]
    marker_handles = [
        plt.Line2D([], [], marker="o", color="#444", ls="", ms=9,
                   markeredgecolor="black", markeredgewidth=0.5,
                   label="oxalate triplet"),
        plt.Line2D([], [], marker="s", color="#444", ls="", ms=9,
                   markeredgecolor="black", markeredgewidth=0.5,
                   label="seed metrics"),
        plt.Line2D([], [], marker="^", color="#444", ls="", ms=9,
                   markeredgecolor="black", markeredgewidth=0.5,
                   label="other traits"),
    ]
    diag_handle = plt.Line2D([], [], color="black", linestyle="--",
                              linewidth=0.9, label="y = x (honest CV)")
    leg1 = ax_leg.legend(
        handles=model_handles, loc="upper center", frameon=False,
        fontsize=9.5, title="Model (panel a colour, panel b bar fill)",
        title_fontsize=10, ncol=min(len(model_handles), 7),
        bbox_to_anchor=(0.5, 1.0),
    )
    ax_leg.add_artist(leg1)
    ax_leg.legend(
        handles=marker_handles + [diag_handle],
        loc="upper center", frameon=False, fontsize=9.5,
        title="Trait group (panel a marker) + panel a reference",
        title_fontsize=10, ncol=4,
        bbox_to_anchor=(0.5, 0.40),
    )

    fig.suptitle(
        "Cross-validation vs final-holdout per-trait winners "
        "(Bellot 2018 small-n pattern)",
        y=0.99, fontsize=12.5,
    )
    fig.savefig(OUT / "figures/fig_cv_vs_holdout.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "figures/fig_cv_vs_holdout.pdf",
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Wrote", OUT / "figures/fig_cv_vs_holdout.png")


if __name__ == "__main__":
    main()
