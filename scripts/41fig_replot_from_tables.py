"""
41fig -- house-style replot of the two SHAP supp figures that can be rebuilt
from saved tables (B26 SHAP-vs-GWAS overlap, B27 ALMT4 rank focus).

The parent script (41_ml_shap_interpretability.py) refits models and recomputes
SHAP on every run, so it is not re-run to restyle. The per-sample SHAP value
matrices needed for the beeswarms (B24/B25) were NOT saved and are therefore
out of scope here; they require a targeted SHAP recompute.

Fixes applied:
  - the intersection glyph rendered as tofu (missing from Helvetica); replaced
    with mathtext $\\cap$
  - trait labels spaced, not underscored
  - _figstyle apply()/publishable_axes()/save_figure()

Inputs  (results/41_ml_shap_interpretability/tables/):
  shap_gwas_overlap.csv  -> B26 bar chart (observed vs chance intersection)
  almt4_focus.csv        -> B27 SHAP rank vs GWAS rank scatter
Outputs (results/41_ml_shap_interpretability/figures/): overwrites fig91/fig92.

Author: Benjamin Narh-Madey
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import (apply, WONG, publishable_axes, save_figure,
                       adjust_labels)
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAB = PROJ / "results/41_ml_shap_interpretability/tables"
FIG = PROJ / "results/41_ml_shap_interpretability/figures"
FIG.mkdir(parents=True, exist_ok=True)


def nice(t: str) -> str:
    return str(t).replace("_", " ")


def fig_overlap(ov: pd.DataFrame) -> None:
    if ov.empty:
        return
    n_shap = int(ov["n_shap_top"].iloc[0])
    n_gwas = int(ov["n_gwas_top"].iloc[0])
    ov = ov.sort_values("n_intersection", ascending=False).reset_index(drop=True)
    xpos = np.arange(len(ov))
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ax.bar(xpos, ov["n_intersection"], width=0.62, color=WONG["vermillion"],
           edgecolor="white", linewidth=0.4,
           label=rf"SHAP top-{n_shap} $\cap$ GWAS top-{n_gwas}")
    ax.plot(xpos, ov["expected_by_chance"], ls="none", marker="_",
            markersize=15, markeredgewidth=1.6, color=WONG["grey"],
            label="Expected by chance")
    ax.set_xticks(xpos)
    ax.set_xticklabels([nice(t) for t in ov["trait"]], rotation=45,
                       ha="right", fontsize=8)
    ax.set_ylabel("Shared marker count")
    ax.set_title(rf"SHAP top-{n_shap} $\cap$ GWAS top-{n_gwas} overlap per "
                 rf"trait (observed vs chance)", loc="left")
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    publishable_axes(ax, grid="y")
    print("[fig]", *save_figure(fig, FIG / "fig91_shap_gwas_overlap_heatmap"))
    plt.close(fig)


def fig_almt4(af: pd.DataFrame) -> None:
    af = af.dropna(subset=["gwas_rank", "shap_rank"]).copy()
    if af.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 6.5), constrained_layout=True)
    ax.scatter(af["gwas_rank"], af["shap_rank"], s=70, color=WONG["blue"],
               edgecolor="black", linewidth=0.5, alpha=0.85, zorder=4)
    texts = [ax.text(r["gwas_rank"], r["shap_rank"], nice(r["trait"]),
                     fontsize=8, color="#222222", zorder=10)
             for _, r in af.iterrows()]
    lim = max(af["gwas_rank"].max(), af["shap_rank"].max()) * 1.05
    ax.plot([1, lim], [1, lim], "k--", linewidth=0.6, alpha=0.5,
            label="Equal-rank line")
    ax.set_xscale("log"); ax.set_yscale("log")
    adjust_labels(texts, ax=ax, expand=(1.3, 1.5), force_text=(0.9, 1.2),
                  only_move={"text": "xy", "static": "xy"})
    ax.set_xlabel("GWAS rank (ALMT4 SNP, Ss10:15.4 Mb)")
    ax.set_ylabel("SHAP rank")
    ax.set_title("ALMT4 SNP (Ss10:15.4 Mb): SHAP rank vs GWAS rank per trait",
                 loc="left")
    ax.legend(fontsize=8.5, framealpha=0.9, edgecolor="none")
    publishable_axes(ax, grid="both", grid_alpha=0.25)
    print("[fig]", *save_figure(fig, FIG / "fig92_almt4_focus"))
    plt.close(fig)


def main() -> None:
    fig_overlap(pd.read_csv(TAB / "shap_gwas_overlap.csv"))
    fig_almt4(pd.read_csv(TAB / "almt4_focus.csv"))


if __name__ == "__main__":
    main()
