"""Standalone re-render of fig86_ml_vs_daetwyler_ceiling from the cached
table. Avoids re-running the full nested-CV ML pipeline in script 40.
Mirrors the same rendering block as script 40 (with the polish edits applied).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _plotstyle import apply, WONG
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAB  = PROJ / "results/40_ml_genomic_prediction/tables"
FIG  = PROJ / "results/40_ml_genomic_prediction/figures"

ceil_df = pd.read_csv(TAB / "ml_vs_daetwyler.csv").sort_values(
    "daetwyler_ceiling").reset_index(drop=True)

fig, ax = plt.subplots(figsize=(11, 6.0), constrained_layout=True)
xpos = np.arange(len(ceil_df))
ax.bar(xpos, ceil_df["daetwyler_ceiling"], width=0.6,
       color="#dddddd", edgecolor="grey", label="Daetwyler ceiling")
ax.errorbar(xpos, ceil_df["best_r"],
            yerr=[ceil_df["best_r"] - ceil_df["best_r_lo"],
                  ceil_df["best_r_hi"] - ceil_df["best_r"]],
            fmt="o", color=WONG["vermillion"],
            markersize=5, elinewidth=0.8, capsize=2,
            label="Best ML model r (95 % CI)")
for i, m in enumerate(ceil_df["best_model"].values):
    bar_top = float(ceil_df["daetwyler_ceiling"].iloc[i])
    ci_top  = float(ceil_df["best_r_hi"].iloc[i])
    y_text  = max(bar_top, ci_top) + 0.025
    ax.text(xpos[i], y_text, m, rotation=90, ha="center", va="bottom",
            fontsize=8, color="grey")
ax.set_xticks(xpos)
ax.set_xticklabels(ceil_df["trait"].values, rotation=35, ha="right",
                    fontsize=9)
ax.set_ylabel("Pearson r")
ax.axhline(0, color="black", linewidth=0.4)
ax.set_title("Best ML model vs Daetwyler ceiling per trait")
ax.legend(loc="upper left", fontsize=8)
ax.grid(True, axis="y", alpha=0.4)
fig.savefig(FIG / "fig86_ml_vs_daetwyler_ceiling.png", dpi=300,
            bbox_inches="tight")
fig.savefig(FIG / "fig86_ml_vs_daetwyler_ceiling.pdf",
            bbox_inches="tight")
plt.close(fig)
print(f"[fig] re-rendered {FIG / 'fig86_ml_vs_daetwyler_ceiling.png'}")
