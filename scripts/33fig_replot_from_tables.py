"""
33fig -- house-style standalone replot of B16, the Stairway Plot 2 N_e(t)
trajectory, from the saved table results/33_stairway/tables/ne_through_time.csv.

B16 is the standalone companion to B4 (the SFS + Stairway two-panel figure).
It carries the same Panel-B semantics as the reworked B4 (77_sfs_stairway_overlay):
log-log axes, the 75 % and 95 % posterior bands, the median trajectory, and the
LD-based contemporary N_e comparator. The parent script (33_stairway.py) runs the
full Stairway Plot 2 Java pipeline on every invocation, so it is not re-run to
restyle; this reads the frozen summary CSV.

Fix carried over from B4: Stairway emits one spurious most-recent grid point at
year ~1e-316 (float underflow from mutation_per_site = 5e-324). Left in, it drags
the log x-axis down to ~1e-316 and squashes the real trajectory. Floor the time
axis to drop only that artifact.

Input : results/33_stairway/tables/ne_through_time.csv
Output: results/33_stairway/figures/fig79_stairway_ne_time.{png,pdf}
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, publishable_axes, save_figure
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
NE = PROJ / "results/33_stairway/tables/ne_through_time.csv"
OUT = PROJ / "results/33_stairway/figures"
OUT.mkdir(parents=True, exist_ok=True)

# LD-derived contemporary N_e reference (results/11_ld_decay), the dashed
# comparator the abstract cites.
LD_NE = 659.0
YEAR_FLOOR = 1e-6


def main() -> None:
    ne = pd.read_csv(NE)
    ne = ne.dropna(subset=["year", "Ne_median"]).copy()
    ne["year"] = ne["year"].astype(float)
    ne["Ne_median"] = ne["Ne_median"].astype(float)
    ne = ne[ne["year"] >= YEAR_FLOOR].drop_duplicates("year").sort_values("year")
    print(f"Stairway grid points after year floor ({YEAR_FLOOR}): {len(ne)}; "
          f"year {ne['year'].min():.4f}..{ne['year'].max():.0f}")

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    if {"Ne_2.5%", "Ne_97.5%"}.issubset(ne.columns):
        ax.fill_between(ne["year"], ne["Ne_2.5%"], ne["Ne_97.5%"],
                        color=WONG["blue"], alpha=0.15, lw=0,
                        label="95 % posterior")
    if {"Ne_12.5%", "Ne_87.5%"}.issubset(ne.columns):
        ax.fill_between(ne["year"], ne["Ne_12.5%"], ne["Ne_87.5%"],
                        color=WONG["blue"], alpha=0.28, lw=0,
                        label="75 % posterior")
    ax.plot(ne["year"], ne["Ne_median"], "-", color=WONG["blue"], lw=1.9,
            label="$N_e$ median (Stairway Plot 2)")
    ax.axhline(LD_NE, color=WONG["vermillion"], lw=1.1, ls="--",
               label=f"LD-based $N_e$ = {LD_NE:.0f}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.9, ne["year"].max() * 1.4)
    ax.set_xlabel("Years before present (log scale, 1 yr per generation)")
    ax.set_ylabel("$N_e$ (effective population size, log scale)")
    ax.set_title("Stairway Plot 2 $N_e(t)$ trajectory "
                 "(African yam bean, n = 95 IITA TSs panel)", loc="left")
    publishable_axes(ax, grid=None)
    ax.grid(True, which="both", color=WONG["grey"], lw=0.4, alpha=0.20)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left")

    fig.tight_layout()
    written = save_figure(fig, OUT / "fig79_stairway_ne_time")
    plt.close(fig)
    print("Wrote", *written)


if __name__ == "__main__":
    main()
