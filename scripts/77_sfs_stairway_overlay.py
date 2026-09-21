"""
77 -- Folded SFS + Stairway demographic-fit overlay.

Two-panel figure that anchors the Paper B bottleneck story:

(A) Observed folded SFS bars at k = 1 .. n/2 alongside the neutral-equilibrium
    expectation E[xi_k] = theta_W * (1/k + 1/(n-k)) / (1 + delta_{k, n-k}).
    Excess at intermediate frequencies = bottleneck signal (Tajima D > 0).

(B) Stairway-inferred N_e(t) trajectory on a log time axis with the median
    and 95 % posterior band.

Saves SFS observed-vs-expected, the Tajima D-style summary statistic, and
the figure.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, publishable_axes, save_figure, panel_label
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
SFS = PROJ / "results/33_stairway/data/sfs.txt"
NE = PROJ / "results/33_stairway/tables/ne_through_time.csv"
OUT = PROJ / "results/77_sfs_stairway_overlay"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

# LD-derived contemporary N_e reference (results/11_ld_decay), overlaid on
# Panel B as the dashed comparator the abstract cites.
LD_NE = 659.0
# Stairway Plot 2 emits one spurious most-recent grid point at year ~1e-316
# (a float underflow from mutation_per_site = 5e-324). It is > 0 so a naive
# `year > 0` filter keeps it and drags the log x-axis down to ~1e-316,
# squashing the real trajectory. Floor the time axis at a physically sane
# value to drop only that artifact.
YEAR_FLOOR = 1e-6


def main() -> None:
    raw = open(SFS).read().strip().split()
    counts = np.array([float(x) for x in raw])
    # Convention: first entry = monomorphic (k = 0); k from 1 .. m where m = n/2
    # in the folded SFS. counts has length ~ n (95 entries here for n = 190 hap).
    # Drop the leading 0 (k = 0) and keep the first m = (len - 1) entries.
    sfs = counts[1:]
    print(f"Folded SFS length = {len(sfs)} (drop k=0 monomorphic)")
    # If the data were already folded with the monomorphic class stripped,
    # this still works -- we are using the segregating-site portion.
    seg = int(sfs.sum())
    print(f"Total segregating sites = {seg}")

    n_chr = 2 * (len(sfs))  # assuming folded SFS for n_chr haploids
    # Note: the actual panel has 95 diploids -> 190 haplotypes. Folded SFS has
    # ceil(190 / 2) = 95 bins for k = 1 .. 95.
    print(f"Inferred haplotype count n_chr = {n_chr} (folded over 1..{n_chr // 2})")

    # Watterson estimator at this n_chr:
    harmonic = sum(1.0 / i for i in range(1, n_chr))
    theta_S = seg / harmonic  # per-locus theta_W (in segregating-site units)
    print(f"Harmonic a_n = {harmonic:.4f}")
    print(f"theta_W (segregating-site units, summed across all loci) = {theta_S:.4f}")

    # Expected folded SFS counts under neutral equilibrium:
    # E[xi_k_folded] = theta * (1/k + 1/(n-k)) / (1 + delta_{k, n-k})
    n = n_chr
    m = len(sfs)  # 95
    exp_counts = np.zeros(m)
    for k in range(1, m + 1):
        if 2 * k == n:
            exp_counts[k - 1] = theta_S * (1.0 / k) * 0.5
        else:
            exp_counts[k - 1] = theta_S * (1.0 / k + 1.0 / (n - k)) * 0.5
    # rescale expected so it matches total observed segregating sites
    exp_counts = exp_counts * (seg / exp_counts.sum())

    df_sfs = pd.DataFrame({
        "k": np.arange(1, m + 1),
        "observed": sfs.astype(int),
        "expected_neutral": exp_counts,
        "obs_minus_exp": sfs - exp_counts,
        "log2_ratio": np.log2((sfs + 0.5) / (exp_counts + 0.5)),
    })
    df_sfs.to_csv(OUT / "tables/observed_vs_neutral_sfs.csv", index=False)

    # Stairway Ne(t)
    ne = pd.read_csv(NE)
    ne = ne.dropna(subset=["year", "Ne_median"]).copy()
    ne["year"] = ne["year"].astype(float)
    ne["Ne_median"] = ne["Ne_median"].astype(float)
    ne = ne[ne["year"] >= YEAR_FLOOR].drop_duplicates("year").sort_values("year")
    print(f"Stairway grid points after year floor ({YEAR_FLOOR}): {len(ne)}; "
          f"year {ne['year'].min():.4f}..{ne['year'].max():.0f}")

    # ---------- FIGURE ----------
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.0))

    # Panel A: SFS observed vs expected
    ax = axes[0]
    k_vals = df_sfs["k"].values
    width = 0.85
    ax.bar(k_vals, df_sfs["observed"], width=width, color=WONG["blue"],
           alpha=0.75, edgecolor="white", lw=0.3, label="Observed (folded SFS)")
    ax.plot(k_vals, df_sfs["expected_neutral"], "o-", color=WONG["vermillion"],
            ms=3, lw=1.3, label=r"Neutral equilibrium $E[\xi_k]$ "
                                r"($\theta_W$ scaled)")
    ax.set_xlabel(f"Minor-allele count $k$ (folded, $n_{{chr}}$ = {n})")
    ax.set_ylabel("SFS bin count")
    ax.set_title("Folded SFS vs neutral-equilibrium expectation", loc="left")
    panel_label(ax, "a")
    ax.legend(loc="upper right")
    publishable_axes(ax, grid="y")

    # Inset: log2 ratio at intermediate frequencies (k = 5..m). Placed in
    # the bottom-right of Panel A so the inset title cannot collide with the
    # parent panel title at the top.
    ax_in = ax.inset_axes([0.50, 0.20, 0.42, 0.32])
    ax_in.bar(k_vals[4:], df_sfs["log2_ratio"].values[4:], width=width,
              color=WONG["green"], alpha=0.85, edgecolor="white", lw=0.3)
    ax_in.axhline(0, color="black", lw=0.6)
    ax_in.set_title(r"log$_2$(obs / expected), $k \geq 5$", fontsize=7.5,
                     pad=2)
    ax_in.tick_params(labelsize=6.5)
    ax_in.set_facecolor("white")
    ax_in.patch.set_alpha(0.92)
    for side in ("top", "right"):
        ax_in.spines[side].set_visible(False)

    # Panel B: Ne(t) Stairway
    ax = axes[1]
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
    # LD-derived contemporary Ne comparator (the curve crosses it ~30 yr BP)
    ax.axhline(LD_NE, color=WONG["vermillion"], lw=1.1, ls="--",
               label=f"LD-based $N_e$ = {LD_NE:.0f}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.9, ne["year"].max() * 1.4)
    ax.set_xlabel("Years before present (log scale, 1 yr per generation)")
    ax.set_ylabel("$N_e$ (effective population size, log scale)")
    ax.set_title("Stairway Plot 2 $N_e(t)$ trajectory", loc="left")
    panel_label(ax, "b")
    publishable_axes(ax, grid=None)
    ax.grid(True, which="both", color=WONG["grey"], lw=0.4, alpha=0.20)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left")

    fig.suptitle("African yam bean folded SFS + Stairway $N_e(t)$ "
                 "(n = 95 IITA TSs panel)", y=1.00)
    fig.tight_layout()
    written = save_figure(fig, OUT / "figures/fig_sfs_stairway_overlay")
    plt.close(fig)
    print("Wrote", *written)

    # quick most-recent Ne summary
    ne_recent = ne.head(50)
    print(f"\nMost-recent Ne_median (closest to present, first 50 grid points): "
          f"median = {float(ne_recent['Ne_median'].median()):.0f}")


if __name__ == "__main__":
    main()
