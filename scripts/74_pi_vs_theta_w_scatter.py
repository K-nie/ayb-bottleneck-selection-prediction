"""
74 -- pi vs theta_W per-window scatter, with chromosome colouring.

Direct visualisation of the bottleneck signature: under neutrality + constant
N_e, pi ~ theta_W; under a recent population bottleneck, theta_W (driven by
the count of segregating sites) drops faster than pi (driven by mean pairwise
divergence at sites that survive). The result is pi > theta_W in many windows,
which maps onto Tajima D > 0.

Source table `results/35_tajima_d/tables/tajima_d_windowed.csv` carries
`pi_per_site` (already normalised to per-site) and `theta_W` (the raw
Watterson estimator, NOT per-site). We normalise theta_W to per-site by
dividing by window length (win_end - win_start) so the two axes are on
the same physical scale -- otherwise the y = x neutrality diagonal is
meaningless and the scatter collapses off-axis.

Plot:
  x = theta_W per site
  y = pi per site
  point colour = chromosome (Wong-palette tab11)
  point size = number of SNPs in window
  reference y = x line (neutrality)
  marginal histograms on each axis
  bottom panel: Tajima D distribution by chromosome (box plot)
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, WONG_CYCLE, publishable_axes, panel_label
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAJ = PROJ / "results/35_tajima_d/tables/tajima_d_windowed.csv"
OUT = PROJ / "results/74_pi_vs_thetaW"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)


def main() -> None:
    df = pd.read_csv(TAJ)
    print(f"Loaded {len(df)} windows from {TAJ.name}.")
    df = df.dropna(subset=["pi_per_site", "theta_W", "tajima_D"]).copy()
    df["pi_per_site"] = df["pi_per_site"].astype(float)
    # theta_W in source is per-window (counts of segregating sites);
    # normalise to per-site so pi and theta_W are on the same scale.
    df["window_bp"] = (df["win_end"] - df["win_start"]).astype(float)
    df["theta_W_per_site"] = df["theta_W"].astype(float) / df["window_bp"]
    df["tajima_D"] = df["tajima_D"].astype(float)

    summary = (
        df.groupby("chr")
        .agg(
            n_windows=("pi_per_site", "size"),
            pi_median=("pi_per_site", "median"),
            theta_w_per_site_median=("theta_W_per_site", "median"),
            tajima_d_median=("tajima_D", "median"),
            frac_d_gt_0=("tajima_D", lambda s: float((s > 0).mean())),
            frac_d_gt_2=("tajima_D", lambda s: float((s > 2).mean())),
        )
        .reset_index()
    )
    summary.to_csv(OUT / "tables/per_chr_summary.csv", index=False)
    print(summary.round(6).to_string(index=False))

    chroms = sorted(df["chr"].unique(),
                    key=lambda s: int(''.join(c for c in s if c.isdigit())))
    # cycle through the Wong palette so the 11 chromosomes are distinguishable
    chr_colour = {c: WONG_CYCLE[i % len(WONG_CYCLE)]
                  for i, c in enumerate(chroms)}

    fig = plt.figure(figsize=(11.5, 8.5), constrained_layout=False)
    gs = fig.add_gridspec(
        4, 3, hspace=0.10, wspace=0.05,
        height_ratios=[0.85, 4.5, 0.18, 2.2], width_ratios=[5.5, 1, 1.3]
    )
    ax = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax)

    # scatter with per-chromosome colour + n_snps-scaled point size
    for c in chroms:
        sub = df[df["chr"] == c]
        ax.scatter(
            sub["theta_W_per_site"], sub["pi_per_site"],
            s=np.clip(sub["n_snps"].values * 1.6, 14, 90),
            color=chr_colour[c], alpha=0.78,
            edgecolor="white", linewidth=0.5,
            label=c, zorder=3,
        )

    # axis limits: pad 5 % around log range
    x_vals = df["theta_W_per_site"].values
    y_vals = df["pi_per_site"].values
    lim_lo = max(min(x_vals.min(), y_vals.min()) * 0.7, 1e-7)
    lim_hi = max(x_vals.max(), y_vals.max()) * 1.4
    # y = x neutrality reference line
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi],
            color="black", linestyle="--", linewidth=0.9,
            label=r"$\pi = \theta_W$ (neutral expectation)", zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_xlabel(r"Watterson $\theta_W$ per site", fontsize=10.5)
    ax.set_ylabel(r"Nucleotide diversity $\pi$ per site", fontsize=10.5)
    publishable_axes(ax, grid="both", grid_alpha=0.20)

    # marginal histograms (top: theta_W per site, right: pi per site)
    bins = np.logspace(np.log10(lim_lo), np.log10(lim_hi), 32)
    ax_top.hist(df["theta_W_per_site"], bins=bins,
                color=WONG["grey"], alpha=0.55, edgecolor="white",
                linewidth=0.3)
    ax_top.set_yticks([])
    for side in ("top", "right", "left"):
        ax_top.spines[side].set_visible(False)
    for tl in ax_top.get_xticklabels():
        tl.set_visible(False)
    ax_top.tick_params(axis="x", which="both", length=0)

    ax_right.hist(df["pi_per_site"], bins=bins, orientation="horizontal",
                  color=WONG["grey"], alpha=0.55, edgecolor="white",
                  linewidth=0.3)
    ax_right.set_xticks([])
    for side in ("top", "right", "bottom"):
        ax_right.spines[side].set_visible(False)
    for tl in ax_right.get_yticklabels():
        tl.set_visible(False)
    ax_right.tick_params(axis="y", which="both", length=0)

    # legend block to the right of the scatter
    ax_leg = fig.add_subplot(gs[0:2, 2])
    ax_leg.axis("off")
    handles, labels = ax.get_legend_handles_labels()
    # put the diagonal-line handle last
    diag_idx = next((i for i, l in enumerate(labels)
                     if "neutral" in l), None)
    if diag_idx is not None:
        handles = handles[:diag_idx] + handles[diag_idx + 1:] + [handles[diag_idx]]
        labels = labels[:diag_idx] + labels[diag_idx + 1:] + [labels[diag_idx]]
    ax_leg.legend(handles, labels, loc="upper left", frameon=False,
                  fontsize=9, title="Chromosome", title_fontsize=10,
                  ncol=1)

    # bottom panel: Tajima D distribution by chromosome
    ax_bd = fig.add_subplot(gs[3, :])
    data_d = [df.loc[df["chr"] == c, "tajima_D"].values for c in chroms]
    bp = ax_bd.boxplot(
        data_d, labels=chroms, patch_artist=True, widths=0.62,
        medianprops={"color": "white", "linewidth": 1.4},
        whiskerprops={"color": "#333", "linewidth": 0.8},
        capprops={"color": "#333", "linewidth": 0.8},
        flierprops={"marker": "o", "markersize": 3,
                    "markerfacecolor": "#888", "markeredgecolor": "none",
                    "alpha": 0.55},
    )
    for patch, c in zip(bp["boxes"], chroms):
        patch.set_facecolor(chr_colour[c])
        patch.set_alpha(0.75)
        patch.set_edgecolor("#333")
        patch.set_linewidth(0.6)
    ax_bd.axhline(0, color="#444", linestyle="--", linewidth=0.8,
                  alpha=0.7, zorder=1, label=r"$D = 0$ (neutral)")
    ax_bd.axhline(2, color=WONG["vermillion"], linestyle=":",
                  linewidth=1.0, alpha=0.85, zorder=1,
                  label=r"$D = 2$ (strong bottleneck threshold)")
    ax_bd.set_ylabel(r"Tajima's $D$", fontsize=10.5)
    ax_bd.set_xlabel("S. stenocarpa pseudo-chromosome", fontsize=10.5)
    ax_bd.legend(loc="upper left", fontsize=9, framealpha=0.92,
                 edgecolor="none")
    publishable_axes(ax_bd, grid="y", grid_alpha=0.20)

    panel_label(ax, "a", x=-0.10, y=1.28)
    panel_label(ax_bd, "b", x=-0.045, y=1.06)

    fig.suptitle(
        r"Bottleneck diagnostic: $\pi$ vs $\theta_W$ per window across "
        r"S. stenocarpa pseudo-chromosomes "
        f"(n = {len(df)} windows; per-chromosome Tajima's $D$ below)",
        y=0.985, fontsize=12,
    )
    fig.savefig(OUT / "figures/fig_pi_vs_thetaW.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "figures/fig_pi_vs_thetaW.pdf",
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\nWrote {OUT / 'figures/fig_pi_vs_thetaW.png'}")


if __name__ == "__main__":
    main()
