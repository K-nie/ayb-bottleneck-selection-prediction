"""
75 -- Per-trait CV-r vs narrow-sense h2 with the Daetwyler ceiling.

Informational-ceiling argument: under the Daetwyler / Goddard genomic-prediction
theory, the maximum attainable prediction correlation is

    r_max = sqrt( h^2 * N_p / (N_p + M_e * h^2) )

where N_p is the training-set size, h^2 is the narrow-sense heritability, and
M_e is the effective number of independent chromosome segments. At our n_work
= 30 -- 75 panel sizes and M_e ~ 1500 -- 3000 (informed by AYB LD decay and
genome size), the ceiling is roughly r_max ~ 0.3 -- 0.6 even at h^2 = 0.7.

Plot: for each trait, observed CV-r of the best estimator vs the REML profile
h^2, with two reference Daetwyler ceilings (M_e = 1500 and M_e = 3000) drawn
as continuous curves at the trait's median n_work. Errorbars on x = h^2 95 %
profile CI; errorbars on y = CV r 95 % CI. Trait labels placed via adjustText
so no overlap; legend pulled out of the data area.

Honest framing: small-n panels with sparse DArTseq markers cannot exceed the
Daetwyler ceiling no matter the estimator. The negative-CV-r traits (Tannin,
Flavonoid, etc.) sit at h^2 = 0.001 (lower boundary) and r_max ~ 0 -- their
"failure" is informational, not algorithmic.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, adjust_labels, publishable_axes
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
H2 = PROJ / "results/24_bootstrap_h2/tables/h2_profile_summary.csv"
CV = PROJ / "results/40_ml_genomic_prediction/tables/cv_predictive_ability.csv"
OUT = PROJ / "results/75_cv_r_vs_h2_daetwyler"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)


def daetwyler_r(h2: float, n_p: float, m_e: float) -> float:
    if h2 <= 0 or n_p <= 0:
        return 0.0
    return float(np.sqrt(h2 * n_p / (n_p + m_e * h2)))


def main() -> None:
    h2 = pd.read_csv(H2)
    cv = pd.read_csv(CV)
    cv_best = (
        cv[cv["model"] != "Dummy"]
        .sort_values(["trait", "r_mean"], ascending=[True, False])
        .drop_duplicates("trait", keep="first")
        .reset_index(drop=True)
    )
    df = cv_best.merge(h2[["trait", "h2_mle", "h2_ci95_low", "h2_ci95_high"]],
                       on="trait", how="left")
    df.to_csv(OUT / "tables/per_trait_cv_vs_h2.csv", index=False)
    print(df[["trait", "model", "n_work", "r_mean", "r_ci_lo", "r_ci_hi",
              "h2_mle", "h2_ci95_low", "h2_ci95_high"]])

    n_median = float(df["n_work"].median())
    print(f"Median trait n_work = {n_median:.0f}")

    # Estimator colour map (Wong-friendly; muted but distinguishable)
    colours = {
        "RR-BLUP":      WONG["blue"],
        "KernelRidge":  WONG["skyblue"],
        "MLP":          WONG["yellow"],
        "ElasticNet":   WONG["purple"],
        "RandomForest": WONG["green"],
        "XGBoost":      WONG["vermillion"],
        "LightGBM":     WONG["orange"],
    }

    # Wider+taller figure so the legend can sit BELOW the data axes
    # in two horizontal rows (ceilings on top, estimators below)
    fig = plt.figure(figsize=(10.5, 8.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[6.0, 1.5], hspace=0.18)
    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[1, 0])
    ax_leg.axis("off")

    # Daetwyler ceiling curves at two M_e values, plus sqrt(h^2) bound
    h2_grid = np.linspace(0.001, 0.999, 200)
    ceiling_handles = []
    for m_e, ls, col, label in [
        (1500, "--", "#117733",
         rf"Daetwyler $r_{{\max}}$ ($M_e$ = 1500, $N$ = {n_median:.0f})"),
        (3000, ":", "#332288",
         rf"Daetwyler $r_{{\max}}$ ($M_e$ = 3000, $N$ = {n_median:.0f})"),
    ]:
        rs = [daetwyler_r(h, n_median, m_e) for h in h2_grid]
        line, = ax.plot(h2_grid, rs, ls, color=col, lw=1.6,
                        label=label, zorder=2)
        ceiling_handles.append(line)
    sqrt_line, = ax.plot(h2_grid, np.sqrt(h2_grid), "-",
                          color="#999", lw=1.0,
                          label=r"$\sqrt{h^2}$ (theoretical upper bound)",
                          zorder=1)
    ceiling_handles.append(sqrt_line)

    # CV winners scatter with x + y error bars and per-trait labels
    trait_label_texts = []
    for _, r in df.iterrows():
        col = colours.get(r["model"], "#666")
        ax.errorbar(
            r["h2_mle"], r["r_mean"],
            yerr=[[r["r_mean"] - r["r_ci_lo"]],
                  [r["r_ci_hi"] - r["r_mean"]]],
            xerr=[[max(r["h2_mle"] - r["h2_ci95_low"], 0)],
                  [max(r["h2_ci95_high"] - r["h2_mle"], 0)]],
            fmt="o", color=col, ecolor=col,
            elinewidth=0.9, capsize=2.5, ms=8,
            alpha=0.92, markeredgecolor="black", markeredgewidth=0.6,
            zorder=4,
        )
        trait_label_texts.append(
            ax.text(r["h2_mle"], r["r_mean"],
                     r["trait"].replace("_", " "),
                     fontsize=9, color="#222", alpha=0.96,
                     zorder=10))
    adjust_labels(trait_label_texts, ax=ax,
                   expand=(1.35, 1.55), force_text=(1.1, 1.4),
                   only_move={"text": "xy", "static": "xy"})

    ax.axhline(0, color="black", lw=0.6, alpha=0.5)
    ax.set_xlim(-0.03, 1.03)
    y_lo = float(min(df["r_ci_lo"].min(), -0.4) - 0.05)
    y_hi = float(max(df["r_ci_hi"].max(), 0.65) + 0.05)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel(r"Narrow-sense heritability $h^2$ (REML profile MLE)",
                  fontsize=11)
    ax.set_ylabel("Observed CV r (best estimator per trait)",
                  fontsize=11)
    ax.set_title(r"Per-trait CV r vs $h^2$ against the Daetwyler informational ceiling",
                  fontsize=11.5)
    publishable_axes(ax, grid="both", grid_alpha=0.20)

    # Bottom legend: ceilings in row 1, CV-winning estimators in row 2
    est_handles = [
        plt.Line2D([], [], marker="o", color=col, ls="", ms=8,
                    markeredgecolor="black", markeredgewidth=0.5,
                    label=name)
        for name, col in colours.items() if name in df["model"].values
    ]
    leg1 = ax_leg.legend(
        ceiling_handles, [h.get_label() for h in ceiling_handles],
        loc="upper center", frameon=False, fontsize=9.5,
        title="Reference curves", title_fontsize=10,
        bbox_to_anchor=(0.5, 1.0), ncol=3,
    )
    ax_leg.add_artist(leg1)
    ax_leg.legend(
        est_handles, [h.get_label() for h in est_handles],
        loc="upper center", frameon=False, fontsize=9.5,
        title="CV-winning estimator (per trait)", title_fontsize=10,
        bbox_to_anchor=(0.5, 0.45),
        ncol=min(len(est_handles), 7),
    )

    fig.tight_layout()
    fig.savefig(OUT / "figures/fig_cv_r_vs_h2.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "figures/fig_cv_r_vs_h2.pdf",
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Wrote", OUT / "figures/fig_cv_r_vs_h2.png")


if __name__ == "__main__":
    main()
