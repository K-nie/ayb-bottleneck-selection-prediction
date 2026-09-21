"""
78 -- SHAP top-5 cross-trait marker overlap.

For each of 13 traits, take the top-5 markers ranked by mean |SHAP| and ask
whether any markers recur across traits (SNP-level pleiotropy under the model
that won CV). At top-5 the chance pairwise overlap is 5 * 5 / 2862 ~ 0.009
markers, so any observed overlap >= 1 would be enriched.

Result: there is none. All 65 top-5 markers are trait-specific -- every
marker appears in exactly one trait's top-5, and all 78 trait pairs share
zero markers. An UpSet plot of this is degenerate (every bar = 1), so we
show the honest picture directly:

  (A) a trait x marker presence "staircase" -- markers grouped by trait,
      coloured by SHAP rank. Perfect block structure = zero column sharing.
  (B) the 13 x 13 pairwise shared-marker heatmap (diagonal = 5, all
      off-diagonal = 0), which quantifies the absence of overlap.
"""
from __future__ import annotations
import sys
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, save_figure, panel_label
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TOP20 = PROJ / "results/41_ml_shap_interpretability/tables/shap_top20_per_trait.csv"
OUT = PROJ / "results/78_shap_upset"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)


TOP_K = 5


def main() -> None:
    df_all = pd.read_csv(TOP20)
    df = df_all[df_all["rank"] <= TOP_K].copy()
    print(f"Loaded {len(df)} top-{TOP_K} SHAP rows "
          f"({df['trait'].nunique()} traits).")

    # Per-marker trait recurrence
    rec = (
        df.groupby("rs")
        .agg(n_traits=("trait", "nunique"),
             traits=("trait", lambda s: ",".join(sorted(set(s)))),
             best_models=("best_model", lambda s: ",".join(sorted(set(s)))))
        .reset_index()
        .sort_values("n_traits", ascending=False)
    )
    rec.to_csv(OUT / "tables/shap_marker_recurrence.csv", index=False)
    print(f"\nMarkers appearing in 2+ trait top-20 sets: "
          f"{(rec['n_traits'] >= 2).sum()}")
    print(rec[rec["n_traits"] >= 2].head(15).to_string(index=False))

    # Build trait -> set(marker) mapping
    traits = sorted(df["trait"].unique())
    trait_sets: dict[str, set[str]] = {
        t: set(df.loc[df["trait"] == t, "rs"]) for t in traits
    }

    # Pairwise overlap matrix
    overlap = pd.DataFrame(0, index=traits, columns=traits, dtype=int)
    for t in traits:
        overlap.loc[t, t] = len(trait_sets[t])
    for t1, t2 in combinations(traits, 2):
        n = len(trait_sets[t1] & trait_sets[t2])
        overlap.loc[t1, t2] = n
        overlap.loc[t2, t1] = n
    overlap.to_csv(OUT / "tables/pairwise_overlap_matrix.csv")
    print(f"\nMean pairwise overlap (off-diagonal, expected "
          f"{TOP_K * TOP_K / 2862:.3f}): "
          f"{overlap.values[np.triu_indices_from(overlap.values, k=1)].mean():.3f}")

    # ---- honest overlap figure ----
    off = overlap.values[np.triu_indices_from(overlap.values, k=1)]
    n_pairs = len(off)
    n_shared_pairs = int((off > 0).sum())
    union_markers = sorted(set().union(*trait_sets.values()))
    n_union = len(union_markers)
    max_recur = int(rec["n_traits"].max()) if len(rec) else 1
    print(f"\nUnion of top-{TOP_K} markers = {n_union}; "
          f"max cross-trait recurrence = {max_recur}; "
          f"pairs sharing >=1 marker = {n_shared_pairs}/{n_pairs}")

    # Build the trait x marker staircase: markers grouped by trait, in rank
    # order. Cell value = SHAP rank (1..TOP_K) rendered so rank 1 is darkest;
    # NaN where the marker is not in that trait's top-K.
    col_markers: list[str] = []
    col_owner: list[str] = []
    for t in traits:
        sub = df.loc[df["trait"] == t].sort_values("rank")
        for m in sub["rs"]:
            col_markers.append(m)
            col_owner.append(t)
    P = np.full((len(traits), len(col_markers)), np.nan)
    rank_by = {(r.trait, r.rs): int(r.rank) for r in df.itertuples()}
    for i, t in enumerate(traits):
        for j, m in enumerate(col_markers):
            if (t, m) in rank_by:
                P[i, j] = TOP_K + 1 - rank_by[(t, m)]  # rank1 -> TOP_K (darkest)

    fig = plt.figure(figsize=(14, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.4, 1.6], wspace=0.22)

    # Panel A: staircase presence heatmap
    axA = fig.add_subplot(gs[0, 0])
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("white")
    axA.imshow(P, aspect="auto", cmap=cmap, vmin=0.5, vmax=TOP_K + 0.5,
               interpolation="nearest")
    # block separators between traits along the marker axis
    cursor = 0
    for t in traits:
        w = sum(1 for o in col_owner if o == t)
        cursor += w
        if cursor < len(col_markers):
            axA.axvline(cursor - 0.5, color=WONG["grey"], lw=0.6, alpha=0.6)
    axA.set_yticks(np.arange(len(traits)))
    axA.set_yticklabels([t.replace("_", " ") for t in traits], fontsize=9)
    axA.set_xticks([])
    axA.set_xlabel(f"{n_union} top-{TOP_K} SHAP markers, grouped by trait "
                   f"(each trait's private block of {TOP_K})")
    axA.set_title(f"Each trait's top-{TOP_K} SHAP markers are its own",
                  loc="left")
    panel_label(axA, "a")
    for s in ("top", "right"):
        axA.spines[s].set_visible(False)

    # Panel B: pairwise shared-marker heatmap (integer counts)
    axB = fig.add_subplot(gs[0, 1])
    O = overlap.loc[traits, traits].values.astype(int)
    Ooff = O.copy().astype(float)
    np.fill_diagonal(Ooff, np.nan)  # colour only off-diagonal; diag is set size
    imB = axB.imshow(Ooff, cmap="Reds", vmin=0, vmax=max(1, off.max()),
                     interpolation="nearest")
    axB.set_xticks(np.arange(len(traits)))
    axB.set_xticklabels([t.replace("_", " ") for t in traits],
                        rotation=90, fontsize=7)
    axB.set_yticks(np.arange(len(traits)))
    axB.set_yticklabels([t.replace("_", " ") for t in traits], fontsize=7)
    for i in range(len(traits)):
        for j in range(len(traits)):
            v = O[i, j]
            axB.text(j, i, str(v), ha="center", va="center", fontsize=6,
                     color=("#999999" if i == j else
                            ("black" if v == 0 else "white")))
    axB.set_title("Shared markers per trait pair\n(diagonal = set size)",
                  loc="left", fontsize=10)
    panel_label(axB, "b")

    fig.suptitle(f"SHAP top-{TOP_K} markers are trait-specific: "
                 f"{n_shared_pairs} of {n_pairs} trait pairs share any marker "
                 f"(max recurrence {max_recur})", y=1.02)
    written = save_figure(fig, OUT / "figures/fig_shap_upset")
    plt.close(fig)
    print("Wrote", *written)


if __name__ == "__main__":
    main()
