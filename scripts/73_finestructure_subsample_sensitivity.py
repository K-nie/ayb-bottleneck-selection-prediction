"""
73 -- K = 24 fineSTRUCTURE subsampling sensitivity.

Reviewer concern: at n = 95 with K = 24, the panel risks over-fitting because
K = n / 4 in a single-population panel.

Pragmatic test: rather than re-running ChromoPainter + fineSTRUCTURE MCMC on
n = 80 / 60 / 40 subsamples (days of wall time), we resample rows + columns of
the existing coancestry chunk-count matrix, then cluster at multiple K values
using hierarchical (Ward) clustering and silhouette scoring on (1 - coancestry)
as the distance metric. This isolates the "is K = 24 robust to sample dropout?"
question to the cluster-stability layer, decoupled from the haplotype-painting
layer.

Two outputs:
  (A) max-silhouette K as a function of subsample size n in {40, 60, 80, 95}.
  (B) Adjusted Rand Index (ARI) of K = 24 cluster membership across pairs of
      bootstrap resamples at each n -- direct stability measure for K = 24.

Honest framing in the README: the silhouette is a coarse approximation of
fineSTRUCTURE's marginal likelihood. The interpretation is "what K does the
coancestry signal natively support at each n", not "what K would fineSTRUCTURE
return at this n".
"""
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score, adjusted_rand_score

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
COANC = PROJ / "results/37_finestructure/tables/coancestry.csv"
CLUSTS = PROJ / "results/37_finestructure/tables/cluster_assignments.csv"
OUT = PROJ / "results/73_finestructure_subsample"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(1492)
N_REPS = 50
SUB_SIZES = [40, 60, 80, 95]
K_GRID = [2, 3, 5, 8, 10, 15, 20, 24, 28]


def hierarchical_labels(D: np.ndarray, k: int) -> np.ndarray:
    """Ward hierarchical clustering on a precomputed distance matrix."""
    # AgglomerativeClustering with ward needs Euclidean features. Use average
    # linkage on the precomputed distance instead -- this is the standard
    # choice for fineSTRUCTURE-style coancestry clustering anyway.
    clf = AgglomerativeClustering(
        n_clusters=k, metric="precomputed", linkage="average"
    )
    return clf.fit_predict(D)


def coancestry_to_distance(M: np.ndarray) -> np.ndarray:
    """Symmetrise and convert coancestry (high = closer) to a distance."""
    S = 0.5 * (M + M.T)
    # max-normalise so distance is in [0, 1]; off-diagonal max excludes diag.
    off = S.copy()
    np.fill_diagonal(off, np.nan)
    s_max = float(np.nanmax(off))
    D = 1.0 - (S / s_max)
    np.fill_diagonal(D, 0.0)
    D = np.clip(D, 0.0, None)
    return D


def silhouette_per_k(D: np.ndarray, k_grid: list[int]) -> dict[int, float]:
    out: dict[int, float] = {}
    n = D.shape[0]
    for k in k_grid:
        if k >= n:
            continue
        labs = hierarchical_labels(D, k)
        # silhouette wants distance matrix; require at least 2 labels with size>=1.
        if len(set(labs)) < 2:
            continue
        try:
            sc = silhouette_score(D, labs, metric="precomputed")
            out[k] = float(sc)
        except Exception:
            continue
    return out


def main() -> None:
    coanc_df = pd.read_csv(COANC, index_col=0)
    samples = coanc_df.index.tolist()
    M_full = coanc_df.values
    assert M_full.shape[0] == M_full.shape[1] == len(samples)
    print(f"Loaded coancestry: {M_full.shape}, {len(samples)} samples.")

    rows = []
    pairwise_ari_rows = []

    for n_sub in SUB_SIZES:
        # ---- per-rep silhouette curves and best-K ----
        sub_labels_at_k24: list[tuple[list[str], np.ndarray]] = []
        for rep in range(N_REPS):
            if n_sub == 95:
                idx = np.arange(95)
            else:
                idx = RNG.choice(95, size=n_sub, replace=False)
            sub_samples = [samples[i] for i in idx]
            M_sub = M_full[np.ix_(idx, idx)]
            D_sub = coancestry_to_distance(M_sub)
            sil = silhouette_per_k(D_sub, K_GRID)
            if not sil:
                continue
            best_k = max(sil, key=lambda k: sil[k])
            for k, s in sil.items():
                rows.append({
                    "n_sub": int(n_sub),
                    "rep": int(rep),
                    "k": int(k),
                    "silhouette": float(s),
                    "is_best_k": bool(k == best_k),
                })
            # store K = 24 labels for ARI calculation
            if 24 in K_GRID and 24 < n_sub:
                labs24 = hierarchical_labels(D_sub, min(24, n_sub - 1))
                sub_labels_at_k24.append((sub_samples, labs24))

        # ---- pairwise ARI at K = 24 across reps that share samples ----
        for i in range(len(sub_labels_at_k24)):
            samples_i, labs_i = sub_labels_at_k24[i]
            for j in range(i + 1, len(sub_labels_at_k24)):
                samples_j, labs_j = sub_labels_at_k24[j]
                shared = set(samples_i) & set(samples_j)
                if len(shared) < 20:
                    continue
                shared_sorted = sorted(shared)
                ai = [labs_i[samples_i.index(s)] for s in shared_sorted]
                aj = [labs_j[samples_j.index(s)] for s in shared_sorted]
                ari = adjusted_rand_score(ai, aj)
                pairwise_ari_rows.append({
                    "n_sub": int(n_sub),
                    "rep_i": int(i),
                    "rep_j": int(j),
                    "n_shared": int(len(shared_sorted)),
                    "ari_k24": float(ari),
                })

    df_sil = pd.DataFrame(rows)
    df_ari = pd.DataFrame(pairwise_ari_rows)
    df_sil.to_csv(OUT / "tables/silhouette_per_n_per_rep.csv", index=False)
    df_ari.to_csv(OUT / "tables/pairwise_ari_k24.csv", index=False)

    # ---- best-K per (n_sub, rep) summary ----
    best_k_df = (
        df_sil[df_sil["is_best_k"]]
        .groupby(["n_sub", "rep"])["k"]
        .first()
        .reset_index()
    )
    best_k_summary = (
        best_k_df.groupby("n_sub")["k"]
        .agg(["mean", "median", "min", "max", lambda s: int(s.mode().iloc[0])])
        .rename(columns={"<lambda_0>": "mode"})
        .reset_index()
    )
    best_k_summary.to_csv(OUT / "tables/best_k_summary.csv", index=False)
    print("\nBest-K (silhouette-maximising) per subsample size:")
    print(best_k_summary)

    # ---- ARI summary ----
    ari_summary = (
        df_ari.groupby("n_sub")["ari_k24"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
    )
    ari_summary.to_csv(OUT / "tables/ari_k24_summary.csv", index=False)
    print("\nPairwise ARI at K = 24 across resampled replicates:")
    print(ari_summary)

    # ---- FIGURE A: silhouette curves by n_sub ----
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    ax = axes[0]
    palette = {40: "#332288", 60: "#117733", 80: "#DDCC77", 95: "#CC6677"}
    for n_sub in SUB_SIZES:
        sub = df_sil[df_sil["n_sub"] == n_sub]
        means = sub.groupby("k")["silhouette"].mean()
        sds = sub.groupby("k")["silhouette"].std()
        ax.plot(means.index, means.values, "-o", color=palette[n_sub],
                label=f"n = {n_sub}", lw=1.8, ms=5)
        ax.fill_between(means.index, (means - sds).values, (means + sds).values,
                        color=palette[n_sub], alpha=0.10)
    ax.set_xlabel("K (number of clusters)")
    ax.set_ylabel("Mean silhouette (Ward, precomputed distance)")
    ax.set_title("Silhouette per K, across subsample sizes")
    ax.text(-0.08, 1.02, "a", transform=ax.transAxes, fontsize=12,
            fontweight="bold", ha="right", va="bottom", zorder=1000)
    ax.axvline(2, color="#888", ls=":", lw=1)
    ax.axvline(3, color="#888", ls=":", lw=1)
    ax.axvline(24, color="#000", ls="--", lw=1)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(True, alpha=0.25)

    # ---- FIGURE B: pairwise ARI at K = 24 ----
    ax = axes[1]
    box_data = [df_ari.loc[df_ari["n_sub"] == n, "ari_k24"].values for n in SUB_SIZES]
    bp = ax.boxplot(box_data, labels=[f"n={n}" for n in SUB_SIZES],
                    patch_artist=True, widths=0.55)
    for patch, n_sub in zip(bp["boxes"], SUB_SIZES):
        patch.set_facecolor(palette[n_sub]); patch.set_alpha(0.55)
    ax.axhline(0.0, color="#666", ls="--", lw=1)
    ax.set_ylabel("Pairwise ARI at K = 24")
    ax.set_title("Cluster-membership ARI across bootstrap pairs")
    ax.text(-0.08, 1.02, "b", transform=ax.transAxes, fontsize=12,
            fontweight="bold", ha="right", va="bottom", zorder=1000)
    ax.grid(True, alpha=0.25)

    fig.suptitle("K = 24 fineSTRUCTURE subsampling sensitivity",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "figures/fig_k24_subsample_sensitivity.png",
                dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "figures/fig_k24_subsample_sensitivity.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print("\nWrote figures to", OUT / "figures")


if __name__ == "__main__":
    main()
