"""
80 -- SHAP top-5 cross-CV-repeat stability via bootstrap proxy.

Per-CV-fold SHAP values were not stored by script 41 (script 41 trains a
single best-CV-winner per trait on the full working set, then computes SHAP
once). The reviewer-grade stability check therefore uses a bootstrap proxy:
for each (trait, best-tree-model) pair we resample the n_work training rows
with replacement, refit the model, recompute SHAP, and record the top-5
markers. Across B = 30 bootstraps we report

  (a) per-marker selection frequency (how many of B reps the marker lands
      in the top-5);
  (b) Spearman rank correlation of full per-marker mean|SHAP| vectors
      between bootstrap pairs;
  (c) Jaccard index of top-5 sets between bootstrap pairs.

We restrict to tree-based winners (Random Forest, LightGBM) because the
Paper B headline-9 result establishes that SHAP recovers GWAS rank only
for tree winners.

Honest framing for the README: bootstrap and k-fold are not equivalent
stability measures. Bootstrap over-samples some rows and drops others
(out-of-bag fraction ~ 36 %); k-fold leaves out one fold at a time. The
bootstrap proxy is a defensible upper bound on rank-stability volatility
because each replicate sees a more-perturbed training set than a single
k-fold leave-out would. The numbers reported here therefore err toward
the conservative side.
"""
from __future__ import annotations
import os
import sys
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
CV40 = ROOT / "results/40_ml_genomic_prediction/tables/cv_predictive_ability.csv"
OUT = ROOT / "results/80_shap_bootstrap_stability"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
HOLDOUT_FRAC = 0.15
HOLDOUT_SEED = 42
TOP_K = 5
N_BOOT = 30

sys.path.insert(0, str(ROOT / "scripts"))
from _pheno import load_phenotypes  # type: ignore
from _figstyle import apply, WONG, publishable_axes, save_figure, panel_label
apply()


def load_data():
    hm = pd.read_csv(HAPMAP, low_memory=False)
    META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sample_cols = [c for c in hm.columns if c not in META]
    ref_alt = hm["alleles"].str.split("/", expand=True)
    ref_alt.columns = ["ref", "alt"]
    dosage = np.full((len(hm), len(sample_cols)), np.nan)
    calls = hm[sample_cols].astype(str)
    for i in range(len(hm)):
        r_, a_ = ref_alt.iloc[i]
        arr = calls.iloc[i].values
        dosage[i, arr == r_ + r_] = 0.0
        dosage[i, (arr == r_ + a_) | (arr == a_ + r_)] = 1.0
        dosage[i, arr == a_ + a_] = 2.0
    dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)
    cr_m = dosage.notna().sum(axis=1) / dosage.shape[1]
    maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2,
                     1 - dosage.mean(axis=1, skipna=True) / 2)
    cr_s = dosage.notna().sum(axis=0) / dosage.shape[0]
    keep_m = ((cr_m >= MARK_CR) & (maf >= MIN_MAF)).values
    keep_s = (cr_s >= SAMP_CR).values
    dosage = dosage.loc[keep_m, keep_s]
    markers = list(dosage.index)
    samples = list(dosage.columns)
    X_full = dosage.T.values.astype(float)

    rng = np.random.default_rng(HOLDOUT_SEED)
    idx_all = np.arange(len(samples))
    rng.shuffle(idx_all)
    n_hold = int(round(len(samples) * HOLDOUT_FRAC))
    work_idx = idx_all[n_hold:]
    work_samples = [samples[i] for i in work_idx]
    X_work = X_full[work_idx]
    mu = np.nanmean(X_work, axis=0)
    X_work_imp = np.where(np.isnan(X_work), mu, X_work)
    scaler = StandardScaler(with_mean=True, with_std=False).fit(X_work_imp)
    X_work_scaled = scaler.transform(X_work_imp)
    return X_work_scaled, work_samples, markers


def fit_tree(name: str, X: np.ndarray, y: np.ndarray):
    if name == "RandomForest":
        return RandomForestRegressor(
            n_estimators=500, max_features=0.1, min_samples_leaf=2,
            random_state=0, n_jobs=-1,
        ).fit(X, y)
    if name == "LightGBM":
        import lightgbm as lgb
        return lgb.LGBMRegressor(
            n_estimators=300, num_leaves=31, learning_rate=0.05,
            feature_fraction=0.3, n_jobs=-1, random_state=0, verbose=-1,
        ).fit(X, y)
    if name == "XGBoost":
        import xgboost as xgb
        return xgb.XGBRegressor(
            n_estimators=300, max_depth=5, learning_rate=0.05, subsample=0.7,
            n_jobs=-1, random_state=0, verbosity=0, tree_method="hist",
        ).fit(X, y)
    raise ValueError(name)


def tree_shap_mean_abs(model, X: np.ndarray) -> np.ndarray:
    import shap
    expl = shap.TreeExplainer(model)
    sv = np.array(expl.shap_values(X))
    if sv.ndim == 3:
        sv = sv[0]
    return np.mean(np.abs(sv), axis=0)


def _short_rs(rs: str) -> str:
    """DArTseq clone id -> compact numeric tag (100035301|F|0-17... -> 100035301)."""
    return str(rs).split("|", 1)[0]


def make_figure(tree_winners: pd.DataFrame, pair_df: pd.DataFrame) -> None:
    """Build the two-panel stability figure from the saved tables only."""
    traits = tree_winners["trait"].tolist()
    models = tree_winners["model"].tolist()
    trait_col = {t: c for t, c in zip(
        traits, [WONG["blue"], WONG["vermillion"], WONG["green"],
                 WONG["orange"], WONG["purple"]])}

    fig = plt.figure(figsize=(13.0, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.30)

    # Panel A: pairwise Spearman + Jaccard boxplots per trait
    axA = fig.add_subplot(gs[0, 0])
    width = 0.35
    x_pos = np.arange(len(traits))
    sp_data = [pair_df.loc[pair_df["trait"] == t, "spearman_full"].values
               for t in traits]
    ja_data = [pair_df.loc[pair_df["trait"] == t, "jaccard_top_k"].values
               for t in traits]
    bp1 = axA.boxplot(sp_data, positions=x_pos - width / 2, widths=width,
                      patch_artist=True, showfliers=False,
                      medianprops={"color": "black", "lw": 1.0})
    bp2 = axA.boxplot(ja_data, positions=x_pos + width / 2, widths=width,
                      patch_artist=True, showfliers=False,
                      medianprops={"color": "black", "lw": 1.0})
    for b in bp1["boxes"]:
        b.set(facecolor=WONG["blue"], alpha=0.65, edgecolor=WONG["blue"])
    for b in bp2["boxes"]:
        b.set(facecolor=WONG["green"], alpha=0.65, edgecolor=WONG["green"])
    axA.set_xticks(x_pos)
    axA.set_xticklabels([f"{t.replace('_', ' ')}\n({m})"
                         for t, m in zip(traits, models)], fontsize=8)
    axA.set_ylabel("Stability metric")
    axA.set_title(f"Bootstrap pairwise stability (B = {N_BOOT})", loc="left")
    panel_label(axA, "a")
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=WONG["blue"], alpha=0.65,
                      label=r"Spearman $\rho$ (full ranks)"),
        plt.Rectangle((0, 0), 1, 1, facecolor=WONG["green"], alpha=0.65,
                      label=f"Jaccard (top-{TOP_K} sets)"),
    ]
    axA.legend(handles=handles, loc="upper left")
    axA.axhline(0, color="black", lw=0.6)
    publishable_axes(axA, grid="y")

    # Panel B: per-trait top-10 marker selection frequency, one row per trait,
    # each bar labelled with its own marker id (markers differ across traits).
    gsB = gs[0, 1].subgridspec(len(traits), 1, hspace=0.65)
    for i, (t, m) in enumerate(zip(traits, models)):
        axb = fig.add_subplot(gsB[i])
        sel = pd.read_csv(OUT / f"tables/{t}_{m}_selection_freq.csv").head(10)
        sel = sel.iloc[::-1]  # largest frequency on top
        y = np.arange(len(sel))
        axb.barh(y, sel["selection_frac"], color=trait_col.get(t, WONG["grey"]),
                 alpha=0.85, edgecolor="white", lw=0.4)
        axb.set_yticks(y)
        axb.set_yticklabels([_short_rs(r) for r in sel["rs"]], fontsize=6)
        axb.set_xlim(0, 1.0)
        axb.tick_params(axis="x", labelsize=7)
        axb.set_title(f"{t.replace('_', ' ')} ({m})", fontsize=8.5, loc="left")
        publishable_axes(axb, grid="x")
        if i == len(traits) - 1:
            axb.set_xlabel(f"Top-{TOP_K} selection frequency across "
                           f"{N_BOOT} bootstraps")
    _bx0 = gs[0, 1].get_position(fig).x0
    fig.text(_bx0 - 0.03, 0.965, "b", fontsize=12, fontweight="bold",
             ha="right", va="bottom")
    fig.text(_bx0 - 0.01, 0.965, "Per-trait top-10 marker selection frequency",
             fontsize=10, ha="left", va="bottom")

    written = save_figure(fig, OUT / "figures/fig_shap_bootstrap_stability")
    plt.close(fig)
    print("Wrote", *written)


def main() -> None:
    cv = pd.read_csv(CV40)
    tree_names = {"RandomForest", "LightGBM", "XGBoost"}
    winners = (
        cv[cv["model"] != "Dummy"]
        .sort_values(["trait", "r_mean"], ascending=[True, False])
        .drop_duplicates("trait", keep="first")
    )
    tree_winners = winners[winners["model"].isin(tree_names)].copy()
    print("\nTree-based CV winners:")
    print(tree_winners[["trait", "model", "r_mean", "n_work"]])

    # Fast path: rebuild the figure from saved tables without the (expensive)
    # 30-bootstrap refit + SHAP recompute. Set FORCE_RECOMPUTE=1 to rerun.
    pair_path = OUT / "tables/pairwise_stability.csv"
    sel_paths = [OUT / f"tables/{r.trait}_{r.model}_selection_freq.csv"
                 for r in tree_winners.itertuples()]
    tables_ready = pair_path.exists() and all(p.exists() for p in sel_paths)
    if tables_ready and os.environ.get("FORCE_RECOMPUTE") != "1":
        print("Saved tables present -> rebuilding figure only "
              "(set FORCE_RECOMPUTE=1 to rerun the bootstrap).")
        make_figure(tree_winners, pd.read_csv(pair_path))
        return

    pheno = load_phenotypes()
    pheno.index = pheno.index.astype(str).str.strip()

    X_work, work_samples, markers = load_data()
    print(f"X_work shape = {X_work.shape}; {len(markers)} markers")

    rng = np.random.default_rng(20260613)

    all_top_k_rows = []
    pairwise_rows = []

    for _, w in tree_winners.iterrows():
        trait = w["trait"]; mname = w["model"]
        y_full = pheno[trait].reindex(work_samples).values.astype(float)
        obs = ~np.isnan(y_full)
        X = X_work[obs]; y = y_full[obs]
        n = X.shape[0]
        print(f"\n[trait={trait} | model={mname}] n_obs = {n}")

        rank_matrix = np.zeros((N_BOOT, X.shape[1]))  # rank per marker per boot
        top_k_sets = []
        for b in range(N_BOOT):
            idx = rng.choice(n, size=n, replace=True)
            Xb = X[idx]; yb = y[idx]
            model = fit_tree(mname, Xb, yb)
            mab = tree_shap_mean_abs(model, Xb)
            # rank (1 = largest |SHAP|)
            ranks = np.empty_like(mab)
            order = np.argsort(-mab)
            ranks[order] = np.arange(1, len(mab) + 1)
            rank_matrix[b] = ranks
            top_k = set(order[:TOP_K])
            top_k_sets.append(top_k)
            for rk, mi in enumerate(order[:TOP_K], start=1):
                all_top_k_rows.append({
                    "trait": trait, "model": mname, "boot": int(b),
                    "rank": int(rk), "rs": markers[int(mi)],
                    "mean_abs_shap": float(mab[int(mi)]),
                })
            if (b + 1) % 5 == 0:
                print(f"  boot {b + 1}/{N_BOOT}")

        # Per-marker selection frequency across bootstraps
        sel_count: dict[int, int] = {}
        for s in top_k_sets:
            for mi in s:
                sel_count[mi] = sel_count.get(mi, 0) + 1
        sel_df = pd.DataFrame([
            {"trait": trait, "model": mname, "rs": markers[mi],
             "selection_freq": int(c), "selection_frac": c / N_BOOT}
            for mi, c in sel_count.items()
        ]).sort_values("selection_freq", ascending=False)
        sel_df.to_csv(
            OUT / f"tables/{trait}_{mname}_selection_freq.csv", index=False
        )

        # Pairwise Spearman correlation of full-marker ranks + Jaccard at top-K
        from scipy.stats import spearmanr
        for i in range(N_BOOT):
            for j in range(i + 1, N_BOOT):
                rho, _ = spearmanr(rank_matrix[i], rank_matrix[j])
                inter = len(top_k_sets[i] & top_k_sets[j])
                union = len(top_k_sets[i] | top_k_sets[j])
                jac = inter / union if union else 0.0
                pairwise_rows.append({
                    "trait": trait, "model": mname,
                    "boot_i": int(i), "boot_j": int(j),
                    "spearman_full": float(rho),
                    "jaccard_top_k": float(jac),
                })

    pd.DataFrame(all_top_k_rows).to_csv(
        OUT / "tables/per_boot_topK_markers.csv", index=False
    )
    pair_df = pd.DataFrame(pairwise_rows)
    pair_df.to_csv(OUT / "tables/pairwise_stability.csv", index=False)

    summary = (
        pair_df.groupby(["trait", "model"])
        [["spearman_full", "jaccard_top_k"]]
        .agg(["mean", "median", "std"])
        .round(3)
    )
    summary.to_csv(OUT / "tables/stability_summary.csv")
    print("\nStability summary:")
    print(summary)

    make_figure(tree_winners, pair_df)


if __name__ == "__main__":
    main()
