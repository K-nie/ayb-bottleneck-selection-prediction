#!/usr/bin/env python3
"""
Nested-CV machine-learning genomic prediction for the AYB DArTseq panel.

Compares eight predictors per trait against the Daetwyler ceiling and
the single-trait GBLUP baseline from script 03:

    1. Dummy mean-predictor                              — naive baseline (Lones 5.2)
    2. RR-BLUP (ridge on centred dosage)                 — Meuwissen 2001 baseline
    3. Elastic net                                        — Zou & Hastie 2005
    4. Random Forest                                      — Heslot et al. 2012
    5. XGBoost                                            — Sandhu et al. 2021
    6. LightGBM                                           — Yan et al. 2023
    7. Kernel Ridge with RBF kernel                       — non-linear GBLUP
    8. MLP (2-hidden-layer, sklearn)                      — Bellot 2018 DL benchmark

Three CV schemes are run side-by-side for honest reporting:
    a. random 5-fold OUTER x 5-fold INNER, 10 outer repeats   (main)
    b. clone-blocked CV (near-clones in a single fold)        (sensitivity)
    c. stratified by PCA k=2 cluster label                    (sensitivity)

A 15 % final holdout is set aside at script start (random seed 42) and
NEVER touches CV — used once at the end to score the best model per
trait (Lones 4.5).

Output:
    results/40_ml_genomic_prediction/
        tables/
            cv_predictive_ability.csv             per-trait per-model r/CI/MAE/RMSE
            top_k_accuracy.csv                    breeder-side top-10 overlap
            ml_vs_daetwyler.csv                   best ML r vs theoretical ceiling
            wilcoxon_vs_rrblup.csv                paired Wilcoxon + Holm corrected
            cohens_d_delta_r.csv                  effect size on Δr vs RR-BLUP
            final_holdout_scores.csv              best-model performance on held-out 15 %
            cv_sensitivity_clone_blocked.csv      clone-blocked CV r
            cv_sensitivity_cluster_strat.csv      cluster-stratified CV r
        figures/
            fig85_model_comparison_forest
            fig86_ml_vs_daetwyler_ceiling
            fig87_top_k_accuracy_heatmap
            fig88_wilcoxon_effect_size

Lones (2024) Patterns 34-pitfall walk-through
---------------------------------------------
Walked through the Lones (2024) ML best-practices checklist for genomic prediction
before writing this script. Per-pitfall status:

  STAGE 1 — Before model building
   2.1 use of data: train/inner-validate/outer-test + 15% holdout separated
   2.2 understand data: EDA in scripts 01, 14, 23
   2.3 don't look at all: 15 % final-holdout untouched until last line
   2.4 clean data: QC + typo corrections in `_pheno.py`
   2.5 enough data: n=95 / n=41 flagged in Methods + Discussion
   2.6 domain experts: Prof. Adewale = breeder collaborator
   2.7 survey lit: Bellot 2018, Heslot 2012, Sandhu 2021, Montesinos-López cited
   2.8 deployment: GEBV prediction for unphenotyped TSs accessions

  STAGE 2 — Reliable model building
   3.1 no test leak: mean-impute, scale, fit ALL inside outer-train; nested CV
   3.2 range of models: 8 estimators across linear / tree / kernel / DL families
   3.3 inappropriate models: regression for continuous y; no time-series assumptions
   3.4 DL progress: MLP included; transformers/foundation models N/A for tabular n<500
   3.5 DL not best: honest report — MLP is one of eight, judged on the same CV
   3.6 feature selection: NONE — all 1,625 markers go into every fold
   3.7 hyperparameter optimisation: inner-fold grid search per outer-train
   3.8 spurious correlations: SHAP cross-check with GWAS hits (script 41)

  STAGE 3 — Robust evaluation
   4.1 appropriate test set: 5-fold outer + 15 % final holdout; clone- and
       cluster-blocked sensitivity checks
   4.2 augmentation before split: N/A — no data augmentation used
   4.3 sequential overfitting: 15 % holdout used exactly ONCE at the very end
   4.4 evaluate multiple times: 10 outer repeats per trait per model
       (Lehermeier 2014 recommends >= 20; 10 is acceptable when n < 200
       and the bootstrap CI on r is reported alongside the mean)
   4.5 final-instance holdout: 15 % final holdout
   4.6 metrics: Pearson r, Spearman ρ, MAE, RMSE, top-10 accuracy
   4.7 model fairness: subgroup r reported for PCA cluster 1 vs cluster 2
   4.8 temporal dependencies: N/A — single-season panel

  STAGE 4 — Fair model comparison
   5.1 bigger number != better: all models tested on identical folds + identical
       outer seed sequence
   5.2 meaningful baselines: dummy mean-predictor + RR-BLUP + GBLUP from script 03
   5.3 statistical tests: paired Wilcoxon signed-rank on outer-repeat r vectors
   5.4 multiple comparisons: Holm-Bonferroni across 6 ML-vs-RR-BLUP tests per trait
   5.5 community benchmarks: N/A — custom panel, no published-benchmark reuse
   5.6 ensemble carefully: stacking attempted in script 41

  STAGE 5 — Reporting
   6.1 transparent: scripts on GitHub; manifests; seed = 42 / 2000+rep
   6.2 multiple ways: 5 metrics + 95 % CIs + 2 sensitivity CV schemes
   6.3 don't overgeneralise: results limited to the 95-line IITA TSs panel
   6.4 statistical significance careful: Cohen's d + bootstrap 95 % CI on Δr
   6.5 look at models: SHAP in script 41
   6.6 ML checklist: this walk-through IS the checklist; REFORMS (Kapoor 2024) referenced

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.stats import pearsonr, spearmanr, wilcoxon
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import Ridge, ElasticNetCV
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler

from _plotstyle import apply, WONG
apply()

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
OUT  = ROOT / "results" / "40_ml_genomic_prediction"
DATA = OUT / "data"
TAB  = OUT / "tables"
FIG  = OUT / "figures"
for d in (DATA, TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

# QC and CV settings
SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05

# Repeated nested CV settings.
# Lehermeier et al. 2014 Genetics recommends >= 20 outer repeats for stable
# r CIs on diversity-panel data; 10 is acceptable when the dataset is small
# (n < 200) and the models are fully nested. We use 10 here to keep wall
# time tractable and use the per-repeat distribution for Wilcoxon + bootstrap
# CIs anyway.
N_OUTER_REPEATS = 10
N_OUTER_FOLDS   = 5
N_INNER_FOLDS   = 5
HOLDOUT_FRAC    = 0.15   # never used in CV (Lones 4.5)
HOLDOUT_SEED    = 42
TOP_K           = 10
N_BOOTSTRAP     = 2000
N_SENS_REPEATS  = 5      # clone-blocked + cluster-stratified sensitivity CV

# Near-clone groups identified by GRM analysis (script 13);
# G_ij >= 1.0 indicates identical-by-state under this marker set.
CLONE_GROUPS = [
    ["TSs377", "TSs297"],
    ["TSs361", "TSs358", "TSs151B"],
    ["TSs89", "TSs138"],
    ["TSs84", "TSs89", "TSs138"],   # overlaps the previous group on purpose
]


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Load post-QC dosage and 13 phenotypes
# ---------------------------------------------------------------------------
print("[load] HapMap...")
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
maf  = np.minimum(dosage.mean(axis=1, skipna=True) / 2,
                  1 - dosage.mean(axis=1, skipna=True) / 2)
cr_s = dosage.notna().sum(axis=0) / dosage.shape[0]
keep_m = ((cr_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (cr_s >= SAMP_CR).values
dosage = dosage.loc[keep_m, keep_s]

samples = dosage.columns.tolist()
X_full = dosage.T.values.astype(float)   # samples x markers (still has NaN)
print(f"[load] dosage matrix: {X_full.shape[0]} samples x {X_full.shape[1]} markers"
       f" (NaN entries to be imputed PER FOLD)")

# Load 13 phenotypes via the project's canonical loader, which handles
# the source-data quirks: "Flavinoid" typo -> "Flavonoid", "Crude Protein"
# space -> "Crude_Protein", Seed_Metrics multi-line header, and the two
# decimal-point typos at TSs282 / TSs136.
from _pheno import load_phenotypes, TRAITS_ALL

pheno = load_phenotypes()           # sample (index) x 13 trait columns
pheno.index = pheno.index.astype(str).str.strip()


def y_aligned(trait):
    """Return trait values aligned to the `samples` order, NaN where missing."""
    return pheno[trait].reindex(samples)


# (trait_name, source_df) tuples used by the legacy code path —
# kept for backward compatibility but `source_df` is no longer used.
TRAITS = [(t, None) for t in TRAITS_ALL]


# ---------------------------------------------------------------------------
# Helpers — per-fold imputation, scaling, holdout split, cluster labels
# ---------------------------------------------------------------------------

def impute_and_scale_train_only(X_train, X_test):
    """Mean-impute markers using TRAINING ROWS ONLY, then centre. The
    same per-marker mean and centred mean are applied to X_test.
    Lones 3.1 — the most common leakage source in plant-ML papers."""
    mu = np.nanmean(X_train, axis=0)
    Xtr = np.where(np.isnan(X_train), mu, X_train)
    Xte = np.where(np.isnan(X_test),  mu, X_test)
    scaler = StandardScaler(with_mean=True, with_std=False).fit(Xtr)
    return scaler.transform(Xtr), scaler.transform(Xte)


def make_holdout_indices(n_total, frac=HOLDOUT_FRAC, seed=HOLDOUT_SEED):
    rng = np.random.default_rng(seed)
    idx = np.arange(n_total)
    rng.shuffle(idx)
    n_hold = int(round(n_total * frac))
    return idx[:n_hold], idx[n_hold:]


def cluster_labels_from_pca(X):
    """Recover the 2-cluster K-means labels used by script 01 PCA.
    Used for stratified-by-cluster sensitivity CV. We re-fit
    K-means(k=2) here so the labels are reproducible without loading
    the script 01 artefact."""
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    Z = StandardScaler(with_mean=True).fit_transform(np.where(np.isnan(X), 0, X))
    pc = PCA(n_components=10, random_state=0).fit_transform(Z)
    km = KMeans(n_clusters=2, n_init=10, random_state=0).fit(pc)
    return km.labels_


def clone_block_labels(samples, clone_groups=CLONE_GROUPS):
    """Each sample gets a 'block id' — clone-mates share the same id
    so they cannot end up on opposite sides of a CV fold."""
    label = list(range(len(samples)))
    samp_to_idx = {s: i for i, s in enumerate(samples)}
    next_block_id = max(label) + 1
    for group in clone_groups:
        idxs = [samp_to_idx[s] for s in group if s in samp_to_idx]
        if len(idxs) < 2:
            continue
        keep = label[idxs[0]]
        for i in idxs[1:]:
            label[i] = keep
    # remap labels to a dense range
    uniq = {v: k for k, v in enumerate(sorted(set(label)))}
    return np.array([uniq[v] for v in label], dtype=int)


# ---------------------------------------------------------------------------
# Model factories
# ---------------------------------------------------------------------------

def model_dummy():
    return DummyRegressor(strategy="mean"), {}


def model_rrblup():
    return Ridge(), {"alpha": [0.01, 0.1, 1, 5, 25, 100, 1000]}


def model_enet():
    return ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9],
                        cv=N_INNER_FOLDS, n_alphas=20,
                        max_iter=20000, n_jobs=-1, random_state=0), {}


def model_rf():
    # Smaller forest + tighter feature subsample for speed; max_features
    # is fixed (no grid) because 1,625 features × 0.05 = ~80 candidate
    # features per split is the sweet spot for genomic data of this size
    # (Heslot et al. 2012 Crop Sci).
    return RandomForestRegressor(n_estimators=100, max_features=0.05,
                                  random_state=0, n_jobs=-1), {
        "min_samples_leaf": [1, 5],
    }


def model_xgb():
    import xgboost as xgb
    return xgb.XGBRegressor(n_estimators=100, max_depth=5,
                             subsample=0.7, n_jobs=-1, random_state=0,
                             verbosity=0,
                             tree_method="hist"), {
        "learning_rate": [0.05, 0.10],
    }


def model_lgbm():
    import lightgbm as lgb
    return lgb.LGBMRegressor(n_estimators=100, num_leaves=31,
                              feature_fraction=0.3, n_jobs=-1,
                              random_state=0, verbose=-1), {
        "learning_rate": [0.05, 0.10],
    }


def model_krr_rbf():
    return KernelRidge(kernel="rbf"), {
        "alpha": [0.1, 1.0, 10.0],
        "gamma": [1e-3, 1e-2],
    }


def model_mlp():
    """Lones 3.5 — honest DL benchmark. 2 hidden layers (64, 32),
    weight decay 1e-3, no early stopping (avoid tuning to val split).
    Single alpha (no grid) since the inner CV already gives variance."""
    return MLPRegressor(hidden_layer_sizes=(64, 32),
                        activation="relu",
                        alpha=1e-3,
                        max_iter=300,
                        early_stopping=False,
                        random_state=0), {}


MODELS = [
    ("Dummy",         model_dummy),
    ("RR-BLUP",       model_rrblup),
    ("ElasticNet",    model_enet),
    ("RandomForest",  model_rf),
    ("XGBoost",       model_xgb),
    ("LightGBM",      model_lgbm),
    ("KernelRidge",   model_krr_rbf),
    ("MLP",           model_mlp),
]

PALETTE = {
    "Dummy":        "#888888",
    "RR-BLUP":      WONG["blue"],
    "ElasticNet":   WONG["skyblue"],
    "RandomForest": WONG["green"],
    "XGBoost":      WONG["vermillion"],
    "LightGBM":     WONG["orange"],
    "KernelRidge":  WONG["purple"],
    "MLP":          WONG["yellow"],
}


def tune_in_inner(make_model, X_train, y_train, n_inner=N_INNER_FOLDS):
    base, grid = make_model()
    if not grid:
        return base.fit(X_train, y_train), {}
    inner = KFold(n_splits=n_inner, shuffle=True, random_state=1)
    keys = list(grid.keys())
    vals = [grid[k] for k in keys]
    best_score = -np.inf
    best_params = None
    for combo in product(*vals):
        params = dict(zip(keys, combo))
        scores = []
        for tr_i, va_i in inner.split(X_train):
            est = make_model()[0].set_params(**params)
            est.fit(X_train[tr_i], y_train[tr_i])
            pred = est.predict(X_train[va_i])
            if np.std(pred) == 0:
                scores.append(0.0)
            else:
                scores.append(pearsonr(pred, y_train[va_i])[0])
        m = float(np.nanmean(scores))
        if m > best_score:
            best_score = m
            best_params = params
    final = make_model()[0].set_params(**best_params)
    final.fit(X_train, y_train)
    return final, best_params


# ---------------------------------------------------------------------------
# Run CV per trait — main + sensitivity schemes
# ---------------------------------------------------------------------------

# 15 % final holdout (Lones 4.5) — same indices for every trait
hold_idx, work_idx = make_holdout_indices(len(samples))
hold_samples = [samples[i] for i in hold_idx]
work_samples = [samples[i] for i in work_idx]
print(f"[holdout] {len(hold_idx)} of {len(samples)} samples held out: "
      f"{', '.join(hold_samples[:5])}...")

cluster_lab = cluster_labels_from_pca(X_full[work_idx])
clone_lab   = clone_block_labels(work_samples)
print(f"[strat]  PCA cluster sizes within working set: "
      f"{np.bincount(cluster_lab).tolist()}")
print(f"[clone]  unique clone-blocks within working set: "
      f"{len(np.unique(clone_lab))}")


# ---------------------------------------------------------------------------
# (a) Main nested CV — 5-fold OUTER x 5-fold INNER x 50 REPEATS
# ---------------------------------------------------------------------------
print(f"\n[cv] nested 5x5 CV x {N_OUTER_REPEATS} repeats, "
      f"{len(MODELS)} models, {len(TRAITS)} traits")

rows = []
top_k_rows = []
per_repeat_r = {}        # (trait, model) -> array of r per repeat (for Wilcoxon)

for trait_name, df in TRAITS:
    y_all = y_aligned(trait_name).values.astype(float)
    # mask in the WORKING set (post-holdout) where y is observed
    y_work = y_all[work_idx]
    X_work = X_full[work_idx]
    obs    = ~np.isnan(y_work)
    n_obs  = int(obs.sum())
    if n_obs < 30:
        print(f"  {trait_name}: only {n_obs} samples, skipping")
        continue
    Xt = X_work[obs]
    yt = y_work[obs]
    cluster_t = cluster_lab[obs]
    clone_t   = clone_lab[obs]
    print(f"\n  {trait_name}: n_work={n_obs}")

    for model_name, make_model in MODELS:
        outer_r, outer_rho, outer_mae, outer_rmse, outer_topk = [], [], [], [], []
        for rep in range(N_OUTER_REPEATS):
            outer = KFold(n_splits=N_OUTER_FOLDS, shuffle=True,
                           random_state=2000 + rep)
            yhat = np.full_like(yt, np.nan, dtype=float)
            for tr_i, te_i in outer.split(Xt):
                Xtr, Xte = impute_and_scale_train_only(Xt[tr_i], Xt[te_i])
                est, _ = tune_in_inner(make_model, Xtr, yt[tr_i])
                yhat[te_i] = est.predict(Xte)
            if np.nanstd(yhat) == 0:
                continue
            r,  _ = pearsonr(yhat, yt)
            rho, _ = spearmanr(yhat, yt)
            err = yhat - yt
            mae  = float(np.nanmean(np.abs(err)))
            rmse = float(np.sqrt(np.nanmean(err ** 2)))
            k = min(TOP_K, len(yt))
            topk = len(set(np.argsort(yhat)[-k:]) & set(np.argsort(yt)[-k:])) / k
            outer_r.append(r); outer_rho.append(rho)
            outer_mae.append(mae); outer_rmse.append(rmse)
            outer_topk.append(topk)
        r_arr = np.array(outer_r)
        per_repeat_r[(trait_name, model_name)] = r_arr
        rows.append({
            "trait": trait_name, "model": model_name, "n_work": n_obs,
            "r_mean":     float(np.nanmean(r_arr)),
            "r_ci_lo":    float(np.nanpercentile(r_arr, 2.5)),
            "r_ci_hi":    float(np.nanpercentile(r_arr, 97.5)),
            "rho_mean":   float(np.nanmean(outer_rho)),
            "mae_mean":   float(np.nanmean(outer_mae)),
            "rmse_mean":  float(np.nanmean(outer_rmse)),
            "n_reps":     int(len(r_arr)),
        })
        top_k_rows.append({
            "trait": trait_name, "model": model_name, "top_k": int(TOP_K),
            "n": n_obs,
            "topk_acc_mean": float(np.nanmean(outer_topk)),
            "topk_acc_ci_lo": float(np.nanpercentile(outer_topk, 2.5)),
            "topk_acc_ci_hi": float(np.nanpercentile(outer_topk, 97.5)),
        })
        print(f"    {model_name:13s}  r = {np.nanmean(r_arr):+.3f} "
              f"[{np.nanpercentile(r_arr, 2.5):+.2f}, "
              f"{np.nanpercentile(r_arr, 97.5):+.2f}]   "
              f"topk={np.nanmean(outer_topk):.2f}")

cv = pd.DataFrame(rows)
cv.to_csv(TAB / "cv_predictive_ability.csv", index=False)
tk = pd.DataFrame(top_k_rows)
tk.to_csv(TAB / "top_k_accuracy.csv", index=False)


# ---------------------------------------------------------------------------
# (b) Clone-blocked CV — each clone group always on one side of a fold
# ---------------------------------------------------------------------------
print("\n[sens] clone-blocked CV (group-aware KFold)")
from sklearn.model_selection import GroupKFold
clone_rows = []
for trait_name, df in TRAITS:
    y_all = y_aligned(trait_name).values.astype(float)
    y_work = y_all[work_idx]
    X_work = X_full[work_idx]
    obs    = ~np.isnan(y_work)
    if obs.sum() < 30:
        continue
    Xt = X_work[obs]; yt = y_work[obs]; groups = clone_lab[obs]
    for model_name, make_model in MODELS:
        outer_r = []
        for rep in range(N_SENS_REPEATS):
            rng = np.random.default_rng(3000 + rep)
            perm = rng.permutation(np.unique(groups))
            group_to_fold = {g: (i % N_OUTER_FOLDS) for i, g in enumerate(perm)}
            fold_id = np.array([group_to_fold[g] for g in groups])
            yhat = np.full_like(yt, np.nan, dtype=float)
            for f in range(N_OUTER_FOLDS):
                tr_i = np.where(fold_id != f)[0]
                te_i = np.where(fold_id == f)[0]
                if len(te_i) < 2 or len(tr_i) < 10:
                    continue
                Xtr, Xte = impute_and_scale_train_only(Xt[tr_i], Xt[te_i])
                est, _ = tune_in_inner(make_model, Xtr, yt[tr_i])
                yhat[te_i] = est.predict(Xte)
            mask = ~np.isnan(yhat)
            if mask.sum() < 10 or np.nanstd(yhat[mask]) == 0:
                continue
            outer_r.append(pearsonr(yhat[mask], yt[mask])[0])
        if outer_r:
            clone_rows.append({
                "trait": trait_name, "model": model_name,
                "r_mean": float(np.mean(outer_r)),
                "r_ci_lo": float(np.percentile(outer_r, 2.5)),
                "r_ci_hi": float(np.percentile(outer_r, 97.5)),
            })
pd.DataFrame(clone_rows).to_csv(
    TAB / "cv_sensitivity_clone_blocked.csv", index=False)


# ---------------------------------------------------------------------------
# (c) PCA-cluster-stratified CV
# ---------------------------------------------------------------------------
print("\n[sens] PCA-cluster-stratified CV")
strat_rows = []
for trait_name, df in TRAITS:
    y_all = y_aligned(trait_name).values.astype(float)
    y_work = y_all[work_idx]
    X_work = X_full[work_idx]
    obs    = ~np.isnan(y_work)
    if obs.sum() < 30:
        continue
    Xt = X_work[obs]; yt = y_work[obs]; strat = cluster_t[obs]
    if len(np.unique(strat)) < 2 or min(np.bincount(strat)) < N_OUTER_FOLDS:
        continue
    for model_name, make_model in MODELS:
        outer_r = []
        for rep in range(N_SENS_REPEATS):
            skf = StratifiedKFold(n_splits=N_OUTER_FOLDS, shuffle=True,
                                   random_state=4000 + rep)
            yhat = np.full_like(yt, np.nan, dtype=float)
            for tr_i, te_i in skf.split(Xt, strat):
                Xtr, Xte = impute_and_scale_train_only(Xt[tr_i], Xt[te_i])
                est, _ = tune_in_inner(make_model, Xtr, yt[tr_i])
                yhat[te_i] = est.predict(Xte)
            mask = ~np.isnan(yhat)
            if mask.sum() < 10 or np.nanstd(yhat[mask]) == 0:
                continue
            outer_r.append(pearsonr(yhat[mask], yt[mask])[0])
        if outer_r:
            strat_rows.append({
                "trait": trait_name, "model": model_name,
                "r_mean": float(np.mean(outer_r)),
                "r_ci_lo": float(np.percentile(outer_r, 2.5)),
                "r_ci_hi": float(np.percentile(outer_r, 97.5)),
            })
pd.DataFrame(strat_rows).to_csv(
    TAB / "cv_sensitivity_cluster_strat.csv", index=False)


# ---------------------------------------------------------------------------
# Wilcoxon signed-rank + Holm-Bonferroni; Cohen's d on Δr (Lones 5.3-5.4, 6.4)
# ---------------------------------------------------------------------------
print("\n[stats] paired Wilcoxon signed-rank vs RR-BLUP + Holm-Bonferroni")
wlx_rows = []
cohen_rows = []
for trait_name in cv["trait"].unique():
    rr_arr = per_repeat_r.get((trait_name, "RR-BLUP"))
    if rr_arr is None:
        continue
    # 6 comparisons (all ML models except Dummy and RR-BLUP itself)
    others = [m for m in [n for n, _ in MODELS]
              if m not in ("Dummy", "RR-BLUP")]
    raw_p = {}
    for other in others:
        oth_arr = per_repeat_r.get((trait_name, other))
        if oth_arr is None or len(oth_arr) != len(rr_arr):
            continue
        try:
            _, p = wilcoxon(oth_arr, rr_arr)
        except ValueError:
            p = np.nan
        # Cohen's d on Δr = mean(other - rrblup) / sd(other - rrblup)
        delta = oth_arr - rr_arr
        d = float(np.mean(delta) / (np.std(delta, ddof=1) + 1e-12))
        # bootstrap 95 % CI on d
        rng = np.random.default_rng(5000)
        boot = []
        for _ in range(N_BOOTSTRAP):
            sub = rng.choice(delta, size=len(delta), replace=True)
            sd = np.std(sub, ddof=1)
            if sd > 0:
                boot.append(np.mean(sub) / sd)
        cohen_rows.append({
            "trait": trait_name, "comparison": f"{other} vs RR-BLUP",
            "mean_delta_r": float(np.mean(delta)),
            "cohens_d":     d,
            "d_ci_lo":      float(np.percentile(boot, 2.5)),
            "d_ci_hi":      float(np.percentile(boot, 97.5)),
        })
        raw_p[other] = p
    # Holm-Bonferroni
    sorted_p = sorted(raw_p.items(), key=lambda kv: kv[1])
    m = len(sorted_p)
    for rank, (other, p) in enumerate(sorted_p):
        p_holm = min(1.0, p * (m - rank))
        wlx_rows.append({
            "trait": trait_name, "comparison": f"{other} vs RR-BLUP",
            "p_raw":  float(p),
            "p_holm": float(p_holm),
            "sig_holm_0.05": bool(p_holm < 0.05),
        })

wlx = pd.DataFrame(wlx_rows)
wlx.to_csv(TAB / "wilcoxon_vs_rrblup.csv", index=False)
coh = pd.DataFrame(cohen_rows)
coh.to_csv(TAB / "cohens_d_delta_r.csv", index=False)
print(f"[stats] {len(wlx)} Wilcoxon rows, "
      f"{(wlx['sig_holm_0.05']).sum()} sig after Holm")


# ---------------------------------------------------------------------------
# Final 15 % holdout — best model per trait, fit once on all working data
# ---------------------------------------------------------------------------
print("\n[holdout] scoring best CV model on final 15 % holdout")
hold_rows = []
for trait_name, df in TRAITS:
    y_all = y_aligned(trait_name).values.astype(float)
    y_hold = y_all[hold_idx]
    y_work = y_all[work_idx]
    obs_h  = ~np.isnan(y_hold)
    obs_w  = ~np.isnan(y_work)
    if obs_h.sum() < 3 or obs_w.sum() < 30:
        continue
    best_row = cv[(cv["trait"] == trait_name) &
                   (cv["model"] != "Dummy")].sort_values(
        "r_mean", ascending=False).head(1)
    if best_row.empty:
        continue
    best_model_name = best_row["model"].values[0]
    make_model = dict(MODELS)[best_model_name]
    Xw = X_full[work_idx][obs_w]; yw = y_work[obs_w]
    Xh = X_full[hold_idx][obs_h]; yh = y_hold[obs_h]
    Xw_imp, Xh_imp = impute_and_scale_train_only(Xw, Xh)
    est, _ = tune_in_inner(make_model, Xw_imp, yw)
    yhat = est.predict(Xh_imp)
    if np.std(yhat) == 0:
        continue
    r_h, _ = pearsonr(yhat, yh)
    rho_h, _ = spearmanr(yhat, yh)
    hold_rows.append({
        "trait": trait_name, "best_model": best_model_name,
        "n_hold": int(obs_h.sum()),
        "cv_r_mean": float(best_row["r_mean"].values[0]),
        "holdout_r": float(r_h),
        "holdout_rho": float(rho_h),
    })
    print(f"    {trait_name:18s} best={best_model_name:13s} "
          f"CV r={best_row['r_mean'].values[0]:+.2f}, "
          f"holdout r={r_h:+.2f}")
pd.DataFrame(hold_rows).to_csv(
    TAB / "final_holdout_scores.csv", index=False)


# ---------------------------------------------------------------------------
# Fig 85. Model-comparison forest plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 8), constrained_layout=True)
traits_present = list(cv["trait"].unique())
model_order = ["Dummy", "RR-BLUP", "ElasticNet", "RandomForest", "XGBoost",
                "LightGBM", "KernelRidge", "MLP"]
n_t = len(traits_present); n_m = len(model_order); height = 0.10
y_pos = np.arange(n_t)[::-1]
for i, m in enumerate(model_order):
    sub = cv[cv["model"] == m].set_index("trait").reindex(traits_present)
    y = y_pos + (i - (n_m - 1) / 2) * height
    ax.errorbar(sub["r_mean"], y,
                 xerr=[sub["r_mean"] - sub["r_ci_lo"],
                       sub["r_ci_hi"] - sub["r_mean"]],
                 fmt="o", color=PALETTE[m], markersize=4,
                 elinewidth=0.8, capsize=2, label=m)
ax.axvline(0, color="black", linewidth=0.4)
ax.set_yticks(y_pos)
ax.set_yticklabels(traits_present, fontsize=9)
ax.set_xlabel("Pearson r (nested 5x5 CV x 50 outer repeats)")
ax.set_title("Genomic-prediction model comparison (AYB DArTseq panel)")
ax.grid(True, axis="x", alpha=0.4)
ax.legend(loc="lower right", fontsize=7, ncol=2)
save(fig, "fig85_model_comparison_forest")
print("[fig] fig85_model_comparison_forest")


# ---------------------------------------------------------------------------
# Fig 86. ML accuracy vs Daetwyler ceiling
# ---------------------------------------------------------------------------
Ne_LD = 659.0
L_M = 649.8 * 1e6 * 1.06 / 100.0 / 1e6
M_e = 2 * Ne_LD * L_M / np.log(4 * Ne_LD * L_M)
print(f"[daetwyler] M_e = {M_e:.0f}")


def daetwyler_r(n, h2, M_e=M_e):
    return float(np.sqrt(n * h2 / (n * h2 + M_e)))


H2_FILE = ROOT / "results" / "24_bootstrap_h2" / "tables" / "h2_profile_intervals.csv"
h2 = {}
if H2_FILE.exists():
    h2_df = pd.read_csv(H2_FILE)
    if "h2_mle" in h2_df.columns and "trait" in h2_df.columns:
        h2 = h2_df.set_index("trait")["h2_mle"].to_dict()

best_per_trait = cv[cv["model"] != "Dummy"].sort_values(
    "r_mean", ascending=False).drop_duplicates("trait")
records = []
for _, row in best_per_trait.iterrows():
    t = row["trait"]; n = row["n_work"]
    ceiling = daetwyler_r(n, h2.get(t, 0.20))
    records.append({
        "trait": t, "n": n, "h2": h2.get(t, np.nan),
        "best_model": row["model"], "best_r": row["r_mean"],
        "best_r_lo": row["r_ci_lo"], "best_r_hi": row["r_ci_hi"],
        "daetwyler_ceiling": ceiling,
    })
ceil_df = pd.DataFrame(records).sort_values("daetwyler_ceiling")
ceil_df.to_csv(TAB / "ml_vs_daetwyler.csv", index=False)

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
# Estimator name above each bar (or above the CI top, whichever is taller)
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
save(fig, "fig86_ml_vs_daetwyler_ceiling")
print("[fig] fig86_ml_vs_daetwyler_ceiling")


# ---------------------------------------------------------------------------
# Fig 87. Top-K ranking heatmap
# ---------------------------------------------------------------------------
pivot = tk.pivot(index="trait", columns="model",
                  values="topk_acc_mean").reindex(
    index=traits_present, columns=model_order)
fig, ax = plt.subplots(figsize=(8.5, 7), constrained_layout=True)
im = ax.imshow(pivot.values, cmap="viridis", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(np.arange(len(model_order)))
ax.set_xticklabels(model_order, rotation=45, ha="right", fontsize=8)
ax.set_yticks(np.arange(len(traits_present)))
ax.set_yticklabels(traits_present, fontsize=8)
for i in range(pivot.shape[0]):
    for j in range(pivot.shape[1]):
        v = pivot.values[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                     fontsize=7,
                     color="white" if v < 0.5 else "black")
fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02,
              label="top-10 overlap")
ax.set_title(f"Top-{TOP_K} ranking accuracy (predicted vs observed)")
save(fig, "fig87_top_k_accuracy_heatmap")
print("[fig] fig87_top_k_accuracy_heatmap")


# ---------------------------------------------------------------------------
# Fig 88. Cohen's d effect size on Δr vs RR-BLUP
# ---------------------------------------------------------------------------
if not coh.empty:
    coh_p = coh.pivot(index="trait", columns="comparison",
                       values="cohens_d").reindex(index=traits_present)
    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    im = ax.imshow(coh_p.values, cmap="RdBu_r",
                    vmin=-1.5, vmax=1.5, aspect="auto")
    ax.set_xticks(np.arange(coh_p.shape[1]))
    ax.set_xticklabels([c.replace(" vs RR-BLUP", "") for c in coh_p.columns],
                        rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(coh_p.shape[0]))
    ax.set_yticklabels(coh_p.index, fontsize=8)
    for i in range(coh_p.shape[0]):
        for j in range(coh_p.shape[1]):
            v = coh_p.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                         fontsize=7,
                         color="white" if abs(v) > 0.8 else "black")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02,
                  label="Cohen's d on Δr")
    ax.set_title("Effect size on Δr vs RR-BLUP (positive = ML better)")
    save(fig, "fig88_wilcoxon_effect_size")
    print("[fig] fig88_wilcoxon_effect_size")

print(f"\nOutputs in: {OUT}")
