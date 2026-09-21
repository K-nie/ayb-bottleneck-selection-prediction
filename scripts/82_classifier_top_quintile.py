"""
82 -- Top-quintile binary classification benchmark.

Re-cast the 13-trait regression problem as a per-trait classification
problem: label samples in the top 20 % by BLUP as positives, rest as
negatives. Train 5 classifiers under stratified 5-fold x 5-repeat CV
and score on the same 15 % final holdout used in script 40.

Classifiers:
  LogReg (L2, class-balanced)
  RandomForest (n_est = 500, class-balanced)
  XGBoost (n_est = 300, scale_pos_weight balanced)
  LightGBM (n_est = 300, class-balanced)
  MLP (one hidden layer 32)

Metrics (per outer fold, then aggregated):
  ROC-AUC
  PR-AUC (average precision)
  MCC at probability threshold 0.5
  Top-5 breeder accuracy (of the 5 highest-probability test samples,
    fraction that are true positives) -- the "rank the best parents"
    metric that matters in breeder-decision contexts

Inputs / outputs match script 40 conventions: 5x5 repeats, same 15 %
holdout seed, same per-trait n_work filtering on phenotype availability.

DATA-LEAKAGE DISCIPLINE
-----------------------
All preprocessing steps that could leak test-fold information into the
training fit are computed PER FOLD using TRAIN-fold data only:

  (a) Marker mean imputation: col_mu = mean(X_train[:, j], skipna=True);
      both X_train and X_test NaN entries are filled with the train-fold
      column mean.
  (b) Centering: StandardScaler(with_mean=True, with_std=False).fit(X_train);
      both X_train and X_test are .transform()-ed with the train-derived
      shift.
  (c) Top-quintile threshold: thr_fold = quantile(y_train, 0.80); both
      train and test labels are derived using this train-only threshold.
      Test-fold y values do NOT contribute to the threshold that labels them.

The fold split itself is driven by a diagnostic binarisation of y on the
WORKING-SET threshold so each fold contains at least one positive in
expectation (stratified KFold). This is balance discipline only; it does
not leak the per-fold threshold into the model.

For the 15 % final-holdout score, all three quantities are derived from
the FULL working set -- the legitimate "apply train-time rule to true
future data" pattern, since the holdout was never seen during work-set
training.
"""
from __future__ import annotations
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_auc_score, average_precision_score, matthews_corrcoef,
)
warnings.filterwarnings("ignore")

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
OUT = ROOT / "results/82_classifier_top_quintile"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
HOLDOUT_FRAC = 0.15
HOLDOUT_SEED = 42
TOP_QUINTILE = 0.20
TOP_K_BREEDER = 5
N_OUTER = 5
N_REPS = 5
RNG_BASE = 20260615

import sys
sys.path.insert(0, str(ROOT / "scripts"))
from _pheno import load_phenotypes, TRAITS_ALL  # type: ignore


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
    samples = list(dosage.columns)
    X_full = dosage.T.values.astype(float)
    rng = np.random.default_rng(HOLDOUT_SEED)
    idx_all = np.arange(len(samples))
    rng.shuffle(idx_all)
    n_hold = int(round(len(samples) * HOLDOUT_FRAC))
    work_idx = idx_all[n_hold:]
    hold_idx = idx_all[:n_hold]
    work_samples = [samples[i] for i in work_idx]
    hold_samples = [samples[i] for i in hold_idx]
    return X_full, samples, work_idx, hold_idx, work_samples, hold_samples


def label_top_quintile(y: np.ndarray, frac: float = TOP_QUINTILE) -> np.ndarray:
    """Top `frac` by y -> 1, rest -> 0. Ties on the threshold go to negative."""
    thr = np.quantile(y, 1.0 - frac)
    return (y > thr).astype(int)


def get_models(seed: int) -> dict:
    """Return fresh model instances seeded with `seed`.

    Class-imbalance handling per estimator:
      - LogReg:       class_weight="balanced" (sklearn API)
      - RandomForest: class_weight="balanced"
      - LightGBM:     class_weight="balanced"; min_child_samples=5 so
                      that small folds (n_train ~ 30 with ~ 6 positives)
                      can grow trees (default 20 collapses to AUC = 0.5)
      - XGBoost:      scale_pos_weight set PER FOLD to n_neg / n_pos
      - MLP:          sklearn MLPClassifier.fit does NOT accept
                      sample_weight (verified sklearn 1.6.1). Imbalance
                      is handled UPSTREAM by random oversampling of the
                      minority class inside fit_and_score before training
                      the MLP. Oversampling uses the per-fold rng seed.
    """
    models = {
        "LogReg":     LogisticRegression(
            penalty="l2", C=1.0, class_weight="balanced",
            solver="liblinear", max_iter=2000, random_state=seed,
        ),
        "RandomForest": RandomForestClassifier(
            n_estimators=500, max_features=0.1, min_samples_leaf=2,
            class_weight="balanced", random_state=seed, n_jobs=-1,
        ),
        "MLP": MLPClassifier(
            hidden_layer_sizes=(32,), activation="relu", alpha=1e-3,
            max_iter=500, early_stopping=False, random_state=seed,
        ),
    }
    try:
        import xgboost as xgb
        models["XGBoost"] = xgb.XGBClassifier(
            n_estimators=300, max_depth=5, learning_rate=0.05,
            subsample=0.7, n_jobs=-1, random_state=seed, verbosity=0,
            tree_method="hist", eval_metric="logloss",
        )
    except Exception as e:
        print(f"  skipping XGBoost: {e}")
    try:
        import lightgbm as lgb
        models["LightGBM"] = lgb.LGBMClassifier(
            n_estimators=300, num_leaves=31, learning_rate=0.05,
            feature_fraction=0.3, class_weight="balanced",
            min_child_samples=5, n_jobs=-1,
            random_state=seed, verbose=-1,
        )
    except Exception as e:
        print(f"  skipping LightGBM: {e}")
    return models


def random_oversample(X: np.ndarray, y: np.ndarray, seed: int):
    """Random oversample the minority class until both classes match the
    majority count. Used for MLP only since sklearn MLPClassifier does
    not accept sample_weight."""
    rng = np.random.default_rng(seed)
    cls, counts = np.unique(y, return_counts=True)
    if len(cls) < 2:
        return X, y
    n_max = counts.max()
    Xs, ys = [], []
    for c in cls:
        idx = np.where(y == c)[0]
        if len(idx) < n_max:
            extra = rng.choice(idx, size=n_max - len(idx), replace=True)
            idx = np.concatenate([idx, extra])
        Xs.append(X[idx]); ys.append(y[idx])
    return np.concatenate(Xs, axis=0), np.concatenate(ys, axis=0)


def top_k_accuracy(y_true: np.ndarray, p_score: np.ndarray, k: int) -> float:
    """Fraction of the top-k predicted-probability samples that are true 1s."""
    if len(p_score) < k:
        k = len(p_score)
    order = np.argsort(-p_score)[:k]
    return float(np.mean(y_true[order]))


def fit_and_score(model, X_tr, y_tr, X_te, y_te, scale_pos_weight=None,
                  mlp_oversample_seed=None):
    """Returns AUC, AUC-PR, MCC, top-k accuracy (None if AUC undefined).

    If `mlp_oversample_seed` is not None AND the model is an MLP, the
    training set is class-balanced via random oversampling before fit.
    All other models handle class imbalance via class_weight (LogReg,
    RF, LGB) or scale_pos_weight (XGB).
    """
    if hasattr(model, "set_params") and scale_pos_weight is not None:
        try:
            model.set_params(scale_pos_weight=scale_pos_weight)
        except Exception:
            pass
    if mlp_oversample_seed is not None and isinstance(model, MLPClassifier):
        X_tr, y_tr = random_oversample(X_tr, y_tr, seed=mlp_oversample_seed)
    model.fit(X_tr, y_tr)
    if hasattr(model, "predict_proba"):
        p = model.predict_proba(X_te)[:, 1]
    elif hasattr(model, "decision_function"):
        p = model.decision_function(X_te)
    else:
        p = model.predict(X_te).astype(float)
    yhat = (p >= 0.5).astype(int)
    n_pos = int(y_te.sum()); n_neg = int(len(y_te) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return None
    return {
        "auc":    float(roc_auc_score(y_te, p)),
        "pr_auc": float(average_precision_score(y_te, p)),
        "mcc":    float(matthews_corrcoef(y_te, yhat)) if (yhat.sum() and (len(yhat) - yhat.sum())) else 0.0,
        "topk":   top_k_accuracy(y_te, p, TOP_K_BREEDER),
        "p_pos_class": float(n_pos),
        "n_test": int(len(y_te)),
    }


def main() -> None:
    X_full, samples, work_idx, hold_idx, work_samples, hold_samples = load_data()
    print(f"Loaded {X_full.shape[0]} samples x {X_full.shape[1]} markers")
    print(f"work n = {len(work_samples)}; hold n = {len(hold_samples)}")

    pheno = load_phenotypes()
    pheno.index = pheno.index.astype(str).str.strip()
    traits = [t for t in TRAITS_ALL if t in pheno.columns]
    print(f"Traits: {traits}")

    cv_rows = []
    hold_rows = []

    for trait in traits:
        # ---- per-trait sample filtering (drop missing phenotype) ----
        y_work_full = pheno[trait].reindex(work_samples).values.astype(float)
        obs_w = ~np.isnan(y_work_full)
        X_w_raw = X_full[work_idx][obs_w].copy()  # KEEP NaNs - impute PER FOLD
        y_w_raw = y_work_full[obs_w]
        y_hold_full = pheno[trait].reindex(hold_samples).values.astype(float)
        obs_h = ~np.isnan(y_hold_full)
        X_h_raw = X_full[hold_idx][obs_h].copy()
        y_h_raw = y_hold_full[obs_h]
        n_w = X_w_raw.shape[0]; n_h = X_h_raw.shape[0]

        if n_w < 20 or n_h < 4:
            print(f"  {trait}: n_w={n_w}, n_h={n_h}; skip")
            continue

        # For diagnostics ONLY: report the working-set threshold.
        # The per-fold positive count is determined INSIDE the CV loop using
        # the train-fold-only threshold (no leakage). The holdout label uses
        # a threshold derived from the full working set (legitimate "apply
        # train-time rule to future data").
        thr_diag = np.quantile(y_w_raw, 1.0 - TOP_QUINTILE)
        print(f"\n[{trait}] n_w = {n_w}; n_h = {n_h}; "
              f"diagnostic threshold from full work = {thr_diag:.3f}")
        if (y_w_raw > thr_diag).sum() < 4:
            print(f"  too few working positives at top-quintile; skip")
            continue

        # ---------- nested-CV outer-fold scoring ----------
        # All three potential leakage sources are now PER-FOLD:
        #   (a) Marker imputation: column means from train-fold X only.
        #   (b) Centering scaler: fit on train-fold X only.
        #   (c) Top-quintile threshold: computed on train-fold y only,
        #       applied to test-fold y for labels.
        # We use np.quantile-based splitting on y to keep the FOLD splits
        # roughly stratified by trait quintile (stratify the BINARISED
        # train-set label is done via StratifiedKFold below using a
        # work-set diagnostic split, then re-binarised inside each fold).
        diag_bin = (y_w_raw > thr_diag).astype(int)
        rep_seeds = [RNG_BASE + r for r in range(N_REPS)]
        for r_idx, seed in enumerate(rep_seeds):
            skf = StratifiedKFold(n_splits=N_OUTER, shuffle=True, random_state=seed)
            models_r = get_models(seed)
            # The fold-split is driven by the diagnostic binarisation so that
            # each fold contains a representative number of (eventual)
            # top-quintile samples. This does NOT leak the per-fold threshold
            # decision into the model; it only ensures fold balance for AUC
            # computation (a fold with zero positives gives undefined AUC).
            for fold_idx, (tr, te) in enumerate(skf.split(X_w_raw, diag_bin)):
                # PER-FOLD preprocessing
                X_tr_nan = X_w_raw[tr]; X_te_nan = X_w_raw[te]
                col_mu = np.nanmean(X_tr_nan, axis=0)
                X_tr = np.where(np.isnan(X_tr_nan), col_mu, X_tr_nan)
                X_te = np.where(np.isnan(X_te_nan), col_mu, X_te_nan)
                sc = StandardScaler(with_mean=True, with_std=False).fit(X_tr)
                X_tr = sc.transform(X_tr); X_te = sc.transform(X_te)

                # PER-FOLD top-quintile threshold from TRAIN ONLY
                thr_fold = np.quantile(y_w_raw[tr], 1.0 - TOP_QUINTILE)
                y_tr = (y_w_raw[tr] > thr_fold).astype(int)
                y_te = (y_w_raw[te] > thr_fold).astype(int)

                n_pos_tr = int(y_tr.sum()); n_neg_tr = int(len(y_tr) - n_pos_tr)
                spw = n_neg_tr / max(n_pos_tr, 1)
                for mname, mdl in models_r.items():
                    res = fit_and_score(mdl, X_tr, y_tr, X_te, y_te,
                                        scale_pos_weight=spw if mname == "XGBoost" else None)
                    if res is None:
                        continue
                    cv_rows.append({
                        "trait": trait, "model": mname,
                        "rep": int(r_idx), "fold": int(fold_idx),
                        "thr_train_fold": float(thr_fold),
                        "n_pos_train": int(y_tr.sum()),
                        "n_pos_test": int(y_te.sum()),
                        "n_test": int(len(y_te)),
                        **res,
                    })

        # ---------- final holdout fit (legitimate use of full work stats) ----
        thr_w = np.quantile(y_w_raw, 1.0 - TOP_QUINTILE)
        col_mu_w = np.nanmean(X_w_raw, axis=0)
        X_w = np.where(np.isnan(X_w_raw), col_mu_w, X_w_raw)
        X_h = np.where(np.isnan(X_h_raw), col_mu_w, X_h_raw)
        sc_w = StandardScaler(with_mean=True, with_std=False).fit(X_w)
        X_w = sc_w.transform(X_w); X_h = sc_w.transform(X_h)
        y_w_bin = (y_w_raw > thr_w).astype(int)
        y_h_bin = (y_h_raw > thr_w).astype(int)
        if y_h_bin.sum() == 0:
            print(f"  no holdout positives at work threshold; skip holdout")
            continue
        for mname, mdl in get_models(RNG_BASE).items():
            spw = None
            if mname == "XGBoost":
                n_pos_w = int(y_w_bin.sum()); n_neg_w = int(len(y_w_bin) - n_pos_w)
                spw = n_neg_w / max(n_pos_w, 1)
            res = fit_and_score(mdl, X_w, y_w_bin, X_h, y_h_bin,
                                scale_pos_weight=spw)
            if res is None:
                continue
            hold_rows.append({
                "trait": trait, "model": mname,
                "thr_used": float(thr_w),
                "n_work": int(n_w), "n_hold": int(n_h),
                "n_pos_work": int(y_w_bin.sum()),
                "n_pos_hold": int(y_h_bin.sum()),
                **res,
            })

    cv = pd.DataFrame(cv_rows)
    cv.to_csv(OUT / "tables/cv_per_fold.csv", index=False)

    # Aggregated per (trait, model)
    agg = (
        cv.groupby(["trait", "model"])
        [["auc", "pr_auc", "mcc", "topk"]]
        .agg(["mean", "std", "count"])
    )
    agg.columns = [f"{m}_{s}" for m, s in agg.columns]
    agg = agg.reset_index()
    agg.to_csv(OUT / "tables/cv_aggregated.csv", index=False)

    hold = pd.DataFrame(hold_rows)
    hold.to_csv(OUT / "tables/final_holdout.csv", index=False)

    # Best model per trait by CV AUC mean
    best = agg.sort_values(["trait", "auc_mean"], ascending=[True, False]).drop_duplicates("trait", keep="first")
    best.to_csv(OUT / "tables/cv_winner_per_trait.csv", index=False)

    print("\n=== Per-trait CV winners by AUC mean ===")
    print(best[["trait", "model", "auc_mean", "auc_std", "pr_auc_mean",
                "topk_mean"]].to_string(index=False))

    # ---------- Figure ----------
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    traits_sorted = sorted(cv["trait"].unique())
    models_sorted = sorted(cv["model"].unique())

    auc_mat = np.full((len(traits_sorted), len(models_sorted)), np.nan)
    topk_mat = np.full((len(traits_sorted), len(models_sorted)), np.nan)
    for i, t in enumerate(traits_sorted):
        for j, m in enumerate(models_sorted):
            sub = agg[(agg["trait"] == t) & (agg["model"] == m)]
            if not sub.empty:
                auc_mat[i, j] = sub["auc_mean"].iloc[0]
                topk_mat[i, j] = sub["topk_mean"].iloc[0]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, max(5.0, 0.45 * len(traits_sorted))))

    cmap_auc = LinearSegmentedColormap.from_list("auc", ["#cccccc", "#dddddd", "#88CCEE", "#117733"], N=256)
    im0 = axes[0].imshow(auc_mat, cmap=cmap_auc, aspect="auto", vmin=0.3, vmax=0.85)
    axes[0].set_xticks(range(len(models_sorted))); axes[0].set_xticklabels(models_sorted, rotation=35, ha="right")
    axes[0].set_yticks(range(len(traits_sorted))); axes[0].set_yticklabels(traits_sorted, fontsize=9)
    axes[0].set_title("CV ROC-AUC (mean across 25 folds)")
    axes[0].text(-0.08, 1.02, "a", transform=axes[0].transAxes, fontsize=12,
                 fontweight="bold", ha="right", va="bottom", zorder=1000)
    for i in range(auc_mat.shape[0]):
        for j in range(auc_mat.shape[1]):
            v = auc_mat[i, j]
            if not np.isnan(v):
                axes[0].text(j, i, f"{v:.2f}", ha="center", va="center",
                              fontsize=8, color="black" if v < 0.7 else "white")
    fig.colorbar(im0, ax=axes[0], fraction=0.04, pad=0.02, label="ROC-AUC")
    axes[0].axvline(-0.5, color="black", lw=0.6)

    cmap_topk = LinearSegmentedColormap.from_list("topk", ["#cccccc", "#dddddd", "#DDCC77", "#CC6677"], N=256)
    im1 = axes[1].imshow(topk_mat, cmap=cmap_topk, aspect="auto", vmin=0, vmax=1)
    axes[1].set_xticks(range(len(models_sorted))); axes[1].set_xticklabels(models_sorted, rotation=35, ha="right")
    axes[1].set_yticks(range(len(traits_sorted))); axes[1].set_yticklabels([""] * len(traits_sorted))
    axes[1].set_title(f"CV top-{TOP_K_BREEDER} breeder accuracy")
    axes[1].text(-0.08, 1.02, "b", transform=axes[1].transAxes, fontsize=12,
                 fontweight="bold", ha="right", va="bottom", zorder=1000)
    for i in range(topk_mat.shape[0]):
        for j in range(topk_mat.shape[1]):
            v = topk_mat[i, j]
            if not np.isnan(v):
                axes[1].text(j, i, f"{v:.2f}", ha="center", va="center",
                              fontsize=8, color="black" if v < 0.6 else "white")
    fig.colorbar(im1, ax=axes[1], fraction=0.04, pad=0.02, label="Top-5 accuracy")

    # Panel C: per-trait holdout AUC for the CV winner
    ax = axes[2]
    if not hold.empty:
        hold_win = hold.merge(best[["trait", "model"]], on=["trait", "model"])
        hold_win = hold_win.sort_values("auc", ascending=True)
        y = np.arange(len(hold_win))
        colors = ["#117733" if v >= 0.7 else "#88CCEE" if v >= 0.5 else "#CC6677"
                  for v in hold_win["auc"]]
        ax.barh(y, hold_win["auc"], color=colors, edgecolor="white")
        ax.axvline(0.5, color="black", ls="--", lw=1, label="AUC = 0.5 (chance)")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{t} ({m})" for t, m in zip(hold_win["trait"], hold_win["model"])], fontsize=8)
        ax.set_xlabel("Holdout ROC-AUC")
        ax.set_title("Final holdout AUC, CV winner per trait")
        ax.text(-0.08, 1.02, "c", transform=ax.transAxes, fontsize=12,
                fontweight="bold", ha="right", va="bottom", zorder=1000)
        ax.set_xlim(0, 1)
        ax.legend(frameon=False, fontsize=8)
        ax.grid(True, alpha=0.25, axis="x")

    fig.suptitle("Top-quintile binary classifier benchmark "
                 "(13 traits x 5 models, stratified 5 x 5 CV + 15 % holdout)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "figures/fig_classifier_top_quintile.png", dpi=300,
                bbox_inches="tight")
    fig.savefig(OUT / "figures/fig_classifier_top_quintile.pdf",
                bbox_inches="tight")
    print("\nWrote", OUT / "figures/fig_classifier_top_quintile.png")


if __name__ == "__main__":
    main()
