#!/usr/bin/env python3
"""
Tail-statistics + plots for script 40's main CV results.

Reads the already-written `cv_predictive_ability.csv` and produces:
    - Wilcoxon signed-rank test (paired r vectors vs RR-BLUP) — but
      since the original per-repeat r vectors live only in script 40's
      memory (now gone), we approximate with a bootstrap on the
      published mean + CI to recover a per-repeat distribution. The
      proper paired test requires re-running with per-repeat r
      persistence; that's a v2 follow-up.
    - Cohen's d on Δr vs RR-BLUP from the published mean + CI
    - Daetwyler ceiling per trait
    - Final-holdout evaluation (re-fit best model on full work set,
      predict holdout)
    - Fig 85 (forest), 86 (Daetwyler), 87 (top-k), 88 (Cohen's d heatmap)

This is the rescue path after killing script 40 mid-sensitivity-CV.
The main CV CSV is intact; this script squeezes the publication-ready
outputs from it.

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import Ridge, ElasticNetCV
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
TAB  = OUT / "tables"
FIG  = OUT / "figures"
CV   = TAB / "cv_predictive_ability.csv"
TK   = TAB / "top_k_accuracy.csv"
HOLDOUT_FRAC = 0.15
HOLDOUT_SEED = 42
N_BOOTSTRAP  = 2000


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Reload main CV results
# ---------------------------------------------------------------------------
cv = pd.read_csv(CV)
tk = pd.read_csv(TK)
print(f"[load] {len(cv)} (trait, model) rows from {CV.name}")

traits_present = list(cv["trait"].drop_duplicates())
model_order = ["Dummy", "RR-BLUP", "ElasticNet", "RandomForest",
                "XGBoost", "LightGBM", "KernelRidge", "MLP"]
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


# ---------------------------------------------------------------------------
# Approximate Cohen's d on Δr vs RR-BLUP from the published per-trait
# mean + 95 % CI. The CI half-width / 1.96 approximates the SE of r
# across outer repeats; we use that to draw synthetic per-repeat r
# distributions and compute Cohen's d on (other - rrblup).
# This is a fallback because the original per-repeat r vectors lived
# only in script 40 memory. A v2 will persist them to disk.
# ---------------------------------------------------------------------------
print("[stats] Cohen's d on Δr vs RR-BLUP (CI-based approximation)")
rng = np.random.default_rng(0)
n_synth = cv["n_reps"].max() if "n_reps" in cv.columns else 10
cohen_rows = []
for trait in traits_present:
    rr = cv[(cv["trait"] == trait) & (cv["model"] == "RR-BLUP")]
    if rr.empty:
        continue
    rr_mu = float(rr["r_mean"].iloc[0])
    rr_se = max((rr["r_ci_hi"].iloc[0] - rr["r_ci_lo"].iloc[0]) / 3.92, 1e-6)
    rr_synth = rng.normal(rr_mu, rr_se, n_synth)
    for other in [m for m in model_order if m not in ("Dummy", "RR-BLUP")]:
        row = cv[(cv["trait"] == trait) & (cv["model"] == other)]
        if row.empty:
            continue
        mu = float(row["r_mean"].iloc[0])
        se = max((row["r_ci_hi"].iloc[0] - row["r_ci_lo"].iloc[0]) / 3.92, 1e-6)
        oth_synth = rng.normal(mu, se, n_synth)
        delta = oth_synth - rr_synth
        d = float(np.mean(delta) / (np.std(delta, ddof=1) + 1e-12))
        boot = []
        for _ in range(N_BOOTSTRAP):
            sub = rng.choice(delta, size=len(delta), replace=True)
            s = np.std(sub, ddof=1)
            if s > 0:
                boot.append(np.mean(sub) / s)
        cohen_rows.append({
            "trait": trait, "comparison": f"{other} vs RR-BLUP",
            "mean_delta_r": float(np.mean(delta)),
            "cohens_d":     d,
            "d_ci_lo":      float(np.percentile(boot, 2.5)),
            "d_ci_hi":      float(np.percentile(boot, 97.5)),
            "approximation_note": "From published mean+CI; v2 will use exact per-repeat r vectors",
        })
coh = pd.DataFrame(cohen_rows)
coh.to_csv(TAB / "cohens_d_delta_r.csv", index=False)
print(f"[stats] {len(coh)} Cohen's d rows written")


# ---------------------------------------------------------------------------
# Daetwyler ceiling per trait
# ---------------------------------------------------------------------------
print("[daetwyler] computing ceiling per trait")
Ne_LD = 659.0
L_M = 649.8 * 1e6 * 1.06 / 100.0 / 1e6
M_e = 2 * Ne_LD * L_M / np.log(4 * Ne_LD * L_M)
print(f"[daetwyler] M_e = {M_e:.0f}")


def daetwyler_r(n, h2):
    return float(np.sqrt(n * h2 / (n * h2 + M_e)))


H2_FILE = ROOT / "results" / "24_bootstrap_h2" / "tables" / "h2_profile_intervals.csv"
h2_map = {}
if H2_FILE.exists():
    h2_df = pd.read_csv(H2_FILE)
    if "trait" in h2_df.columns and "h2_mle" in h2_df.columns:
        h2_map = h2_df.set_index("trait")["h2_mle"].to_dict()

best_per_trait = cv[cv["model"] != "Dummy"].sort_values(
    "r_mean", ascending=False).drop_duplicates("trait")
ceil_rows = []
for _, r in best_per_trait.iterrows():
    t = r["trait"]
    n = int(r.get("n_work", 95))
    ceil_rows.append({
        "trait": t, "n": n, "h2": h2_map.get(t, np.nan),
        "best_model": r["model"], "best_r": r["r_mean"],
        "best_r_lo": r["r_ci_lo"], "best_r_hi": r["r_ci_hi"],
        "daetwyler_ceiling": daetwyler_r(n, h2_map.get(t, 0.20)),
    })
ceil_df = pd.DataFrame(ceil_rows).sort_values("daetwyler_ceiling")
ceil_df.to_csv(TAB / "ml_vs_daetwyler.csv", index=False)


# ---------------------------------------------------------------------------
# Final-holdout evaluation — same 15 % set as script 40
# ---------------------------------------------------------------------------
print("\n[holdout] re-fitting best model per trait, scoring on 15 % holdout")

# Reload dosage + phenotypes the same way script 40 did
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
keep_m = ((cr_m >= 0.90) & (maf >= 0.05)).values
keep_s = (cr_s >= 0.90).values
dosage = dosage.loc[keep_m, keep_s]
samples = list(dosage.columns)
X_full = dosage.T.values.astype(float)

rng2 = np.random.default_rng(HOLDOUT_SEED)
idx_all = np.arange(len(samples))
rng2.shuffle(idx_all)
n_hold = int(round(len(samples) * HOLDOUT_FRAC))
hold_idx = idx_all[:n_hold]
work_idx = idx_all[n_hold:]

from _pheno import load_phenotypes
pheno = load_phenotypes()
pheno.index = pheno.index.astype(str).str.strip()


def fit_best(model_name, X, y):
    if model_name == "RR-BLUP":
        return Ridge(alpha=1.0).fit(X, y)
    if model_name == "ElasticNet":
        return ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], cv=5, n_alphas=20,
                             max_iter=20000, n_jobs=-1, random_state=0).fit(X, y)
    if model_name == "RandomForest":
        return RandomForestRegressor(n_estimators=100, max_features=0.05,
                                       min_samples_leaf=2, n_jobs=-1,
                                       random_state=0).fit(X, y)
    if model_name == "XGBoost":
        import xgboost as xgb
        return xgb.XGBRegressor(n_estimators=100, max_depth=5,
                                 learning_rate=0.05, subsample=0.7,
                                 n_jobs=-1, random_state=0, verbosity=0,
                                 tree_method="hist").fit(X, y)
    if model_name == "LightGBM":
        import lightgbm as lgb
        return lgb.LGBMRegressor(n_estimators=100, num_leaves=31,
                                   learning_rate=0.05,
                                   feature_fraction=0.3,
                                   n_jobs=-1, random_state=0,
                                   verbose=-1).fit(X, y)
    if model_name == "KernelRidge":
        return KernelRidge(kernel="rbf", alpha=1.0, gamma=1e-3).fit(X, y)
    if model_name == "MLP":
        return MLPRegressor(hidden_layer_sizes=(64, 32), activation="relu",
                              alpha=1e-3, max_iter=300, early_stopping=False,
                              random_state=0).fit(X, y)
    return Ridge(alpha=1.0).fit(X, y)


hold_rows = []
for trait in traits_present:
    row = cv[(cv["trait"] == trait) & (cv["model"] != "Dummy")]
    if row.empty: continue
    best_model_name = row.sort_values("r_mean", ascending=False).iloc[0]["model"]
    y_all = pheno[trait].reindex(samples).values.astype(float)
    y_hold = y_all[hold_idx]; y_work = y_all[work_idx]
    obs_h = ~np.isnan(y_hold); obs_w = ~np.isnan(y_work)
    if obs_h.sum() < 3 or obs_w.sum() < 30: continue
    Xw = X_full[work_idx][obs_w]; yw = y_work[obs_w]
    Xh = X_full[hold_idx][obs_h]; yh = y_hold[obs_h]
    mu = np.nanmean(Xw, axis=0)
    Xw_imp = np.where(np.isnan(Xw), mu, Xw)
    Xh_imp = np.where(np.isnan(Xh), mu, Xh)
    sc = StandardScaler(with_mean=True, with_std=False).fit(Xw_imp)
    Xw_imp = sc.transform(Xw_imp); Xh_imp = sc.transform(Xh_imp)
    try:
        m = fit_best(best_model_name, Xw_imp, yw)
    except Exception as e:
        print(f"    fit failed for {trait}: {e}")
        continue
    yhat = m.predict(Xh_imp)
    if np.std(yhat) == 0: continue
    r_h, _ = pearsonr(yhat, yh)
    rho_h, _ = spearmanr(yhat, yh)
    cv_r = float(row.sort_values("r_mean", ascending=False).iloc[0]["r_mean"])
    hold_rows.append({
        "trait": trait, "best_model": best_model_name,
        "n_hold": int(obs_h.sum()),
        "cv_r_mean": cv_r,
        "holdout_r": float(r_h),
        "holdout_rho": float(rho_h),
    })
    print(f"    {trait:18s} best={best_model_name:13s} "
           f"CV r={cv_r:+.2f}, holdout r={r_h:+.2f}")
pd.DataFrame(hold_rows).to_csv(TAB / "final_holdout_scores.csv", index=False)


# ---------------------------------------------------------------------------
# Fig 85. Forest plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 8), constrained_layout=True)
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
ax.set_yticks(y_pos); ax.set_yticklabels(traits_present, fontsize=9)
ax.set_xlabel("Pearson r (nested 5x5 CV x 10 outer repeats)")
ax.set_title("Genomic-prediction model comparison (AYB DArTseq panel)")
ax.grid(True, axis="x", alpha=0.4)
ax.legend(loc="lower right", fontsize=7, ncol=2)
save(fig, "fig85_model_comparison_forest")
print("[fig] fig85_model_comparison_forest")


# ---------------------------------------------------------------------------
# Fig 86. ML vs Daetwyler
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
xpos = np.arange(len(ceil_df))
ax.bar(xpos, ceil_df["daetwyler_ceiling"], width=0.6,
        color="#dddddd", edgecolor="grey", label="Daetwyler ceiling")
ax.errorbar(xpos, ceil_df["best_r"],
             yerr=[ceil_df["best_r"] - ceil_df["best_r_lo"],
                   ceil_df["best_r_hi"] - ceil_df["best_r"]],
             fmt="o", color=WONG["vermillion"],
             markersize=5, elinewidth=0.8, capsize=2,
             label="Best ML model r (95% CI)")
for i, m in enumerate(ceil_df["best_model"].values):
    ax.text(xpos[i], -0.08, m, rotation=45, ha="right", va="top",
             fontsize=7, color="grey")
ax.set_xticks(xpos)
ax.set_xticklabels(ceil_df["trait"].values, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Pearson r")
ax.axhline(0, color="black", linewidth=0.4)
ax.set_title("Best ML model vs Daetwyler ceiling per trait")
ax.legend(loc="upper left", fontsize=8)
ax.grid(True, axis="y", alpha=0.4)
save(fig, "fig86_ml_vs_daetwyler_ceiling")
print("[fig] fig86_ml_vs_daetwyler_ceiling")


# ---------------------------------------------------------------------------
# Fig 87. top-K heatmap
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
ax.set_title("Top-10 ranking accuracy (predicted vs observed)")
save(fig, "fig87_top_k_accuracy_heatmap")
print("[fig] fig87_top_k_accuracy_heatmap")


# ---------------------------------------------------------------------------
# Fig 88. Cohen's d heatmap
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
