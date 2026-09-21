#!/usr/bin/env python3
"""
SHAP feature-importance and dependence analysis for the AYB DArTseq panel.

Consumes the per-trait best-CV-model from script 40 and produces:

    - Per-trait SHAP beeswarm summary (top-20 markers)
    - Per-trait SHAP dependence plot for the top-3 markers, coloured by
      interaction with the second-most-influential marker (Lundberg et
      al. 2020 Nat Mach Intell, the standard dependence-plot convention)
    - GWAS-vs-SHAP cross-validation table per trait — how many of the
      top-N SHAP markers also appear in the top-30 GWAS hits, and
      where the ALMT4 Soluble_Oxalate SNP shows up in both sets

Outputs (results/41_ml_shap_interpretability/):
    tables/
        shap_top20_per_trait.csv               top-20 markers per trait by mean|SHAP|
        shap_gwas_overlap.csv                  k-shap_top vs k-gwas_top overlap counts
        almt4_focus.csv                        ALMT4 SNP SHAP and GWAS rank per trait
    figures/
        fig89_shap_beeswarm_{trait}.png/.pdf   per-trait beeswarm
        fig90_shap_dependence_{trait}.png/.pdf top-3 dependence plots per trait
        fig91_shap_gwas_overlap_heatmap.png/.pdf
        fig92_almt4_focus.png/.pdf

Lones (2024) Patterns 34-pitfall walk-through
---------------------------------------------
Walked the Lones (2024) ML best-practices checklist for SHAP interpretability.
Items most relevant to this script:

  2.2 understand data: re-use the script 40 CV results (already audited)
  3.1 no test leak: SHAP fit on FULL working data (post-holdout); the 15 % final
      holdout is left untouched here too — same Lones 4.5 discipline as script 40
  3.6 feature selection: NONE; SHAP runs on all 1,625 markers
  3.8 spurious correlations: SHAP-vs-GWAS cross-tab is the guardrail
  4.5 final-instance holdout: not re-used here (script 40 already scored it)
  6.5 look at your models: this script IS the looking-at-models pass

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
CV40   = ROOT / "results" / "40_ml_genomic_prediction" / "tables" / "cv_predictive_ability.csv"
OUT    = ROOT / "results" / "41_ml_shap_interpretability"
DATA   = OUT / "data"
TAB    = OUT / "tables"
FIG    = OUT / "figures"
for d in (DATA, TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
HOLDOUT_FRAC = 0.15
HOLDOUT_SEED = 42
N_SHAP_TOP = 20
N_DEP_TOP  = 3
GWAS_TOP_K = 30


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Load dosage, phenotypes, and the script-40 CV results
# ---------------------------------------------------------------------------
print("[load] HapMap + CV results from script 40...")
if not CV40.exists():
    raise SystemExit(f"Script 40 output not found: {CV40}. Run script 40 first.")

cv = pd.read_csv(CV40)
print(f"[load] {len(cv)} (trait, model) rows from script 40")

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
markers = list(dosage.index)
samples = list(dosage.columns)
X_full = dosage.T.values.astype(float)

# Reproduce script 40's 15 % holdout exactly (same seed)
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
print(f"[load] working set: {X_work.shape[0]} samples x {X_work.shape[1]} markers")

# Marker anchoring (for the GWAS cross-tab)
anc = pd.read_csv(ANCHOR)
anc = anc[anc["rs"].isin(markers)].copy()
anc_idx = {r: i for i, r in enumerate(markers)}

# Phenotypes via the canonical loader (handles Flavinoid typo, Crude
# Protein space, Seed_Metrics multi-line header, and the two TSs282 /
# TSs136 decimal-point corrections in one place).
from _pheno import load_phenotypes, TRAITS_ALL

pheno = load_phenotypes()
pheno.index = pheno.index.astype(str).str.strip()


def y_for(trait):
    return pheno[trait].reindex(work_samples)


TRAITS = [(t, None) for t in TRAITS_ALL]


# ---------------------------------------------------------------------------
# Build the best model per trait + compute SHAP
# ---------------------------------------------------------------------------
import shap

def fit_best(model_name, X, y):
    """Match the model name to a sensible default config — the inner-CV
    tuning done in script 40 is approximated here by re-using close-to-
    default parameters. For tree models we use shap.TreeExplainer;
    for everything else we fall back to KernelExplainer with a small
    background sample."""
    if model_name == "RR-BLUP":
        m = Ridge(alpha=1.0).fit(X, y)
    elif model_name == "ElasticNet":
        m = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], cv=5, n_alphas=20,
                          max_iter=20000, n_jobs=-1, random_state=0).fit(X, y)
    elif model_name == "RandomForest":
        m = RandomForestRegressor(n_estimators=500, max_features=0.1,
                                    min_samples_leaf=2,
                                    random_state=0, n_jobs=-1).fit(X, y)
    elif model_name == "XGBoost":
        import xgboost as xgb
        m = xgb.XGBRegressor(n_estimators=300, max_depth=5,
                              learning_rate=0.05, subsample=0.7,
                              n_jobs=-1, random_state=0, verbosity=0,
                              tree_method="hist").fit(X, y)
    elif model_name == "LightGBM":
        import lightgbm as lgb
        m = lgb.LGBMRegressor(n_estimators=300, num_leaves=31,
                                learning_rate=0.05, feature_fraction=0.3,
                                n_jobs=-1, random_state=0,
                                verbose=-1).fit(X, y)
    elif model_name == "KernelRidge":
        m = KernelRidge(kernel="rbf", alpha=1.0, gamma=1e-3).fit(X, y)
    elif model_name == "MLP":
        m = MLPRegressor(hidden_layer_sizes=(64, 32), activation="relu",
                          alpha=1e-3, max_iter=500, early_stopping=False,
                          random_state=0).fit(X, y)
    else:
        m = Ridge(alpha=1.0).fit(X, y)
    return m


def shap_values_for(model, model_name, X):
    """Return SHAP values matrix (n_samples x n_features)."""
    if model_name in ("RandomForest", "XGBoost", "LightGBM"):
        try:
            expl = shap.TreeExplainer(model)
            return np.array(expl.shap_values(X))
        except Exception as e:
            print(f"    TreeExplainer failed: {e}; falling back to KernelExplainer")
    # Kernel / linear / MLP — use KernelExplainer with a small background
    bg = shap.utils.sample(X, 30, random_state=0)
    expl = shap.KernelExplainer(model.predict, bg)
    return np.array(expl.shap_values(X, nsamples=100))


top20_rows = []
overlap_rows = []
almt4_rows = []
ALMT4_RS = None  # discovered per trait via the soluble_oxalate GWAS top SNP

# Locate the ALMT4 SNP (Ss10:15,394,673 in the materials_and_methods)
# We try to find the marker row whose anchored position is closest.
if "snp_pos_ayb" in anc.columns and "chr_ayb" in anc.columns:
    near = anc[(anc["chr_ayb"] == "Ss10") &
                (anc["snp_pos_ayb"].between(15300000, 15500000))]
    if not near.empty:
        ALMT4_RS = near.iloc[
            (near["snp_pos_ayb"] - 15394673).abs().argsort()[:1]
        ]["rs"].iloc[0]
        print(f"[focus] ALMT4 SNP = {ALMT4_RS} at Ss10:"
               f"{near.loc[near['rs'] == ALMT4_RS, 'snp_pos_ayb'].iloc[0]:,}")


for trait_name, df in TRAITS:
    y = y_for(trait_name).values.astype(float)
    obs = ~np.isnan(y)
    if obs.sum() < 30:
        print(f"  {trait_name}: only {obs.sum()} samples — skip")
        continue
    Xt = X_work_scaled[obs]
    yt = y[obs]

    # best non-dummy model per trait by CV r
    row = cv[(cv["trait"] == trait_name) & (cv["model"] != "Dummy")]
    if row.empty:
        continue
    best_model_name = row.sort_values("r_mean", ascending=False).iloc[0]["model"]
    print(f"\n  {trait_name}: n={obs.sum()}; best model = {best_model_name}")

    try:
        m = fit_best(best_model_name, Xt, yt)
    except Exception as e:
        print(f"    fit failed: {e}")
        continue

    try:
        sv = shap_values_for(m, best_model_name, Xt)
    except Exception as e:
        print(f"    SHAP failed: {e}")
        continue
    if sv.ndim == 3:   # xgboost can return (n_outputs, n_samples, n_features)
        sv = sv[0]

    mean_abs = np.abs(sv).mean(axis=0)
    top_idx  = np.argsort(mean_abs)[::-1][:N_SHAP_TOP]
    for rank, mi in enumerate(top_idx, start=1):
        rs = markers[mi]
        top20_rows.append({
            "trait": trait_name, "rank": rank, "rs": rs,
            "mean_abs_shap": float(mean_abs[mi]),
            "best_model": best_model_name,
        })

    # GWAS top-30 per trait
    g_paths = [
        ROOT / "results" / "03_gwas_mlm" / "tables" / f"gwas_{trait_name}_M1_K.csv",
        ROOT / "results" / "15_gwas_new_traits" / "tables" / f"gwas_{trait_name}_M1_K.csv",
    ]
    gwas_top = set()
    for gp in g_paths:
        if gp.exists():
            gdf = pd.read_csv(gp).sort_values("p").head(GWAS_TOP_K)
            gwas_top = set(gdf["rs"].astype(str))
            break
    shap_top = {markers[i] for i in top_idx}
    inter = shap_top & gwas_top
    overlap_rows.append({
        "trait": trait_name,
        "best_model": best_model_name,
        "n_shap_top": len(shap_top),
        "n_gwas_top": len(gwas_top),
        "n_intersection": len(inter),
        "expected_by_chance": len(shap_top) * len(gwas_top) / len(markers),
    })

    # ALMT4 focus row
    if ALMT4_RS is not None:
        if ALMT4_RS in markers:
            ami = anc_idx[ALMT4_RS]
            shap_rank = int(np.where(np.argsort(mean_abs)[::-1] == ami)[0][0]) + 1
            gwas_rank = None
            if g_paths[0].exists():
                gdf = pd.read_csv(g_paths[0]).sort_values("p").reset_index(drop=True)
                hit = gdf[gdf["rs"] == ALMT4_RS]
                if not hit.empty:
                    gwas_rank = int(hit.index[0]) + 1
            elif g_paths[1].exists():
                gdf = pd.read_csv(g_paths[1]).sort_values("p").reset_index(drop=True)
                hit = gdf[gdf["rs"] == ALMT4_RS]
                if not hit.empty:
                    gwas_rank = int(hit.index[0]) + 1
            almt4_rows.append({
                "trait": trait_name, "best_model": best_model_name,
                "shap_rank":  shap_rank,
                "gwas_rank":  gwas_rank,
                "shap_mean_abs": float(mean_abs[ami]),
            })

    # --- Fig 89. Beeswarm summary ---
    try:
        fig = plt.figure(figsize=(7.5, 6))
        ev = shap.Explanation(
            values=sv,
            data=Xt,
            feature_names=[markers[i] for i in range(Xt.shape[1])],
        )
        shap.plots.beeswarm(ev, max_display=N_SHAP_TOP, show=False)
        fig.suptitle(f"SHAP beeswarm — {trait_name} (best: {best_model_name})",
                      fontsize=10)
        save(fig, f"fig89_shap_beeswarm_{trait_name}")
        print(f"    [fig] fig89_shap_beeswarm_{trait_name}")
    except Exception as e:
        print(f"    beeswarm failed: {e}")

    # --- Fig 90. Dependence plots for top-3 markers ---
    try:
        top3 = top_idx[:N_DEP_TOP]
        interact_idx = top_idx[1]   # colour by the 2nd most influential
        fig, axes = plt.subplots(1, N_DEP_TOP, figsize=(13, 4),
                                  constrained_layout=True)
        if N_DEP_TOP == 1:
            axes = [axes]
        for k, mi in enumerate(top3):
            ax = axes[k]
            colour_vals = Xt[:, interact_idx]
            sc = ax.scatter(Xt[:, mi], sv[:, mi],
                             c=colour_vals, cmap="viridis", s=14,
                             alpha=0.85, edgecolor="none")
            ax.axhline(0, color="grey", linewidth=0.4)
            ax.set_xlabel(f"dosage at {markers[mi]}", fontsize=8)
            ax.set_ylabel(f"SHAP value", fontsize=8)
            ax.set_title(f"top-{k+1}: {markers[mi]}", fontsize=9)
            cb = fig.colorbar(sc, ax=ax, fraction=0.06, pad=0.02)
            cb.set_label(f"interact: {markers[interact_idx]}",
                          fontsize=7)
        fig.suptitle(f"SHAP dependence — {trait_name} (best: {best_model_name})",
                      fontsize=10)
        save(fig, f"fig90_shap_dependence_{trait_name}")
        print(f"    [fig] fig90_shap_dependence_{trait_name}")
    except Exception as e:
        print(f"    dependence-plot failed: {e}")


# ---------------------------------------------------------------------------
# Aggregate tables
# ---------------------------------------------------------------------------
pd.DataFrame(top20_rows).to_csv(TAB / "shap_top20_per_trait.csv", index=False)
overlap = pd.DataFrame(overlap_rows)
overlap.to_csv(TAB / "shap_gwas_overlap.csv", index=False)
pd.DataFrame(almt4_rows).to_csv(TAB / "almt4_focus.csv", index=False)
print(f"\n[done] tables in {TAB}")


# ---------------------------------------------------------------------------
# Fig 91. SHAP-vs-GWAS overlap heatmap (one column, sorted by intersection)
# ---------------------------------------------------------------------------
if not overlap.empty:
    overlap_sorted = overlap.sort_values("n_intersection", ascending=False)
    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    xpos = np.arange(len(overlap_sorted))
    ax.bar(xpos, overlap_sorted["n_intersection"], width=0.6,
            color=WONG["vermillion"],
            label=f"shap-top-{N_SHAP_TOP} ∩ gwas-top-{GWAS_TOP_K}")
    ax.errorbar(xpos, overlap_sorted["expected_by_chance"],
                 yerr=0, fmt="_", color="grey",
                 markersize=14,
                 label="expected by chance")
    ax.set_xticks(xpos)
    ax.set_xticklabels(overlap_sorted["trait"].values,
                        rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("marker count")
    ax.set_title(f"SHAP top-{N_SHAP_TOP} ∩ GWAS top-{GWAS_TOP_K} per trait")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, axis="y", alpha=0.4)
    save(fig, "fig91_shap_gwas_overlap_heatmap")
    print("[fig] fig91_shap_gwas_overlap_heatmap")


# ---------------------------------------------------------------------------
# Fig 92. ALMT4 SNP across traits — SHAP rank vs GWAS rank
# ---------------------------------------------------------------------------
if almt4_rows:
    af = pd.DataFrame(almt4_rows).dropna()
    if not af.empty:
        from _figstyle import adjust_labels, publishable_axes
        fig, ax = plt.subplots(figsize=(8, 6.5), constrained_layout=True)
        ax.scatter(af["gwas_rank"], af["shap_rank"], s=70,
                    color=WONG["blue"], edgecolor="black", alpha=0.85,
                    zorder=4)
        # Collect trait labels and let adjustText push them apart so the
        # cluster near (rank 1, rank ~1000) stays readable.
        texts = []
        for _, row in af.iterrows():
            texts.append(ax.text(row["gwas_rank"], row["shap_rank"],
                                  row["trait"].replace("_", " "),
                                  fontsize=8, color="#222222", zorder=10))
        lim = max(af["gwas_rank"].max(), af["shap_rank"].max()) * 1.05
        ax.plot([0, lim], [0, lim], "k--", linewidth=0.6, alpha=0.5,
                 label="equal-rank line")
        ax.set_xscale("log"); ax.set_yscale("log")
        adjust_labels(texts, ax=ax,
                       expand=(1.3, 1.5), force_text=(0.9, 1.2),
                       only_move={"text": "xy", "static": "xy"})
        ax.set_xlabel(f"GWAS rank ({ALMT4_RS})")
        ax.set_ylabel(f"SHAP rank")
        ax.set_title(f"ALMT4 SNP (Ss10:15.4 Mb) -- cross-method evidence")
        ax.legend(fontsize=8.5, framealpha=0.9, edgecolor="none")
        publishable_axes(ax, grid="both", grid_alpha=0.25)
        save(fig, "fig92_almt4_focus")
        print("[fig] fig92_almt4_focus")

print(f"\nOutputs in: {OUT}")
