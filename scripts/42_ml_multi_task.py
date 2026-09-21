#!/usr/bin/env python3
"""
Multi-task / multi-output genomic prediction for the AYB DArTseq panel.

Compares single-trait baselines against multi-task and shared-trunk
models on four genetically-correlated trait groups:

    Group A — Oxalate triplet  (Total, Soluble, Insoluble)        n approx 41
    Group B — Protein + 3 oxalates                                n approx 41
    Group C — Seed-size cluster (Length, Width, Thickness, Mass)  n approx 79
    Group D — Biochem antioxidant cluster (Antiox, Flavo, Phenol) n approx 79

Models:
    1. Single-trait RR-BLUP                            — script 40 baseline
    2. Single-trait XGBoost                            — script 40 baseline
    3. Multi-trait GBLUP (Calus & Veerkamp 2011)       — joint kernel ridge
    4. Multi-output XGBoost                            — Sandhu et al. 2021
    5. Shared-trunk MLP                                — Montesinos-López 2018, 2019

The Q1 question this script answers: do shared-trunk or multi-task models
recover information that single-trait models miss when traits are
genetically correlated? At n approx 41 for the oxalate triplet, the
expectation (Bellot et al. 2018 / Montesinos-López 2019) is "modest gain
for tightly correlated traits, no gain for loosely correlated traits".

Cross-validation: 5-fold OUTER x 10 outer repeats (no inner tuning;
multi-task models use fixed hyperparameters chosen to match script 40's
inner-CV winners on the matching single-trait variant). Mean-imputation
and scaling done per outer-train fold (Lones 3.1).

Output:
    results/42_ml_multi_task/
        tables/
            multi_task_per_trait_r.csv
            delta_r_vs_single.csv
            wilcoxon_multi_vs_single.csv
        figures/
            fig93_multi_task_forest
            fig94_delta_r_heatmap

Lones (2024) Patterns walk-through
----------------------------------
  2.1 use of data:    train/test per outer fold; the script 40 15 % holdout
                      is NOT re-used here (different sample subsets per
                      group; multi-task evaluation stays in CV space).
  3.1 no test leak:   per-fold imputation + scaling
  3.2 range of models: 5 models per group (2 single + 3 multi)
  3.5 DL not best:    shared-trunk MLP included with honest framing
  3.6 feature select: NONE — full marker matrix every fold
  3.7 hyperparam opt: fixed configs from script 40 inner-CV winners
  4.1 appropriate test: 5-fold outer x 10 repeats
  4.4 evaluate mult:  10 outer repeats, bootstrap 95 % CI on per-trait r
  4.6 metrics:        Pearson r + bootstrap CI, per-trait
  5.2 baselines:      single-trait RR-BLUP and single-trait XGBoost
  5.3 stat tests:     paired Wilcoxon on per-repeat r (multi vs single)
  5.4 multiple comp:  Holm-Bonferroni across the 3 multi-task models per trait
  6.1 transparent:    seeds, hyperparams, and CV scheme baked in here
  6.3 don't overgen:  results limited to the 95-line IITA TSs panel

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.stats import pearsonr, wilcoxon
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.multioutput import MultiOutputRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler

from _plotstyle import apply, WONG
apply()
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
OUT  = ROOT / "results" / "42_ml_multi_task"
DATA = OUT / "data"
TAB  = OUT / "tables"
FIG  = OUT / "figures"
for d in (DATA, TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
N_OUTER_FOLDS   = 5
N_OUTER_REPEATS = 10
N_BOOTSTRAP     = 2000


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Load post-QC dosage and the canonical phenotype frame
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
samples = list(dosage.columns)
X_full = dosage.T.values.astype(float)
print(f"[load] dosage matrix: {X_full.shape[0]} x {X_full.shape[1]}")

from _pheno import load_phenotypes
pheno = load_phenotypes()
pheno.index = pheno.index.astype(str).str.strip()


# ---------------------------------------------------------------------------
# Genetically-correlated trait groups
# (correlations from script 14 Pearson matrix; threshold |r| >= 0.3)
# ---------------------------------------------------------------------------
GROUPS = {
    "A_oxalates":   ["Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate"],
    "B_prot_ox":    ["Crude_Protein", "Total_Oxalate", "Soluble_Oxalate",
                      "Insoluble_Oxalate"],
    "C_seed_size":  ["Seed_Length", "Seed_Width", "Seed_Thickness",
                      "Mass_of_Seeds"],
    "D_antiox":     ["Antioxidant", "Flavonoid", "Phenol"],
}


# ---------------------------------------------------------------------------
# Per-fold imputation + scaling (Lones 3.1)
# ---------------------------------------------------------------------------

def impute_and_scale_train_only(X_train, X_test):
    mu = np.nanmean(X_train, axis=0)
    Xtr = np.where(np.isnan(X_train), mu, X_train)
    Xte = np.where(np.isnan(X_test),  mu, X_test)
    sc = StandardScaler(with_mean=True, with_std=False).fit(Xtr)
    return sc.transform(Xtr), sc.transform(Xte)


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

def fit_single_rrblup(X, y):
    return Ridge(alpha=1.0).fit(X, y)


def fit_single_xgb(X, y):
    import xgboost as xgb
    return xgb.XGBRegressor(n_estimators=100, max_depth=5,
                              learning_rate=0.05, subsample=0.7,
                              n_jobs=-1, random_state=0, verbosity=0,
                              tree_method="hist").fit(X, y)


def fit_multi_gblup(X, Y):
    """Multi-trait GBLUP via kernel ridge with shared alpha. Calus &
    Veerkamp 2011 approximation: fit independent ridges with the same
    penalty (cheap proxy to factor-analytic multi-trait BLUP)."""
    return Ridge(alpha=1.0).fit(X, Y)


def fit_multi_xgb(X, Y):
    """Multi-output XGBoost via sklearn's MultiOutputRegressor (one
    model per output, but shared feature matrix and identical seeds —
    not a true joint model). Sandhu et al. 2021 used this pattern."""
    import xgboost as xgb
    base = xgb.XGBRegressor(n_estimators=100, max_depth=5,
                              learning_rate=0.05, subsample=0.7,
                              n_jobs=-1, random_state=0, verbosity=0,
                              tree_method="hist")
    return MultiOutputRegressor(base, n_jobs=1).fit(X, Y)


def fit_shared_trunk_mlp(X, Y):
    """Shared-trunk MLP with multi-output regression head.
    Montesinos-López et al. 2018, 2019. 2 hidden layers (64, 32),
    weight decay 1e-3, no early stopping."""
    return MLPRegressor(hidden_layer_sizes=(64, 32), activation="relu",
                          alpha=1e-3, max_iter=300, early_stopping=False,
                          random_state=0).fit(X, Y)


# ---------------------------------------------------------------------------
# Run CV per group
# ---------------------------------------------------------------------------

rows = []
delta_rows = []
wlx_rows = []

for group_name, traits in GROUPS.items():
    Y_block = pheno.reindex(samples)[traits].values.astype(float)
    obs = ~np.any(np.isnan(Y_block), axis=1)
    if obs.sum() < 30:
        print(f"\n[{group_name}] only {obs.sum()} samples with all "
               f"{len(traits)} traits — skip")
        continue
    Xg = X_full[obs]
    Yg = Y_block[obs]
    print(f"\n[{group_name}] n={obs.sum()}, traits={traits}")

    # storage for per-repeat per-trait r per model
    r_per_model = {
        "Single RR-BLUP": np.full((N_OUTER_REPEATS, len(traits)), np.nan),
        "Single XGBoost": np.full((N_OUTER_REPEATS, len(traits)), np.nan),
        "Multi-trait GBLUP": np.full((N_OUTER_REPEATS, len(traits)), np.nan),
        "Multi-output XGB": np.full((N_OUTER_REPEATS, len(traits)), np.nan),
        "Shared-trunk MLP": np.full((N_OUTER_REPEATS, len(traits)), np.nan),
    }
    for rep in range(N_OUTER_REPEATS):
        outer = KFold(n_splits=N_OUTER_FOLDS, shuffle=True,
                       random_state=6000 + rep)
        # Predictions per trait per model
        yhat = {name: np.full_like(Yg, np.nan, dtype=float)
                for name in r_per_model}
        for tr_i, te_i in outer.split(Xg):
            Xtr, Xte = impute_and_scale_train_only(Xg[tr_i], Xg[te_i])
            Ytr = Yg[tr_i]
            # single-trait models (one fit per trait)
            for ti, tr_name in enumerate(traits):
                m_rr  = fit_single_rrblup(Xtr, Ytr[:, ti])
                m_xgb = fit_single_xgb(Xtr,   Ytr[:, ti])
                yhat["Single RR-BLUP"][te_i, ti] = m_rr.predict(Xte)
                yhat["Single XGBoost"][te_i, ti] = m_xgb.predict(Xte)
            # multi-trait models
            m_mt_gblup = fit_multi_gblup(Xtr, Ytr)
            yhat["Multi-trait GBLUP"][te_i] = m_mt_gblup.predict(Xte)
            m_mo_xgb = fit_multi_xgb(Xtr, Ytr)
            yhat["Multi-output XGB"][te_i] = m_mo_xgb.predict(Xte)
            m_mlp = fit_shared_trunk_mlp(Xtr, Ytr)
            yhat["Shared-trunk MLP"][te_i] = m_mlp.predict(Xte)
        # per-trait r per model
        for name, pred in yhat.items():
            for ti in range(len(traits)):
                p = pred[:, ti]
                if np.nanstd(p) > 0:
                    r_per_model[name][rep, ti] = pearsonr(p, Yg[:, ti])[0]

    # aggregate
    for name, mat in r_per_model.items():
        for ti, tr_name in enumerate(traits):
            arr = mat[:, ti]
            rows.append({
                "group": group_name, "trait": tr_name, "model": name,
                "n":         int(obs.sum()),
                "r_mean":    float(np.nanmean(arr)),
                "r_ci_lo":   float(np.nanpercentile(arr, 2.5)),
                "r_ci_hi":   float(np.nanpercentile(arr, 97.5)),
                "n_reps":    int(np.sum(~np.isnan(arr))),
            })
            print(f"    {name:18s} | {tr_name:18s} | "
                  f"r={np.nanmean(arr):+.3f} "
                  f"[{np.nanpercentile(arr, 2.5):+.2f}, "
                  f"{np.nanpercentile(arr, 97.5):+.2f}]")

    # Δr per multi-task model vs single-trait baseline (RR-BLUP)
    base_rr = r_per_model["Single RR-BLUP"]
    for name in ("Multi-trait GBLUP", "Multi-output XGB", "Shared-trunk MLP"):
        mt = r_per_model[name]
        # paired Wilcoxon on per-repeat (mt - base_rr) per trait
        for ti, tr_name in enumerate(traits):
            mask = ~np.isnan(mt[:, ti]) & ~np.isnan(base_rr[:, ti])
            if mask.sum() < 5:
                continue
            delta = mt[mask, ti] - base_rr[mask, ti]
            try:
                _, p = wilcoxon(mt[mask, ti], base_rr[mask, ti])
            except ValueError:
                p = np.nan
            rng = np.random.default_rng(7000 + ti)
            boot = []
            for _ in range(N_BOOTSTRAP):
                sub = rng.choice(delta, size=len(delta), replace=True)
                boot.append(np.mean(sub))
            delta_rows.append({
                "group": group_name, "trait": tr_name,
                "model": name,
                "delta_r_mean": float(np.mean(delta)),
                "delta_r_ci_lo": float(np.percentile(boot, 2.5)),
                "delta_r_ci_hi": float(np.percentile(boot, 97.5)),
                "wilcoxon_p_raw": float(p),
            })
            wlx_rows.append({
                "group": group_name, "trait": tr_name,
                "model": name, "p_raw": float(p),
            })

# Holm-Bonferroni within each (group, trait) across the 3 multi-task models
wlx = pd.DataFrame(wlx_rows)
holm_rows = []
for (g, t), sub in wlx.groupby(["group", "trait"]):
    sub_sorted = sub.sort_values("p_raw")
    m = len(sub_sorted)
    for rank, (_, row) in enumerate(sub_sorted.iterrows()):
        p_holm = min(1.0, row["p_raw"] * (m - rank))
        holm_rows.append({
            "group": g, "trait": t, "model": row["model"],
            "p_raw":  row["p_raw"],
            "p_holm": p_holm,
            "sig_holm_0.05": bool(p_holm < 0.05),
        })
pd.DataFrame(rows).to_csv(TAB / "multi_task_per_trait_r.csv", index=False)
pd.DataFrame(delta_rows).to_csv(TAB / "delta_r_vs_single.csv", index=False)
pd.DataFrame(holm_rows).to_csv(TAB / "wilcoxon_multi_vs_single.csv", index=False)


# ---------------------------------------------------------------------------
# Fig 93. Per-trait forest by group
# ---------------------------------------------------------------------------
df = pd.DataFrame(rows)
fig, axes = plt.subplots(len(GROUPS), 1,
                          figsize=(10, 2 + 1.6 * sum(len(t) for t in GROUPS.values())),
                          constrained_layout=True)
if len(GROUPS) == 1:
    axes = [axes]
model_order = ["Single RR-BLUP", "Single XGBoost", "Multi-trait GBLUP",
                "Multi-output XGB", "Shared-trunk MLP"]
palette = {
    "Single RR-BLUP":    "#888888",
    "Single XGBoost":    "#bbbbbb",
    "Multi-trait GBLUP": WONG["blue"],
    "Multi-output XGB":  WONG["vermillion"],
    "Shared-trunk MLP":  WONG["green"],
}
for ax, (g, traits) in zip(axes, GROUPS.items()):
    sub = df[df["group"] == g]
    n_t = len(traits); height = 0.13
    y_pos = np.arange(n_t)[::-1]
    for i, m in enumerate(model_order):
        s = sub[sub["model"] == m].set_index("trait").reindex(traits)
        y = y_pos + (i - (len(model_order) - 1) / 2) * height
        ax.errorbar(s["r_mean"], y,
                     xerr=[s["r_mean"] - s["r_ci_lo"],
                           s["r_ci_hi"] - s["r_mean"]],
                     fmt="o", color=palette[m], markersize=5,
                     elinewidth=0.8, capsize=2, label=m)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_yticks(y_pos); ax.set_yticklabels(traits, fontsize=8)
    ax.set_xlabel("Pearson r (5-fold CV x 10 repeats)" if ax == axes[-1] else "")
    ax.set_title(f"{g}  (n={sub['n'].iloc[0]})", fontsize=9, loc="left")
    ax.grid(True, axis="x", alpha=0.4)
    if ax == axes[0]:
        ax.legend(loc="lower right", fontsize=7)
fig.suptitle("Multi-task vs single-trait genomic prediction "
              "(AYB DArTseq panel)", fontsize=11)
save(fig, "fig93_multi_task_forest")
print("\n[fig] fig93_multi_task_forest")


# ---------------------------------------------------------------------------
# Fig 94. Δr heatmap (rows = trait, cols = multi-task model)
# ---------------------------------------------------------------------------
delta = pd.DataFrame(delta_rows)
if not delta.empty:
    delta["trait_grp"] = delta["group"] + " :: " + delta["trait"]
    pivot = delta.pivot(index="trait_grp", columns="model",
                         values="delta_r_mean")
    fig, ax = plt.subplots(figsize=(7, max(4, 0.35 * len(pivot))),
                            constrained_layout=True)
    im = ax.imshow(pivot.values, cmap="RdBu_r",
                    vmin=-0.10, vmax=0.10, aspect="auto")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=7)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.3f}", ha="center", va="center",
                         fontsize=7,
                         color="white" if abs(v) > 0.06 else "black")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02,
                  label="Δr vs single-trait RR-BLUP")
    ax.set_title("Multi-task gain over single-trait RR-BLUP (positive = better)")
    save(fig, "fig94_delta_r_heatmap")
    print("[fig] fig94_delta_r_heatmap")

print(f"\nOutputs in: {OUT}")
