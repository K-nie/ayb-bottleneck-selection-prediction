"""
41b -- targeted SHAP recompute + house-style beeswarm for B24 (Seed_Length).

The Seed_Length SHAP beeswarm is the one SHAP figure in Paper B that must be
recomputed rather than restyled from a saved table: script 41 never wrote the
per-sample SHAP value matrices, only the top-20 mean|SHAP| ranks. This script
recomputes SHAP for Seed_Length, freezes the matrix to
results/41_ml_shap_interpretability/data/, and renders the beeswarm in the
house (_figstyle) style so it matches the rest of the Paper B supplement.

Why Seed_Length:
  Seed_Length -- the success case. Best CV model RandomForest (r = 0.214); 13 of
                 its top-20 SHAP markers also fall in the top-30 GWAS hits, the
                 strongest SHAP-vs-GWAS concordance in the panel.

The companion Soluble_Oxalate x RR-BLUP beeswarm (former B25) was dropped from
the supplement: its best model is a linear Ridge, and script 41 attributed it
with KernelExplainer. With n = 37 samples and p = 1,625 markers KernelExplainer
is a stochastic coalition solve whose exact draw depends on the (unrecorded)
Jun-2022 shap/sklearn versions. Script 41 never saved the per-sample matrix, so
the Jun-2022 top-20 cannot be reproduced -- a fresh run in any stable environment
returns physically sane (~1-2.5) but different top markers. Rather than ship a
beeswarm whose markers contradict the frozen shap_top20_per_trait.csv that the
manuscript cites, the figure was removed. The success case (Seed_Length, below)
uses TreeExplainer, which is exact and deterministic and reproduces the frozen
table byte-for-byte.

Outputs:
  data/shap_matrix_Seed_Length.npz     sv (n x p), X (n x p), feature ids, samples
  figures/fig89_shap_beeswarm_Seed_Length.{png,pdf}   house-style beeswarm
"""
from __future__ import annotations
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, save_figure
apply()

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
OUT = ROOT / "results" / "41_ml_shap_interpretability"
DATA = OUT / "data"
FIG = OUT / "figures"
for d in (DATA, FIG):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
HOLDOUT_FRAC = 0.15
HOLDOUT_SEED = 42
N_SHAP_TOP = 20

# (trait, model) to recompute -- the Seed_Length success-case winner from
# script 40's CV. The Soluble_Oxalate x RR-BLUP case (former B25) was dropped;
# see the module docstring for why it is not reproducible.
TARGETS = [("Seed_Length", "RandomForest")]

# SHAP diverging colour map, Wong blue (low dosage) -> grey -> vermillion (high).
SHAP_CMAP = LinearSegmentedColormap.from_list(
    "wong_div", [WONG["blue"], "#d9d9d9", WONG["vermillion"]], N=256)


def load_working_set():
    """Reproduce script 41's dosage matrix + 15 % holdout working set."""
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
    X_work_scaled = StandardScaler(with_mean=True, with_std=False).fit_transform(
        X_work_imp)
    print(f"[load] working set {X_work_scaled.shape[0]} samples "
          f"x {X_work_scaled.shape[1]} markers")
    return X_work_scaled, markers, work_samples


def anchor_labels(markers) -> dict:
    """rs -> 'Ss10:15.4 Mb' when anchored, else a shortened rs."""
    lab = {}
    if ANCHOR.exists():
        anc = pd.read_csv(ANCHOR)
        anc = anc[anc["rs"].isin(markers)]
        for _, r in anc.iterrows():
            if pd.notna(r.get("chr_ayb")) and pd.notna(r.get("snp_pos_ayb")):
                lab[r["rs"]] = f"{r['chr_ayb']}:{r['snp_pos_ayb'] / 1e6:.1f} Mb"
    for m in markers:
        lab.setdefault(m, str(m).split("|")[0])
    return lab


def compute_shap(model_name, X, y):
    import shap
    # Seed_Length winner is a RandomForest -> TreeExplainer is exact and
    # deterministic, so it reproduces script 41's frozen top-20 exactly.
    m = RandomForestRegressor(n_estimators=500, max_features=0.1,
                              min_samples_leaf=2, random_state=0,
                              n_jobs=-1).fit(X, y)
    sv = np.array(shap.TreeExplainer(m).shap_values(X))
    if sv.ndim == 3:
        sv = sv[0]
    return sv


def beeswarm(sv, X, labels, trait, model_name, out_stub):
    """House-style SHAP beeswarm: top-N markers by mean|SHAP|, density-jittered
    points coloured by (standardised) dosage."""
    mean_abs = np.abs(sv).mean(axis=0)
    order = np.argsort(mean_abs)[::-1][:N_SHAP_TOP][::-1]  # bottom = most important
    n = len(order)
    fig, ax = plt.subplots(figsize=(7.6, 0.34 * n + 1.6))

    for row, mi in enumerate(order):
        s = sv[:, mi]
        fv = X[:, mi]
        # colour by dosage rank (robust to the mean-centred scaling)
        rank = pd.Series(fv).rank(pct=True).values if np.ptp(fv) > 0 \
            else np.full_like(fv, 0.5)
        # density-based vertical jitter (histogram of SHAP values -> spread)
        counts, edges = np.histogram(s, bins=max(6, len(s) // 4))
        bin_idx = np.clip(np.digitize(s, edges[1:-1]), 0, len(counts) - 1)
        jitter = np.zeros(len(s))
        for b in np.unique(bin_idx):
            sel = np.where(bin_idx == b)[0]
            k = len(sel)
            spread = 0.34 * min(1.0, k / 6.0)
            jitter[sel] = np.linspace(-spread, spread, k) if k > 1 else 0.0
        ax.scatter(s, np.full(len(s), row) + jitter, c=rank, cmap=SHAP_CMAP,
                   s=16, alpha=0.85, edgecolor="none", vmin=0, vmax=1,
                   rasterized=True, zorder=3)

    ax.axvline(0, color=WONG["grey"], lw=0.7, zorder=1)
    ax.set_yticks(range(n))
    ax.set_yticklabels([labels.get(_MARKERS[mi], str(_MARKERS[mi]))
                        for mi in order], fontsize=8)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("SHAP value (impact on predicted trait)")
    ax.set_title(f"SHAP beeswarm -- {trait.replace('_', ' ')} "
                 f"(best model: {model_name})", loc="left")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(True, axis="x", alpha=0.20)
    ax.set_axisbelow(True)

    sm = plt.cm.ScalarMappable(cmap=SHAP_CMAP, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cb.set_ticks([0, 1])
    cb.set_ticklabels(["low", "high"])
    cb.set_label("Marker dosage")
    fig.tight_layout()
    written = save_figure(fig, out_stub)
    plt.close(fig)
    return written


# marker list carried globally so the beeswarm can map column index -> rs id
_MARKERS: list = []


def main() -> None:
    global _MARKERS
    from _pheno import load_phenotypes
    X, markers, work_samples = load_working_set()
    _MARKERS = markers
    labels = anchor_labels(markers)
    pheno = load_phenotypes()
    pheno.index = pheno.index.astype(str).str.strip()

    for trait, model_name in TARGETS:
        y = pheno[trait].reindex(work_samples).values.astype(float)
        obs = ~np.isnan(y)
        Xt, yt = X[obs], y[obs]
        print(f"\n[{trait}] n={obs.sum()}; model={model_name}")
        sv = compute_shap(model_name, Xt, yt)

        mean_abs = np.abs(sv).mean(axis=0)
        top = np.argsort(mean_abs)[::-1][:N_SHAP_TOP]
        np.savez_compressed(
            DATA / f"shap_matrix_{trait}.npz",
            sv=sv, X=Xt, markers=np.array(markers, dtype=object),
            samples=np.array([s for s, o in zip(work_samples, obs) if o],
                             dtype=object),
            model=model_name)
        print(f"  saved SHAP matrix {sv.shape}; top marker "
              f"{labels.get(markers[top[0]], markers[top[0]])} "
              f"(mean|SHAP| {mean_abs[top[0]]:.4g})")

        written = beeswarm(sv, Xt, labels, trait, model_name,
                           FIG / f"fig89_shap_beeswarm_{trait}")
        print("  wrote", *written)


if __name__ == "__main__":
    main()
