"""
83 -- Shiu Lab ML-Pipeline classifier wrapper.

Re-runs the per-trait top-quintile binary classifier benchmark from
script 82, but using the Shiu Lab pipeline
(github.com/ShiuLab/ML-Pipeline). Shiu Lab handles class imbalance
natively via balanced downsampling x N replicates (default n = 100),
which is the gold-standard plant-genomics approach -- avoids the
per-model class_weight / scale_pos_weight / oversample dispatch table
we wrote ourselves in script 82.

What this wrapper does
----------------------
1. Build a 95-sample (ID, Class, dosage_features) TSV per trait.
   - Class = 1 if top-quintile by BLUP, 0 otherwise.
   - The top-quintile threshold is computed on the WORKING SET (n=81)
     to avoid letting the 14 holdout samples' BLUP values define their
     own positive boundary (data-leakage discipline).
2. Write a test_instances.txt holding the SAME 14 holdout sample IDs
   used by script 40 / 41 / 82 (HOLDOUT_SEED = 42).
3. For each (trait, algorithm) in {RF, GB, LogReg, SVMrbf}, invoke
   ML_classification.py with:
     -cl_train 1,0     (binary, top-quintile vs rest)
     -pos 1            (top-quintile is the positive class)
     -test test_instances.txt   (the 14 holdout samples)
     -n  20            (20 balanced down-sampled replicates per fit)
     -cv_num 5         (5-fold CV, matches script 40 / 82 convention)
     -gs f             (no grid search; use sensible defaults for speed)
     -x_norm f         (DArTseq dosage is already {0, 1, 2}, no need)
     -threshold_test F1
   Each run appends one row to RESULTS.txt in the run directory.
4. Aggregate RESULTS.txt into a per-(trait, alg) AUC table + a 3-panel
   figure that's directly comparable to script 82's output.

Honest notes
------------
- Shiu Lab balanced downsampling discards majority-class samples per
  replicate (n=20 replicates -> 20 different downsample draws).
  Aggregated AUC is the mean across the 20 reps. This is a stronger
  imbalance discipline than the class_weight / scale_pos_weight knobs
  used in script 82.
- All preprocessing (NaN imputation, scaling) is applied to the FULL
  input TSV before any train/test split, but the holdout split is
  defined first and Shiu Lab's CV loop only sees working samples in
  training. The dosage matrix has NO per-fold imputation here; we
  pre-impute the full matrix using column means computed on WORKING
  SAMPLES ONLY, which preserves the leakage discipline from script 82.
- Grid search OFF for compute (each grid-on run is ~10x slower); the
  default RF / GB / SVMrbf / LogReg hyperparameters are robust enough
  for the small-n single-pass comparison.
"""
from __future__ import annotations
from pathlib import Path
import subprocess
import shutil
import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
SHIU = Path("/tmp/shiu_ml_pipeline")
OUT = ROOT / "results/83_shiu_lab_classifier"
(OUT / "data").mkdir(parents=True, exist_ok=True)
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
(OUT / "runs").mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
HOLDOUT_FRAC = 0.15
HOLDOUT_SEED = 42
TOP_QUINTILE = 0.20
N_BALANCED = 20
CV_NUM = 5
ALGOS = ["RF", "GB", "LogReg", "SVMrbf"]

import sys
sys.path.insert(0, str(ROOT / "scripts"))
from _pheno import load_phenotypes, TRAITS_ALL  # type: ignore


def load_dosage():
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
    return dosage


def split_work_hold(samples: list[str]) -> tuple[list[str], list[str]]:
    rng = np.random.default_rng(HOLDOUT_SEED)
    idx = np.arange(len(samples))
    rng.shuffle(idx)
    n_hold = int(round(len(samples) * HOLDOUT_FRAC))
    hold = [samples[i] for i in idx[:n_hold]]
    work = [samples[i] for i in idx[n_hold:]]
    return work, hold


def build_trait_inputs(dosage: pd.DataFrame, pheno: pd.DataFrame,
                       work: list[str], hold: list[str]) -> dict[str, dict]:
    """Returns trait -> dict(input_tsv_path, test_ids_path, n_pos_work,
    n_pos_hold, threshold)."""
    samples = dosage.columns.tolist()
    X = dosage.T  # samples x markers
    # column-mean imputation using WORKING samples only (no leakage)
    col_mu = X.loc[work].mean(axis=0, skipna=True)
    X_imp = X.fillna(col_mu)

    out = {}
    for trait in TRAITS_ALL:
        if trait not in pheno.columns:
            continue
        y_work = pheno[trait].reindex(work)
        y_hold = pheno[trait].reindex(hold)
        n_obs_work = y_work.notna().sum()
        n_obs_hold = y_hold.notna().sum()
        if n_obs_work < 20 or n_obs_hold < 4:
            print(f"  {trait}: n_work_obs={n_obs_work}, n_hold_obs={n_obs_hold}; skip")
            continue
        thr = float(y_work.dropna().quantile(1.0 - TOP_QUINTILE))
        # binary labels for all samples
        label = pheno[trait].reindex(samples).apply(
            lambda v: 1 if (pd.notna(v) and v > thr) else
                      (0 if pd.notna(v) else np.nan)
        )
        # drop samples with no phenotype (they cannot get a label)
        present = label.dropna().index.tolist()
        n_pos_work_eff = int(((pheno[trait].reindex(work) > thr).fillna(False)).sum())
        n_pos_hold_eff = int(((pheno[trait].reindex(hold) > thr).fillna(False)).sum())
        if n_pos_work_eff < 4 or n_pos_hold_eff == 0:
            print(f"  {trait}: pos_work={n_pos_work_eff}, "
                  f"pos_hold={n_pos_hold_eff}; skip")
            continue
        # Build the Shiu-Lab TSV: ID, Class, f1, f2, ...
        feats = X_imp.loc[present]
        df = pd.DataFrame(
            feats.values, index=feats.index,
            columns=[f"f{i}" for i in range(feats.shape[1])],
        )
        df.insert(0, "Class", label.loc[present].astype(int).values)
        df.index.name = "ID"
        # Shiu Lab expects header with ID column (default y_name = "Class")
        tsv_path = OUT / "data" / f"{trait}_input.tsv"
        df.to_csv(tsv_path, sep="\t", index=True)
        # Test instance IDs: the subset of holdout samples that have a label
        held_with_label = [s for s in hold if s in present]
        test_ids_path = OUT / "data" / f"{trait}_test_ids.txt"
        with open(test_ids_path, "w") as f:
            for s in held_with_label:
                f.write(f"{s}\n")
        out[trait] = {
            "tsv": tsv_path, "test_ids": test_ids_path,
            "threshold": thr,
            "n_present": int(len(present)),
            "n_test": int(len(held_with_label)),
            "n_pos_work": int(n_pos_work_eff),
            "n_pos_hold": int(n_pos_hold_eff),
        }
        print(f"  {trait}: n_present={len(present)}, "
              f"n_test={len(held_with_label)}, "
              f"pos_work={n_pos_work_eff}, pos_hold={n_pos_hold_eff}")
    return out


def run_shiu(trait: str, alg: str, info: dict) -> tuple[int, str]:
    """Invoke ML_classification.py for one (trait, alg). Returns
    (returncode, save_prefix)."""
    save_prefix = f"{trait}_{alg}"
    cwd = OUT / "runs" / f"{trait}_{alg}"
    cwd.mkdir(parents=True, exist_ok=True)
    # Copy input + test ids to the run dir (so RESULTS.txt is per-(trait,alg))
    tsv_local = cwd / info["tsv"].name
    test_local = cwd / info["test_ids"].name
    shutil.copy(info["tsv"], tsv_local)
    shutil.copy(info["test_ids"], test_local)
    cmd = [
        "/Users/black_einstein/miniconda3/envs/yeast-viz/bin/python",
        str(SHIU / "ML_classification.py"),
        "-df", str(tsv_local.name),
        "-test", str(test_local.name),
        "-alg", alg,
        "-cl_train", "1,0",
        "-pos", "1",
        "-n", str(N_BALANCED),
        "-cv_num", str(CV_NUM),
        "-gs", "f",
        "-x_norm", "f",
        "-drop_na", "f",
        "-threshold_test", "F1",
        "-save", save_prefix,
        "-tag", trait,
    ]
    print(f"  RUN {trait} x {alg}: cwd = {cwd}")
    proc = subprocess.run(
        cmd, cwd=cwd, capture_output=True, text=True, timeout=3600,
    )
    log = cwd / f"{save_prefix}.log"
    with open(log, "w") as f:
        f.write("STDOUT\n======\n" + proc.stdout
                + "\n\nSTDERR\n======\n" + proc.stderr)
    return proc.returncode, save_prefix


def parse_results(run_dir: Path) -> dict | None:
    """Read RESULTS.txt from one run dir; return the LAST data row as a
    dict, or None if no rows."""
    rf = run_dir / "RESULTS.txt"
    if not rf.exists():
        return None
    with open(rf) as f:
        lines = [l.strip() for l in f if l.strip()]
    if len(lines) < 2:
        return None
    header = lines[0].split("\t")
    last = lines[-1].split("\t")
    return dict(zip(header, last))


def main() -> None:
    print("[load] dosage matrix...")
    dosage = load_dosage()
    samples = dosage.columns.tolist()
    print(f"  {dosage.shape}: {len(samples)} samples x {dosage.shape[0]} markers")

    work, hold = split_work_hold(samples)
    print(f"[split] work n = {len(work)}; hold n = {len(hold)}")

    pheno = load_phenotypes()
    pheno.index = pheno.index.astype(str).str.strip()
    print(f"[pheno] {pheno.shape[0]} sample rows, "
          f"{pheno.shape[1]} traits")

    print("\n[prep] writing per-trait input TSVs...")
    trait_info = build_trait_inputs(dosage, pheno, work, hold)
    print(f"  prepared {len(trait_info)} traits")

    print(f"\n[run] {len(trait_info)} traits x {len(ALGOS)} algos "
          f"= {len(trait_info) * len(ALGOS)} runs...")
    summary_rows = []
    for trait, info in trait_info.items():
        for alg in ALGOS:
            try:
                rc, prefix = run_shiu(trait, alg, info)
            except subprocess.TimeoutExpired:
                print(f"    TIMEOUT: {trait} x {alg}")
                continue
            run_dir = OUT / "runs" / f"{trait}_{alg}"
            row = parse_results(run_dir)
            if row is None:
                print(f"    FAIL: {trait} x {alg}; returncode = {rc}; "
                      f"no RESULTS.txt rows")
                summary_rows.append({
                    "trait": trait, "alg": alg, "status": "fail",
                    "returncode": rc,
                })
                continue
            # cast numeric columns
            for k in ("AUCROC_val", "AUCROC_val_sd", "AUCPRc_val",
                      "F1_val", "AUCROC_test", "AUCPRc_test", "F1_test"):
                if k in row:
                    try:
                        row[k] = float(row[k])
                    except Exception:
                        pass
            row["trait"] = trait
            row["alg"] = alg
            row["status"] = "ok"
            summary_rows.append(row)
            print(f"    OK: {trait} x {alg}  "
                  f"AUCROC_val = {row.get('AUCROC_val', 'NA')}  "
                  f"AUCROC_test = {row.get('AUCROC_test', 'NA')}")

    df = pd.DataFrame(summary_rows)
    df.to_csv(OUT / "tables/per_trait_per_alg_results.csv", index=False)

    ok = df[df["status"] == "ok"].copy()
    pivot_auc = ok.pivot(index="trait", columns="alg",
                          values="AUCROC_val") if not ok.empty else None
    if pivot_auc is not None:
        pivot_auc.to_csv(OUT / "tables/auc_val_pivot.csv")
        print("\n=== Shiu Lab CV ROC-AUC val (mean across balanced reps) ===")
        print(pivot_auc.round(3).to_string())


if __name__ == "__main__":
    main()
