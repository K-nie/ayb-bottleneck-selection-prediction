"""
83b -- Rerun only the (trait, algo) pairs that didn't produce a
valid RESULTS.txt in the first 83 sweep. In practice that's all
GB + SVMrbf runs (the first sweep crashed on np.str and on
GradientBoostingClassifier(loss='deviance') / SVMrbf kernel
UnboundLocalError, all of which are now patched in
/tmp/shiu_ml_pipeline/).

After this rerun, every (trait, algo) cell should have a RESULTS.txt
and we can aggregate the 11 traits x 4 algos = 44 cells into the
final AUC table.
"""
from __future__ import annotations
import shutil
from pathlib import Path
import sys
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
sys.path.insert(0, str(ROOT / "scripts"))
from importlib import util
spec = util.spec_from_file_location(
    "runner", ROOT / "scripts" / "83_shiu_lab_classifier_runner.py"
)
m = util.module_from_spec(spec); spec.loader.exec_module(m)  # type: ignore

OUT = m.OUT
RUNS = OUT / "runs"
ALGOS_TO_RERUN = ["GB", "SVMrbf"]


def needs_rerun(trait: str, alg: str) -> bool:
    rdir = RUNS / f"{trait}_{alg}"
    rf = rdir / "RESULTS.txt"
    if not rf.exists():
        return True
    try:
        lines = [l.strip() for l in open(rf) if l.strip()]
    except Exception:
        return True
    if len(lines) < 2:
        return True
    # Verify the data row parses
    header = lines[0].split("\t")
    last = lines[-1].split("\t")
    row = dict(zip(header, last))
    try:
        float(row.get("AUCROC_val", "nan"))
        return False
    except Exception:
        return True


def main() -> None:
    print("[load] dosage + pheno...")
    dosage = m.load_dosage()
    samples = dosage.columns.tolist()
    work, hold = m.split_work_hold(samples)
    pheno = m.load_phenotypes()
    pheno.index = pheno.index.astype(str).str.strip()

    print("[prep] re-using existing input TSVs if present; otherwise rebuild...")
    trait_info = m.build_trait_inputs(dosage, pheno, work, hold)
    print(f"  trait_info has {len(trait_info)} traits")

    to_run = []
    for trait in trait_info:
        for alg in ALGOS_TO_RERUN:
            if needs_rerun(trait, alg):
                to_run.append((trait, alg))
    print(f"\n[plan] {len(to_run)} (trait, algo) cells to rerun:")
    for t, a in to_run:
        print(f"  - {t} x {a}")

    rows = []
    for trait, alg in to_run:
        # Clean the run dir to avoid stale RESULTS.txt entries
        rdir = RUNS / f"{trait}_{alg}"
        if rdir.exists():
            shutil.rmtree(rdir)
        rdir.mkdir(parents=True, exist_ok=True)
        try:
            rc, prefix = m.run_shiu(trait, alg, trait_info[trait])
        except Exception as e:
            print(f"    EXC: {trait} x {alg}: {e}")
            continue
        row = m.parse_results(rdir)
        if row is None:
            print(f"    FAIL: {trait} x {alg}; rc = {rc}")
            rows.append({
                "trait": trait, "alg": alg, "status": "fail",
                "returncode": rc,
            })
            continue
        for k in ("AUCROC_val", "AUCROC_val_sd", "AUCPRc_val", "F1_val",
                  "AUCROC_test", "AUCPRc_test", "F1_test"):
            if k in row:
                try:
                    row[k] = float(row[k])
                except Exception:
                    pass
        row["trait"] = trait
        row["alg"] = alg
        row["status"] = "ok"
        rows.append(row)
        print(f"    OK: {trait} x {alg}  "
              f"AUCROC_val = {row.get('AUCROC_val', 'NA')}  "
              f"AUCROC_test = {row.get('AUCROC_test', 'NA')}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tables/rerun_gb_svmrbf_results.csv", index=False)
    print("\nWrote", OUT / "tables/rerun_gb_svmrbf_results.csv")


if __name__ == "__main__":
    main()
