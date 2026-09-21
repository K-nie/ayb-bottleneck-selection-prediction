"""Smoke test: prep inputs only + run ONE (trait, alg) to validate the
Shiu Lab pipeline end-to-end on our data before launching all 52 runs."""
import sys
from pathlib import Path
ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
sys.path.insert(0, str(ROOT / "scripts"))
# Import functions from the runner
from importlib import util
spec = util.spec_from_file_location(
    "runner", ROOT / "scripts" / "83_shiu_lab_classifier_runner.py"
)
m = util.module_from_spec(spec); spec.loader.exec_module(m)  # type: ignore

print("[load] dosage...")
dosage = m.load_dosage()
samples = dosage.columns.tolist()
work, hold = m.split_work_hold(samples)
print(f"[split] work {len(work)}, hold {len(hold)}")
pheno = m.load_phenotypes()
pheno.index = pheno.index.astype(str).str.strip()
print("[prep] inputs...")
ti = m.build_trait_inputs(dosage, pheno, work, hold)
print(f"  prepared {len(ti)} traits")

# Run only Soluble_Oxalate x RF as smoke test
target = "Soluble_Oxalate"
if target in ti:
    rc, prefix = m.run_shiu(target, "RF", ti[target])
    print(f"  smoke run rc = {rc}")
    rdir = m.OUT / "runs" / f"{target}_RF"
    row = m.parse_results(rdir)
    print("  RESULTS.txt row:", row)
else:
    print(f"  {target} not in trait_info; available: {list(ti.keys())}")
