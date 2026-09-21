"""
83c -- Aggregate all RESULTS.txt under results/83_shiu_lab_classifier/runs/
into a single per-(trait, alg) AUC pivot table, then save figures + CSVs.

Runs after 83 + 83b. Reports the final CV AUC matrix and holdout AUC for
the per-trait winners.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
OUT = PROJ / "results/83_shiu_lab_classifier"


def parse(rf: Path) -> dict | None:
    if not rf.exists():
        return None
    lines = [l.strip() for l in open(rf) if l.strip()]
    if len(lines) < 2:
        return None
    header = lines[0].split("\t")
    last = lines[-1].split("\t")
    row = dict(zip(header, last))
    for k in ("AUCROC_val", "AUCROC_val_sd", "AUCPRc_val", "F1_val",
              "AUCROC_test", "AUCPRc_test", "F1_test", "FeatureNum",
              "BalancedSize", "CVfold", "BalancedRuns"):
        if k in row:
            try:
                row[k] = float(row[k])
            except Exception:
                pass
    return row


def main() -> None:
    rows = []
    for rdir in sorted((OUT / "runs").iterdir()):
        if not rdir.is_dir():
            continue
        name = rdir.name
        parsed = parse(rdir / "RESULTS.txt")
        if parsed is None:
            rows.append({"run": name, "status": "FAIL"})
            continue
        # split run name into trait + algo by the LAST underscore prefix
        # that matches one of the algo tags
        algo = None
        for a in ("RF", "GB", "LogReg", "SVMrbf"):
            if name.endswith("_" + a):
                algo = a
                trait = name[: -(len(a) + 1)]
                break
        if algo is None:
            continue
        parsed["trait"] = trait
        parsed["algo"] = algo
        parsed["status"] = "OK"
        rows.append(parsed)

    df = pd.DataFrame(rows)
    ok = df[df["status"] == "OK"].copy()
    print(f"OK rows: {len(ok)};  FAIL rows: {(df['status'] == 'FAIL').sum()}")

    auc_val = ok.pivot(index="trait", columns="algo", values="AUCROC_val")
    auc_val_sd = ok.pivot(index="trait", columns="algo", values="AUCROC_val_sd")
    auc_test = ok.pivot(index="trait", columns="algo", values="AUCROC_test")
    f1_val = ok.pivot(index="trait", columns="algo", values="F1_val")
    pr_val = ok.pivot(index="trait", columns="algo", values="AUCPRc_val")

    auc_val.to_csv(OUT / "tables/auc_val_pivot.csv")
    auc_val_sd.to_csv(OUT / "tables/auc_val_sd_pivot.csv")
    auc_test.to_csv(OUT / "tables/auc_test_pivot.csv")
    pr_val.to_csv(OUT / "tables/prc_val_pivot.csv")
    f1_val.to_csv(OUT / "tables/f1_val_pivot.csv")
    ok.to_csv(OUT / "tables/all_runs_long.csv", index=False)

    print("\n=== Shiu Lab CV ROC-AUC (val, mean across 20 balanced reps) ===")
    print(auc_val.round(3).to_string())
    print("\n=== Shiu Lab CV ROC-AUC sd ===")
    print(auc_val_sd.round(3).to_string())
    print("\n=== Holdout ROC-AUC (test) ===")
    print(auc_test.round(3).to_string())

    # Per-trait CV winner
    winners = ok.sort_values(["trait", "AUCROC_val"],
                              ascending=[True, False]).drop_duplicates("trait")
    winners.to_csv(OUT / "tables/cv_winners.csv", index=False)
    print("\n=== Per-trait CV winners ===")
    print(winners[["trait", "algo", "AUCROC_val", "AUCROC_val_sd",
                    "AUCPRc_val", "AUCROC_test"]].round(3).to_string(index=False))

    # ---- Figure ----
    traits = sorted(auc_val.index)
    algos = ["LogReg", "RF", "GB", "SVMrbf"]
    algos = [a for a in algos if a in auc_val.columns]
    mat_auc = auc_val.reindex(index=traits, columns=algos).values
    mat_test = auc_test.reindex(index=traits, columns=algos).values

    fig, axes = plt.subplots(1, 2, figsize=(13, max(5.2, 0.42 * len(traits))))

    cmap = LinearSegmentedColormap.from_list(
        "auc", ["#cccccc", "#dddddd", "#88CCEE", "#117733", "#332288"], N=256,
    )
    im0 = axes[0].imshow(mat_auc, cmap=cmap, aspect="auto", vmin=0.3, vmax=0.85)
    axes[0].set_xticks(range(len(algos))); axes[0].set_xticklabels(algos)
    axes[0].set_yticks(range(len(traits))); axes[0].set_yticklabels(traits, fontsize=9)
    axes[0].set_title("A. Shiu Lab CV ROC-AUC (val, mean across 20 balanced reps)")
    for i in range(mat_auc.shape[0]):
        for j in range(mat_auc.shape[1]):
            v = mat_auc[i, j]
            if not np.isnan(v):
                axes[0].text(j, i, f"{v:.2f}", ha="center", va="center",
                              fontsize=8, color="black" if v < 0.7 else "white")
    fig.colorbar(im0, ax=axes[0], fraction=0.04, pad=0.02)

    im1 = axes[1].imshow(mat_test, cmap=cmap, aspect="auto", vmin=0.3, vmax=1.0)
    axes[1].set_xticks(range(len(algos))); axes[1].set_xticklabels(algos)
    axes[1].set_yticks(range(len(traits))); axes[1].set_yticklabels([""] * len(traits))
    axes[1].set_title("B. Holdout ROC-AUC (test)")
    for i in range(mat_test.shape[0]):
        for j in range(mat_test.shape[1]):
            v = mat_test[i, j]
            if not np.isnan(v):
                axes[1].text(j, i, f"{v:.2f}", ha="center", va="center",
                              fontsize=8, color="black" if v < 0.8 else "white")
    fig.colorbar(im1, ax=axes[1], fraction=0.04, pad=0.02)

    fig.suptitle("Top-quintile binary classifier benchmark via Shiu Lab "
                 "ML-Pipeline (balanced downsampling x 20 reps)", fontsize=12,
                 y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "figures/fig_shiu_classifier_auc.png", dpi=300,
                bbox_inches="tight")
    fig.savefig(OUT / "figures/fig_shiu_classifier_auc.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print("\nWrote figures + tables.")


if __name__ == "__main__":
    main()
