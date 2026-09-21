# B24 SHAP beeswarm recompute (Seed_Length) — and why B25 was dropped

## What this documents and why

Paper B carries one per-trait SHAP beeswarm supplementary figure:

- **B24 — Seed_Length** (success case): the trait whose top-20 SHAP markers most
  strongly overlap the GWAS top-30 hits (13/20), and whose best cross-validated
  genomic-prediction model is a RandomForest.

A second per-trait beeswarm was planned — **B25 — Soluble_Oxalate** (the linear
RR-BLUP contrast case) — but it was **dropped** from the supplement. The reason
is documented below; it is a reproducibility limitation, not a cosmetic choice.

The parent script `scripts/41_ml_shap_interpretability.py` renders these
beeswarms but **does not save the per-sample SHAP value matrices** — only the
top-20 mean|SHAP| ranks (`tables/shap_top20_per_trait.csv`). The beeswarms
therefore cannot be restyled from a frozen table like the other Paper B supp
figures. `scripts/41b_shap_beeswarm_two_traits.py` recomputes SHAP for
Seed_Length, freezes the matrix here, and renders the house-style
(`_figstyle`) beeswarm.

## Inputs (provenance)

- Genotypes: `AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv`
  (DArTseq HapMap, biallelic dosage encoded ref/ref=0, het=1, alt/alt=2).
- Marker anchoring: `refs/ayb_genome/ayb_marker_anchoring.csv`
  (rs -> `chr_ayb`, `snp_pos_ayb`; used only for beeswarm y-axis labels).
- Phenotypes: canonical `_pheno.load_phenotypes()` loader.
- Best-model-per-trait: `results/40_ml_genomic_prediction/tables/cv_predictive_ability.csv`
  (Seed_Length winner = RandomForest, r_mean = 0.214).

## Reproduce

```
conda activate yeast-viz
python scripts/41b_shap_beeswarm_two_traits.py
```

## Parameters / thresholds (identical to script 41, no re-tuning)

- Marker call rate >= 0.90; sample call rate >= 0.90; MAF >= 0.05
  -> 81 working samples x 1,625 markers after filtering.
- 15 % final holdout removed with `HOLDOUT_SEED = 42` (same shuffle as scripts
  40 and 41), so SHAP is fit on the **working set only** — no test leak
  (Lones 2024 item 4.5).
- Dosage mean-centred (`StandardScaler(with_mean=True, with_std=False)`), NaNs
  imputed to the working-set marker mean, exactly as script 41.
- RandomForest: 500 trees, `max_features=0.1`, `min_samples_leaf=2`,
  `random_state=0`.
- Beeswarm shows top-20 markers by mean|SHAP|; N observed for Seed_Length = 81.

## Explainer choice — B24 (Seed_Length)

- **RandomForest -> `shap.TreeExplainer`** (exact, deterministic).

Because TreeExplainer is exact and the RandomForest is seeded
(`random_state=0`), the recompute reproduces script 41's frozen ranks
byte-for-byte. Verified 2026-09-19: the recomputed top-5 markers and mean|SHAP|
values match `tables/shap_top20_per_trait.csv` exactly —
100056034 (0.00864), 100041853 (0.00706), 100013936 (0.00699),
100020050 (0.00658), 100014483 (0.00648).

## Why B25 (Soluble_Oxalate) was dropped

Soluble_Oxalate's best CV model is the linear RR-BLUP (Ridge), which script 41
attributed with `shap.KernelExplainer`. With n = 37 observed samples and
p = 1,625 markers, KernelExplainer is an underdetermined, **stochastic**
coalition regression whose exact solution depends on the shap/sklearn versions
in use (background sampling, the internal `l1_reg` LASSO path, and RNG state).
Two independent problems make the Jun-2022 figure non-reproducible:

1. **The per-sample matrix was never saved** — only the top-20 ranks in
   `tables/shap_top20_per_trait.csv`. A beeswarm needs the full matrix.
2. **The exact Jun-2022 shap/sklearn versions were never recorded.** They were a
   pre-1.4 sklearn (the frozen values are physically sane, ~1–2.5), but the exact
   build is lost.

Verified 2026-09-19, three re-computations give three different top-5 marker
sets, none matching the frozen table:

| Source | top-5 markers (mean\|SHAP\|) |
|---|---|
| Frozen Jun-2022 table (cited by the manuscript) | 100045081 (2.56), 100030229 (2.34), 100016088 (1.85), 100041505 (1.72), 100021842 (1.61) |
| KernelExplainer, pinned sklearn 1.3.2 + shap 0.49.1 (`shap-repro` env) | 100015014 (2.26), 100030054 (1.68), 100040047 (1.60), 100014023 (1.45), 100015287 (1.45) |
| LinearExplainer (exact analytic phi = coef·(x−E[x])) | 100013478 (1.47), 100016704 (1.31), 100020485 (1.22), 100041762 (1.21), 100021519 (1.20) |

The current sklearn 1.6.1 in `yeast-viz` is worse still: its `LassoLarsIC`
(used by KernelExplainer's `l1_reg="auto"` when n < p) degenerates and produces
artifact attributions of order ~1e13–1e15, which is physically meaningless
(soluble oxalate is on a scale of tens of mg/100 g). This is a sklearn >=1.4
regression, not a data problem — hence the pinned `shap-repro` env for the
probe above.

**Decision (2026-09-19, B. Narh-Madey):** because no re-computation can
reproduce the frozen `shap_top20_per_trait.csv` markers that the manuscript
cites, and shipping a beeswarm whose y-axis markers contradict that table would
introduce a manuscript inconsistency, the Soluble_Oxalate beeswarm (former
Figure S10 / B25) was removed rather than replaced. The frozen top-20 table and
the SHAP-vs-GWAS overlap analysis (script 41) are unaffected and still stand in
the text; the linear-shrinkage contrast is still carried by the overlap heatmap
(B26 / Fig S10) and the ALMT4 rank scatter (B27 / Fig S11). Downstream SI
figures were renumbered (old S11–S16 -> S10–S15).

## Outputs

- `data/shap_matrix_Seed_Length.npz` — `sv` (81 x 1625), `X` (81 x 1625),
  `markers`, `samples`, `model`.
- `figures/fig89_shap_beeswarm_Seed_Length.{png,pdf}` — B24.
- Share copy: `manuscript/paper_B/share/figures/supp/B24_*`.

## Caveats a reviewer may ask about

- The beeswarm point colour is the **within-marker dosage rank** (percentile),
  not the mean-centred value, so the low/high colour scale is comparable across
  markers regardless of allele frequency.
- SHAP here is descriptive of the fitted model on the working set; it is not a
  causal or out-of-sample importance measure. The SHAP-vs-GWAS overlap table
  (script 41) is the independent-evidence guardrail against reading too much
  into a single model's attributions.

## Software versions (pinned)

- Recompute (B24, deterministic): shap 0.49.1, scikit-learn 1.6.1, numpy 2.0.2,
  Python 3.9.23 (conda env `yeast-viz`).
- B25 reproducibility probe only: scikit-learn 1.3.2 + shap 0.49.1, Python 3.9
  (conda env `shap-repro`).
