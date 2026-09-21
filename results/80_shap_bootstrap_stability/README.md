# 80 -- SHAP top-5 cross-CV stability via bootstrap proxy

Reviewer-grade stability check on SHAP top-K rankings. Script 41 trained
the per-trait CV winner on the full working set once and computed SHAP
once -- the 10 nested-CV repeat-level SHAP values were not stored. To
test ranking stability at the resolution of cross-CV-repeat variance, we
substitute a **bootstrap proxy**: for each tree-based CV winner we
resample the n_work training rows with replacement, refit, recompute
SHAP, and record the top-5 markers across B = 30 bootstraps.

## Tree-based winners only

Paper B finding #9 establishes that SHAP recovers GWAS rank ONLY when the
CV winner is a tree-based estimator. The bootstrap-stability check is
therefore restricted to the three tree winners:

  - **Insoluble_Oxalate x RandomForest** (n_work = 37, CV r = 0.134)
  - **Seed_Length x RandomForest** (n_work = 81, CV r = 0.209)
  - **Seed_Thickness x LightGBM** (n_work = 81, CV r = 0.160)

## Three stability metrics

For each bootstrap pair (i, j) we compute:

(a) **Spearman rho on the full-marker mean|SHAP| rank vector.** Measures
    whether the overall importance ordering is stable.

(b) **Jaccard index on the top-5 sets.** Measures whether the same five
    markers come back at the top each time.

(c) **Per-marker selection frequency** across the 30 bootstraps -- the
    fraction of bootstraps in which a marker lands in the top-5.

## Findings

| Trait              | Model | Spearman_rho mean | Jaccard top-5 mean | Jaccard top-5 median |
|--------------------|-------|------------------:|-------------------:|---------------------:|
| Insoluble_Oxalate  | RF    |             0.275 |              ~0.01 |                  0.0 |
| Seed_Length        | RF    |             0.305 |              ~0.02 |                  0.0 |
| Seed_Thickness     | LGB   |             0.705 |              ~0.03 |                  0.0 |

The full-rank Spearman is moderate for RandomForest (~ 0.28 -- 0.31) and
substantially better for LightGBM (~ 0.71). LightGBM's feature-fraction
sampling per tree produces more globally stable feature importance than
RF's max_features-driven tree-level random subsetting at this n.

**However, the top-5 Jaccard is ~ 0 across all three traits.** Median
Jaccard at top-5 across 435 bootstrap pairs per trait is exactly zero --
the top-5 marker sets are essentially disjoint across bootstrap replicates.

Per-marker selection frequencies at top-5 across 30 bootstraps:

**Seed_Length x RF** (top 5 most-frequently selected):
  - `100020050|F|0-34:T>G-34:T>G` -- in 10 / 30 bootstraps (33 %)
  - `100013031|F|0-8:C>T-8:C>T`   -- in 6 / 30 bootstraps (20 %)
  - `100033134|F|0-15:C>T-15:C>T` -- in 5 / 30 bootstraps (17 %)
  - `100041853|F|0-27:G>A-27:G>A` -- in 5 / 30 bootstraps (17 %)
  - `100047059|F|0-34:G>T-34:G>T` -- in 4 / 30 bootstraps (13 %)

**Seed_Thickness x LGB** (top 5):
  - `100035958|F|0-25:C>T-25:C>T` -- in 10 / 30 bootstraps (33 %)
  - `100041762|F|0-5:T>C-5:T>C`   -- in 10 / 30 bootstraps (33 %)
  - `100032753|F|0-47:A>G-47:A>G` -- in 9  / 30 bootstraps (30 %)
  - `100021017|F|0-55:C>T-55:C>T` -- in 7  / 30 bootstraps (23 %)
  - `100016587|F|0-41:C>T-41:C>T` -- in 6  / 30 bootstraps (20 %)

**Insoluble_Oxalate x RF** (top 5):
  - `100035301|F|0-17:C>T-17:C>T` -- in 6 / 30 bootstraps (20 %)
  - `100035901|F|0-29:T>A-29:T>A` -- in 5 / 30 bootstraps (17 %)
  - `100035171|F|0-38:C>T-38:C>T` -- in 5 / 30 bootstraps (17 %)
  - `100029472|F|0-28:G>A-28:G>A` -- in 4 / 30 bootstraps (13 %)
  - `100045385|F|0-50:G>A-50:G>A` -- in 4 / 30 bootstraps (13 %)

## Headline takeaways

1. **Full-rank importance is stable; top-K is not.** The Spearman of
   the full-marker importance ranking is 0.7 for LGB and 0.3 for RF
   across bootstraps. But the top-5 set is almost completely re-chosen
   each bootstrap. SHAP top-K interpretation is therefore fragile at our
   n = 37 -- 81 panel sizes.

2. **Selection-frequency framing replaces top-K framing.** Even the
   most-stable marker per trait shows up in only 33 % of bootstraps.
   Reporting "this marker is in the SHAP top-5 of trait X" is
   methodologically thin at this n. The honest framing is "this marker
   appears in the SHAP top-5 of trait X in 1/3 of bootstrap replicates"
   -- which still rules out the chance baseline (5 / 1,625 = 0.3 %) by
   ~ 100 fold, but cuts the over-claim about a "single best marker".

3. **Methodological caveat for the Paper B SHAP discussion.** The text
   should report SHAP top-5 markers with the per-marker bootstrap
   selection frequency rather than as a singleton list.

## Bootstrap vs k-fold cross-validation stability

Bootstrap perturbs the training set more aggressively than k-fold -- each
bootstrap replicate sees ~ 36 % of rows duplicated and ~ 64 % of unique
rows. K-fold leaves out one fold at a time. The bootstrap-based stability
metric is therefore a **conservative lower bound** on top-K stability at
k-fold resolution; the true Jaccard across 10 CV repeats would likely be
slightly higher than the values reported here, but the qualitative
"top-K is unstable" finding will hold.

## Outputs

- `tables/per_boot_topK_markers.csv` -- per-bootstrap, per-trait
  top-5 marker rows.
- `tables/Insoluble_Oxalate_RandomForest_selection_freq.csv` -- per-marker
  bootstrap selection frequency for the trait.
- `tables/Seed_Length_RandomForest_selection_freq.csv` -- same for
  Seed_Length.
- `tables/Seed_Thickness_LightGBM_selection_freq.csv` -- same for
  Seed_Thickness.
- `tables/pairwise_stability.csv` -- 1,305 rows: per-(trait, model)
  pairwise Spearman + Jaccard top-5 per bootstrap pair.
- `tables/stability_summary.csv` -- aggregated stability per trait.
- `figures/fig_shap_bootstrap_stability.png` / `.pdf` -- two-panel
  figure: per-trait stability boxplot + selection-frequency bar chart.

## Caveats

- **Bootstrap is not k-fold.** The proxy gives conservative bounds.
- **B = 30 is modest.** A defensible production run would use B = 100;
  the qualitative conclusion (top-K unstable, full-rank stable for LGB)
  holds at B = 30.
- **Per-bootstrap shap.TreeExplainer.shap_values(X_boot) is computed on
  the bootstrap sample, not held-out.** This is the same convention as
  script 41. A held-out version would use the out-of-bag rows; this is
  worth doing as a follow-up if reviewers push back on the bootstrap-
  in-bootstrap-out framing.

## Script

`scripts/80_shap_bootstrap_stability.py`. Runtime: ~ 4 minutes on the
local workstation (B = 30 x 3 traits with RF n_estimators = 500 and LGB
n_estimators = 300).
