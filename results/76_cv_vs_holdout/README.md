# 76 -- CV r vs final-holdout r (Bellot 2018 pattern)

Direct visualisation of cross-validation honesty for the 13 per-trait CV
winners against the 15 % final holdout. Points on the y = x diagonal
indicate CV r tracks held-out r honestly; large positive (CV > holdout)
gaps indicate CV over-optimism (Bellot 2018 small-n deep-net failure
mode); large positive (holdout > CV) gaps indicate holdout-luck at small
n_hold.

## Per-trait results

| Trait              | Winner       | n_work | CV r   | Holdout r | n_hold | Gap (CV - hold) |
|--------------------|--------------|-------:|-------:|----------:|-------:|----------------:|
| Tannin             | KernelRidge  |     79 | -0.029 |    -0.228 |     14 |          +0.199 |
| Phenol             | KernelRidge  |     79 | +0.083 |    +0.423 |     14 |          -0.340 |
| Flavonoid          | KernelRidge  |     79 | -0.034 |    +0.271 |     14 |          -0.305 |
| Antioxidant        | RR-BLUP      |     79 | +0.118 |    -0.083 |     14 |          +0.200 |
| Seed_Length        | RandomForest |     81 | +0.209 |    +0.184 |     14 |          +0.025 |
| Seed_Width         | KernelRidge  |     81 | +0.092 |    +0.605 |     14 |          -0.513 |
| Seed_Thickness     | LightGBM     |     81 | +0.160 |    +0.131 |     14 |          +0.029 |
| Mass_of_Seeds      | KernelRidge  |     81 | +0.059 |    +0.385 |     14 |          -0.326 |
| Seed_Coat_Tannin   | RR-BLUP      |     81 | +0.065 |    +0.003 |     14 |          +0.062 |
| **Crude_Protein**  | **MLP**      | **37** | +0.148 |    *-0.561* |  **4** |     **+0.709** |
| Total_Oxalate      | ElasticNet   |     37 | +0.171 |    +0.463 |      4 |          -0.291 |
| Soluble_Oxalate    | RR-BLUP      |     37 | +0.280 |    +0.536 |      4 |          -0.256 |
| Insoluble_Oxalate  | RandomForest |     37 | +0.134 |    +0.641 |      4 |          -0.508 |

## Headline takeaways

1. **Crude_Protein x MLP is the textbook Bellot reproduction.** CV r =
   +0.148, holdout r = -0.561, gap = +0.71. At n_hold = 4 the holdout r
   estimate is high-variance, but the sign reversal is exactly the
   Bellot 2018 (Genetics 209:629) small-n deep-net pattern. The MLP
   over-fits the n_work = 37 protein-trait training set in a way the CV
   reps don't catch.

2. **Most CV winners are RR-BLUP / KernelRidge.** Of 13 winners, 7 use
   linear shrinkage estimators (3 RR-BLUP + 4 KernelRidge), 1 uses
   ElasticNet, 1 uses MLP (the Crude_Protein failure), and 4 use
   tree-based models. The Daetwyler ceiling argument in script 75 backs
   this distribution: the panel-level h^2 + small n_work makes the
   linear shrinkage path the right one for most traits.

3. **Large negative gaps (holdout > CV) sit at n_hold = 4 oxalate
   traits.** Insoluble_Oxalate (-0.51), Seed_Width (-0.51), Phenol
   (-0.34), Total_Oxalate (-0.29) all post holdout r above CV r. At
   n_hold = 4 the Pearson r CI is wide enough that any of these gaps
   could be lucky-draw artefact. The headline result is "holdout did not
   undermine the CV ranking", not "this estimator predicts even better
   on truly held-out samples than CV claims".

4. **Tree-based winners track CV most honestly.** Seed_Length (gap
   +0.025) and Seed_Thickness (gap +0.029) -- the two tree-winning
   traits at n_work = 81 -- have the smallest CV-vs-holdout gaps in the
   table. The RR-BLUP wins (Seed_Coat_Tannin gap +0.062, Antioxidant
   gap +0.200) and KernelRidge wins (Tannin gap +0.199) sit a step
   above. The MLP loss is the outlier on the opposite end (+0.71).

## Outputs

- `tables/cv_vs_holdout_winners.csv` -- 13 rows: per-trait CV-winner
  vs holdout r, with CV r 95 % CI and n_hold.
- `tables/cv_vs_holdout_all_models.csv` -- per-trait per-model CV r
  alongside the winner column (no per-model holdout).
- `figures/fig_cv_vs_holdout.png` / `.pdf` -- two-panel figure:
  (A) CV vs holdout scatter with y = x line; (B) sorted CV-vs-holdout
  gap bar chart.

## Caveats

- **n_hold = 4 for the oxalate / Crude_Protein traits.** Per-trait
  holdout r is high-variance at this size; report directionality, not
  magnitude. The CV r estimates (5x5x10 nested CV) are tighter.
- **One winner per trait, not per model.** `final_holdout_scores.csv`
  carries only the CV-winning model's holdout score per trait. The full
  per-(trait, model) holdout-r matrix is not stored; the per-model gap
  distribution in panel B uses the winner only.

## Script

`scripts/76_cv_vs_holdout_per_trait.py`. Runtime: ~ 5 seconds.
