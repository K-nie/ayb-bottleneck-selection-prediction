# 75 -- Per-trait CV r vs narrow-sense h^2 with the Daetwyler ceiling

The Daetwyler / Goddard genomic-prediction theory sets an informational
ceiling for prediction accuracy as a function of three quantities:

    r_max = sqrt( h^2 * N_p / (N_p + M_e * h^2) )

where N_p is the training set size, h^2 is narrow-sense heritability, and
M_e is the effective number of independent chromosome segments. At our
n_work = 30 -- 75 panel sizes with M_e in the 1,500 -- 3,000 range
(informed by AYB LD decay r^2 ~ 0.1 by ~ 200 kb and a ~ 600 Mb genome),
the ceiling sits at r_max ~ 0.3 -- 0.5 even when h^2 reaches 0.7.

## Plot

x axis: REML-profile narrow-sense h^2 with one-sided CI (from script 24).
y axis: observed CV r of the best non-Dummy estimator per trait, with
95 % paired-rep CI from script 40. Trait labels annotated. Two reference
Daetwyler ceiling curves drawn at the median n_work = 79 -- one at M_e =
1,500, the other at M_e = 3,000. The `sqrt(h^2)` upper bound is also
plotted as a soft reference.

## Headline takeaways

1. **All 13 traits land at or below the Daetwyler r_max curve.** No
   estimator beats the informational ceiling. Soluble_Oxalate (h^2 MLE
   pinned at the upper edge with one-sided CI 0.14 -- 0.999, RR-BLUP CV r
   = 0.280) is the only trait where the ceiling is comfortably above 0.3.

2. **Negative-CV-r traits sit at the h^2 = 0.001 lower boundary.** Tannin
   (KernelRidge, CV r = -0.029), Flavonoid (CV r = -0.034), Antioxidant
   (CV r = 0.118 -- but lower bound -0.10) all have REML MLE pinned at
   0.001 with upper CI under 0.31. Their poor predictive performance is
   informational, not algorithmic -- there's no extractable signal at
   this n.

3. **Seed-size traits show a Daetwyler-consistent pattern.** Seed_Length
   (RF, CV r = 0.209) and Seed_Thickness (LGB, CV r = 0.160) sit just
   below the M_e = 1,500 ceiling. Their h^2 estimates are middling
   (0.55 -- 0.66 upper CI) but the small N_p = 81 keeps the ceiling near
   r_max ~ 0.4 -- 0.5.

4. **The estimator legend confirms the linear / kernel-shrinkage win.**
   RR-BLUP wins 3 traits; KernelRidge wins 5. Tree-based estimators win
   only Seed_Length / Seed_Thickness / Insoluble_Oxalate (3 traits).
   ElasticNet wins Total_Oxalate. MLP wins Crude_Protein (with the Bellot
   holdout collapse documented in script 76). At this panel size, the
   Daetwyler ceiling argues that the linear / kernel-shrinkage estimators
   absorb the available signal -- there's no algorithmic ceiling left
   for non-linear models to break.

## Outputs

- `tables/per_trait_cv_vs_h2.csv` -- 13 rows: best-CV-model per trait
  joined with h^2 MLE + CI.
- `figures/fig_cv_r_vs_h2.png` / `.pdf` -- scatter with two Daetwyler
  ceiling curves and per-trait CIs.

## Caveats

- **M_e is approximate.** AYB has no published M_e estimate at the
  DArTseq scale. M_e = 1,500 corresponds to one independent block per
  ~ 400 kb (AYB LD-decay regime); M_e = 3,000 corresponds to one block
  per ~ 200 kb (tighter LD). Either value bounds the ceiling well below
  perfection at our small n.
- **h^2 lower-CI floor is the grid limit, not the true posterior**. The
  one-sided CI for many traits hits 0.001 at the grid floor -- which is
  consistent with no detectable h^2 at this n, but the true posterior
  could fall slightly above zero.
- **N_p in the formula is per-trait n_work**, not the panel size (95).
  Crude_Protein, Total_Oxalate, Soluble_Oxalate, Insoluble_Oxalate sit
  at n_work = 37 (missing phenotypes drop the trait coverage). The
  Daetwyler ceiling curves in the plot use the median across all 13
  traits (n = 79).

## Script

`scripts/75_cv_r_vs_h2_daetwyler.py`. Runtime: ~ 5 seconds.
