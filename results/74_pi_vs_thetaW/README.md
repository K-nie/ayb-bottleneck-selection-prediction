# 74 -- pi vs theta_W per-window scatter

Bottleneck visualisation. Under neutrality with constant N_e, pi (mean
pairwise diversity at sites that segregate) tracks theta_W (Watterson
estimator from segregating-site count). Under a recent bottleneck, theta_W
drops faster than pi because rare alleles are lost first; the result is
windows with pi > theta_W and positive Tajima D.

## Method

Read the 255 windowed-Tajima-D rows from `results/35_tajima_d/tables/
tajima_d_windowed.csv`. Each row carries per-window pi-per-site, theta_W,
and Tajima D from the AYB-anchored marker set.

Plot pi_per_site against theta_W on log-log axes; point colour codes
chromosome; point size codes n_snps; the y = x diagonal marks neutrality.
Marginal histograms on each axis. Below, a per-chromosome boxplot of
Tajima D.

## Findings

Across 255 windows, the cloud sits above the y = x line at intermediate
diversity values -- the bottleneck signature. Per-chromosome Tajima D
medians (255 windows total):

| chr   | n_win | pi_median | theta_W_median | D_median | frac D > 0 | frac D > 2 |
|-------|------:|----------:|---------------:|---------:|-----------:|-----------:|
| Ss01  |    25 |  7e-6     |    -           |     2.30 |       0.92 |       0.52 |
| Ss02  |    27 |  6e-6     |    -           |     1.74 |       0.89 |       0.41 |
| Ss03  |    35 |  6e-6     |    -           |     1.22 |       0.80 |       0.43 |
| Ss04  |    40 |  7e-6     |    -           |     2.13 |       1.00 |       0.55 |
| Ss05  |    22 |  6e-6     |    -           |     1.30 |       0.82 |       0.27 |
| Ss06  |    18 |  6e-6     |    -           |     1.68 |       1.00 |       0.39 |
| Ss07  |    11 |  10e-6    |    -           |     1.63 |       1.00 |       0.45 |
| Ss08  |    27 |  8e-6     |    -           |     1.87 |       0.96 |       0.48 |
| Ss09  |    11 |  10e-6    |    -           |    *2.77*|       1.00 |     **0.73** |
| Ss10  |    31 |  7e-6     |    -           |     1.95 |       0.97 |       0.48 |
| Ss11  |     8 |  7e-6     |    -           |     1.87 |       1.00 |       0.50 |

The bottleneck signature is genome-wide (every chromosome's median D > 0),
but its strength varies. Ss09 is the outlier: median D = 2.77 and 73 % of
windows exceed D = 2 -- the strongest bottleneck signature in the panel.
Ss04 carries the most windows (n = 40) and the second-strongest D
distribution (median 2.13, 55 % > 2). Ss05 sits at the bottom of the
chromosome-level D ranking (median 1.30, 27 % > 2).

## Headline takeaways

1. **Bottleneck is genome-wide.** Every chromosome has median Tajima D > 0
   and frac(D > 0) >= 0.80. The recent-bottleneck signature recovered by
   the genome-wide median (D = +1.90, panel-wide; Paper B finding #2) is
   not driven by a small subset of chromosomes.

2. **Ss09 is the strongest-bottleneck chromosome.** Median D = 2.77, 73 %
   of windows D > 2. Worth flagging in the Discussion if any candidate-gene
   work touches Ss09 in the future -- the demographic distortion is
   particularly pronounced there.

3. **Ss05 is the weakest-bottleneck chromosome.** Median D = 1.30, 27 % of
   windows D > 2. Still positive (consistent with a bottleneck) but well
   below the panel mean -- relevant context for Paper A1's Cluster-2
   attribution work that names 13 of 20 Cluster-2-driven sites as Ss05
   markers.

## Outputs

- `tables/per_chr_summary.csv` -- per-chromosome window count, pi/theta_W
  medians, Tajima D medians, frac(D > 0) and frac(D > 2).
- `figures/fig_pi_vs_thetaW.png` / `.pdf` -- pi-vs-theta_W scatter with
  marginal histograms and per-chromosome Tajima D boxplot.

## Script

`scripts/74_pi_vs_theta_w_scatter.py`. Runtime: ~ 5 seconds.
