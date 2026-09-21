# 77 -- Folded SFS + Stairway model fit overlay

Two-panel figure that visualises the Paper B demographic story:

## Panel A. Folded SFS vs neutral-equilibrium expectation

For each minor-allele count k = 1 .. 94 the folded SFS bin count is plotted
as a bar; the neutral-equilibrium expectation
E[xi_k] = theta * (1/k + 1/(n-k)) / (1 + delta_{k, n-k}) is drawn on top
with theta_W rescaled so the integrated expected count matches the
1,473-segregating-site total.

n_chr = 188 (95 diploids -> 190 haplotypes -> folded over 1 -- 94 with
the half-bin singletons-vs-doubletons edge).

Inset: log2(observed / expected) for k >= 5, showing the bottleneck
signature -- the bars are biased upward (excess at intermediate frequencies)
relative to the neutral 1/k tail.

## Panel B. Stairway Plot 2 N_e(t) inferred trajectory

The Stairway-inferred N_e median curve with 75 % and 95 % posterior
shading on log-log axes. Most-recent N_e median = 108 (across the 50
nearest-to-present grid points).

## Headline takeaways

1. **The folded SFS is enriched at intermediate frequencies relative to the
   neutral expectation.** This is the SFS signature that explains the
   genome-wide Tajima D = +1.90 median: a bottleneck leaves common alleles
   over-represented relative to a constant-N_e expectation.

2. **N_e collapses into the most recent 1,000 years.** The Stairway
   trajectory rises from ~ 100 at present into the 1,000 -- 200,000 range
   in the mid-Holocene, then varies above that into deeper time. The
   recovery profile reads as "African yam bean has been demographically
   stressed in the recent past at a scale that the LD-based N_e estimate
   (659 from script 11) is consistent with at the most-recent end".

3. **The Stairway curve and Tajima D = +1.90 tell the same story.** This
   panel is what Paper B finding #1 and #2 cite when claiming "consistent
   demographic-history evidence from two independent SFS-based methods".

## Outputs

- `tables/observed_vs_neutral_sfs.csv` -- 94 rows: k, observed,
  expected_neutral, obs_minus_exp, log2_ratio.
- `figures/fig_sfs_stairway_overlay.png` / `.pdf` -- two-panel figure.

## Caveats

- **DArTseq ascertainment**. The site-frequency spectrum is computed on
  DArTseq markers selected against a complexity-reduced fraction of the
  AYB genome. Ascertainment biases the SFS toward intermediate frequencies
  compared to whole-genome SNP calling. Stairway's bottleneck inference is
  robust to mild ascertainment (Liu and Fu 2020) but a strong ascertainment
  bias toward common alleles could inflate the apparent contraction signal.
- **Folded vs unfolded**. We use folded SFS because AYB has no high-quality
  outgroup reference at the SNP level. Folding halves the resolution of
  the demographic signal at high k.
- **n_chr = 188 vs the canonical 190**. Two samples drop out of one or
  more chromosomes due to the QC sample-call-rate threshold.

## Script

`scripts/77_sfs_stairway_overlay.py`. Runtime: ~ 5 seconds.
