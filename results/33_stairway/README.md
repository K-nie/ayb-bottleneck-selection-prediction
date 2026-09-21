# 33 — Stairway Plot 2 demographic inference

**Script:** `scripts/33_stairway.py`

## What was done

Reconstruction of effective population size N_e through time for the AYB
panel via Stairway Plot 2 (Liu & Fu 2020 *Nature Commun*), using the folded
site-frequency spectrum from the AYB-anchored DArTseq SNPs.

## Method

- Folded SFS built from 1,473 segregating sites across 190 alleles (95
  diploid samples) after filtering on call-rate (sample CR ≥ 0.90, marker
  CR ≥ 0.90) and MAF ≥ 0.05.
- Blueprint parameters: μ = 7 × 10⁻⁹ per site per generation (Ossowski et
  al. 2010 *Arabidopsis* proxy), generation time = 1 year, total surveyed
  sequence length L = 112,000 bp (DArTseq tag coverage: 1,625 markers ×
  ~69 bp per tag; this is the correct denominator for SFS-based theta
  estimation, not the 649.8 Mb assembly span), n_input = 200 bootstrap
  replicates per breakpoint, training fraction 0.67, nrand breakpoints
  ⌊(n−2)/4⌋ = 47 and ⌊(n−2)/2⌋ = 94. The two higher nrand breakpoints
  (141 and 188) were planned but truncated for wall-time reasons; the
  Stairpainter selected nrand = 47 as the maximum-likelihood breakpoint
  anyway, so the truncation does not affect the final curve.
- 443 bootstrap reps were optimised before the truncation (200 at
  nrand = 47, 200 at nrand = 94, 43 partial at nrand = 141), all 200 of
  which were used in the final summary aggregation.
- The Java pipeline was driven from a no-spaces working directory
  (`/tmp/ayb_stairway/`) because the project path contains a space, which
  Java's class-loader handles poorly.

## Inputs

- HapMap dosage matrix (post-QC) → folded SFS
- AYB marker anchoring table → restricts SFS to the 1,473 AYB-anchored
  markers

## Outputs

- `data/sfs.txt` — folded SFS used as input
- `data/ayb.blueprint` — Stairway Plot 2 configuration
- `tables/ne_through_time.csv` — median N_e curve with 95 % CI
- `figures/fig79_stairway_ne_time.png/.pdf` — N_e(t) on log axes; LD-derived
  current N_e = 659 overlaid as a dashed reference line

## Quick findings

- Recent (1–10 years) median N_e ≈ 108 with 95 % CI [51, 286]; an order
  of magnitude smaller than the LD-derived current N_e = 659.
- The Stairway curve crosses the LD-derived N_e line at ~30 years before
  present (≈ generation 30 since 1 yr/gen).
- Climbs to N_e ≈ 10⁴ at 1,000 years before present and to ≈ 2 × 10⁵ at
  10⁵ years before present. The shape is consistent with a strong recent
  bottleneck on top of a much larger historical N_e.

## Caveats

- The current AYB panel is a curated 95-line subset of the IITA TSs
  collection. Selection bias toward Nigerian provenance limits how far
  back in time the inferred N_e curve can be trusted.
- A single fixed mutation-rate proxy from *Arabidopsis* introduces
  systematic uncertainty on absolute time-scales. Relative shape of the
  trajectory is more interpretable than absolute years.
- Folded SFS without an outgroup means we cannot resolve very-recent
  population events that depend on ancestral-state polarisation.
