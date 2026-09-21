# 81 -- Per-chromosome triangular LD r^2 heatmap (Ss05, Ss10, Ss04)

Direct LD-architecture visualisation for the three chromosomes that carry
the paper's load-bearing signals.

## Method

For each chromosome we filter the AYB-anchored marker set to that
chromosome, mean-impute remaining missingness in the dosage matrix, and
compute pairwise r^2 across the chromosome's marker positions. Plot
shows a rotated triangular layout: the x axis carries genomic position
in Mb; the y axis carries pairwise half-distance (pos_j - pos_i) / 2;
cell colour is r^2.

This is the canonical Haploview-style LD heatmap, rotated 45 deg so
distant pairs sit higher on the figure.

## Chromosomes

- **Ss05** -- carries the Paper A1 Cluster-2 attribution-driven sites
  (13 of 20). Examining LD architecture here addresses "are those 13
  sites in tight LD with one another, or are they distributed across
  independent LD blocks on Ss05?"
- **Ss10** -- carries ALMT4_2 at 15,394,673 (Paper A2 four-method
  consensus locus). LD structure around this position bears on the
  per-locus interpretation of the GWAS signal.
- **Ss04** -- carries the suggestive MATE-cluster oxalate locus at
  67,792,286 and is the chromosome with the strongest Tajima D
  distribution at chromosome scale (median 2.13, 55 % of windows D > 2).

## Findings

| Chr  | n_anchored | n_pairs | r^2 median | r^2 95th pct | span_Mb |
|------|-----------:|--------:|-----------:|-------------:|--------:|
| Ss05 |        107 |   5,671 |    0.0125  |       0.169  |    56.3 |
| Ss10 |        134 |   8,911 |    0.0131  |       0.162  |    18.8 |
| Ss04 |        210 |  21,945 |    0.0121  |       0.122  |    76.0 |

## Headline takeaways

1. **LD is short-range across all three chromosomes.** Pair-wise median
   r^2 sits at ~ 0.012 (the background-noise level for an n = 95 sample
   at the MAF distribution observed). The 95th-percentile r^2 ranges
   from 0.12 (Ss04) to 0.17 (Ss05) -- a handful of pairs in tight LD on
   each chromosome.

2. **Ss10 is the densest-marker chromosome at the shortest span.** 134
   markers across 18.8 Mb gives a marker every ~ 140 kb, the highest
   density among the three. ALMT4_2 at 15.4 Mb sits in the middle of
   the chromosome; the triangular heatmap shows the local-LD
   neighbourhood around it.

3. **The ALMT4_2 region (Ss10:15.4 Mb) does not show a deep red triangle
   spanning Mb-scale distances.** Local LD around the ALMT4_2 site is
   short-range. This is internally consistent with the Paper B iHS
   cross-overlay result (script 72): ALMT4_2 is not in a sweep, the
   trait-associated allele tracks standing diversity rather than a
   recent positive-selection event that would generate extended
   haplotype block structure.

4. **Ss04 has the lowest 95th-pct r^2** (0.12 vs 0.17 on Ss05). At 210
   markers across 76 Mb, the chromosome carries the most independent
   LD blocks per Mb -- consistent with the Tajima D distribution being
   the most bottleneck-distorted at chromosome scale (median D = 2.13),
   but each individual block remaining short-range.

## Outputs

- `tables/ld_pairs_Ss05.csv` -- 5,671 pairs: rs_i, rs_j, pos_i, pos_j, r^2.
- `tables/ld_pairs_Ss10.csv` -- 8,911 pairs.
- `tables/ld_pairs_Ss04.csv` -- 21,945 pairs.
- `tables/per_chr_summary.csv` -- 3 rows: per-chr n_markers, n_pairs,
  r^2 median, r^2 95th pct, span_Mb.
- `figures/fig_per_chr_ld.png` / `.pdf` -- 3-panel triangular heatmap
  with ALMT4_2 and Ss04 MATE-cluster focal-position annotations.

## Caveats

- **DArTseq marker density is ~ one per 100 -- 700 kb across these
  chromosomes** -- insufficient to resolve LD blocks at the few-kb
  scale where most plant LD breakdown occurs. The triangle plots
  capture the medium-range (50 kb -- 5 Mb) structure but cannot speak
  to within-block LD architecture.
- **Mean-imputation for missing genotypes** shrinks r^2 toward zero for
  markers with substantial missingness. The reported values are
  conservative.
- **No phasing.** r^2 is computed on diploid dosage, which approximates
  haplotype-level r^2 well in autogamous selfers (low-heterozygosity
  regime) but slightly understates haplotype-level LD in outbred
  individuals.

## Script

`scripts/81_per_chr_ld_heatmap.py`. Runtime: ~ 60 seconds.
