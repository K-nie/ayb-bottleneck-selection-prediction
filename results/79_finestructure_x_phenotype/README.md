# 79 -- fineSTRUCTURE coancestry x phenotype side-annotation

Cross-paper bridge B <-> A2: do the K = 24 haplotype-resolution clusters
recovered by ChromoPainter + fineSTRUCTURE stratify any of the 13 phenotype
BLUPs?

## Method

Heatmap: 95 x 95 coancestry chunk-count matrix, samples reordered to put
K = 24 cluster members adjacent. Right-hand annotations: 13 columns, one
per trait, coloured by per-sample BLUP z-score. Top bar: cluster identity.

Statistics: per-trait Kruskal-Wallis H test of BLUP ~ K = 24 cluster
membership, with BH-FDR across the 13 traits.

## Findings

Kruskal-Wallis trait ~ K = 24 cluster (only clusters with n >= 2 in
non-missing trait coverage):

| Trait              | n_used | n_groups | H      | p_raw  | p_BH   |
|--------------------|-------:|---------:|-------:|-------:|-------:|
| Seed_Thickness     |     94 |       23 |  37.20 | 0.0225 | 0.2214 |
| Seed_Length        |     94 |       23 |  35.53 | 0.0341 | 0.2214 |
| Mass_of_Seeds      |     94 |       23 |  33.83 | 0.0511 | 0.2214 |
| Seed_Width         |     94 |       23 |  29.37 | 0.1347 | 0.4377 |
| Soluble_Oxalate    |     36 |       13 |  14.54 | 0.2673 | 0.6951 |
| Tannin             |     92 |       23 |  24.38 | 0.3275 | 0.7096 |
| Phenol             |     92 |       23 |  21.58 | 0.4852 | 0.9011 |
| Antioxidant        |     92 |       23 |  20.26 | 0.5669 | 0.9212 |
| Total_Oxalate      |     36 |       13 |   6.63 | 0.8809 | 0.9809 |
| Insoluble_Oxalate  |     36 |       13 |   6.10 | 0.9110 | 0.9809 |
| Seed_Coat_Tannin   |     94 |       23 |  13.68 | 0.9126 | 0.9809 |
| Crude_Protein      |     91 |       23 |  11.91 | 0.9592 | 0.9809 |
| Flavonoid          |     92 |       23 |  10.52 | 0.9809 | 0.9809 |

## Headline takeaway (honest negative result)

**No trait survives BH-FDR correction.** Seed-size traits (Seed_Thickness,
Seed_Length, Mass_of_Seeds) sit at p_raw = 0.02 -- 0.05 but the lowest
BH-corrected q is 0.22 -- well above any conventional FDR threshold.
The K = 24 fineSTRUCTURE clusters do NOT stratify phenotype variance at
levels that survive multiple-testing correction.

This is the expected pattern at the panel size: with n = 95 split across
24 clusters, most clusters carry 1 -- 6 samples. The H statistic at K = 24
is a high-degrees-of-freedom test and the per-cluster n is too small to
detect cluster-level trait mean differences except where the underlying
effect is large.

The directionality is suggestive for seed-size traits, consistent with
the K = 24 clusters reflecting breeding pedigrees within IITA's TSs
collection that may stratify the seed-trait variance, but the panel is
under-powered to claim phenotype association at this resolution.

## What the figure does

(i) Visualises the K = 24 block-diagonal coancestry structure that
fineSTRUCTURE recovers. Cluster boundaries are drawn as faint white lines
on the heatmap. The 95 x 95 matrix is sharply block-diagonal under the
cluster reorder -- a clean visual that K = 24 is reading something real
in the haplotype-painting signal.

(ii) Visualises trait BLUP z-scores per sample alongside the cluster
identity. Visual inspection of the seed-size column block shows weak
horizontal banding within clusters but no clean run-of-positive or
run-of-negative cluster blocks -- which matches the statistical result.

(iii) Documents the honest negative result: the K = 24 partition is
genealogically real (high coancestry within cluster) but phenotypically
modest at the BLUP level.

## Outputs

- `tables/kruskal_trait_vs_K24.csv` -- 13 rows: per-trait Kruskal-Wallis
  H, p_raw, p_BH.
- `tables/per_cluster_trait_means.csv` -- per-(trait, cluster) BLUP mean
  and n.
- `figures/fig_coancestry_x_phenotype.png` / `.pdf` -- coancestry
  heatmap with trait BLUP side-annotation.

## Caveats

- **K = 24 / n = 95 inflates degrees of freedom.** Kruskal-Wallis at 22
  df under H_0 has expected H = 22. The seed-size H values (33 -- 37)
  are above expectation by a margin that doesn't survive BH-FDR but
  would warrant follow-up at a larger n.
- **Selfing-rate uncertainty.** AYB is a presumed autogamous selfer but
  the F_IS / F_ROH paradox (Paper B finding #5) suggests the panel
  contains accessions with mixed mating histories. K = 24 clusters may
  partly reflect that mating-history axis rather than pedigree-line
  structure alone.

## Script

`scripts/79_finestructure_phenotype_overlay.py`. Runtime: ~ 15 seconds.
