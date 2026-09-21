# 72 — iHS × A2 GWAS-top-30 cross-overlay

Quantifies the relationship between Paper A2's per-trait GWAS top-30 SNP sets and Paper B's selscan iHS recent-positive-selection scan. The headline finding **ALMT4_2 is not in an iHS sweep** is the cross-paper integration that ties the trait-genetics story of A2 to the demographic-history story of B.

## Method

1. **Per trait** — look up each of the 30 top GWAS SNPs (from `results/21_candidate_genes_ayb/tables/snp_top30_anchored_<trait>.csv`) in the iHS table. Record `norm_iHS`, `|norm_iHS|`, the binary selscan `crit` flag, and the genomic distance to the nearest |iHS| > 2 SNP on the same chromosome (the GWAS top-30 set is a subset of the AYB-anchored 2,862 markers; the iHS table runs on 314 SNPs that survived Beagle phasing, so most GWAS SNPs do not have an exact-position iHS score and we report the nearest-neighbour).

2. **Per A2 candidate-gene focal region** — compute windowed iHS density at ± 500 kb, ± 1 Mb, ± 2 Mb around the focal SNP. Two focal regions: Soluble_Oxalate × ALMT4_2 at Ss10:15,394,673 (the A2 four-method consensus locus); Insoluble_Oxalate × Ss04 MATE cluster at Ss04:67,792,286 (the A2 suggestive partition-side hit).

3. **Per trait overlap-vs-panel-rate test** — under the panel-wide |iHS| > 2 rate (17 of 314 = 5.4 %), the expected number of |iHS| > 2 SNPs in a random 30-SNP sample is 1.62. We test each trait's observed count against this baseline under a one-sided binomial tail.

## Findings

**A2 focal regions × iHS density:**

| A2 focal region | Chr | Position | |iHS|>2 within ±500 kb | within ±1 Mb | within ±2 Mb | Nearest |iHS|>2 |
|---|---|---|---|---|---|---|
| **ALMT4_2** | Ss10 | 15,394,673 | **0** | **0** | **0** | **NONE on Ss10** |
| Ss04 MATE cluster | Ss04 | 67,792,286 | 0 | 0 | 1 | 1.29 Mb away at norm_iHS = +2.26 |

**Per-trait GWAS top-30 × |iHS| > 2 overlap (binomial test against panel-wide 5.4 % rate):**

| Trait | n top-30 | with iHS at exact SNP | at |iHS| > 2 | Expected | Binom p (1-sided) |
|---|---|---|---|---|---|
| Total_Oxalate | 30 | 12 | 2 | 1.62 | 0.49 |
| Flavonoid | 30 | 11 | 1 | 1.62 | 0.81 |
| Insoluble_Oxalate | 30 | 8 | 1 | 1.62 | 0.81 |
| Phenol | 30 | 4 | 1 | 1.62 | 0.81 |
| Seed_Width | 30 | 10 | 1 | 1.62 | 0.81 |
| Tannin | 30 | 8 | 1 | 1.62 | 0.81 |
| Antioxidant | 30 | 7 | 0 | 1.62 | 1.00 |
| Crude_Protein | 30 | 3 | 0 | 1.62 | 1.00 |
| Mass_of_Seeds | 30 | 8 | 0 | 1.62 | 1.00 |
| Seed_Coat_Tannin | 30 | 9 | 0 | 1.62 | 1.00 |
| Seed_Length | 30 | 8 | 0 | 1.62 | 1.00 |
| Seed_Thickness | 30 | 10 | 0 | 1.62 | 1.00 |
| **Soluble_Oxalate** | 30 | 5 | **0** | 1.62 | 1.00 |

**Three takeaways:**

1. **ALMT4_2 is not in an iHS sweep.** Ss10 carries zero |iHS| > 2 outliers anywhere; the closest is on a different chromosome. The Soluble_Oxalate × ALMT4 trait-genetic signal from A2 tracks **standing diversity** in the panel rather than an **ongoing recent positive selection** event. This is internally consistent with AYB's underutilised-crop status — the panel has not been intensively selected at modern timescales, so the trait-associated alleles segregating at the ALMT4 locus reflect the standing-variation pool the bottleneck preserved rather than a sweep on a single advantageous allele.

2. **The Ss04 MATE cluster sits ~ 1.3 Mb from one iHS outlier** (Ss04:66,505,361, norm_iHS = +2.26). At this n the 2 Mb resolution is too coarse to call the MATE-cluster region "in a sweep" — the outlier could be unrelated and 1.29 Mb of genomic distance leaves room for many independent LD blocks. The Insoluble_Oxalate mechanism call rests on the A2 evidence (single-marker p = 0.18 + tandem MATE-family gene proximity), not on iHS support.

3. **No trait shows GWAS-top-30 × iHS-outlier enrichment over the panel rate.** Total_Oxalate has 2 of 30 SNPs at |iHS| > 2 (vs expected 1.62; one-sided binomial p = 0.49); all other traits have 0 or 1 hits. The panel's recent-positive-selection footprint is small (17 outliers across 314 SNPs) and does not co-locate with the trait-associated GWAS regions at this n. The interpretation matches the §3 framing of the Paper B Discussion: the bottleneck signature (Tajima D median = +1.90; 47 % of windows > 2) is genome-wide, not locus-specific, and the per-locus trait signals come from standing variation rather than from recent sweeps.

## Outputs

- `tables/per_trait_top30_ihs_lookup.csv` — 390 rows: per-trait per-top-30-SNP iHS lookup (13 traits × 30 SNPs).
- `tables/a2_focal_region_density.csv` — 2 rows: windowed iHS density at the ALMT4_2 and Ss04 MATE-cluster focal regions.
- `tables/per_trait_overlap_summary.csv` — 13 rows: per-trait overlap statistic with the binomial one-sided p-value.
- `figures/fig_ihs_at_a2_callouts.png` / `.pdf` — iHS Manhattan with the two A2 candidate-gene focal regions highlighted (vermillion bands).

## Caveats

- **iHS coverage is sparse.** Only 314 SNPs in the iHS table vs 2,862 AYB-anchored markers. selscan requires phased haplotypes and minor allele frequency filters that drop a substantial fraction. The "nearest-iHS" distance for most GWAS SNPs is therefore large; the per-trait lookup table reports this explicitly.
- **n = 95 limits iHS power.** The expected number of |iHS| > 2 SNPs under a strong sweep at this panel size is small to begin with. The absence of an ALMT4-region sweep signal is consistent both with "no sweep" and with "weak sweep we cannot detect at this n" — we report the former as the parsimonious reading but flag the alternative.
- **Cross-paper integration framing.** The result is a B finding, not an A2 finding. The A2 trait-genetic case for ALMT4 stands on the four-method GWAS concordance + h² CI + protein domain + phylogenetic clade. Paper B contributes the demographic context: the trait-associated allele segregates in standing variation, not in a sweep.

## Script

`scripts/72_ihs_vs_a2_gwas_overlay.py`. Runtime: ~ 5 seconds.
