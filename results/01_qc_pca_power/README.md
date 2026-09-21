# 01 — QC, PCA, and GWAS power analysis

**Script:** `scripts/01_qc_pca_power.py`

## What was done

Marker and sample quality control on the raw DArTseq DAf18-2580 HapMap file (3,204 SNPs × 105 accessions), principal component analysis on the QC-passing dosage matrix, K-means clustering for population structure, and a per-(MAF × h²_QTL) GWAS power grid for an n = 95 panel at Bonferroni α / m_eff.

## Method

HapMap calls were converted to 0/1/2 alternate-allele dosage. Per-marker call rate, MAF, observed heterozygosity, and PIC were computed; per-sample call rate and heterozygosity were also tabulated. Filters: sample call rate ≥ **0.90**, marker call rate ≥ **0.90**, MAF ≥ **0.05**. Missing genotypes were mean-imputed per marker. PCA used scikit-learn 1.4 on the standardised dosage; K-means *k* selected by silhouette score across k = 2..7. Power was computed as χ² CDF on the non-centrality parameter `NCP = n h² / (1 − h²) × 2 MAF (1 − MAF)` against the m_eff Bonferroni threshold.

## Quick findings

- Retained **95 samples × 1,625 markers** post-QC.
- Median MAF 0.20, median PIC 0.27.
- Best K-means *k* = **2** (silhouette = 0.42). Cluster sizes: large = 84, small = 11.
- PC1 + PC2 explain **14.1 %** of marker variance.
- Power analysis: at n = 95 with m_eff = 1,625, detection requires h²_QTL ≥ 0.20 at MAF = 0.20 for >80 % power.

## Caveats

- Cluster labels are an artefact of unsupervised k-means on the SNP-derived PC space; they are not independent of the marker structure.
- Power grid uses *raw m* = 1,625 not the Li & Ji m_eff = 93 (which is reported in `12_power_meff`). The relaxed m_eff threshold is the one actually used in the GWAS scans (scripts 03, 15).
- No outside provenance metadata to validate the cluster interpretation; the 11-line outlier subgroup could be a genebank-collection-history artefact rather than a population-genetics signal.
