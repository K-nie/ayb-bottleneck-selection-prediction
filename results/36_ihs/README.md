# 36 — Beagle 5 phasing + selscan iHS scan

**Script:** `scripts/36_ihs.py`

## What was done

Beagle 5 statistical phasing followed by selscan integrated-haplotype-
homozygosity (iHS) scan for recent positive selection (Voight et al. 2006
*PLoS Biol*; Browning & Browning 2018 *Am J Hum Genet*; Szpiech & Hernandez
2014 *Mol Biol Evol*).

## Method

- Built a diploid VCF from the AYB-anchored dosage matrix (95 samples ×
  1,473 markers; 19 duplicate-position markers dropped → 1,454).
- Phased with Beagle 5.5 using default parameters, 4 threads.
- Per-chromosome selscan iHS with a uniform 1.06 cM/Mb genetic map
  (Lonardi 2019 cowpea proxy). DArTseq marker spacing (~500 kb median) is
  much sparser than the WGS density selscan was designed for, so
  `--max-gap 5,000,000` and `--max-extend 20,000,000` were used to allow
  EHH to decay before selscan aborts the calculation.
- iHS scores standardised within MAF bins across chromosomes via
  selscan's `norm` tool.

## Inputs

- HapMap dosage matrix (post-QC, MAF ≥ 0.05)
- AYB marker anchoring table

## Outputs

- `data/ayb.vcf.gz` — unphased input VCF
- `data/ayb.phased.vcf.gz` — Beagle 5 output (reused by scripts 37 and 39)
- `data/Ss*.ihs.ihs.out` — per-chromosome raw iHS
- `data/Ss*.ihs.ihs.out.100bins.norm` — standardised iHS
- `tables/ihs_combined.csv` — combined per-SNP normalised iHS
- `figures/fig80_ihs_manhattan.png/.pdf` — Manhattan plot of |iHS| with the
  |iHS| = 2 reference line

## Quick findings

- 1,454 phased markers across Ss01–Ss11.
- 314 SNPs received a finite normalised iHS score (the rest fell within
  marker-spacing gaps where EHH could not decay reliably).
- 17 SNPs with |iHS| > 2 — candidate positive-selection signals.

## Caveats

- iHS is calibrated for WGS density (1–5 kb median spacing). At DArTseq
  density the test loses power and the |iHS| > 2 set should be treated as
  a *suggestive* outlier list, not a formal scan.
- The uniform-rate genetic map ignores fine-scale recombination variation;
  the recombination heterogeneity inferred by LDhelmet in script 39
  partially addresses this for downstream re-analysis.
- Phasing accuracy is limited by panel size (95 samples); long-range
  switch errors are expected at the chromosome ends.
