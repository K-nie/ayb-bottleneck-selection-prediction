# 37 — ChromoPainter + fineSTRUCTURE coancestry

**Script:** `scripts/37_finestructure.py`

## What was done

Haplotype-based population structure inference with ChromoPainter +
fineSTRUCTURE in linked mode (Lawson et al. 2012 *PLoS Genet*).

## Method

- Used the Beagle-phased VCF from script 36 (95 samples × 1,454 markers).
- Converted to per-chromosome ChromoPainter `.phase` files (190 haplotype
  rows + a position line) and uniform-rate `.recombfile` blocks at
  1.06 cM/Mb (cowpea proxy).
- Ran the full fs pipeline: parameter EM (stage 1), full painting and
  coancestry-matrix combination (stage 2), fineSTRUCTURE MCMC (stage 3),
  tree-building MCMC (stage 4). Pipeline driven from `/tmp/ayb_fs/` to
  avoid the path-with-spaces issue that breaks fs's stage-1 log routing.

## Inputs

- `results/36_ihs/data/ayb.phased.vcf.gz`

## Outputs

- `data/Ss*.phase`, `Ss*.recombfile`, `ayb.ids` — fs input
- `data/ayb_fs_linked.chunkcounts.out` — combined 95 × 95 coancestry
  (chunk-count) matrix from chromocombine
- `data/ayb_fs_linked_mcmc.xml` — finestructure MCMC trace
- `data/ayb_fs_linked_tree.xml` — MAP partition + dendrogram
- `data/ayb_fs_linked.newick` — Newick of the fs coancestry tree
- `tables/coancestry.csv` — coancestry matrix in CSV form
- `tables/cluster_assignments.csv` — sample → fs cluster (K = 24 MAP)
- `figures/fig81_coancestry_heatmap.png/.pdf` — log(1 + chunkcount)
- `figures/fig82_finestructure_tree.png/.pdf` — sample dendrogram

## Quick findings

- MAP partition: **K = 24 fineSTRUCTURE clusters** across the 95 samples
  — much finer resolution than the K = 4 admixture solution in script 09,
  consistent with the higher information content of haplotype-block-level
  inference vs. unlinked-genotype STRUCTURE.
- Coancestry matrix shows a block-diagonal pattern aligned with the fs
  tree, indicating real sub-structure rather than smooth admixture.

## Caveats

- DArTseq sparse markers (~500 kb median spacing) reduce the linkage
  signal that drives fineSTRUCTURE's haplotype-level resolution. Compare
  K against the K from STRUCTURE / DAPC for stability rather than
  treating either as authoritative.
- The uniform 1.06 cM/Mb recombination map is a coarse approximation.
  The LDhelmet result in script 39 will be used to validate whether
  recombination heterogeneity affects the fs clustering.
- Phasing accuracy is the upstream limit; long-range switch errors will
  shift haplotype-block boundaries and shrink apparent coancestry between
  truly-related samples.
