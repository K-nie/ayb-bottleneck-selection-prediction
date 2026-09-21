# 88 — Orientation figure (paper B, Figure B1): SNP-PCA + ADMIXTURE K=3

## What this does and why
Builds the reader-orientation figure for paper B: panel-level population
structure of the 95-line *Sphenostylis stenocarpa* DArTseq panel. It answers
"how is this panel structured?" in one view before the demographic and
selection analyses that follow. Two panels:

- **A** — SNP-PCA scatter (PC1 vs PC2), each accession coloured by its majority
  ADMIXTURE K=3 ancestry component. Ties the K=2 split the panel reads on PC1 to
  the K=3 ADMIXTURE partition.
- **B** — ADMIXTURE K=3 stacked-ancestry bar plot, one bar per accession,
  ordered by majority component then dominant ancestry fraction.

## Key decisions / caveats (reviewer-facing)
- **K=3, not K=4.** The B1 manifest row originally specified a K=4 STRUCTURE
  bar plot. The ADMIXTURE cross-validation error is minimised at **K=3**
  (`results/09_admixture/tables/admixture_cv_error.csv`: CV = 0.5143 at K=3 vs
  0.52122 at K=4, 0.52143 at K=2), and the manuscript body reports the panel
  structure as "K=2 SNP-PCA / K=3 ADMIXTURE". A K=4 bar plot would contradict
  both the CV result and the manuscript, so the figure uses the CV-optimal K=3.
- **No sample map.** The manifest row also specified a geographic sample map.
  This is a single IITA TSs genebank collection with **no per-accession
  latitude/longitude** in any project data file (Genesys-PGR holds no
  *Sphenostylis* records). No coordinates were fabricated; the PCA panel serves
  as the spatial orientation instead.
- Majority-ancestry group sizes at K=3: cluster 1 = 29, cluster 2 = 12,
  cluster 3 = 54 accessions. Cluster 2 is the PC1-negative outlier set.

## Inputs (provenance)
- `results/01_qc_pca_power/tables/pca_coords.csv` — 95 accessions × PC1–PC10,
  from the QC + PCA step on the 95 × 1,625 QC-filtered marker matrix.
- `results/01_qc_pca_power/tables/pca_variance.csv` — per-PC variance explained
  (PC1 = 8.48%, PC2 = 5.57%).
- `results/09_admixture/tables/admixture_Q_K3.csv` — ADMIXTURE K=3 Q-matrix
  (95 accessions), from `scripts/09_admixture.py` (ADMIXTURE on the same panel).
- `results/09_admixture/tables/admixture_cv_error.csv` — K-selection CV error.

## Reproduce
```bash
conda run -n yeast-viz python scripts/88_orientation_pca_admixture.py
```
No parameters or thresholds beyond majority-component assignment
(argmax over the three Q columns).

## Outputs
- `figures/fig_orientation_pca_admixture.{png,pdf}` — the composite (300 dpi).
- `tables/sample_ancestry_assignment.csv` — per-accession PC1/PC2, Q1–Q3, and
  majority-cluster assignment.

Manuscript placement: **Figure B1 (supplementary)**, per
`manuscript/paper_B/figure_manifest_B.csv`. Share copy synced to
`manuscript/paper_B/share/figures/supp/`.

## Software versions
- ADMIXTURE Q-matrix produced upstream by `scripts/09_admixture.py`.
- Figure: matplotlib via the `yeast-viz` conda env; shared style
  `scripts/_figstyle.py` (Wong colour-vision-friendly palette, 300 dpi PNG+PDF).
