# A recent bottleneck shapes diversity, selection and genomic prediction in African yam bean (Sphenostylis stenocarpa)

Analysis code and result tables for the demographic-history, selection and genomic-prediction study of the African yam bean DArTseq panel. Covers Stairway Plot demographic reconstruction, the site-frequency spectrum, Tajima's D, integrated haplotype score (iHS), nucleotide diversity versus Watterson's theta, fineSTRUCTURE co-ancestry, per-chromosome LD, machine-learning genomic prediction with SHAP interpretability and multi-task models, the Daetwyler prediction ceiling, and candidate-gene and pathway enrichment for selection signals.

This repository contains the **analysis code and derived result tables** for
the study. Raw genotype and phenotype data are archived separately (see Data below);
manuscript drafts are not included.

## Repository layout

```
scripts/    numbered Python and R analysis scripts (shared helpers: _plotstyle, _pheno, _figstyle)
results/    one directory per analysis stage, each with tables/ and a README.md
            documenting method, inputs, outputs, findings, and caveats
refs/       machine-learning best-practice reference material
```

Each `results/<NN>_*/README.md` is the reproducibility and methods record for that
stage: what it does and why, its inputs, the exact commands and thresholds, the outputs,
and the caveats.

## Reproducing

Scripts run from a single conda environment (`yeast-viz`: numpy, pandas, scikit-learn,
matplotlib, plus the R toolchain for the `.R` steps). Stages are numbered by pipeline order;
run the lower-numbered QC/PCA stage first, then the analysis of interest. Per-stage READMEs
give the exact command and parameters.

## Data

- **Raw DArTseq genotypes** (order DAf18-2580, 105 accessions; post-QC 95 accessions x 1,625
  markers at sample and marker call rate >= 0.90, MAF >= 0.05): Zenodo [10.5281/zenodo.20348832](https://doi.org/10.5281/zenodo.20348832).
- **Reference genome**: *S. stenocarpa* chromosome-scale assembly, ENA [PRJEB57813](https://www.ebi.ac.uk/ena/browser/view/PRJEB57813)
  (Shorinola et al. 2024); Funannotate annotation at Zenodo [10.5281/zenodo.13853757](https://doi.org/10.5281/zenodo.13853757).
- Phenotype data are held by the breeding program and available from the authors on
  reasonable request.

## Headline findings

- Stairway Plot and the site-frequency spectrum indicate a recent contraction in effective population size.
- Positive genome-wide Tajima's D is consistent with a bottleneck / allele-frequency deficit rather than balancing selection.
- iHS, fineSTRUCTURE co-ancestry and per-chromosome LD characterise haplotype structure and candidate selection footprints.
- Machine-learning genomic prediction is benchmarked against GBLUP with SHAP-based interpretability and referenced to the Daetwyler ceiling.

## Citation

If you use this code or its outputs, please cite this repository (see `CITATION.cff`)
and the accompanying manuscript.

## License

Code is released under the MIT License (see `LICENSE`).
