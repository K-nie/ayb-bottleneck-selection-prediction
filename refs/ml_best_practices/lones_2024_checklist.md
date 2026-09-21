# Lones 2024 — Avoiding common machine learning pitfalls

Structured checklist of the 34 pitfalls in Lones, M. A. (2024). *Avoiding
common machine learning pitfalls*. **Patterns** (Cell Press), DOI
10.1016/j.patter.2024.101046. arXiv v5: 2108.02497.

This file is the canonical reference for any ML-related agent in this
environment. Both `ml-genomics` and `ml-breeding` agents are required to
consult it before designing a model-comparison study.

The PDF lives alongside this file at `lones_2024_pitfalls.pdf`.

---

## How to use this checklist

Before writing or running any ML pipeline, walk through the 34 pitfalls
in order and write a short justification for each one — either "addressed
by X" or "N/A because Y". The justification should appear in the script
header docstring AND in the manuscript Methods section. Reviewers at Q1
journals increasingly cite this paper directly; matching its language
shortens the review cycle.

The pitfalls are organised by the five stages of an ML project.

---

## Stage 1 — Before you start to build models

### 2.1 Do think about how and where you will use data
Plan up front for **training**, **validation** (hyperparameter selection),
and **test** (final generalisation estimate) — these are typically three
disjoint subsets, not one train/test split.

### 2.2 Do take the time to understand your data
Treat exploratory data analysis as a precondition for model building.
Document data provenance, collection protocol, and any known limitations.
Garbage in → garbage out.

### 2.3 Don't look at all your data
Hold a portion of the data unseen until the very end. Looking at every
row during EDA risks untestable assumptions leaking into the model. The
test set is for testing, not for inspection.

### 2.4 Do clean your data
Inspect for duplicates, missing values, outliers, and inconsistent
records. Imputation must be done with care to avoid leakage (see 3.1).
Drop records that are clearly errors; keep records that are extreme but
plausible.

### 2.5 Do make sure you have enough data
Working out whether n is sufficient is hard before you build a model.
Options when n is small: data augmentation (after splitting — see 4.2),
transfer learning, or limit model complexity. Be explicit about the n
constraint in the paper.

### 2.6 Do talk to domain experts
Domain experts help select problems worth solving, the right feature
set, the right model, and the right journal. In plant breeding this
means working with the breeder, not just the bioinformatician.

### 2.7 Do survey the literature
Read the prior work on the same problem in the same domain. Cite it. Use
it as the baseline.

### 2.8 Do think about how your model will be deployed
For breeding: deployment = predicting GEBVs for new accessions in future
selection cycles. Model complexity, interpretability, and update
cadence must serve that workflow.

---

## Stage 2 — How to reliably build models

### 3.1 Don't allow test data to leak into the training process
The most common ML mistake in published papers. Forms of leakage:
- Scaling fit on the full dataset (do it per training fold).
- Imputation fit on the full dataset (do it per training fold).
- Feature selection on the full dataset (do it per training fold).
- Hyperparameter tuning on the test fold (use nested CV — the inner
  loop tunes, the outer loop tests).
- Test data appearing in training under any preprocessing or
  augmentation step.
- Time-series leakage (look-ahead bias — see 4.8).

### 3.2 Do try out a range of different models
No Free Lunch theorem: no single best model. Try a span — linear, tree
ensemble, kernel, neural — to find what fits the data's inductive biases.

### 3.3 Don't use inappropriate models
Match model to data type. Don't apply a classification model where
regression is the right framing. Don't apply a model that assumes
independence to data with structure (e.g., families, environments,
time).

### 3.4 Do keep up with progress in deep learning (and its pitfalls)
If you use DL: be aware of transformers, transfer learning, foundation
models, and the privacy / memorisation issues that come with them.

### 3.5 Don't assume deep learning will be the best approach
For tabular data with limited n, tree-based models (XGBoost, LightGBM)
typically outperform deep nets. Grinsztajn et al. 2022 NeurIPS is the
canonical citation for "tree-based models still beat DL on tabular".
**For plant genomic prediction at n < 500, DL has not been shown to beat
GBLUP on additive traits** (Bellot, de los Campos & Pérez-Enciso 2018
Genetics). State this explicitly when reporting DL results.

### 3.6 Do be careful where and how you do feature selection
Feature selection must happen **inside** the training fold of cross-
validation, never on the full dataset before splitting. The same rule
applies to dimensionality reduction (PCA loadings on training only).

### 3.7 Do optimise your model's hyperparameters
Use a real strategy — grid search, random search, Bayesian optimisation,
or AutoML. Always tune INSIDE an inner CV loop, never on the test fold.

### 3.8 Do avoid learning spurious correlations
Features correlated with the target but with no semantic meaning will
fool you. Inspect what the model actually learnt — SHAP, LIME,
permutation importance, ablation studies. Domain expertise is the
guardrail against spurious correlations.

---

## Stage 3 — How to robustly evaluate models

### 4.1 Do use an appropriate test set
The test set must not overlap with training, must be representative of
the deployment population, and (for related-sample data) must not share
identity with training samples. **For pedigreed plant panels, use
leave-one-family-out or kinship-blocked CV in addition to random
k-fold.**

### 4.2 Don't do data augmentation before splitting your data
Augment only the training set. If you augment the full set before
splitting, augmented copies of train samples can end up in the test
set — catastrophic leakage.

### 4.3 Do avoid sequential overfitting
Re-using the same test set across many model variants gradually
overfits the test set. Use a held-out validation set for model
selection, and reserve the test set for the final model.

### 4.4 Do evaluate a model multiple times
Train and evaluate stochastic models multiple times (different random
seeds, different fold partitions). Report mean ± CI, not point
estimates. Cross-validation with repeated random splits (≥ 20 repeats,
Lehermeier et al. 2014) gives stable variance estimates.

### 4.5 Do save some data to evaluate your final model instance
Cross-validation tells you about the model **family**. Performance of
the **specific instance** that will ship needs its own held-out test
set. Don't pick the best CV fold and report that as performance — it's
biased upward.

### 4.6 Do choose metrics carefully
- **Classification with class imbalance**: accuracy is misleading. Use
  F1, Cohen's κ, MCC.
- **Regression**: report Pearson r AND Spearman ρ AND RMSE AND MAE,
  not just one. Different metrics surface different model behaviours.
- **Breeder ranking**: top-k accuracy (overlap between predicted top-k
  and observed top-k) — Schopp et al. 2017 G3.
- **Genomic prediction specifically**: compare against the Daetwyler
  ceiling = sqrt(n × h² / (n × h² + M_e)) — any r above the ceiling
  is suspect (Daetwyler et al. 2008 PLoS ONE; Goddard 2009 Genetica).

### 4.7 Do consider model fairness
For models deployed on human subjects (medical, hiring, financial,
criminal-justice ML), test for performance parity across protected
groups. Less directly relevant for plant ML, but the analogue is
performance parity across **subpopulations** (clusters, geographic
provenance, breeding programmes).

### 4.8 Don't ignore temporal dependencies in time series data
For longitudinal phenotype data, year-blocked CV; never random k-fold
on multi-year trials. Look-ahead bias is catastrophic.

---

## Stage 4 — How to compare models fairly

### 5.1 Don't assume a bigger number means a better model
A higher headline r doesn't make a model better if (a) the test sets
differ, (b) the hyperparameter optimisation effort differs, or (c) the
data preprocessing differs. Cross-paper comparisons require
re-running competitor methods on the same data.

### 5.2 Do use meaningful baselines
- A **dummy** baseline that predicts the mean (regression) or the
  majority class (classification) — your model should beat this trivially.
- A **simple** baseline (linear regression, logistic regression,
  RR-BLUP for genomic prediction).
- A **state-of-the-art** baseline (the published method you claim to
  beat).
If your model only beats the dummy and not the simple baseline, that's
the result — report it.

### 5.3 Do use statistical tests when comparing models
For paired predictions on the same test instances: McNemar's test
(classification) or paired t-test / Wilcoxon signed-rank (regression).
For independent samples: Mann-Whitney U if non-normal. Carrasco et al.
2020 surveys recent recommendations.

### 5.4 Do correct for multiple comparisons
Every pairwise model comparison consumes one test from your family-wise
error budget. Use Bonferroni, Holm-Bonferroni, or Benjamini-Hochberg
FDR. Don't report 20 uncorrected pairwise p < 0.05 values as "significant".

### 5.5 Don't always believe results from community benchmarks
Repeated use of the same public benchmark by the community gradually
overfits the benchmark. Foundation models may have been trained on the
benchmark's test data.

### 5.6 Do combine models (carefully)
Ensembles (bagging, boosting, stacking) often beat single models.
Watch for leakage: the test data used to evaluate the ensemble must
not have been used to train any of its components.

---

## Stage 5 — How to report your results

### 6.1 Do be transparent
Share the scripts. Share the data (or document why you can't). Share
the hyperparameter grids. Share the random seeds. Reproducibility
checklists (REFORMS, Kapoor et al. 2024) are increasingly mandatory.

### 6.2 Do report performance in multiple ways
Multiple datasets, multiple metrics, multiple slices of the test set.
A single headline number obscures trade-offs.

### 6.3 Don't generalise beyond the data
Claims must be supported by the data. "Our model improves genomic
prediction in AYB" is supported. "Our model improves genomic
prediction in legumes" requires multi-species testing.

### 6.4 Do be careful when reporting statistical significance
- Don't worship the p < 0.05 threshold. Report exact p AND effect size.
- Cohen's d (parametric) or Cliff's δ / Kolmogorov-Smirnov (robust).
- Bayesian alternatives where natural.
- Significance ≠ importance.

### 6.5 Do look at your models
Don't ship a model without inspecting what it learnt. SHAP, LIME,
permutation importance, partial-dependence plots, ablation studies.
For trees, look at split conditions. For linear models, look at
coefficients.

### 6.6 Do use a machine learning checklist
The REFORMS checklist (Kapoor, Lones et al. 2024 Science Advances) is
the consensus-based version for ML-in-science. Other domain checklists
exist; consult them for medical / clinical / regulatory work.

---

## Quick lookup — pitfalls especially relevant for plant-breeding ML

For genomic prediction in a low-n plant breeding panel (n < 500
typically), these are the high-impact items:

1. **3.1 — Test-data leakage via per-fold scaling, imputation, and
   feature selection.** Almost every published plant-ML paper has at
   least one leakage issue.
2. **3.5 — Don't assume DL will win.** At n < 500, GBLUP / RR-BLUP
   usually wins or ties for additive traits. State this honestly.
3. **4.1 — Use a kinship-aware CV scheme.** Random k-fold on a
   pedigreed panel leaks family information across the fold boundary.
   Leave-one-family-out or clone-blocked CV is the honest comparison.
4. **4.5 — Hold out a final test set.** Don't report the best CV fold;
   report a held-out instance evaluation.
5. **4.6 — Daetwyler ceiling as the upper bound.** Any model that
   "beats" the ceiling is overfit or has leaked. Always compare ML
   models against this theoretical anchor.
6. **5.2 — RR-BLUP as the baseline.** This is the field's default
   genomic-prediction model since Meuwissen, Hayes & Goddard 2001
   Genetics. Don't introduce a new ML method without comparing
   against it.
7. **5.3 + 5.4 — Wilcoxon signed-rank on paired CV r with
   Holm-Bonferroni correction.** Every reviewer at Plant Genome or
   TAG now expects this.
8. **6.4 — Effect size on Δr.** A 0.02 improvement in r over GBLUP
   may be statistically significant but breeder-irrelevant. Report
   Cohen's d or the bootstrap 95 % CI on Δr.

---

## Plant-breeding-specific extensions (not in Lones)

Beyond the 34 Lones pitfalls, these are field-specific best practices
the agent should add when working on plant-breeding ML:

- **Population structure**: report performance separately for each
  fineSTRUCTURE / ADMIXTURE cluster if the population is structured.
- **G × E interaction**: if multi-environment data, factor-analytic
  mixed models (Smith, Cullis & Thompson 2005) outperform single-
  environment ML in most cases. Don't replace; complement.
- **Pedigree-aware models**: ssGBLUP (Misztal et al. 2009) when both
  genotyped and ungenotyped accessions exist.
- **Heritability bounds**: ML accuracy cannot exceed sqrt(h²)
  asymptotically. Report each trait's h² CI from REML profile
  likelihood (Wilks 1938).
- **Multi-trait models**: when traits are genetically correlated
  (off-diagonals in the genetic-correlation matrix far from 0),
  multi-trait GBLUP / multi-task ML can help. Test honestly: the
  improvement should be larger than the within-trait gain from
  ensembling.
- **Breeder direction**: report top-k accuracy weighted by the
  breeder's selection intensity (top-1 %, top-5 %, top-10 %),
  not just default top-10.
- **GEBV calibration**: predicted GEBVs should be regressed on
  observed phenotypes; the slope should be near 1 if the model is
  calibrated (Legarra & Reverter 2018 GSE).

---

## References (consensus must-cites for plant-ML papers)

- Lones MA (2024). Avoiding common machine learning pitfalls. *Patterns* 5, 101046.
- Kapoor S, Cantrell EM, ..., Lones MA, ..., Narayanan A (2024). REFORMS: Consensus-based recommendations for machine-learning-based science. *Science Advances* 10, eadk3452.
- Meuwissen THE, Hayes BJ, Goddard ME (2001). Prediction of total genetic value using genome-wide dense marker maps. *Genetics* 157, 1819–1829.
- Habier D, Fernando RL, Dekkers JCM (2007). The impact of genetic relationship information on genome-assisted breeding values. *Genetics* 177, 2389–2397.
- VanRaden PM (2008). Efficient methods to compute genomic predictions. *J Dairy Sci* 91, 4414–4423.
- Daetwyler HD, Villanueva B, Woolliams JA (2008). Accuracy of predicting the genetic risk of disease using a genome-wide approach. *PLoS ONE* 3, e3395.
- Goddard ME (2009). Genomic selection: prediction of accuracy and maximisation of long term response. *Genetica* 136, 245–257.
- Heslot N, Yang H-P, Sorrells ME, Jannink J-L (2012). Genomic selection in plant breeding: a comparison of models. *Crop Sci* 52, 146–160.
- Lehermeier C, Schön C-C, de los Campos G (2014). Assessment of genetic heterogeneity in structured plant populations using multivariate whole-genome regression models. *Genetics* 198, 3–16.
- Bellot P, de los Campos G, Pérez-Enciso M (2018). Can deep learning improve genomic prediction of complex human traits? *Genetics* 210, 809–819.
- Montesinos-López OA, Montesinos-López A, Crossa J et al. (2018, 2019). Multi-trait, multi-environment deep learning modelling for genomic-enabled prediction of plant traits. *G3*, *The Plant Genome*, *Heredity*.
- Sandhu KS, Mihalyov PD, Lewien MJ, Pumphrey MO, Carter AH (2021). Multi-trait machine and deep learning models for genomic selection using spectral information in a wheat breeding program. *The Plant Genome* 14, e20119.
- Wang K, Abid MA, Rasheed A, Crossa J, Hearne S, Li H (2022). DNNGP, a deep neural network-based method for genomic prediction. *Plant Methods* 18, 24.
- Lundberg SM, Lee S-I (2017). A unified approach to interpreting model predictions. *NeurIPS* 30.
- Lundberg SM, Erion G, Chen H et al. (2020). From local explanations to global understanding with explainable AI for trees. *Nat Mach Intell* 2, 56–67.
- Grinsztajn L, Oyallon E, Varoquaux G (2022). Why do tree-based models still outperform deep learning on typical tabular data? *NeurIPS* 35.
- Legarra A, Reverter A (2018). Semi-parametric estimates of population accuracy and bias of predictions of breeding values and future phenotypes using the LR method. *Genet Sel Evol* 50, 53.
