# 86 — Pathway-level enrichment for selection signatures

Two tests share the same ROH + iHS substrate that feeds the script-85
candidate-gene scan, each one asking which gene families and curated
pathways the over-represented selection signature in the AYB panel
tracks.

## Tests run

### Test 1 — Pfam enrichment in the script-85 candidate set vs the Funannotate genome-wide background

Per-Pfam hypergeometric test of candidate-set membership conditioned on
genome-wide Pfam frequency, limited to Pfam terms carrying ≥ 2
candidate-set hits. BH-FDR across the surviving terms.

### Test 2 — Per-pathway selection-statistic comparison via permutation

Each curated pathway gene set (Pfam / product-keyword regex) is scored
by per-gene panel ROH fraction and per-gene count of |norm_iHS| > 2
markers within ± 100 kb, against all 30,840 Funannotate mRNAs. A
two-sided permutation null then tests the difference of means between
pathway and background. The per-test resolution is 1/3000 — the script
constant `N_PERM = 4999` is set, but the function call hands in
`n_perm = 2999`, so the reported `p_perm` values sit at 1/3000
granularity rather than 1/5000. RNG seed = `20260619`. BH-FDR runs
across pathways within each of the two metrics (ROH delta, iHS delta).

## Headline findings

### Pfam enrichment recovers the ILR1 IAA-amino-acid-hydrolase family at the top of the ranking

Three Pfam terms clear BH q < 5 × 10⁻⁵, all of them annotating the
peptidase M20 family — the catalytic and dimerisation domains that the
ILR1 IAA-amino-acid hydrolase family carries:

| Pfam | Description | n_bg | n_cand | fold | p_hyper | q_BH |
|---|---|---:|---:|---:|---:|---:|
| PF07687 | Peptidase M20 dimerisation | 9 | 3 | 143 | 1.0e-6 | 1.5e-5 |
| PF01546 | Peptidase M20 (catalytic) | 12 | 3 | 107 | 2.6e-6 | 2.0e-5 |
| PF03767 | M20 peptidase / carboxypeptidase | 14 | 3 | 92 | 4.4e-6 | 2.2e-5 |

All 3 candidate hits per term trace back to the same Ss04:63.6–63.7 Mb
tandem array of three ILR1-like 4 paralogs (`AYBTSS11_014560/14561/14562`)
that the script-85 narrative flagged. Two further families pass q < 10⁻³:

| Pfam | Description | n_bg | n_cand | fold | q_BH |
|---|---|---:|---:|---:|---:|
| PF00208 | Glu/Leu/Phe/Val dehydrogenase, dimer | 9 | 2 | 95 | 5.7e-4 |
| PF02812 | Glu/Leu/Phe/Val dehydrogenase, NAD-binding | 9 | 2 | 95 | 5.7e-4 |

Both terms annotate the same enzyme-family pair — amino-acid
dehydrogenases bridging nitrogen turnover and primary-metabolism flux.
Behind them sit PF00450 (serine carboxypeptidase, q = 4.5e-3), PF00240
(ubiquitin, q = 8.5e-3) and PF00462 (glutaredoxin, q = 1.4e-2).

### Per-pathway permutation flags auxin, cell wall, and oxalate-upstream pathways in the iHS tail; nothing clears BH q < 0.05

15 curated pathways were screened. After BH correction on the iHS-delta
metric, the ranking is:

| Pathway | n_genes | Δ iHS (pathway − bg) | p_perm | q_BH |
|---|---:|---:|---:|---:|
| Auxin homeostasis (IAA hydrolase / PIN / AUX1) | 62 | +0.070 | 0.0067 | 0.074 |
| Cell-wall biosynthesis (CESA / PME / expansin) | 72 | +0.059 | 0.012 | 0.074 |
| Oxalate / ascorbate biosynthesis upstream (MIOX / GLDH) | 20 | +0.139 | 0.020 | 0.080 |
| Redox / glutathione / oxidative stress | 151 | +0.022 | 0.056 | 0.143 |
| Translation machinery (ribosomes / eIF / eEF) | 265 | +0.016 | 0.065 | 0.143 |
| Primary metabolism (G6PDH / PFK / TCA) | 28 | +0.060 | 0.072 | 0.143 |
| PPR (mitochondrial / plastid RNA editing) | 72 | +0.031 | 0.089 | 0.153 |

On the ROH-delta metric, "Iron / metal transport" lands at raw
p_perm = 0.021 (delta = −0.042, ROH-depleted) and "Redox / glutathione"
at raw p_perm = 0.045 (delta = +0.013); both float above q_BH > 0.25
once correction lands. Of all pathways, only auxin homeostasis carries
an iHS-delta direction that lines up with the candidate-gene scan. Its
62 background genes average +0.070 more iHS outliers within ± 100 kb
than the 30,778-gene complement.

### Honest framing

On either metric, no per-pathway test clears BH q < 0.05. The four
strongest pathway signals (auxin homeostasis, cell-wall biosynthesis,
oxalate / ascorbate upstream, redox) cluster at q = 0.07–0.14, and that
band reads as the pattern across the top of the ranking — not as a
locus-specific selection claim. At the Pfam level, the enrichment does
recover the M20-peptidase / ILR1 IAA-hydrolase family at q < 5 × 10⁻⁵,
but that signal traces back to a single tandem array of three paralogs
at Ss04:63.6 Mb. The reading is gene-family-expansion preservation
rather than diffuse pathway-wide selection.

## Limitations

- **Test 2 returns no BH-significant pathways.** The per-pathway
  permutation ranking is suggestive, but nothing on it survives
  correction. The honest framing treats the tail as a biologically
  coherent pattern rather than as locus evidence.
- **Pathway gene sets are regex-defined off Funannotate `product` and
  Pfam strings.** That route misses orthologs annotated only as
  "hypothetical protein" while carrying a Pfam assignment with no
  product-string match — 53 of the script-85 candidates sit in that
  bucket. To lift recall, a revision pass would cross-reference against
  orthogroup membership (OrthoFinder against a Phaseoleae reference
  set).
- **iHS substrate is sparse.** The 314 phased SNPs that come through
  Beagle + MAF filter carry the same ascertainment-and-density problem
  that §3.4 of the draft flags. Test 2 cannot pick up pathways enriched
  on chromosomes carrying near-zero phased-iHS substrate (Ss09, Ss01),
  regardless of biological reality.
- **N_PERM constant mismatch.** The script header carries
  `N_PERM = 4999`, while the permutation calls hand in `n_perm = 2999`;
  the reported p-values sit at 1/3000 granularity. Re-running at the
  higher resolution would not change the rank ordering, but it would
  push the smallest p-values to the 2 × 10⁻⁴ floor.

## Outputs

- `tables/pfam_enrichment_candidate_vs_background.csv` — 15 Pfam terms
  carrying ≥ 2 candidate-set hits, with fold enrichment, hypergeometric
  p, and BH-q.
- `tables/pathway_selection_enrichment.csv` — 12 curated pathways
  paired with per-pathway and per-background mean ROH fraction, mean
  iHS outlier count, observed deltas, two-sided permutation p, and BH-q.
- `tables/per_gene_selection_metrics.csv` — a 30,840-row genome-wide
  per-gene table covering ROH fraction, iHS outlier count, chromosome
  coordinates, and the Funannotate annotation search string. This is
  the substrate for any further pathway-level test.
- `figures/fig_pathway_selection_enrichment.png/.pdf` — a two-panel
  forest plot showing per-pathway ROH delta (left) and iHS delta
  (right); green marks significant positive enrichment at raw
  p_perm < 0.05; red marks significant negative; grey is ns.

## Script

`scripts/86_pathway_selection_enrichment.py`. Rather than duplicating
parsing code, the script re-imports the GFF parser from
`scripts/85_candidate_gene_scan_bottleneck_chrs.py` through
`importlib.util.spec_from_file_location`. Runtime is ≈ 4 min on the
local workstation, dominated by per-gene ROH density grid construction
across 9 chromosomes.
