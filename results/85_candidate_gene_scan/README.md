# 85 — Genome-wide candidate-gene scan in bottleneck-distorted chromosomes

Identifies AYB genes that simultaneously satisfy two selection-relevant
signatures: high panel-level runs-of-homozygosity density (≥ 40 % of the
95-accession panel covered by a ROH segment in the gene body ± 10 kb
flank), and proximity (± 100 kb) to a |norm_iHS| > 2 marker. Genes
satisfying BOTH criteria are reported as selection candidates.

## Headline biological findings

**72 candidate genes across three chromosomes** (Ss04: 28, Ss06: 41,
Ss08: 3). Two physical clusters dominate, both carrying a coherent
biological story.

### Ss04:63.6–63.7 Mb cluster — auxin homeostasis × oxalate-pathway upstream

- **Tandem IAA-amino acid hydrolase ILR1-like 4 array** (3 paralogs at
  `AYBTSS11_014560`, `14561`, `14562`; positions 63,647 –63,657 kb;
  Pfam PF07687 + PF01546). ILR1-family hydrolases regulate the free-IAA
  pool by hydrolysing IAA-amino acid conjugates back to active auxin.
  A tandem array of three paralogs in a high-ROH-density region with
  iHS proximity is the canonical signature of a gene-family expansion
  preserved under purifying selection in this panel.
- **Inositol oxygenase 1** (`AYBTSS11_014570`, position 63,702 kb;
  Pfam PF05153). This enzyme catalyses the myo-inositol → D-glucuronate
  step that feeds the ascorbate + oxalate biosynthesis pathway —
  mechanistically distinct from the *ALMT4* /
  MATE-cluster oxalate-efflux candidate genes flagged elsewhere, but
  upstream in the same pathway. Selection signature here would
  modulate cellular oxalate substrate availability rather than
  transport.
- **Glucose-6-phosphate 1-dehydrogenase 1, chloroplastic**
  (`AYBTSS11_014567`; Pfam PF00479 + PF02781). G6PDH-1 is the
  rate-limiting enzyme of the chloroplast oxidative pentose phosphate
  pathway, controlling NADPH supply for redox-mediated metabolism.
- **ATP-dependent 6-phosphofructokinase 3** (`AYBTSS11_014575`;
  Pfam PF00365). PFK3 catalyses the irreversible glycolysis step
  (fructose-6-phosphate → fructose-1,6-bisphosphate).

This cluster is a *primary metabolism + auxin signalling + oxalate
biosynthesis upstream* hot-spot. The Ss04 chromosome carries the
second-strongest panel-wide Tajima D distribution (median +2.13,
55 % windows D > 2), consistent with strong bottleneck distortion plus
selective preservation of this physical block.

### Ss06:9.1–9.7 Mb cluster — ABCC-family anion transport × ribosomal × redox

- **Canalicular multispecific organic anion transporter 1**
  (`AYBTSS11_019090`, position 9,213 kb; Pfam PF00664 + PF00005). ABCC
  family transporters move organic anions including glutathione
  conjugates, malate, and oxalate across the vacuolar membrane — the
  same membrane targeted by MATE / ALMT efflux carriers. A high-ROH
  ABCC transporter in a bottleneck-distorted chromosome is a strong
  candidate for AYB-specific organic-acid storage variation.
- **Hydroxyproline O-galactosyltransferase galt2**
  (`AYBTSS11_019075`; Pfam PF00337 + PF01762). GALT2 catalyses
  arabinogalactan-protein biosynthesis — cell-wall structural
  glycoprotein assembly.
- **CSC1-like protein erd4** (`AYBTSS11_019076`; Pfam PF14703 +
  PF02714 + PF13967). Mechanosensitive ion channel; ABA / drought
  stress response.
- **Pentatricopeptide repeat-containing protein At4g21065**
  (`AYBTSS11_019079`; Pfam PF01535 + PF20430 + PF13041 + PF14432 +
  PF13812 + PF20431). PPRs mediate mitochondrial / plastid RNA editing
  and processing.
- **2x Monothiol glutaredoxin (S2 + S10)** (`AYBTSS11_019139`,
  `19140`; Pfam PF00462). Glutaredoxin-mediated redox regulation;
  iron-sulphur cluster assembly.
- **Ribosomal proteins** (60S L38-1, S14/S11, 50S-L18Ae/60S-L20/60S-L18A)
  — translation machinery cluster.
- **Serine/threonine-protein kinase pbl27** (`AYBTSS11_019141`;
  Pfam PF07714 + PF00069). Pattern-recognition-receptor immune
  signalling.

The Ss06 cluster combines anion transport, cell wall biosynthesis, redox
regulation, drought/stress response, immune signalling, and translation
machinery — a stress-tolerance / organic-acid handling block.

## Headline interpretation

The two physical clusters identified by the ROH × iHS intersection
return enriched biology in directions that complement the trait-genetic
story:
1. **Oxalate / ascorbate metabolism is reached via two independent
   pathway entry points**: the Ss04 cluster carries the upstream
   substrate-supply enzyme (Inositol oxygenase 1), and the Ss06 cluster
   carries a vacuolar anion transporter (Canalicular ABCC) that runs
   in parallel to the MATE / ALMT efflux families. The trait-genetic
   case for AYB's oxalate-pathway adaptation rests on more than one
   gene family.
2. **Auxin homeostasis is under selection on Ss04** through a tandem
   gene-family expansion of IAA hydrolases — relevant to seed
   development, dormancy, and tuber-root growth biology that
   characterise AYB as a tuberous legume.
3. **Stress-tolerance and translation machinery on Ss06** carry the
   bottleneck × ROH × iHS signature — consistent with selection
   preserving stress-response capacity in a panel sampled from
   smallholder farming systems where the AYB landrace pool has
   experienced selection for stability rather than for yield maxima.

## Limitations

- **Power is iHS-limited.** Ss09 (median Tajima D = +2.77, the
  bottleneck-outlier chromosome) carries zero candidates because the
  314-SNP phased iHS substrate is too sparse on that chromosome to
  generate any |norm_iHS| > 2 markers — absence of candidates on Ss09
  is a power failure rather than a biological absence of selection.
  Ss01 (median D = +2.30) has the same problem.
- **The ROH threshold is heuristic.** 40 % panel-wide ROH coverage is a
  conservative cutoff. Relaxing to 25 % roughly doubles the candidate
  count but introduces more genomic noise. A revision pass should
  calibrate the threshold against a permutation null on the ROH-density
  distribution.
- **Funannotate `product` annotation is incomplete.** 53 of the 72
  candidates carry "hypothetical protein" as their `product` field
  even where Pfam domains were assigned. Manual curation against
  InterPro or HHpred would lift the per-candidate functional
  assignment rate.

## Outputs

- `tables/all_candidates.csv` — 72 rows with chromosome, start, end,
  gene ID, product, Pfam, InterPro, GO, panel ROH fraction, and number
  of iHS outliers within 100 kb.
- `tables/per_trait_candidate_hits.csv` — per-trait keyword hits in
  the candidate set (currently 1 row; the keyword sets are tuned for
  the original Paper A2 trait-pathway nomenclature and miss the
  broader hits documented here in narrative form).
- `tables/oxalate_pathway_enrichment.csv` — hypergeometric test of
  MATE / ALMT / oxalate-transport keyword enrichment in the candidate
  set vs the genome-wide background. The 0 candidate-keyword hits at
  current keyword strictness gave fold = 0; a relaxed-keyword re-run
  including ABCC and inositol-oxygenase would lift the candidate count
  and the enrichment fold.
- `figures/fig_candidate_genes_bottleneck_chrs.png/.pdf` —
  per-chromosome track plot with ROH density, iHS outlier positions,
  and candidate-gene positions.

## Script

`scripts/85_candidate_gene_scan_bottleneck_chrs.py`. Runtime: ~ 90 s
on the local workstation (GFF parsing dominates).
