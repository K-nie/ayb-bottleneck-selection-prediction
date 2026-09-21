# 35 — Windowed Tajima's D, π, and Watterson's θ

**Script:** `scripts/35_tajima_d.py`

## What was done

Windowed scan of three classical neutrality / diversity statistics across
the AYB-anchored DArTseq markers: nucleotide diversity π (Nei & Li 1979),
Watterson's θ_W (Watterson 1975), and Tajima's D (Tajima 1989).

## Method

- 250 kb windows with 100 kb step across Ss01–Ss11.
- Per window with ≥ 4 segregating SNPs:
  - π_window = Σ over sites of (2 × c_alt × (n − c_alt)) / (n × (n − 1)),
    where c_alt is the alt-allele count and n = 2 × n_samples = 190.
  - θ_W = S / a₁ with a₁ = Σ_{i=1}^{n−1} 1/i.
  - Tajima's D = (π − θ_W) / √Var, using Tajima 1989's c₁, c₂, e₁, e₂
    variance expansion.

## Inputs

- HapMap dosage matrix (post-QC, MAF ≥ 0.05)
- AYB marker anchoring table → chromosome + position per marker

## Outputs

- `tables/tajima_d_windowed.csv` — chrom, win_start, win_end, mid_bp,
  n_snps, π_total, π_per_site, θ_W, tajima_D
- `figures/fig77_tajima_d_manhattan.png/.pdf` — Manhattan-style scan with
  the ±2 reference lines
- `figures/fig78_pi_vs_theta.png/.pdf` — per-window π vs θ_W (the geometry
  underlying Tajima's D)

## Quick findings

- 255 windows passed the SNP threshold.
- Median D = +1.895 across windows.
- 119 windows have D > +2 (consistent with balancing selection, population
  contraction, or both); 0 windows with D < −2.

## Caveats

- DArTseq is ascertained on common-variant SNPs, which inflates π relative
  to θ_W and pushes Tajima's D upward in a way that does not reflect
  evolutionary process. The genome-wide positive D should not be over-
  interpreted on absolute scale; the relative distribution across windows
  is more meaningful.
- Marker density (~500 kb mean spacing) limits resolution to chromosome-
  scale signal. Sub-windowed peaks are not interpretable at this density.
