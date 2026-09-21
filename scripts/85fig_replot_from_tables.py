"""
85fig -- house-style replot of B29, the per-chromosome candidate-gene track plot
across the bottleneck-distorted chromosomes, from saved tables.

The parent script (85_candidate_gene_scan_bottleneck_chrs.py) parses the ~30 s
Funannotate GFF to call the candidate genes, then writes all_candidates.csv. The
figure only needs the three panel layers -- panel ROH density, iHS outlier
positions, and the candidate-gene call-outs -- all of which are recoverable from
the frozen CSVs without re-parsing the GFF:

  panel ROH density  recomputed per chromosome from roh_calls.csv (cheap)
  iHS outlier lines  |norm_iHS| > 2 markers from ihs_combined.csv (cheap)
  candidate genes    gene-body midpoints from all_candidates.csv

One stacked panel per bottleneck chromosome: filled ROH-density track (Wong
blue), the 0.40 panel-fraction threshold line, iHS outlier ticks along the top,
and candidate-gene markers.

Inputs:
  results/34_roh/tables/roh_calls.csv
  results/36_ihs/tables/ihs_combined.csv
  results/85_candidate_gene_scan/tables/all_candidates.csv
Output:
  results/85_candidate_gene_scan/figures/fig_candidate_genes_bottleneck_chrs.{png,pdf}
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, save_figure
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
ROH = PROJ / "results/34_roh/tables/roh_calls.csv"
IHS = PROJ / "results/36_ihs/tables/ihs_combined.csv"
CAND = PROJ / "results/85_candidate_gene_scan/tables/all_candidates.csv"
OUT = PROJ / "results/85_candidate_gene_scan/figures"
OUT.mkdir(parents=True, exist_ok=True)

BOTTLENECK_CHRS = ["Ss09", "Ss04", "Ss01", "Ss10", "Ss06", "Ss08", "Ss02"]
ROH_PANEL_FRAC_THRESHOLD = 0.40
IHS_THRESHOLD = 2.0
FLANK_KB = 100
BIN_KB = 50


def roh_density(roh: pd.DataFrame, chrom: str, chr_max_bp: int,
                bin_kb: int = BIN_KB) -> pd.Series:
    n_samples = roh["sample"].nunique()
    sub = roh[roh["chr"] == chrom]
    bin_bp = bin_kb * 1000
    n_bins = chr_max_bp // bin_bp + 1
    seen = {b: set() for b in range(n_bins)}
    for _, row in sub.iterrows():
        s_bin, e_bin = int(row["start_bp"]) // bin_bp, int(row["end_bp"]) // bin_bp
        for b in range(s_bin, min(e_bin + 1, n_bins)):
            seen[b].add(row["sample"])
    counts = np.array([len(seen[b]) for b in range(n_bins)])
    return pd.Series(counts / n_samples, index=np.arange(n_bins) * bin_bp,
                     name="roh_panel_frac")


def main() -> None:
    roh = pd.read_csv(ROH)
    ihs = pd.read_csv(IHS)
    cand = pd.read_csv(CAND)
    print(f"ROH calls {len(roh)}; iHS SNPs {len(ihs)}; candidates {len(cand)}")

    chrs = [c for c in BOTTLENECK_CHRS
            if (roh["chr"] == c).any() or (ihs["chr"] == c).any()]
    n_chr = len(chrs)
    fig, axes = plt.subplots(n_chr, 1, figsize=(11.0, 1.9 * n_chr), sharex=False)
    if n_chr == 1:
        axes = [axes]

    for ax, chrom in zip(axes, chrs):
        ihs_chr = ihs[ihs["chr"] == chrom]
        roh_chr = roh[roh["chr"] == chrom]
        cand_chr = cand[cand["chr"] == chrom]
        chr_max_bp = int(max(
            roh_chr["end_bp"].max() if len(roh_chr) else 0,
            ihs_chr["pos"].max() if len(ihs_chr) else 0,
            cand_chr["end"].max() if len(cand_chr) else 0,
        )) + 500_000
        dens = roh_density(roh, chrom, chr_max_bp)
        pos_mb = dens.index / 1e6

        ax.fill_between(pos_mb, 0, dens.values, color=WONG["blue"], alpha=0.35,
                        lw=0, label="Panel ROH density")
        ax.axhline(ROH_PANEL_FRAC_THRESHOLD, color=WONG["blue"], lw=0.8, ls=":",
                   alpha=0.8)

        outlier_pos = ihs_chr.loc[ihs_chr["norm_iHS"].abs() > IHS_THRESHOLD, "pos"]
        for p in outlier_pos:
            ax.axvline(p / 1e6, color=WONG["vermillion"], lw=0.5, alpha=0.35,
                       zorder=1)

        mids = (cand_chr["start"] + cand_chr["end"]) / 2e6
        ax.scatter(mids, np.full(len(mids), 0.93), s=26, marker="v",
                   color=WONG["green"], edgecolor="black", lw=0.4, zorder=5)

        ax.set_ylim(0, 1.0)
        ax.set_xlim(0, pos_mb.max())
        ax.set_ylabel("ROH frac.", fontsize=9)
        ax.text(0.01, 0.90, f"{chrom}  ({len(cand_chr)} candidate"
                f"{'s' if len(cand_chr) != 1 else ''})",
                transform=ax.transAxes, fontsize=9.5, fontweight="bold",
                va="top", ha="left")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.grid(True, axis="y", alpha=0.20)
        ax.set_axisbelow(True)

    axes[-1].set_xlabel(r"Position along $\it{S.\ stenocarpa}$ "
                        "pseudo-chromosome (Mb)")

    handles = [
        plt.Line2D([], [], color=WONG["blue"], lw=6, alpha=0.35,
                   label="Panel ROH density"),
        plt.Line2D([], [], color=WONG["blue"], lw=0.8, ls=":",
                   label=f"ROH threshold = {ROH_PANEL_FRAC_THRESHOLD:.2f}"),
        plt.Line2D([], [], color=WONG["vermillion"], lw=0.8, alpha=0.5,
                   label=f"|iHS| > {IHS_THRESHOLD:.0f} marker"),
        plt.Line2D([], [], color=WONG["green"], marker="v", ls="none",
                   mec="black", mew=0.4, ms=8,
                   label=f"Selection candidate (ROH + iHS within {FLANK_KB} kb)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.015))
    fig.suptitle("Candidate-gene scan in bottleneck-distorted chromosomes "
                 "(ROH density x iHS outliers x Funannotate genes)",
                 fontsize=12, y=1.0)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    written = save_figure(fig, OUT / "fig_candidate_genes_bottleneck_chrs")
    plt.close(fig)
    print("Wrote", *written)


if __name__ == "__main__":
    main()
