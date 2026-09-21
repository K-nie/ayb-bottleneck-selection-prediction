"""
85 -- Genome-wide candidate-gene scan in bottleneck-distorted chromosomes.

For each chromosome with strong bottleneck distortion (Tajima D > 2 windows),
identify genes that sit at the intersection of two selection-relevant
signatures:
  (a) high panel-level runs-of-homozygosity (ROH) density -- fraction of
      the 95-accession panel carrying a ROH segment covering the gene body
      or its immediate ±10 kb flanks. Elevated ROH density at a locus is
      consistent with a recent selective sweep or with a region preserved
      under purifying selection in the AYB genebank panel.
  (b) iHS outlier proximity -- gene body or ±100 kb flanks containing at
      least one |norm_iHS| > 2 marker. Captures genes near loci with
      extended haplotype homozygosity signatures.

Genes that satisfy BOTH criteria simultaneously are reported as "selection
candidates"; we cross-reference Pfam / InterPro / Funannotate `product`
annotations to assign biological function. Per-trait relevance is assigned
by keyword match (phenylpropanoid pathway -> biochem traits, MATE / ALMT /
oxalate metabolism -> oxalate fractions, seed-storage proteins ->
Crude_Protein, seed-development regulators -> seed metrics).

Two figures: (i) per-chromosome track plot showing ROH density across the
chromosome + iHS Manhattan + candidate-gene call-out positions; (ii)
per-trait candidate-gene biology summary table.
"""
from __future__ import annotations
import re
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
GFF = PROJ / "refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3"
ROH = PROJ / "results/34_roh/tables/roh_calls.csv"
IHS = PROJ / "results/36_ihs/tables/ihs_combined.csv"
TAJ = PROJ / "results/35_tajima_d/tables/tajima_d_windowed.csv"
OUT = PROJ / "results/85_candidate_gene_scan"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

# Chromosomes ranked by Tajima D distortion (median D from results/74).
# Ss09 has the strongest bottleneck signature (median D = +2.77, 73 % windows
# D > 2); Ss04 second strongest (D = +2.13, 55 %); Ss01 (D = +2.30, 52 %);
# Ss10 (the *ALMT4_2* chromosome; D = +1.95, 48 %).
BOTTLENECK_CHRS = ["Ss09", "Ss04", "Ss01", "Ss10", "Ss06", "Ss08", "Ss02"]
FLANK_BP = 100_000  # iHS outlier window
ROH_PANEL_FRAC_THRESHOLD = 0.40  # fraction of 95 samples with ROH at locus
IHS_THRESHOLD = 2.0  # |norm_iHS| > 2

TRAIT_KEYWORD_SETS = {
    "Crude_Protein":
        ["legumin", "vicilin", "convicilin", "albumin", "globulin",
         "phaseolin", "PR-1", "seed storage", "11S", "7S", "cupin"],
    "Soluble_Oxalate":
        ["malate transporter", "ALMT", "MATE family", "MATE efflux",
         "DETOXIFICATION", "oxalate", "vacuolar", "anion channel"],
    "Total_Oxalate":
        ["oxalate", "MATE", "vacuolar", "anion channel", "ascorbate",
         "glyoxylate"],
    "Insoluble_Oxalate":
        ["MATE", "DETOXIFICATION", "vacuolar", "oxalate", "calcium binding",
         "calcium transport"],
    "Antioxidant":
        ["peroxidase", "superoxide dismutase", "catalase", "glutathione",
         "thioredoxin", "ascorbate", "tocopherol"],
    "Flavonoid":
        ["chalcone synthase", "CHS", "chalcone isomerase", "flavonol",
         "flavanone", "anthocyanin", "MYB", "DFR", "leucoanthocyanidin"],
    "Phenol":
        ["phenylalanine ammonia-lyase", "PAL", "cinnamate", "phenylpropanoid",
         "lignin", "C4H", "4CL", "HCT"],
    "Tannin":
        ["proanthocyanidin", "tannin", "BANYULS", "anthocyanidin reductase",
         "LAR", "leucoanthocyanidin reductase", "tannase"],
    "Seed_Coat_Tannin":
        ["proanthocyanidin", "BANYULS", "TT2", "TT8", "TT12", "MATE", "tannin"],
    "Seed_Length":
        ["AINTEGUMENTA", "ANT", "AP2", "TTG1", "STK", "DA1", "DA2",
         "BIG SEEDS", "DELLA", "BS1", "TRANSPARENT TESTA"],
    "Seed_Width":
        ["pectin methylesterase", "PME", "AINTEGUMENTA", "expansin",
         "cellulose synthase"],
    "Seed_Thickness":
        ["AINTEGUMENTA", "expansin", "cellulose", "xylem", "ANT"],
    "Mass_of_Seeds":
        ["AINTEGUMENTA", "TTG2", "GW2", "GS3", "GW8", "GS5", "phytochrome",
         "PROTEINASE INHIBITOR", "trypsin inhibitor"],
}


def parse_gff(path: Path) -> pd.DataFrame:
    """Parse the Funannotate GFF and return a per-mRNA DataFrame with
    chr, start, end, strand, gene_id, product, PFAM list, InterPro list,
    GO terms, and the concatenated keyword-searchable annotation string."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]), int(f[4]), f[6], f[8]
            attr_d = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_d[k.strip()] = v.strip()
            gid = attr_d.get("Parent", attr_d.get("ID", ""))
            name = attr_d.get("Name", "")
            product = attr_d.get("product", "")
            dbx = attr_d.get("Dbxref", "")
            pfam = ",".join(re.findall(r"PFAM:(PF\d+)", dbx))
            ipr = ",".join(re.findall(r"InterPro:(IPR\d+)", dbx))
            go = attr_d.get("Ontology_term", "")
            note = attr_d.get("note", "")
            rows.append({
                "chr": chrom, "start": start, "end": end, "strand": strand,
                "gene_id": gid, "name": name, "product": product,
                "pfam": pfam, "interpro": ipr, "go": go, "note": note,
            })
    df = pd.DataFrame(rows)
    df["annot_search"] = (df["product"].fillna("") + " | " +
                          df["note"].fillna("") + " | " +
                          df["pfam"].fillna("") + " | " +
                          df["interpro"].fillna("") + " | " +
                          df["name"].fillna(""))
    return df


def per_locus_roh_density(roh: pd.DataFrame, chrom: str, bin_kb: int = 50,
                           chr_max_bp: int = 100_000_000) -> pd.Series:
    """Per-bin fraction of the 95-accession panel carrying a ROH segment
    overlapping each bin on this chromosome."""
    n_samples = roh["sample"].nunique()
    sub = roh[roh["chr"] == chrom]
    bin_bp = bin_kb * 1000
    n_bins = chr_max_bp // bin_bp + 1
    counts = np.zeros(n_bins, dtype=int)
    seen = {b: set() for b in range(n_bins)}
    for _, row in sub.iterrows():
        s_bin = int(row["start_bp"]) // bin_bp
        e_bin = int(row["end_bp"]) // bin_bp
        for b in range(s_bin, e_bin + 1):
            if b < n_bins:
                seen[b].add(row["sample"])
    counts = np.array([len(seen[b]) for b in range(n_bins)])
    return pd.Series(counts / n_samples, index=np.arange(n_bins) * bin_bp,
                     name="roh_panel_frac")


def gene_roh_frac(gene_start: int, gene_end: int, roh_panel_frac: pd.Series,
                   bin_bp: int, flank_bp: int = 10_000) -> float:
    """Mean panel ROH fraction in the gene-body window plus ± flank_bp."""
    s = max(0, gene_start - flank_bp)
    e = gene_end + flank_bp
    s_bin, e_bin = s // bin_bp, e // bin_bp
    bins = roh_panel_frac.iloc[s_bin: e_bin + 1]
    return float(bins.mean()) if len(bins) else 0.0


def main() -> None:
    print("[load] Funannotate GFF (this takes ~ 30 s)")
    genes = parse_gff(GFF)
    print(f"  parsed {len(genes)} mRNAs across {genes['chr'].nunique()} sequences")

    roh = pd.read_csv(ROH)
    ihs = pd.read_csv(IHS)
    taj = pd.read_csv(TAJ)
    print(f"  ROH calls: {len(roh)}; iHS SNPs: {len(ihs)}; "
          f"Tajima windows: {len(taj)}")

    bin_kb = 50; bin_bp = bin_kb * 1000

    all_candidates = []
    chrom_panels = {}

    for chrom in BOTTLENECK_CHRS:
        sub_genes = genes[genes["chr"] == chrom].copy()
        if sub_genes.empty:
            print(f"  {chrom}: no genes in GFF; skip")
            continue
        chr_max_bp = int(sub_genes["end"].max()) + 1_000_000
        roh_density = per_locus_roh_density(roh, chrom, bin_kb=bin_kb,
                                              chr_max_bp=chr_max_bp)
        ihs_chr = ihs[ihs["chr"] == chrom]
        outlier_pos = ihs_chr.loc[ihs_chr["norm_iHS"].abs() > IHS_THRESHOLD,
                                   "pos"].values
        sub_genes["roh_panel_frac"] = sub_genes.apply(
            lambda r: gene_roh_frac(int(r["start"]), int(r["end"]),
                                     roh_density, bin_bp), axis=1
        )
        if len(outlier_pos):
            def ihs_within(start, end):
                lo, hi = start - FLANK_BP, end + FLANK_BP
                return int(((outlier_pos >= lo) & (outlier_pos <= hi)).sum())
            sub_genes["n_ihs_outliers_within_100kb"] = sub_genes.apply(
                lambda r: ihs_within(int(r["start"]), int(r["end"])), axis=1
            )
        else:
            sub_genes["n_ihs_outliers_within_100kb"] = 0

        sub_genes["passes_roh_threshold"] = (
            sub_genes["roh_panel_frac"] >= ROH_PANEL_FRAC_THRESHOLD
        )
        sub_genes["passes_ihs_threshold"] = (
            sub_genes["n_ihs_outliers_within_100kb"] > 0
        )
        sub_genes["selection_candidate"] = (
            sub_genes["passes_roh_threshold"]
            & sub_genes["passes_ihs_threshold"]
        )
        n_cand = sub_genes["selection_candidate"].sum()
        print(f"  {chrom}: {len(sub_genes)} mRNAs; "
               f"ROH ≥ {ROH_PANEL_FRAC_THRESHOLD}: "
               f"{sub_genes['passes_roh_threshold'].sum()}; "
               f"iHS within 100 kb: {sub_genes['passes_ihs_threshold'].sum()}; "
               f"BOTH: {n_cand}")

        candidates = sub_genes[sub_genes["selection_candidate"]].copy()
        all_candidates.append(candidates)
        chrom_panels[chrom] = {
            "roh_density": roh_density, "ihs_chr": ihs_chr,
            "outlier_pos": outlier_pos, "candidates": candidates,
            "all_genes": sub_genes,
        }

    if not all_candidates:
        print("No candidates anywhere; abort")
        return
    cand = pd.concat(all_candidates, ignore_index=True)
    print(f"\nTotal candidates across {len(BOTTLENECK_CHRS)} chromosomes: "
          f"{len(cand)}")

    # ---- per-trait keyword mapping ----
    trait_rows = []
    for tname, keywords in TRAIT_KEYWORD_SETS.items():
        pat = re.compile("|".join(re.escape(k) for k in keywords),
                          flags=re.IGNORECASE)
        hits = cand[cand["annot_search"].str.contains(pat, regex=True, na=False)]
        for _, h in hits.iterrows():
            trait_rows.append({
                "trait": tname, "chr": h["chr"], "start": h["start"],
                "end": h["end"], "gene_id": h["gene_id"],
                "product": h["product"], "pfam": h["pfam"],
                "interpro": h["interpro"],
                "roh_panel_frac": h["roh_panel_frac"],
                "n_ihs_outliers_within_100kb": h["n_ihs_outliers_within_100kb"],
            })
    trait_hits = pd.DataFrame(trait_rows)
    trait_hits.to_csv(OUT / "tables/per_trait_candidate_hits.csv", index=False)
    print(f"Per-trait keyword hits: {len(trait_hits)}")
    if len(trait_hits):
        print(trait_hits.groupby("trait").size().sort_values(ascending=False))

    cand.to_csv(OUT / "tables/all_candidates.csv", index=False)

    # ---- enrichment test: are candidates enriched for *Soluble_Oxalate*
    # MATE/ALMT keywords vs the genome-wide background? ----
    oxalate_pat = re.compile(
        r"(?i)(ALMT|aluminum-activated|MATE|DETOXIFICATION|vacuolar.*anion|"
        r"oxalate|malate transporter)"
    )
    bg_hits = int(genes["annot_search"].str.contains(oxalate_pat, na=False).sum())
    cand_hits = int(cand["annot_search"].str.contains(oxalate_pat, na=False).sum())
    bg_total = len(genes)
    cand_total = len(cand)
    print(f"\nMATE/ALMT/oxalate background hits: {bg_hits} / {bg_total} "
           f"({100 * bg_hits / bg_total:.2f} %)")
    print(f"MATE/ALMT/oxalate candidate hits:   {cand_hits} / {cand_total} "
           f"({100 * cand_hits / max(cand_total, 1):.2f} %)")
    from scipy.stats import hypergeom
    p_enrich = hypergeom.sf(cand_hits - 1, bg_total, bg_hits, cand_total)
    fold = (cand_hits / max(cand_total, 1)) / (bg_hits / bg_total)
    print(f"Hypergeometric enrichment p = {p_enrich:.3g}; fold = {fold:.2f}")
    summary = pd.DataFrame([{
        "test": "MATE/ALMT/oxalate enrichment in bottleneck-chr candidates",
        "n_candidates": cand_total, "n_candidates_with_keyword": cand_hits,
        "n_background": bg_total, "n_background_with_keyword": bg_hits,
        "fold_enrichment": fold, "p_hypergeometric": p_enrich,
    }])
    summary.to_csv(OUT / "tables/oxalate_pathway_enrichment.csv", index=False)

    # ---- per-chromosome track plot ----
    n_chr = len(chrom_panels)
    fig, axes = plt.subplots(n_chr, 1, figsize=(11.5, 2.4 * n_chr),
                               sharex=False)
    if n_chr == 1:
        axes = [axes]
    for ax, (chrom, panel) in zip(axes, chrom_panels.items()):
        roh_density = panel["roh_density"]
        positions_mb = roh_density.index / 1e6
        ax.fill_between(positions_mb, 0, roh_density.values, color="#332288",
                          alpha=0.45, label="Panel ROH density")
        ax.axhline(ROH_PANEL_FRAC_THRESHOLD, color="#332288", lw=0.8, ls=":",
                    alpha=0.7)
        if len(panel["outlier_pos"]):
            for p in panel["outlier_pos"]:
                ax.axvline(p / 1e6, color="#CC6677", lw=0.7, alpha=0.6)
        for _, c_row in panel["candidates"].iterrows():
            ax.scatter([(c_row["start"] + c_row["end"]) / 2e6], [0.92], s=22,
                        marker="v", color="#117733",
                        edgecolor="black", lw=0.4, zorder=5)
        ax.set_ylim(0, 1.0)
        ax.set_ylabel("ROH density")
        ax.set_title(
            f"{chrom}  --  {len(panel['candidates'])} selection candidates "
            f"(ROH-density >= {ROH_PANEL_FRAC_THRESHOLD} AND iHS within {FLANK_BP // 1000} kb)",
            fontsize=10,
        )
        ax.set_xlabel(f"{chrom} position (Mb)")
        ax.grid(True, alpha=0.20)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "figures/fig_candidate_genes_bottleneck_chrs.png",
                 dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "figures/fig_candidate_genes_bottleneck_chrs.pdf",
                 bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote figure: {OUT / 'figures/fig_candidate_genes_bottleneck_chrs.png'}")


if __name__ == "__main__":
    main()
