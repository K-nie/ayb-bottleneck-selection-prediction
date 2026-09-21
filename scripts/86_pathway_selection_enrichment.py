"""
86 -- Pathway-level selection-statistic enrichment.

Two complementary tests on the same underlying ROH and iHS data, asking
biology-forward questions about what is under selection in AYB:

  Test 1 -- GO / Pfam / product-keyword enrichment in the ROH x iHS
  selection-candidate set (script 85's 72-gene call) vs the genome-wide
  Funannotate background. Hypergeometric per term, BH-FDR across terms.

  Test 2 -- Per-pathway selection-statistic comparison: for each curated
  pathway gene set defined by Pfam / product keyword, compare per-gene
  mean panel ROH density and per-gene mean number of iHS outliers within
  +/- 100 kb against the genome-wide background of all 30,840 mRNAs.
  Two-sided permutation p (n_perm = 4999) on the difference of means.

This is the formal biology-forward complement to script 85: it asks
"which pathway gene families are over-represented in the bottleneck-and-
selection signature?" rather than just listing candidate genes.
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import hypergeom

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
sys.path.insert(0, str(PROJ / "scripts"))

# Re-use the GFF parser from script 85 to avoid double-coding
SCR85 = PROJ / "scripts/85_candidate_gene_scan_bottleneck_chrs.py"
from importlib.util import spec_from_file_location, module_from_spec
spec = spec_from_file_location("s85", SCR85)
s85 = module_from_spec(spec); spec.loader.exec_module(s85)  # type: ignore

OUT = PROJ / "results/86_pathway_selection_enrichment"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

CAND_CSV = PROJ / "results/85_candidate_gene_scan/tables/all_candidates.csv"
ROH = PROJ / "results/34_roh/tables/roh_calls.csv"
IHS = PROJ / "results/36_ihs/tables/ihs_combined.csv"

RNG = np.random.default_rng(20260619)
N_PERM = 4999
MIN_GENES_PER_TERM = 5
IHS_THRESHOLD = 2.0
FLANK_BP = 100_000

# Curated biology-forward pathway gene sets. Each set is a Pfam / product
# keyword regex; gene sets are defined off the full Funannotate background.
PATHWAYS = {
    "Auxin homeostasis (IAA hydrolase / PIN / AUX1)":
        r"(?i)(IAA.amino.acid.hydrolase|auxin|AUX1|PIN[1-9]|YUCCA|TAR2)",
    "Phenylpropanoid pathway (PAL / C4H / 4CL / HCT / DFR)":
        r"(?i)(phenylalanine.ammonia.lyase|cinnamate.4|4.coumarate|HCT|"
        r"chalcone|flavonol|flavanone|leucoanthocyanidin|anthocyanin)",
    "Oxalate / organic-acid efflux (ALMT / MATE / ABCC)":
        r"(?i)(ALMT|aluminum.activated.malate.transporter|MATE.efflux|"
        r"DETOXIFICATION|canalicular|multispecific.organic.anion|"
        r"vacuolar.anion|tonoplast)",
    "Oxalate / ascorbate biosynthesis upstream (MIOX / GLDH / GalUR)":
        r"(?i)(inositol.oxygenase|MIOX|GDP.mannose|galacturonate|"
        r"L.galactose|glyoxylate)",
    "Seed-storage proteins (legumin / vicilin / convicilin / phaseolin)":
        r"(?i)(legumin|vicilin|convicilin|phaseolin|11S.globulin|7S.globulin|"
        r"cupin)",
    "Seed-development regulators (AINTEGUMENTA / TTG / DA1 / GW2)":
        r"(?i)(AINTEGUMENTA|AP2.transcription|TTG[1-3]|DA1|GW2|GS3|GW8|"
        r"BIG.SEEDS|STK)",
    "Cell-wall biosynthesis (CESA / PME / expansin / galactosyltransferase)":
        r"(?i)(cellulose.synthase|CESA|pectin.methylesterase|expansin|"
        r"hydroxyproline.O.galactosyltransferase|arabinogalactan)",
    "Drought / ABA / stress response":
        r"(?i)(abscisic.acid|ABA.responsive|ABF|AREB|RD22|RD29|drought|"
        r"erd4|dehydrin|LEA)",
    "Redox / glutathione / oxidative stress":
        r"(?i)(glutathione|glutaredoxin|peroxidase|superoxide.dismutase|"
        r"catalase|thioredoxin|ascorbate.peroxidase)",
    "Iron / metal transport":
        r"(?i)(iron.transporter|ferritin|FRO|IRT|ferric.reductase|"
        r"NRAMP|zinc.transporter)",
    "Primary metabolism (G6PDH / PFK / TCA / glycolysis)":
        r"(?i)(glucose.6.phosphate|phosphofructokinase|pyruvate.kinase|"
        r"aconitase|citrate.synthase|malate.dehydrogenase)",
    "Translation machinery (ribosomal proteins, eIF / eEF)":
        r"(?i)(ribosomal.protein|elongation.factor|translation.initiation)",
    "PPR (mitochondrial / plastid RNA editing)":
        r"(?i)(pentatricopeptide.repeat|PPR|RNA.editing)",
    "Immune signalling (kinases, LRR-RLK, NB-LRR)":
        r"(?i)(LRR.receptor.kinase|NB.LRR|LRR.NB|pattern.recognition|"
        r"PR.protein|chitin.receptor|pbl27|RPS5)",
}


def main() -> None:
    print("[load] GFF + candidate set + selection data")
    genes = s85.parse_gff(PROJ / "refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3")
    cand = pd.read_csv(CAND_CSV)
    cand_ids = set(cand["gene_id"])
    print(f"  background: {len(genes)} mRNAs; candidate set: {len(cand)} mRNAs")

    # ---- Test 1: Pfam enrichment in candidate vs background ----
    print("\n[test 1] Pfam enrichment in 72-gene candidate vs full background")
    def explode_pfam(df):
        rows = []
        for _, r in df.iterrows():
            for p in str(r.get("pfam", "")).split(","):
                p = p.strip()
                if p:
                    rows.append({"gene_id": r["gene_id"], "pfam": p})
        return pd.DataFrame(rows)

    bg_pfam = explode_pfam(genes)
    cand_pfam = explode_pfam(cand)
    bg_total = genes["gene_id"].nunique()
    cand_total = cand["gene_id"].nunique()
    pfam_counts = bg_pfam.groupby("pfam")["gene_id"].nunique().reset_index(
        name="n_background")
    pfam_cand = cand_pfam.groupby("pfam")["gene_id"].nunique().reset_index(
        name="n_candidate")
    enr = pfam_counts.merge(pfam_cand, on="pfam", how="left").fillna(0)
    enr["n_candidate"] = enr["n_candidate"].astype(int)
    enr = enr[enr["n_candidate"] >= 2].copy()  # focus on Pfam terms with at least 2 candidates
    enr["fold"] = (enr["n_candidate"] / cand_total) / (enr["n_background"] / bg_total)
    enr["p_hyper"] = enr.apply(lambda r: hypergeom.sf(
        r["n_candidate"] - 1, bg_total, r["n_background"], cand_total), axis=1)
    from statsmodels.stats.multitest import multipletests
    enr["p_bh"] = multipletests(enr["p_hyper"], method="fdr_bh")[1]
    enr = enr.sort_values("p_hyper").reset_index(drop=True)
    enr.to_csv(OUT / "tables/pfam_enrichment_candidate_vs_background.csv",
                index=False)
    print(f"  Pfam terms with >= 2 candidates: {len(enr)}")
    print(enr.head(20).to_string(index=False))

    # ---- Test 2: Per-pathway selection-statistic comparison ----
    print("\n[test 2] Per-pathway selection-statistic enrichment")
    # Need per-gene ROH density and per-gene n_ihs_outliers. Use the same
    # computation as script 85 but for ALL genes, not just bottleneck chrs.

    roh = pd.read_csv(ROH)
    ihs = pd.read_csv(IHS)

    # per-chrom ROH-density grid
    print("  pre-computing per-chromosome ROH density (50 kb bins)")
    density_by_chr = {}
    for chrom in sorted(genes["chr"].unique()):
        sub_genes = genes[genes["chr"] == chrom]
        chr_max = int(sub_genes["end"].max()) + 1_000_000
        density_by_chr[chrom] = s85.per_locus_roh_density(
            roh, chrom, bin_kb=50, chr_max_bp=chr_max
        )

    print("  per-gene ROH fraction + iHS outlier count")
    gene_metrics = []
    for _, gene in genes.iterrows():
        chrom = gene["chr"]
        if chrom not in density_by_chr:
            continue
        roh_panel_frac = s85.gene_roh_frac(
            int(gene["start"]), int(gene["end"]), density_by_chr[chrom],
            bin_bp=50_000,
        )
        outlier_pos = ihs.loc[(ihs["chr"] == chrom)
                              & (ihs["norm_iHS"].abs() > IHS_THRESHOLD),
                              "pos"].values
        if len(outlier_pos):
            lo = int(gene["start"]) - FLANK_BP
            hi = int(gene["end"]) + FLANK_BP
            n_ihs = int(((outlier_pos >= lo) & (outlier_pos <= hi)).sum())
        else:
            n_ihs = 0
        gene_metrics.append({
            "gene_id": gene["gene_id"], "chr": chrom,
            "start": gene["start"], "end": gene["end"],
            "roh_panel_frac": roh_panel_frac, "n_ihs_within_100kb": n_ihs,
            "annot_search": gene["annot_search"],
        })
    gm = pd.DataFrame(gene_metrics)
    gm.to_csv(OUT / "tables/per_gene_selection_metrics.csv", index=False)

    print(f"  ROH-density background mean = {gm['roh_panel_frac'].mean():.4f}")
    print(f"  iHS-count background mean = {gm['n_ihs_within_100kb'].mean():.4f}")

    # Per-pathway test
    pathway_rows = []
    for name, pat in PATHWAYS.items():
        mask = gm["annot_search"].str.contains(pat, regex=True, na=False)
        n_pw = int(mask.sum())
        if n_pw < MIN_GENES_PER_TERM:
            print(f"  SKIP {name}: only {n_pw} genes (< {MIN_GENES_PER_TERM})")
            continue
        roh_pw = gm.loc[mask, "roh_panel_frac"].values
        ihs_pw = gm.loc[mask, "n_ihs_within_100kb"].values
        roh_bg = gm.loc[~mask, "roh_panel_frac"].values
        ihs_bg = gm.loc[~mask, "n_ihs_within_100kb"].values

        # Two-sample permutation test on the difference of means
        def perm_p(a, b, n_perm=N_PERM):
            obs = a.mean() - b.mean()
            pool = np.concatenate([a, b])
            n_a = len(a)
            cnt = 1
            for _ in range(n_perm):
                p = RNG.permutation(pool)
                diff = p[:n_a].mean() - p[n_a:].mean()
                if abs(diff) >= abs(obs):
                    cnt += 1
            return float(obs), cnt / (n_perm + 1)

        obs_roh, p_roh = perm_p(roh_pw, roh_bg, n_perm=2999)
        obs_ihs, p_ihs = perm_p(ihs_pw, ihs_bg, n_perm=2999)
        pathway_rows.append({
            "pathway": name, "n_genes": n_pw,
            "roh_mean_pathway": float(roh_pw.mean()),
            "roh_mean_background": float(roh_bg.mean()),
            "roh_delta": obs_roh, "p_perm_roh": p_roh,
            "ihs_mean_pathway": float(ihs_pw.mean()),
            "ihs_mean_background": float(ihs_bg.mean()),
            "ihs_delta": obs_ihs, "p_perm_ihs": p_ihs,
        })
    pdf = pd.DataFrame(pathway_rows)
    if len(pdf):
        from statsmodels.stats.multitest import multipletests
        for col in ("p_perm_roh", "p_perm_ihs"):
            mask = pdf[col].notna()
            pdf.loc[mask, col + "_bh"] = multipletests(
                pdf.loc[mask, col], method="fdr_bh"
            )[1]
        pdf = pdf.sort_values("p_perm_roh").reset_index(drop=True)
        pdf.to_csv(OUT / "tables/pathway_selection_enrichment.csv",
                   index=False)
        print("\n[pathway enrichment results]")
        print(pdf[["pathway", "n_genes",
                    "roh_mean_pathway", "roh_mean_background", "roh_delta",
                    "p_perm_roh", "p_perm_roh_bh",
                    "ihs_mean_pathway", "ihs_delta",
                    "p_perm_ihs", "p_perm_ihs_bh"]].round(4).to_string(index=False))

        # ---- Figure: per-pathway forest of ROH delta and iHS delta ----
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
        pdf_plot = pdf.sort_values("roh_delta")

        ax = axes[0]
        y = np.arange(len(pdf_plot))
        colors = ["#117733" if (p < 0.05 and d > 0) else
                  "#CC6677" if (p < 0.05 and d < 0) else "#999999"
                  for p, d in zip(pdf_plot["p_perm_roh"], pdf_plot["roh_delta"])]
        ax.barh(y, pdf_plot["roh_delta"], color=colors, edgecolor="white")
        ax.axvline(0, color="black", lw=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{n}\n(n={ng})" for n, ng in
                            zip(pdf_plot["pathway"], pdf_plot["n_genes"])],
                            fontsize=8)
        ax.set_xlabel("Delta panel ROH fraction (pathway − background)")
        ax.set_title("Per-pathway ROH-density delta vs background")
        ax.text(-0.08, 1.02, "a", transform=ax.transAxes, fontsize=12,
                fontweight="bold", ha="right", va="bottom", zorder=1000)
        ax.grid(True, alpha=0.25, axis="x")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

        ax = axes[1]
        pdf_plot2 = pdf.sort_values("ihs_delta")
        y2 = np.arange(len(pdf_plot2))
        colors2 = ["#117733" if (p < 0.05 and d > 0) else
                   "#CC6677" if (p < 0.05 and d < 0) else "#999999"
                   for p, d in zip(pdf_plot2["p_perm_ihs"], pdf_plot2["ihs_delta"])]
        ax.barh(y2, pdf_plot2["ihs_delta"], color=colors2, edgecolor="white")
        ax.axvline(0, color="black", lw=0.7)
        ax.set_yticks(y2)
        ax.set_yticklabels([f"{n}\n(n={ng})" for n, ng in
                             zip(pdf_plot2["pathway"], pdf_plot2["n_genes"])],
                             fontsize=8)
        ax.set_xlabel("Delta number of iHS outliers within ±100 kb")
        ax.set_title("Per-pathway iHS-outlier-proximity delta vs background")
        ax.text(-0.08, 1.02, "b", transform=ax.transAxes, fontsize=12,
                fontweight="bold", ha="right", va="bottom", zorder=1000)
        ax.grid(True, alpha=0.25, axis="x")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

        fig.suptitle("Per-pathway selection-statistic enrichment in the AYB DArTseq panel\n"
                     "Green = significant positive enrichment (p_perm < 0.05); "
                     "red = significant negative; grey = ns",
                     fontsize=11, y=1.03)
        fig.tight_layout()
        fig.savefig(OUT / "figures/fig_pathway_selection_enrichment.png",
                    dpi=300, bbox_inches="tight")
        fig.savefig(OUT / "figures/fig_pathway_selection_enrichment.pdf",
                    bbox_inches="tight")
        plt.close(fig)
        print(f"\nWrote figure: {OUT / 'figures/fig_pathway_selection_enrichment.png'}")
    else:
        print("\n[pathway enrichment] no pathways met the n_genes threshold")


if __name__ == "__main__":
    main()
