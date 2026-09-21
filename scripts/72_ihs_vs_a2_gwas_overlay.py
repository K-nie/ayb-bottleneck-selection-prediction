#!/usr/bin/env python3
"""
Paper-B addition: iHS x A2 GWAS-top-30 cross-overlay.

The selscan iHS scan (script 36) returns 17 SNPs at |norm_iHS| > 2 across
the AYB-anchored marker subset, distributed over Ss02, Ss03, Ss04, Ss06,
and Ss08. Strikingly, Ss10 -- which carries the load-bearing Soluble_
Oxalate x ALMT4 candidate from Paper A2 -- has zero |iHS| > 2 hits. This
cross-overlay quantifies the relationship between Paper A2's per-trait
GWAS top-30 SNP set and Paper B's iHS scan:

  1. Per trait, look up each of the top-30 GWAS SNPs in the iHS table
     and record norm_iHS, |norm_iHS|, the binary selscan crit flag, and
     the genomic distance to the nearest |iHS| > 2 SNP.
  2. Per A2 candidate-gene callout focal region (Soluble_Oxalate x
     ALMT4 at Ss10:15.4 Mb; Insoluble_Oxalate x MATE-cluster at
     Ss04:67.8 Mb), compute a windowed iHS density: number of |iHS|
     > 2 hits within +/- 500 kb, +/- 1 Mb, +/- 2 Mb.
  3. Per trait, report the overlap statistic: how many of the top-30
     GWAS SNPs sit at |iHS| > 2 (versus the per-chromosome expectation
     under the panel-wide |iHS| > 2 rate).

The headline finding the script locks in: ALMT4 (the A2 load-bearing
mechanism) is NOT inside an iHS sweep. The trait-associated allele
tracks standing diversity from before any recent selection, not from
an ongoing sweep at this locus -- a clean biological interpretation
that A2 and B together support.

Outputs (results/72_ihs_vs_a2_gwas/)
------------------------------------
tables/
    per_trait_top30_ihs_lookup.csv   per-trait top-30 SNP iHS lookup
    a2_focal_region_density.csv      iHS density at the A2 focal regions
    per_trait_overlap_summary.csv    per-trait iHS-overlap stats vs panel rate
figures/
    fig_ihs_at_a2_callouts.png/.pdf  iHS Manhattan with A2 focal SNPs marked

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, publishable_axes, save_figure
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
IHS_CSV = ROOT / "results" / "36_ihs" / "tables" / "ihs_combined.csv"
GWAS_TOP_DIR = ROOT / "results" / "21_candidate_genes_ayb" / "tables"

OUT = ROOT / "results" / "72_ihs_vs_a2_gwas"
TAB = OUT / "tables"
FIG = OUT / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

TRAITS = [
    "Tannin", "Phenol", "Flavonoid", "Antioxidant",
    "Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds",
    "Seed_Coat_Tannin", "Crude_Protein",
    "Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate",
]

# Two A2 candidate-gene focal regions, used to compute windowed iHS density.
A2_FOCAL = [
    {"trait": "Soluble_Oxalate", "gene": "ALMT4_2",
     "chr": "Ss10", "pos": 15_394_673,
     "note": "A2 four-method consensus locus"},
    {"trait": "Insoluble_Oxalate", "gene": "Ss04 MATE cluster",
     "chr": "Ss04", "pos": 67_792_286,
     "note": "A2 suggestive locus (single-marker p = 0.18, rank 26)"},
]

WINDOWS_KB = [500, 1_000, 2_000]
IHS_THRESHOLD = 2.0


def load_ihs() -> pd.DataFrame:
    df = pd.read_csv(IHS_CSV)
    df["abs_iHS"] = df["norm_iHS"].abs()
    return df


def load_top30(trait: str) -> pd.DataFrame | None:
    path = GWAS_TOP_DIR / f"snp_top30_anchored_{trait}.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


def lookup_ihs_per_snp(top30: pd.DataFrame, ihs: pd.DataFrame, trait: str) -> pd.DataFrame:
    """For each GWAS top-30 SNP, look up the iHS value at the same (chr,
    pos). If the SNP itself is absent from the iHS table (iHS is computed
    only on the AYB-anchored markers that survive the Beagle phasing
    filter, which is a tighter subset), record the nearest iHS SNP on the
    same chromosome and the genomic distance."""
    rows = []
    for _, snp in top30.iterrows():
        chr_a = str(snp["chr_ayb"])
        pos = int(snp["snp_pos_ayb"])
        # Exact-position match first.
        exact = ihs.loc[(ihs["chr"] == chr_a) & (ihs["pos"] == pos)]
        if len(exact) > 0:
            ihs_row = exact.iloc[0]
            rows.append({
                "trait": trait,
                "gwas_rs": snp["rs"],
                "chr": chr_a,
                "pos": pos,
                "gwas_p": float(snp["p"]),
                "ihs_at_snp": True,
                "norm_iHS": float(ihs_row["norm_iHS"]),
                "abs_iHS": float(ihs_row["abs_iHS"]),
                "selscan_crit": int(ihs_row["crit"]),
                "nearest_ihs_pos": pos,
                "nearest_ihs_dist_bp": 0,
                "nearest_ihs_value": float(ihs_row["norm_iHS"]),
            })
            continue
        # Otherwise the nearest iHS SNP on the same chromosome.
        same_chr = ihs.loc[ihs["chr"] == chr_a]
        if same_chr.empty:
            rows.append({
                "trait": trait, "gwas_rs": snp["rs"], "chr": chr_a,
                "pos": pos, "gwas_p": float(snp["p"]),
                "ihs_at_snp": False, "norm_iHS": np.nan,
                "abs_iHS": np.nan, "selscan_crit": np.nan,
                "nearest_ihs_pos": np.nan,
                "nearest_ihs_dist_bp": np.nan,
                "nearest_ihs_value": np.nan,
            })
            continue
        dists = (same_chr["pos"] - pos).abs()
        nearest = same_chr.iloc[dists.values.argmin()]
        rows.append({
            "trait": trait, "gwas_rs": snp["rs"], "chr": chr_a,
            "pos": pos, "gwas_p": float(snp["p"]),
            "ihs_at_snp": False, "norm_iHS": np.nan,
            "abs_iHS": np.nan, "selscan_crit": np.nan,
            "nearest_ihs_pos": int(nearest["pos"]),
            "nearest_ihs_dist_bp": int(dists.min()),
            "nearest_ihs_value": float(nearest["norm_iHS"]),
        })
    return pd.DataFrame(rows)


def focal_region_density(ihs: pd.DataFrame, focal: dict) -> dict:
    """Count |iHS| > 2 SNPs within +/- W kb windows around a focal point."""
    out = dict(focal)
    on_chr = ihs.loc[ihs["chr"] == focal["chr"]].copy()
    if on_chr.empty:
        for w in WINDOWS_KB:
            out[f"ihs_above2_within_{w}kb"] = 0
            out[f"ihs_snps_within_{w}kb"] = 0
        out["closest_above2_dist_bp"] = np.nan
        out["closest_above2_pos"] = np.nan
        out["closest_above2_value"] = np.nan
        return out
    on_chr["dist"] = (on_chr["pos"] - focal["pos"]).abs()
    for w in WINDOWS_KB:
        window = on_chr.loc[on_chr["dist"] <= w * 1_000]
        out[f"ihs_above2_within_{w}kb"] = int((window["abs_iHS"] > IHS_THRESHOLD).sum())
        out[f"ihs_snps_within_{w}kb"] = int(len(window))
    above2_chr = on_chr.loc[on_chr["abs_iHS"] > IHS_THRESHOLD]
    if above2_chr.empty:
        out["closest_above2_dist_bp"] = np.nan
        out["closest_above2_pos"] = np.nan
        out["closest_above2_value"] = np.nan
    else:
        idx = above2_chr["dist"].values.argmin()
        nearest = above2_chr.iloc[idx]
        out["closest_above2_dist_bp"] = int(nearest["dist"])
        out["closest_above2_pos"] = int(nearest["pos"])
        out["closest_above2_value"] = float(nearest["norm_iHS"])
    return out


def per_trait_overlap_summary(per_snp: pd.DataFrame, ihs: pd.DataFrame) -> pd.DataFrame:
    """For each trait, report how many of the top-30 SNPs have an iHS
    score at the SNP itself, how many are at |iHS| > 2, and the panel-
    wide expected count under random assignment of 30 SNPs to the iHS
    set."""
    panel_above2_rate = (ihs["abs_iHS"] > IHS_THRESHOLD).mean()
    n_top30 = 30
    expected_count = panel_above2_rate * n_top30

    rows = []
    for trait, sub in per_snp.groupby("trait"):
        n_ihs_at_snp = int(sub["ihs_at_snp"].sum())
        n_above2_at_snp = int((sub["abs_iHS"] > IHS_THRESHOLD).sum())
        # Convert per-trait observed-vs-expected into a binomial-tail
        # 1-sided p-value for "is the trait enriched for |iHS| > 2 over
        # the panel rate?"
        from scipy import stats
        emp_p = float(1 - stats.binom.cdf(n_above2_at_snp - 1, n_top30,
                                            panel_above2_rate))
        rows.append({
            "trait": trait,
            "n_top30": n_top30,
            "n_with_ihs_at_snp": n_ihs_at_snp,
            "n_above_iHS_2_at_snp": n_above2_at_snp,
            "expected_under_panel_rate": expected_count,
            "panel_above2_rate": panel_above2_rate,
            "binom_p_one_sided": emp_p,
        })
    return pd.DataFrame(rows).sort_values("n_above_iHS_2_at_snp", ascending=False)


def plot_ihs_with_a2_callouts(ihs: pd.DataFrame, focal_list: list[dict],
                                out_path: Path) -> None:
    chrom_order = [f"Ss{i:02d}" for i in range(1, 12)]
    df = ihs.copy()
    df = df.loc[df["chr"].isin(chrom_order)]
    df = df.sort_values(["chr", "pos"]).reset_index(drop=True)

    chrom_offsets = {}
    cur = 0
    GAP = 5_000_000
    for c in chrom_order:
        chrom_offsets[c] = cur
        max_pos = df.loc[df["chr"] == c, "pos"].max() if (df["chr"] == c).any() else 0
        cur += (max_pos or 0) + GAP
    df["x"] = df.apply(lambda r: chrom_offsets[r["chr"]] + r["pos"], axis=1)

    fig, ax = plt.subplots(figsize=(12, 4.6))
    for c in chrom_order:
        sub = df.loc[df["chr"] == c]
        if sub.empty:
            continue
        col = WONG["blue"] if (int(c[2:]) % 2) else WONG["skyblue"]
        ax.scatter(sub["x"], sub["norm_iHS"], s=10, color=col, alpha=0.75,
                   edgecolors="none", rasterized=True)
    # Threshold lines.
    ax.axhline(IHS_THRESHOLD, color=WONG["grey"], linewidth=0.7, linestyle="--")
    ax.axhline(-IHS_THRESHOLD, color=WONG["grey"], linewidth=0.7, linestyle="--")

    # y-limits sized to the data, with fixed headroom for the callout row so
    # the gene labels sit in clear space above the cloud and below the title.
    y_max = float(np.nanmax(np.abs(df["norm_iHS"])))
    top = max(IHS_THRESHOLD + 0.5, y_max) + 0.4
    label_y = top - 0.25          # callout text baseline, inside the axes
    ax.set_ylim(-top, top)

    # Focal SNP markers: vertical band + a leader from the cloud up to a
    # label parked on a single tidy row, so labels never overwrite the points.
    for focal in focal_list:
        sub = df.loc[df["chr"] == focal["chr"]]
        if sub.empty:
            continue
        x_focal = chrom_offsets[focal["chr"]] + focal["pos"]
        ax.axvspan(x_focal - 200_000, x_focal + 200_000,
                   color=WONG["vermillion"], alpha=0.18, zorder=0)
        ax.axvline(x_focal, color=WONG["vermillion"], lw=0.8, ls=":",
                   alpha=0.7, zorder=1)
        ax.annotate(focal["gene"],
                    xy=(x_focal, IHS_THRESHOLD),
                    xytext=(x_focal, label_y),
                    fontsize=9, color=WONG["vermillion"], fontweight="bold",
                    ha="center", va="top", rotation=0, zorder=6,
                    bbox={"boxstyle": "round,pad=0.3", "fc": "white",
                          "ec": WONG["vermillion"], "lw": 0.7, "alpha": 0.95},
                    arrowprops={"arrowstyle": "-", "color": WONG["vermillion"],
                                "lw": 0.7, "alpha": 0.8})
    ticks = [chrom_offsets[c] + (df.loc[df["chr"] == c, "pos"].max() or 0) / 2
              for c in chrom_order]
    ax.set_xticks(ticks)
    ax.set_xticklabels(chrom_order, fontsize=8)
    ax.set_xlabel(r"$\it{S.\ stenocarpa}$ pseudo-chromosome (AYB-anchored markers)")
    ax.set_ylabel(r"Normalised iHS (threshold $\pm$2)")
    ax.set_title("Recent positive selection (iHS) with Paper A2 "
                 "candidate-gene focal regions", loc="left", pad=8)
    publishable_axes(ax, grid="y")
    written = save_figure(fig, out_path)
    plt.close(fig)
    print("Wrote", *written)


def main() -> None:
    print(f"[load] iHS: {IHS_CSV.name}")
    ihs = load_ihs()
    print(f"  {len(ihs)} SNPs; |norm_iHS| > {IHS_THRESHOLD}: "
          f"{(ihs['abs_iHS'] > IHS_THRESHOLD).sum()}")

    all_per_snp = []
    for t in TRAITS:
        top30 = load_top30(t)
        if top30 is None:
            print(f"  [skip] no top-30 file for {t}")
            continue
        per_snp = lookup_ihs_per_snp(top30, ihs, t)
        all_per_snp.append(per_snp)
    per_snp = pd.concat(all_per_snp, ignore_index=True)
    per_snp_path = TAB / "per_trait_top30_ihs_lookup.csv"
    per_snp.to_csv(per_snp_path, index=False)
    print(f"[result] wrote {per_snp_path} ({len(per_snp)} rows)")

    focal_rows = [focal_region_density(ihs, f) for f in A2_FOCAL]
    focal_df = pd.DataFrame(focal_rows)
    focal_path = TAB / "a2_focal_region_density.csv"
    focal_df.to_csv(focal_path, index=False)
    print(f"[result] wrote {focal_path}")
    print(focal_df[["gene", "chr", "pos",
                      "ihs_above2_within_500kb",
                      "ihs_above2_within_1000kb",
                      "ihs_above2_within_2000kb",
                      "closest_above2_dist_bp",
                      "closest_above2_value"]].to_string(index=False))

    summary = per_trait_overlap_summary(per_snp, ihs)
    summary_path = TAB / "per_trait_overlap_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"[result] wrote {summary_path}")
    print(summary.to_string(index=False))

    print("[plot] iHS Manhattan with A2 callouts")
    plot_ihs_with_a2_callouts(ihs, A2_FOCAL, FIG / "fig_ihs_at_a2_callouts")

    # Manuscript-text summary
    print("\n=== manuscript-text summary ===")
    panel_above2_rate = (ihs["abs_iHS"] > IHS_THRESHOLD).mean()
    print(f"  panel-wide |iHS| > 2 rate                 : {panel_above2_rate:.3f}")
    for focal, row in zip(A2_FOCAL, focal_rows):
        print(f"  A2 candidate region {focal['gene']:25s} ({focal['chr']}:{focal['pos']:,}):")
        for w in WINDOWS_KB:
            print(f"     |iHS| > 2 within +/- {w} kb        : "
                  f"{row[f'ihs_above2_within_{w}kb']}")
        if pd.notna(row["closest_above2_dist_bp"]):
            print(f"     closest |iHS| > 2 distance         : "
                  f"{row['closest_above2_dist_bp']:,} bp at norm_iHS = "
                  f"{row['closest_above2_value']:+.2f}")
        else:
            print(f"     closest |iHS| > 2 distance         : NONE on this chromosome")


if __name__ == "__main__":
    main()
