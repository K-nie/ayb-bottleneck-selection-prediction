"""Shared publication-quality figure style for AYB DArTseq manuscripts.

Layered on top of the original `_plotstyle.py` palette + rc settings. Adds:

  * `apply()`            -- journal-grade matplotlib rcParams (300 dpi PDF + PNG)
  * `publishable_axes()` -- per-axes spine + tick + grid cleanup
  * `save_figure()`      -- writes PNG and PDF side-by-side at 300 dpi
  * `adjust_labels()`    -- adjustText wrapper with sensible repulsion defaults
  * `label_top_n()`      -- pick top-N points by |y| (or arbitrary score) and
                            place non-overlapping labels via adjustText
  * `manhattan_axes()`   -- assemble a Manhattan x-axis with alternating chrom
                            shading + centered chromosome ticks
  * `manhattan_with_genes()` -- plot per-SNP -log10(p) + gene callouts at
                                q < threshold (uses adjustText)
  * `WONG` / `WONG_CYCLE` / `CLUSTER_PAL` / `TRAIT_PAL` / `TRAIT_PAL_FULL`
    -- colour-vision-friendly palettes

adjustText is treated as an optional dep -- helpers fall back to plain
`ax.text()` if the package is missing. The yeast-viz conda env ships
adjustText 1.4.0; activate that env when re-rendering the collision-flagged
figures.

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from adjustText import adjust_text as _adjust_text
    HAVE_ADJUSTTEXT = True
except ImportError:
    HAVE_ADJUSTTEXT = False


# ---------------------------------------------------------------------------
# Palettes (Okabe-Ito / Wong colour-vision-friendly)
# ---------------------------------------------------------------------------
WONG = {
    "black": "#000000",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#999999",
}

WONG_CYCLE = [
    WONG["blue"], WONG["vermillion"], WONG["green"], WONG["orange"],
    WONG["purple"], WONG["skyblue"], WONG["yellow"], WONG["black"],
]

CLUSTER_PAL = [WONG["blue"], WONG["vermillion"], WONG["green"], WONG["orange"]]

TRAIT_PAL = {
    "Tannin": WONG["vermillion"],
    "Phenol": WONG["blue"],
    "Flavonoid": WONG["green"],
    "Antioxidant": WONG["purple"],
}

TRAIT_PAL_FULL = {
    "Tannin": WONG["vermillion"],
    "Phenol": WONG["blue"],
    "Flavonoid": WONG["green"],
    "Antioxidant": WONG["purple"],
    "Crude_Protein": WONG["orange"],
    "Total_Oxalate": "#7E1E1E",
    "Soluble_Oxalate": WONG["vermillion"],
    "Insoluble_Oxalate": "#9C6E00",
    "Seed_Coat_Tannin": "#5C3A21",
    "Seed_Length": WONG["skyblue"],
    "Seed_Width": "#2E5E8C",
    "Seed_Thickness": "#0E3550",
    "Mass_of_Seeds": WONG["black"],
}

# Pleiotropy categories used in selection-index / OCS plots
SCENARIO_PAL = {
    "Balanced": WONG["blue"],
    "Nutritional": WONG["blue"],
    "Soluble_Oxalate priority": WONG["vermillion"],
    "Anti Nutritional": WONG["vermillion"],
    "Anti-Nutritional": WONG["vermillion"],
    "Crude_Protein priority": WONG["green"],
}


# ---------------------------------------------------------------------------
# rcParams + per-axes cleanup
# ---------------------------------------------------------------------------
def apply() -> None:
    """Set matplotlib rcParams for journal-grade output (300 dpi, no
    top / right spines, sans-serif, vector-safe PDF fonts)."""
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "font.family": "sans-serif",
        # Arial / DejaVu ship a *separate* bold face; macOS "Helvetica" is a
        # single .ttc whose bold cannot be selected by matplotlib, so a
        # fontweight="bold" request silently renders regular. Keep a
        # true-bold family first so panel letters and bold titles are bold.
        "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica Neue", "Helvetica"],
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "figure.titlesize": 12,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.9,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "grid.linewidth": 0.5,
        "grid.alpha": 0.35,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def panel_label(ax: plt.Axes, letter: str, *,
                 x: float = -0.08, y: float = 1.02,
                 fontsize: float = 12.0,
                 va: str = "bottom", ha: str = "right") -> plt.Text:
    """Stamp a Nature-style panel letter at the top-left of `ax`.

    Panel letters are **lowercase and bold** (a, b, c ...) in the sans-serif
    figure font, placed just outside the top-left corner of the axes in axes-
    fraction coordinates. Any uppercase letter passed in is lowered so legacy
    'A'/'B' callers convert automatically. Returns the Text handle.
    """
    return ax.text(x, y, str(letter).lower(), transform=ax.transAxes,
                    fontsize=fontsize, fontweight="bold",
                    va=va, ha=ha, zorder=1000)


def panel_labels(axes: Iterable[plt.Axes], *,
                  start: int = 0, **kwargs) -> list[plt.Text]:
    """Apply :func:`panel_label` to an ordered iterable of axes, lettering
    them a, b, c ... from ``start`` (0 -> 'a'). Returns the Text handles."""
    out = []
    for i, ax in enumerate(axes):
        out.append(panel_label(ax, chr(ord("a") + start + i), **kwargs))
    return out


def publishable_axes(ax: plt.Axes, grid: str | None = "y",
                       grid_alpha: float = 0.25) -> None:
    """Strip top/right spines, lighten ticks, and (optionally) draw a faint
    grid on the requested axis. Call once per Axes when working with
    `subplots`."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(direction="out", width=0.9, length=3.5)
    if grid is None:
        ax.grid(False)
    else:
        ax.grid(axis=grid, color=WONG["grey"], linewidth=0.5,
                alpha=grid_alpha, zorder=0)
        ax.set_axisbelow(True)


# ---------------------------------------------------------------------------
# Saving
# ---------------------------------------------------------------------------
def save_figure(fig: plt.Figure, path_stub: str | Path,
                  formats: Sequence[str] = ("png", "pdf"),
                  dpi: int = 300) -> list[Path]:
    """Save `fig` to `<path_stub>.<ext>` for each ext in `formats`.

    `path_stub` may be passed without an extension. Returns the list of
    written paths.
    """
    stub = Path(path_stub)
    if stub.suffix.lower() in {".png", ".pdf", ".svg"}:
        stub = stub.with_suffix("")
    stub.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for ext in formats:
        out = stub.with_suffix(f".{ext}")
        fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
        written.append(out)
    return written


# ---------------------------------------------------------------------------
# adjustText label helpers
# ---------------------------------------------------------------------------
def adjust_labels(texts: list, ax: plt.Axes | None = None,
                    arrow_color: str = "#444444",
                    arrow_lw: float = 0.6,
                    expand: tuple[float, float] = (1.20, 1.40),
                    force_text: tuple[float, float] = (0.6, 0.8),
                    force_points: tuple[float, float] = (0.4, 0.6),
                    only_move: dict | None = None,
                    iter_lim: int = 250) -> None:
    """Wrap `adjustText.adjust_text` with sensible defaults for journal
    Manhattan / scatter / forest plots. Silently no-ops if adjustText is not
    importable -- callers should ensure the yeast-viz env is active.
    """
    if not HAVE_ADJUSTTEXT or not texts:
        return
    _adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color=arrow_color, lw=arrow_lw,
                          alpha=0.85, shrinkA=2, shrinkB=2),
        expand=expand,
        force_text=force_text,
        force_static=force_points,
        only_move=only_move or {"text": "xy", "static": "xy"},
        max_move=18,
        time_lim=8.0,
    )


def label_top_n(ax: plt.Axes, x: Sequence[float], y: Sequence[float],
                labels: Sequence[str], n: int = 15,
                score: Sequence[float] | None = None,
                fontsize: float = 7.5,
                color: str = "#222222",
                **adjust_kwargs) -> list:
    """Label the top-`n` points of (x, y) by |score| (defaults to |y|) and
    place the labels with adjustText so they do not collide. Returns the
    list of `Text` objects placed (so the caller can style further).
    """
    x = np.asarray(x)
    y = np.asarray(y)
    labels = np.asarray(labels)
    score = np.abs(np.asarray(score if score is not None else y))
    if n >= len(x):
        idx = np.arange(len(x))
    else:
        idx = np.argpartition(-score, n - 1)[:n]
    texts = []
    for i in idx:
        if pd.isna(labels[i]) or str(labels[i]).strip() == "":
            continue
        texts.append(ax.text(x[i], y[i], str(labels[i]),
                              fontsize=fontsize, color=color, zorder=10))
    adjust_labels(texts, ax=ax, **adjust_kwargs)
    return texts


# ---------------------------------------------------------------------------
# Manhattan helpers
# ---------------------------------------------------------------------------
def _chrom_offsets(df: pd.DataFrame, chrom_col: str, pos_col: str,
                    gap_bp: float = 5e6) -> tuple[pd.Series, list[str], list[float]]:
    """Compute cumulative x positions across pseudo-chromosomes."""
    chroms = sorted(df[chrom_col].dropna().unique().tolist())
    offsets = {}
    centers = []
    cursor = 0.0
    for c in chroms:
        sub = df.loc[df[chrom_col] == c, pos_col]
        if sub.empty:
            offsets[c] = cursor
            continue
        offsets[c] = cursor
        centers.append(cursor + sub.max() / 2.0)
        cursor += sub.max() + gap_bp
    x_global = df.apply(lambda r: offsets.get(r[chrom_col], np.nan) + r[pos_col],
                         axis=1)
    return x_global, chroms, centers


def manhattan_axes(ax: plt.Axes, df: pd.DataFrame, *,
                     chrom_col: str = "chrom",
                     pos_col: str = "bp_pos",
                     band_pal: tuple[str, str] = ("#222B45", "#5896C9"),
                     band_alpha: float = 1.0,
                     gap_bp: float = 5e6) -> pd.Series:
    """Lay out a Manhattan x-axis: cumulative coordinate across chromosomes,
    alternating per-chromosome colours, tick labels centred on each chrom.
    Returns the global x coordinate as a Series aligned with `df`.

    The returned x is also stored as `df["_x_global"]` (mutates input) so the
    caller can use it to position scatter / callout labels without recomputing.
    """
    x_global, chroms, centers = _chrom_offsets(df, chrom_col, pos_col,
                                                  gap_bp=gap_bp)
    df["_x_global"] = x_global
    df["_chrom_color"] = [band_pal[i % 2] for i in
                            [chroms.index(c) for c in df[chrom_col]]]
    ax.set_xticks(centers)
    ax.set_xticklabels(chroms, fontsize=8.5)
    ax.set_xlim(x_global.min() - gap_bp * 0.5,
                 x_global.max() + gap_bp * 0.5)
    publishable_axes(ax, grid="y")
    return x_global


def manhattan_with_genes(ax: plt.Axes, df: pd.DataFrame, *,
                            chrom_col: str = "chrom",
                            pos_col: str = "bp_pos",
                            p_col: str = "p",
                            q_col: str | None = "q",
                            gene_col: str = "gene_label",
                            q_threshold: float = 0.10,
                            point_size: float = 8.0,
                            point_alpha: float = 0.55,
                            callout_size: float = 7.5,
                            callout_color: str = WONG["vermillion"],
                            band_pal: tuple[str, str] = ("#222B45", "#5896C9"),
                            highlight_color: str = WONG["vermillion"],
                            ylabel: str | None = None) -> dict:
    """Plot a publication-grade Manhattan with non-overlapping gene callouts.

    `df` must carry per-SNP rows with at least `chrom_col`, `pos_col`, `p_col`.
    SNPs whose `q_col` value is <= `q_threshold` AND whose `gene_col` is not
    empty are labelled with their gene symbol using adjustText so the labels
    never sit on top of each other or on the data.

    Returns a dict with handles to the scatter, the highlighted scatter, and
    the list of Text objects placed.
    """
    df = df.copy()
    df["_neglog10p"] = -np.log10(df[p_col].clip(lower=1e-300))
    x = manhattan_axes(ax, df, chrom_col=chrom_col, pos_col=pos_col,
                          band_pal=band_pal)

    is_hit = (df[q_col] <= q_threshold) if q_col and q_col in df.columns \
        else (df[p_col] <= 1e-4)
    base = df.loc[~is_hit]
    hit = df.loc[is_hit]

    sc = ax.scatter(base["_x_global"], base["_neglog10p"],
                     c=base["_chrom_color"], s=point_size, alpha=point_alpha,
                     linewidths=0, zorder=2)
    sc_hit = ax.scatter(hit["_x_global"], hit["_neglog10p"],
                          c=highlight_color, s=point_size * 2.0,
                          edgecolors="black", linewidths=0.4, zorder=4)

    # significance line at q_threshold (BH bound via observed -log10p at
    # boundary) -- if no hits, just skip
    if not hit.empty:
        cutoff = -np.log10(hit[p_col].max())
        ax.axhline(cutoff, color=highlight_color, lw=0.7, ls="--",
                     alpha=0.6, zorder=1)

    # gene callouts -- only label rows that pass the q threshold AND have a
    # non-empty gene label
    callout = hit[hit[gene_col].notna() & (hit[gene_col].astype(str).str.strip() != "")]
    texts = []
    for _, r in callout.iterrows():
        texts.append(ax.text(r["_x_global"], r["_neglog10p"],
                              str(r[gene_col]),
                              fontsize=callout_size, color=callout_color,
                              fontstyle="italic", zorder=10))
    adjust_labels(texts, ax=ax)

    ax.set_xlabel("S. stenocarpa pseudo-chromosome", fontsize=10)
    ax.set_ylabel(ylabel or r"$-\log_{10}\,p$", fontsize=10)
    return dict(scatter=sc, scatter_hit=sc_hit, texts=texts)


# ---------------------------------------------------------------------------
# Misc helpers
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Funannotate nearest-gene lookup (cached)
# ---------------------------------------------------------------------------
_GFF_CACHE: dict[str, "pd.DataFrame"] = {}


def _load_funannotate_gff(gff_path: str | Path) -> "pd.DataFrame":
    """Parse the Funannotate GFF once per path, cache the result. Returns a
    DataFrame with columns [chr, start, end, name, product] for genes on
    Ss01..Ss11."""
    import re as _re
    key = str(gff_path)
    if key in _GFF_CACHE:
        return _GFF_CACHE[key]
    attr_re = _re.compile(r"(\w+)=([^;]+)")
    chrom_re = _re.compile(r"^Ss\d+$")
    genes = []
    mrna_product: dict[str, str] = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not chrom_re.match(f[0]):
                continue
            attrs = dict(attr_re.findall(f[8]))
            if f[2] == "gene":
                genes.append({"chr": f[0], "start": int(f[3]),
                                "end": int(f[4]),
                                "gene_id": attrs.get("ID", ""),
                                "name": attrs.get("Name", "")})
            elif f[2] == "mRNA":
                pid = attrs.get("Parent", "")
                if pid and pid not in mrna_product:
                    mrna_product[pid] = attrs.get("product", "")
    gdf = pd.DataFrame(genes)
    gdf["product"] = gdf["gene_id"].map(mrna_product).fillna("")
    gdf["mid"] = (gdf["start"] + gdf["end"]) / 2.0
    _GFF_CACHE[key] = gdf
    return gdf


def nearest_named_gene(chrom: str, pos: float, gff_path: str | Path,
                         window_bp: int = 500_000) -> str | None:
    """Return the nearest *named* (non-Spste.*) Funannotate gene to
    (chrom, pos) within +/- `window_bp`. Falls back to a compact Spste tag
    (Spste.TSs11.05G175210.1 -> Ss05g175210) when no named gene sits
    inside the window.
    """
    import re as _re
    try:
        gdf = _load_funannotate_gff(gff_path)
    except FileNotFoundError:
        return None
    sub = gdf[gdf["chr"] == chrom]
    if sub.empty:
        return None
    dist = (sub["mid"] - pos).abs()
    window = sub.loc[dist <= window_bp].copy()
    window["dist"] = dist.loc[window.index]
    named = window.loc[~window["name"].str.startswith("Spste.")]
    if not named.empty:
        return str(named.sort_values("dist").iloc[0]["name"]).strip()
    i = int(dist.idxmin())
    nm = str(sub.loc[i, "name"]).strip()
    m = _re.match(r"Spste\.\w+\.(\d+)G(\d+)", nm)
    if m:
        return f"Ss{int(m.group(1)):02d}g{m.group(2)}"
    return nm or None


def wrap_xticks(ax: plt.Axes, rotation: float = 45, ha: str = "right",
                  fontsize: float = 8.0) -> None:
    """Rotate x-tick labels to defuse crowding on the 95-accession axis."""
    for lab in ax.get_xticklabels():
        lab.set_rotation(rotation)
        lab.set_ha(ha)
        lab.set_fontsize(fontsize)


def thin_xticks(ax: plt.Axes, every: int = 5) -> None:
    """Show only every-Nth x-tick label (useful when the axis has >50
    categorical labels). Pairs with `wrap_xticks` for the F-per-sample and
    F_ROH per-sample bar charts."""
    labels = ax.get_xticklabels()
    for i, lab in enumerate(labels):
        if i % every != 0:
            lab.set_visible(False)


__all__ = [
    "apply", "publishable_axes", "save_figure",
    "panel_label", "panel_labels",
    "adjust_labels", "label_top_n",
    "manhattan_axes", "manhattan_with_genes",
    "wrap_xticks", "thin_xticks",
    "WONG", "WONG_CYCLE", "CLUSTER_PAL", "TRAIT_PAL", "TRAIT_PAL_FULL",
    "SCENARIO_PAL",
    "HAVE_ADJUSTTEXT",
    "nearest_named_gene",
]
