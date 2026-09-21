"""Shared publication-quality matplotlib settings for AYB figures.

Designed for journal output: 300 dpi, sans-serif body, modest axis weight,
no top/right spines, colour-vision-friendly palettes.

Author: Benjamin Narh-Madey
"""

import matplotlib as mpl
import matplotlib.pyplot as plt


def apply():
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
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
    })


# Okabe-Ito colour-vision-friendly palette
WONG = {
    "black": "#000000",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
}

CLUSTER_PAL = [WONG["blue"], WONG["vermillion"], WONG["green"], WONG["orange"]]
TRAIT_PAL = {
    "Tannin": WONG["vermillion"],
    "Phenol": WONG["blue"],
    "Flavonoid": WONG["green"],
    "Antioxidant": WONG["purple"],
}
