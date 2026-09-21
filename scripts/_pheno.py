"""Unified AYB phenotype loader (v2 — BLUP-enabled).

Three on-disk sources, harmonised into a single wide DataFrame indexed
by sample (genotype) ID. As of 2026-05-18 the wet-chemistry file was
upgraded to Updated_Wet_Chemistry_Data.xlsx, which carries
replicate-level data for 9 of the 13 main traits plus a new
moisture-content trait. This loader exposes both the replicate-level
data and the per-genotype BLUPs from a REML mixed model, with means
as a fallback for traits that still lack replicates.

Data sources
------------
1. Updated_Wet_Chemistry_Data.xlsx (2026-05-18 onward)
     Means              105 x 4 means (Tannin, Phenol, Flavinoid, Antioxidant)
                         — fallback only; BLUPs preferred when reps exist.
     Raw_ANF             3 reps x 105 genotypes for the 4 wet-chemistry traits.
     Crude_Protein       2 reps x 106 genotypes (note: was 46-genotype means
                         in Just_46_Accessions.xlsx before).
     Moisture_Content    2-4 reps x 104 genotypes (NEW; used as a fixed
                         covariate in the Crude_Protein LMM, not as a
                         headline trait).
     Seed_Metrics        10 reps x 105 genotypes for Length / Width /
                         Thickness / Weight (renamed to Seed_Length /
                         Seed_Width / Seed_Thickness / Mass_of_Seeds).

2. Seed_Metrics.xlsx (legacy, still on disk)
     Data_on_Seeds       105 means for Seed_Length / Seed_Width /
                         Seed_Thickness / Mass_of_Seeds plus the only
                         source for Seed_Coat_Tannin (no replicates
                         available for that one trait).

3. Just_46_Accessions.xlsx (legacy)
     Sheet1              46 means for Crude_Protein, Total_Oxalate,
                         Soluble_Oxalate, Insoluble_Oxalate. The
                         Crude_Protein column here is superseded by
                         Updated_Wet_Chemistry_Data.xlsx's 106-genotype
                         replicate data; the three oxalate fractions
                         remain means-only.

Data-quality decisions baked in
-------------------------------
  * "Flavinoid" -> "Flavonoid" (typo in source spreadsheets).
  * "Crude Protein" (space) -> "Crude_Protein".
  * Seed-metrics column rename: Length / Width / Thickness / Weight ->
    Seed_Length / Seed_Width / Seed_Thickness / Mass_of_Seeds.
  * TYPO_CORRECTIONS table (two decimal-point typos) applied to the
    final wide BLUP / means table the same way as before.

Public functions
----------------
load_replicates(trait=None)   long-format replicate-level data
load_blups(...)               wide BLUPs table from REML LMM
load_means(...)               legacy means table (no LMM)
load_phenotypes(...)          DEFAULT = load_blups; means fallback per trait
trait_n_table(df)             per-trait sample size and descriptive stats
heritability_table()          per-trait Holland H2 from the BLUP fit

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import numpy as np

# statsmodels' optimisation layer still references the removed `np.Inf`
# alias. Restore it before any statsmodels import so the LMM solver
# does not crash on NumPy 2.x.
if not hasattr(np, "Inf"):
    np.Inf = np.inf  # type: ignore[attr-defined]

import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
WET  = ROOT / "Updated_Wet_Chemistry_Data.xlsx"
SEED = ROOT / "Seed_Metrics.xlsx"
J46  = ROOT / "Just_46_Accessions.xlsx"

TRAITS_BIOCHEM = ["Tannin", "Phenol", "Flavonoid", "Antioxidant"]
TRAITS_SEED    = ["Seed_Length", "Seed_Width", "Seed_Thickness",
                  "Mass_of_Seeds", "Seed_Coat_Tannin"]
TRAITS_PROT_OX = ["Crude_Protein", "Total_Oxalate",
                  "Soluble_Oxalate", "Insoluble_Oxalate"]
TRAITS_ALL = TRAITS_BIOCHEM + TRAITS_SEED + TRAITS_PROT_OX

# Replicated traits (eligible for BLUP) vs means-only traits.
TRAITS_REPLICATED = ["Tannin", "Phenol", "Flavonoid", "Antioxidant",
                     "Seed_Length", "Seed_Width", "Seed_Thickness",
                     "Mass_of_Seeds", "Crude_Protein"]
TRAITS_MEANS_ONLY = ["Seed_Coat_Tannin", "Total_Oxalate",
                     "Soluble_Oxalate", "Insoluble_Oxalate"]

TYPO_CORRECTIONS = [
    {"sample": "TSs282", "trait": "Seed_Thickness",
     "raw_value": 63.747, "corrected_value": 6.3747,
     "reason": "Decimal-point typo in the legacy means file: 63.747 mm is "
               "10x larger than every other line. The replicate-level data "
               "shows the source is TSs282 Rep 2 = 583.00 (see "
               "REPLICATE_TYPO_CORRECTIONS); the mean of 583 plus nine ~6 mm "
               "reps yields 63.747."},
    {"sample": "TSs136", "trait": "Seed_Coat_Tannin",
     "raw_value": 7063.722747, "corrected_value": 706.3722747,
     "reason": "Decimal-point typo: 7063.7 is ~5x the next-highest value "
               "(1381.4, TSs44C; 99th pct = 1377.2). Corrected to 706.37, "
               "which slots near the panel median (805.8). Seed_Coat_Tannin "
               "has no replicate data; correction applied at the means "
               "level only."},
]

REPLICATE_TYPO_CORRECTIONS = [
    {"sample": "TSs282", "rep": 2, "trait": "Seed_Thickness",
     "raw_value": 583.00, "corrected_value": 5.83,
     "reason": "100x decimal-point typo. Other 9 reps for TSs282 sit at "
               "~5.8-6.5 mm. Correcting this row brings the genotype mean "
               "from 63.747 mm (the documented legacy outlier) into the "
               "normal range and lets the LMM fit cleanly."},
    {"sample": "TSs49", "rep": 5, "trait": "Seed_Thickness",
     "raw_value": 6024.00, "corrected_value": 6.024,
     "reason": "1000x decimal-point typo. Other 9 reps for TSs49 sit at "
               "~6 mm. The legacy means file already smoothed this one out "
               "via aggregation, but it must be corrected here before the "
               "LMM fit otherwise the genotype variance estimate blows up."},
]


# ---------------------------------------------------------------------------
# Long-format replicate loaders (private)
# ---------------------------------------------------------------------------

def _wet_reps() -> pd.DataFrame:
    """3 reps x 105 genotypes for the four wet-chemistry traits."""
    df = pd.read_excel(WET, sheet_name="Raw_ANF").rename(
        columns={"Genotypes": "sample", "Flavinoid": "Flavonoid"})
    return df.melt(id_vars=["sample", "Rep"],
                    value_vars=["Tannin", "Phenol", "Flavonoid", "Antioxidant"],
                    var_name="trait", value_name="value")


def _crude_protein_reps() -> pd.DataFrame:
    """2 reps x 106 genotypes."""
    df = pd.read_excel(WET, sheet_name="Crude_Protein").rename(
        columns={"Genotypes": "sample"})
    df["trait"] = "Crude_Protein"
    return df[["sample", "Rep", "trait", "Crude_Protein"]].rename(
        columns={"Crude_Protein": "value"})


def _moisture_reps() -> pd.DataFrame:
    """2-4 reps x 104 genotypes. NEW; covariate for Crude_Protein LMM."""
    df = pd.read_excel(WET, sheet_name="Moisture_Content").rename(
        columns={"Genotypes": "sample"})
    df["trait"] = "Moisture_Content"
    return df[["sample", "Rep", "trait", "MC"]].rename(columns={"MC": "value"})


def _seed_metrics_reps() -> pd.DataFrame:
    """10 reps x 105 genotypes for Length / Width / Thickness / Weight."""
    df = pd.read_excel(WET, sheet_name="Seed_Metrics").rename(
        columns={"Genotypes": "sample", "Reps": "Rep",
                 "Length": "Seed_Length", "Width": "Seed_Width",
                 "Thickness": "Seed_Thickness", "Weight": "Mass_of_Seeds"})
    return df.melt(id_vars=["sample", "Rep"],
                    value_vars=["Seed_Length", "Seed_Width",
                                "Seed_Thickness", "Mass_of_Seeds"],
                    var_name="trait", value_name="value")


def load_replicates(trait: Optional[str] = None,
                     apply_corrections: bool = True) -> pd.DataFrame:
    """Return long-format replicate-level data.

    Columns: sample, Rep, trait, value.
    If `trait` is given, return only rows for that trait.
    Means-only traits (Seed_Coat_Tannin + 3 oxalate fractions) are not
    in this table — they have no replicates anywhere.

    When `apply_corrections=True` (default), the REPLICATE_TYPO_CORRECTIONS
    table is applied — see that constant for the audit trail.
    """
    parts = [_wet_reps(), _crude_protein_reps(),
             _moisture_reps(), _seed_metrics_reps()]
    out = pd.concat(parts, ignore_index=True)
    out["sample"] = out["sample"].astype(str).str.strip()
    if apply_corrections:
        for c in REPLICATE_TYPO_CORRECTIONS:
            mask = ((out["sample"] == c["sample"]) &
                    (out["Rep"] == c["rep"]) &
                    (out["trait"] == c["trait"]))
            out.loc[mask, "value"] = c["corrected_value"]
    out = out.dropna(subset=["value"])
    if trait is not None:
        out = out[out["trait"] == trait].copy()
    return out


# ---------------------------------------------------------------------------
# Means loaders (private, legacy)
# ---------------------------------------------------------------------------

def _wet_means() -> pd.DataFrame:
    df = pd.read_excel(WET, sheet_name="Means").rename(
        columns={"Genotypes": "sample", "Flavinoid": "Flavonoid"})
    return df[["sample"] + TRAITS_BIOCHEM]


def _seed_metrics_means() -> pd.DataFrame:
    """Legacy Seed_Metrics.xlsx — still the only source for Seed_Coat_Tannin."""
    df = pd.read_excel(SEED, sheet_name="Data_on_Seeds", header=1).rename(
        columns={"Genotypes": "sample", "Seed Width": "Seed_Width"})
    return df[["sample"] + TRAITS_SEED]


def _j46_means() -> pd.DataFrame:
    df = pd.read_excel(J46).rename(columns={"Genotypes": "sample",
                                              "Crude Protein": "Crude_Protein"})
    return df[["sample"] + TRAITS_PROT_OX]


def load_means(apply_corrections: bool = True) -> pd.DataFrame:
    """Legacy means table — sample x 13 traits.

    Uses the Means sheet of the new wet-chemistry file (same as before),
    the legacy Seed_Metrics.xlsx (still the only source for
    Seed_Coat_Tannin), and Just_46_Accessions.xlsx (for the four
    protein-and-oxalate traits, n=46).

    Use load_blups() instead when you want REML-corrected genotype values.
    """
    wet = _wet_means()
    seed = _seed_metrics_means()
    j46 = _j46_means()
    samples = sorted(set(wet["sample"]) | set(seed["sample"]) | set(j46["sample"]))
    out = pd.DataFrame({"sample": samples})
    out = (out.merge(wet, on="sample", how="left")
              .merge(seed, on="sample", how="left")
              .merge(j46, on="sample", how="left"))
    out["sample"] = out["sample"].astype(str).str.strip()
    out = out.set_index("sample")
    if apply_corrections:
        for m in TYPO_CORRECTIONS:
            if m["sample"] in out.index and m["trait"] in out.columns:
                out.at[m["sample"], m["trait"]] = m["corrected_value"]
    return out[TRAITS_ALL]


# ---------------------------------------------------------------------------
# BLUP fitting (REML LMM via statsmodels)
# ---------------------------------------------------------------------------

def _fit_blup_one_trait(reps_df: pd.DataFrame, trait: str,
                         covariate_column: Optional[str] = None,
                         covariate_per_sample: Optional[pd.Series] = None,
                         rep_fixed_effect: bool = False,
                         ) -> tuple[pd.Series, dict]:
    """REML LMM for one trait. Returns (BLUP per sample, variance components).

    Model: value_ij = mu + (covariate_ij)? + (rep_j)? + g_i + e_ij

    Covariate handling:
      * If `covariate_column` is given, that column must already exist
        on `reps_df` (this is the rep-matched path).
      * If `covariate_per_sample` is given, the covariate is broadcast
        to every rep of a sample (sample-mean fallback path). The
        rep-matched path is preferred when MC differs across reps.

    Variance components reported:
      * H2_holland   line-mean broad-sense H2 = var_g / (var_g + var_e / n_rep_median).
                     Standard but biased under unbalanced n_rep.
      * H2_cullis    generalised line-mean H2 from average prediction-error
                     variance (Cullis et al. 2006). For sample i with n_i reps
                     the per-sample reliability is lambda_i = var_g /
                     (var_g + var_e / n_i); H2_cullis = mean(lambda_i).
      * LRT_Rep_p    likelihood-ratio p-value for the Rep main effect
                     (sensitivity for confounding day / plate / position).
                     Only populated when `rep_fixed_effect=True` is also
                     run separately and supplied via the caller.
    """
    import statsmodels.formula.api as smf
    from scipy.stats import chi2

    sub = reps_df[reps_df["trait"] == trait].copy()
    if covariate_column is not None and covariate_column in sub.columns:
        sub = sub.dropna(subset=[covariate_column])
        cov_terms = covariate_column
    elif covariate_per_sample is not None:
        sub["_cov"] = sub["sample"].map(covariate_per_sample)
        sub = sub.dropna(subset=["_cov"])
        cov_terms = "_cov"
    else:
        cov_terms = None

    parts = []
    if cov_terms is not None:
        parts.append(cov_terms)
    if rep_fixed_effect:
        parts.append("C(Rep)")
    formula = "value ~ " + (" + ".join(parts) if parts else "1")

    try:
        m = smf.mixedlm(formula, data=sub, groups=sub["sample"]).fit(reml=True)
    except Exception as e:
        return (sub.groupby("sample")["value"].mean(),
                {"trait": trait, "var_g": np.nan, "var_e": np.nan,
                 "n_rep_median": np.nan, "H2_holland": np.nan,
                 "H2_cullis": np.nan, "converged": False,
                 "on_boundary": False,
                 "lmm_status": f"failed: {e}"})

    # Evaluate the fixed part of the model at the population-mean
    # covariate value (when present) so BLUPs are on the same scale
    # as raw genotype means.
    fixed_offset = float(m.fe_params.iloc[0])
    if cov_terms is not None and cov_terms in m.fe_params.index:
        cov_mean = float(sub[cov_terms].mean())
        fixed_offset += float(m.fe_params[cov_terms]) * cov_mean
    if rep_fixed_effect:
        # average the Rep dummies; the BLUP is then evaluated at the
        # population-average rep effect, which is the conventional choice
        rep_terms = [c for c in m.fe_params.index if c.startswith("C(Rep)[")]
        if rep_terms:
            fixed_offset += float(np.mean([float(m.fe_params[c]) for c in rep_terms]))

    re = m.random_effects
    blups = pd.Series({g: float(rf.iloc[0]) for g, rf in re.items()})
    blup_per_sample = blups + fixed_offset

    var_g = float(m.cov_re.iloc[0, 0])
    var_e = float(m.scale)
    n_per_sample = sub.groupby("sample").size()
    n_rep_median = int(n_per_sample.median())

    if (var_g + var_e) > 0:
        H2_holland = var_g / (var_g + var_e / n_rep_median)
        # Cullis H2 via per-sample reliability lambda_i
        lambda_i = var_g / (var_g + var_e / n_per_sample)
        H2_cullis = float(lambda_i.mean())
        reliability = lambda_i  # for downstream de-regression
    else:
        H2_holland = np.nan
        H2_cullis = np.nan
        reliability = pd.Series(np.nan, index=n_per_sample.index)

    # REML boundary diagnostic: var_e at the lower numerical boundary
    # indicates the optimiser has hit the floor of the parameter space
    # and the point estimate of H2 should not be interpreted past the
    # second decimal place.
    on_boundary = bool(var_e < 1e-8 * (var_g + 1e-12) or var_g < 1e-8)

    info = {
        "trait": trait,
        "var_g": var_g, "var_e": var_e,
        "n_genotypes": int(n_per_sample.shape[0]),
        "n_rep_median": n_rep_median,
        "n_rep_min": int(n_per_sample.min()),
        "n_rep_max": int(n_per_sample.max()),
        "H2_holland": float(H2_holland),
        "H2_cullis": float(H2_cullis),
        "converged": bool(getattr(m, "converged", True)),
        "on_boundary": on_boundary,
        "lmm_status": "ok",
    }
    # attach reliability so the caller can save it per-sample
    info["_reliability"] = reliability
    return blup_per_sample, info


def _lrt_rep_effect(reps_df: pd.DataFrame, trait: str,
                     covariate_column: Optional[str] = None,
                     covariate_per_sample: Optional[pd.Series] = None
                     ) -> dict:
    """Likelihood-ratio test for a Rep fixed effect in the BLUP LMM.

    Fits the same LMM with and without `C(Rep)` and returns
    (LR-statistic, df, p-value). Used as a sensitivity check that
    omitting Rep does not absorb systematic measurement-run variation
    into the genotype random effect.
    """
    import statsmodels.formula.api as smf
    from scipy.stats import chi2

    sub = reps_df[reps_df["trait"] == trait].copy()
    if covariate_column is not None and covariate_column in sub.columns:
        sub = sub.dropna(subset=[covariate_column])
        cov_terms = covariate_column
    elif covariate_per_sample is not None:
        sub["_cov"] = sub["sample"].map(covariate_per_sample)
        sub = sub.dropna(subset=["_cov"])
        cov_terms = "_cov"
    else:
        cov_terms = None

    base_terms = [cov_terms] if cov_terms is not None else []
    f0 = "value ~ " + (" + ".join(base_terms) if base_terms else "1")
    f1 = "value ~ " + (" + ".join(base_terms + ["C(Rep)"]))

    try:
        # LRT is only well-defined under ML, not REML, when comparing
        # nested fixed-effect structures. Refit both with reml=False.
        m0 = smf.mixedlm(f0, data=sub, groups=sub["sample"]).fit(reml=False)
        m1 = smf.mixedlm(f1, data=sub, groups=sub["sample"]).fit(reml=False)
        n_levels = sub["Rep"].nunique()
        df = max(n_levels - 1, 1)
        lr = 2.0 * (m1.llf - m0.llf)
        p = float(chi2.sf(lr, df))
        return {"trait": trait, "LRT_Rep_chi2": float(lr),
                "LRT_Rep_df": int(df), "LRT_Rep_p": p}
    except Exception as e:
        return {"trait": trait, "LRT_Rep_chi2": np.nan,
                "LRT_Rep_df": np.nan, "LRT_Rep_p": np.nan,
                "LRT_Rep_status": f"failed: {e}"}


def _attach_rep_matched_mc(reps: pd.DataFrame) -> pd.DataFrame:
    """Add a per-row `MC` column to the replicate table by joining on
    (sample, Rep). Reps without a matching MC value get NaN; the
    Crude_Protein LMM drops those rows. 208/212 (98%) of Crude_Protein
    observations have a rep-matched MC value in the panel.
    """
    mc_rows = reps[reps["trait"] == "Moisture_Content"][["sample", "Rep", "value"]]
    mc_rows = mc_rows.rename(columns={"value": "MC"})
    return reps.merge(mc_rows, on=["sample", "Rep"], how="left")


def load_blups(use_mc_covariate_for_protein: bool = True,
                apply_corrections: bool = True,
                return_variance_components: bool = False,
                return_reliability: bool = False,
                run_rep_lrt: bool = False,
                ) -> pd.DataFrame:
    """Return one wide DataFrame: sample x 13 traits.

    For the 9 replicated traits, values are REML BLUPs from the LMM
    `value ~ 1 | sample` (or `value ~ MC + 1 | sample` for Crude_Protein
    when `use_mc_covariate_for_protein` is True; the MC covariate is
    rep-matched to each Crude_Protein observation, not averaged across
    reps). For the 4 means-only traits (Seed_Coat_Tannin + 3 oxalates),
    the means table is passed through.

    Returns:
        BLUP table (always)
        + per-trait variance-component table when
          return_variance_components=True
        + per-(sample, trait) reliability table when
          return_reliability=True (used by Garrick de-regression for
          downstream GWAS / GBLUP work)

    If `run_rep_lrt=True`, the variance-component table also carries a
    LRT_Rep_p column from a likelihood-ratio test for a Rep main
    effect. This is a sensitivity check that omitting Rep does not
    absorb measurement-run variation into the genotype random effect.
    """
    reps = load_replicates()
    means = load_means(apply_corrections=False)
    samples = sorted(set(reps["sample"].unique()) | set(means.index))
    blup_df = pd.DataFrame(index=samples)
    var_rows = []
    reliab_rows = []

    # Rep-matched MC (preferred over sample-mean MC; the legacy path
    # collapsed across-rep variation in MC, which is the very thing the
    # covariate is supposed to absorb).
    if use_mc_covariate_for_protein:
        reps_with_mc = _attach_rep_matched_mc(reps)
    else:
        reps_with_mc = reps

    for t in TRAITS_REPLICATED:
        if t == "Crude_Protein" and use_mc_covariate_for_protein:
            blup, vc = _fit_blup_one_trait(reps_with_mc, t,
                                             covariate_column="MC")
            if run_rep_lrt:
                vc.update(_lrt_rep_effect(reps_with_mc, t,
                                            covariate_column="MC"))
        else:
            blup, vc = _fit_blup_one_trait(reps, t)
            if run_rep_lrt:
                vc.update(_lrt_rep_effect(reps, t))
        blup_df[t] = blup.reindex(samples)

        # peel reliability series out of the info dict before saving
        rel = vc.pop("_reliability", None)
        if rel is not None and isinstance(rel, pd.Series):
            for s, lam in rel.items():
                reliab_rows.append({"sample": s, "trait": t,
                                     "n_rep": int(reps[(reps.sample == s) &
                                                       (reps.trait == t)].shape[0]),
                                     "reliability": float(lam)})
        var_rows.append(vc)

    for t in TRAITS_MEANS_ONLY:
        blup_df[t] = means[t].reindex(samples)
        var_rows.append({
            "trait": t, "var_g": np.nan, "var_e": np.nan,
            "n_genotypes": int(means[t].notna().sum()),
            "n_rep_median": 1, "n_rep_min": 1, "n_rep_max": 1,
            "H2_holland": np.nan, "H2_cullis": np.nan,
            "converged": True, "on_boundary": False,
            "lmm_status": "means_only_no_replicates",
        })

    blup_df.index.name = "sample"
    blup_df.index = blup_df.index.astype(str).str.strip()
    if apply_corrections:
        for m in TYPO_CORRECTIONS:
            if m["sample"] in blup_df.index and m["trait"] in blup_df.columns:
                blup_df.at[m["sample"], m["trait"]] = m["corrected_value"]
    blup_df = blup_df[TRAITS_ALL]

    out: list = [blup_df]
    if return_variance_components:
        out.append(pd.DataFrame(var_rows))
    if return_reliability:
        out.append(pd.DataFrame(reliab_rows))
    if len(out) == 1:
        return out[0]
    return tuple(out)  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# Canonical loader — defaults to BLUPs (the new behaviour). Old callers
# get an automatic upgrade unless they explicitly ask for means.
# ---------------------------------------------------------------------------

def load_phenotypes(use_blups: bool = True,
                     apply_corrections: bool = True) -> pd.DataFrame:
    """Default entry point. Returns BLUP-corrected genotype values where
    replicates exist; means as a fallback for the four means-only
    traits.

    Set `use_blups=False` for the legacy means-only behaviour.
    """
    if use_blups:
        return load_blups(use_mc_covariate_for_protein=True,
                           apply_corrections=apply_corrections)
    return load_means(apply_corrections=apply_corrections)


def heritability_table() -> pd.DataFrame:
    """Per-trait Holland H2 from the REML LMM, plus variance components.

    Means-only traits return NaN for H2 (no replicates to estimate the
    error variance).
    """
    _, vc = load_blups(return_variance_components=True)
    return vc


def trait_n_table(df: pd.DataFrame) -> pd.DataFrame:
    """Per-trait sample size and descriptive stats; useful for figure captions."""
    rows = []
    for t in TRAITS_ALL:
        if t not in df.columns:
            continue
        s = df[t].dropna()
        rows.append({
            "trait": t,
            "n": len(s),
            "mean": float(s.mean()),
            "sd": float(s.std(ddof=1)),
            "min": float(s.min()),
            "max": float(s.max()),
        })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print("=== BLUP fit (default load_phenotypes) ===")
    df = load_phenotypes()
    print(f"loaded {df.shape[0]} samples x {df.shape[1]} traits")
    print(trait_n_table(df).to_string(index=False))
    print()
    print("=== variance components + Holland H2 ===")
    print(heritability_table().to_string(index=False))
