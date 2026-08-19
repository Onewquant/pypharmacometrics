"""
Infliximab pharmacogenomic association analysis
================================================

Purpose
-------
Analyze the association between candidate genetic variants and infliximab
clearance (CL) or anti-drug antibody (ADA), reflecting the following rules:

1. Analyze infliximab only.
2. Analyze IND, MAINT, and ALL separately.
3. Analyze CL and ADA as separate endpoints.
4. For each PHASE x END_POINT x GENETIC_MODEL stratum, apply Benjamini-
   Hochberg FDR correction across the rsIDs tested in that stratum.
5. CL is analyzed after natural-log transformation:
       log(CL) ~ genotype + SEX + WEIGHT + ALBUMIN + ADA
   The genotype effect is back-transformed and reported as a geometric mean
   ratio (GMR).
6. ADA is binary and is analyzed using logistic regression:
       ADA ~ genotype + SEX + WEIGHT + ALBUMIN
   The genotype effect is reported as an odds ratio (OR).
7. Three prespecified genetic models are produced as separate analysis
   families:
       additive  : dosage 0/1/2, effect per variant allele
       dominant  : carrier (1/2) vs non-carrier (0)
       recessive : variant homozygote (2) vs others (0/1)

Important interpretation of ALL
-------------------------------
ALL is a patient-level whole-period summary, following the structure of the
previous scripts: mean CL and mean continuous covariates across IND/MAINT,
and maximum ADA across IND/MAINT. It is not a repeated-measures interaction
model. IND and MAINT remain the phase-specific primary analyses.

The script automatically finds the input CSV files in the same folder as this
.py file, including filenames with suffixes such as '(3)' or '(4)'.
"""

from __future__ import annotations

import json
import re
import warnings
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import shapiro
from statsmodels.stats.multitest import multipletests


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

DRUG = "infliximab"
PHASES = ("IND", "MAINT", "ALL")
ENDPOINTS = ("CL", "ADA")
GENETIC_MODELS = ("additive", "dominant", "recessive")

MIN_BINARY_GROUP_N = 8
MIN_ADDITIVE_LEVEL_N = 8
MIN_TOTAL_N = 16
FDR_ALPHA = 0.05

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "infliximab_stratified_fdr_results"

RSID_GENE_DICT = {
    "rs9828223": "CD96",
    "rs2097432": "HLA-DQA1",
    "rs396991": "FCGR3A",
    "rs1800629": "TNF (TNF-alpha)",
    "rs4149570": "TNFRSF1A",
    "rs3397": "TNFRSF1B",
    "rs1061624": "TNFRSF1B",
    "rs5030728": "TLR4",
    "rs3804099": "TLR2",
    "rs10499563": "IL6",
    "rs2275913": "IL17A",
    "rs1800872": "IL10",
    "rs3024505": "IL10",
    "rs361525": "TNF (TNF-alpha)",
    "rs767455": "TNF (TNF-alpha)",
    "rs1061622": "TNF (TNF-alpha)",
    "rs765249238": "SLCO2A1",
    "rs776813259": "SLCO2A1",
}

PHASE_ORDER = {"IND": 0, "MAINT": 1, "ALL": 2}
ENDPOINT_ORDER = {"CL": 0, "ADA": 1}
MODEL_ORDER = {"additive": 0, "dominant": 1, "recessive": 2}


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def find_input_file(exact_name: str, glob_pattern: str) -> Path:
    """Find an exact filename first, then a suffixed variant in SCRIPT_DIR."""
    exact_path = SCRIPT_DIR / exact_name
    if exact_path.exists():
        return exact_path

    candidates = sorted(
        SCRIPT_DIR.glob(glob_pattern),
        key=lambda p: (p.stat().st_mtime, p.name),
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            f"Could not find '{exact_name}' or a file matching "
            f"'{glob_pattern}' in {SCRIPT_DIR}"
        )
    return candidates[0]


def normalize_uid(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
    )


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def base_rsid(genotype_column: str) -> str:
    match = re.match(r"^(rs\d+)", genotype_column)
    return match.group(1) if match else genotype_column


def parse_alleles(genotype_column: str) -> tuple[str, str]:
    match = re.search(r"\(0=([^,]+),\s*1=([^\)]+)\)", genotype_column)
    if not match:
        return "", ""
    return match.group(1).strip(), match.group(2).strip()


def safe_float(value: Any) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return np.nan
    return value if np.isfinite(value) else np.nan


def safe_exp(value: Any) -> float:
    value = safe_float(value)
    if not np.isfinite(value) or value > 709 or value < -745:
        return np.nan
    return float(np.exp(value))


def geometric_mean_ci(values: Iterable[float]) -> str:
    x = pd.Series(values, dtype="float64").dropna()
    x = x[x > 0]
    n = len(x)
    if n == 0:
        return "NA"

    lx = np.log(x)
    mean_log = lx.mean()
    gmean = np.exp(mean_log)

    if n == 1:
        return f"{gmean:.4g}, n={n}"

    se = lx.std(ddof=1) / np.sqrt(n)
    lower = np.exp(mean_log - 1.96 * se)
    upper = np.exp(mean_log + 1.96 * se)
    return f"{gmean:.4g} ({lower:.4g}-{upper:.4g}), n={n}"


def binary_count(values: Iterable[float]) -> str:
    x = pd.Series(values).dropna()
    n = len(x)
    if n == 0:
        return "NA"
    events = int((x == 1).sum())
    return f"{events}/{n} ({100 * events / n:.1f}%)"


def format_group_estimates(
    df: pd.DataFrame,
    endpoint: str,
    predictor: str,
    genetic_model: str,
) -> str:
    """Return JSON-formatted unadjusted group summaries."""
    summaries: dict[str, str] = {}

    if genetic_model == "additive":
        groups = [(0, "dosage_0"), (1, "dosage_1"), (2, "dosage_2")]
    else:
        groups = [(0, "group_0"), (1, "group_1")]

    for value, label in groups:
        endpoint_values = df.loc[df[predictor] == value, endpoint]
        if endpoint == "CL":
            summaries[label] = geometric_mean_ci(endpoint_values)
        else:
            summaries[label] = binary_count(endpoint_values)

    return json.dumps(summaries, ensure_ascii=False)


def build_genetic_predictor(
    genotype_df: pd.DataFrame,
    genetic_model: str,
) -> tuple[pd.DataFrame, str, str, str]:
    """Create the predictor and return labels/effect description."""
    out = genotype_df.copy()
    dosage = numeric(out["GENOTYPE_DOSAGE"])
    dosage = dosage.where(dosage.isin([0, 1, 2]))
    out["GENOTYPE_DOSAGE"] = dosage

    if genetic_model == "additive":
        predictor = "GENOTYPE_EFFECT"
        out[predictor] = dosage
        group_0_label = "dosage 0"
        group_1_label = "per variant allele"
    elif genetic_model == "dominant":
        predictor = "GENOTYPE_EFFECT"
        out[predictor] = np.where(
            dosage.isin([1, 2]),
            1.0,
            np.where(dosage == 0, 0.0, np.nan),
        )
        group_0_label = "non-carrier (dosage 0)"
        group_1_label = "carrier (dosage 1/2)"
    elif genetic_model == "recessive":
        predictor = "GENOTYPE_EFFECT"
        out[predictor] = np.where(
            dosage == 2,
            1.0,
            np.where(dosage.isin([0, 1]), 0.0, np.nan),
        )
        group_0_label = "others (dosage 0/1)"
        group_1_label = "variant homozygote (dosage 2)"
    else:
        raise ValueError(f"Unknown genetic model: {genetic_model}")

    return out, predictor, group_0_label, group_1_label


def model_eligibility(
    df: pd.DataFrame,
    predictor: str,
    genetic_model: str,
) -> tuple[bool, str, dict[str, int]]:
    counts = (
        df[predictor]
        .dropna()
        .value_counts()
        .sort_index()
        .astype(int)
        .to_dict()
    )
    counts_str = {str(k): int(v) for k, v in counts.items()}

    if len(df) < MIN_TOTAL_N:
        return False, f"TOTAL_N_LT_{MIN_TOTAL_N}", counts_str

    if genetic_model in {"dominant", "recessive"}:
        n0 = int(counts.get(0.0, counts.get(0, 0)))
        n1 = int(counts.get(1.0, counts.get(1, 0)))
        if n0 < MIN_BINARY_GROUP_N or n1 < MIN_BINARY_GROUP_N:
            return (
                False,
                f"GROUP_N_LT_{MIN_BINARY_GROUP_N}",
                counts_str,
            )
        return True, "OK", counts_str

    observed_counts = [int(v) for v in counts.values()]
    if len(observed_counts) < 2:
        return False, "FEWER_THAN_2_DOSAGE_LEVELS", counts_str
    if min(observed_counts) < MIN_ADDITIVE_LEVEL_N:
        return (
            False,
            f"OBSERVED_DOSAGE_LEVEL_N_LT_{MIN_ADDITIVE_LEVEL_N}",
            counts_str,
        )
    return True, "OK", counts_str


def extract_covariate_effects(
    fit: Any,
    covariates: list[str],
    endpoint: str,
) -> str:
    effects: dict[str, dict[str, float | str]] = {}
    confidence = fit.conf_int()

    for covariate in covariates:
        if covariate not in fit.params.index:
            continue
        p_value = safe_float(fit.pvalues.get(covariate, np.nan))
        if not np.isfinite(p_value) or p_value >= 0.05:
            continue

        beta = safe_float(fit.params[covariate])
        lower_beta = safe_float(confidence.loc[covariate, 0])
        upper_beta = safe_float(confidence.loc[covariate, 1])

        if endpoint == "CL":
            effect_name = "GMR"
            effect = safe_exp(beta)
            lower = safe_exp(lower_beta)
            upper = safe_exp(upper_beta)
        else:
            effect_name = "OR"
            effect = safe_exp(beta)
            lower = safe_exp(lower_beta)
            upper = safe_exp(upper_beta)

        effects[covariate] = {
            "effect_name": effect_name,
            "effect": round(safe_float(effect), 5),
            "ci_lower": round(safe_float(lower), 5),
            "ci_upper": round(safe_float(upper), 5),
            "p_value": round(p_value, 6),
        }

    return json.dumps(effects, ensure_ascii=False)


def fit_association_model(
    model_df: pd.DataFrame,
    endpoint: str,
    predictor: str,
    covariates: list[str],
) -> dict[str, Any]:
    x_columns = [predictor] + covariates
    X = sm.add_constant(model_df[x_columns], has_constant="add")

    result: dict[str, Any] = {
        "STATUS": "MODEL_ERROR",
        "MODEL_METHOD": "",
        "BETA": np.nan,
        "SE": np.nan,
        "EFFECT_NAME": "",
        "EFFECT": np.nan,
        "PERCENT_CHANGE": np.nan,
        "CI_LOWER": np.nan,
        "CI_UPPER": np.nan,
        "P_VALUE": np.nan,
        "P_VALUE_MODEL_RAW": np.nan,
        "ADJ_R2": np.nan,
        "RESIDUAL_SHAPIRO_P": np.nan,
        "CONVERGED": np.nan,
        "EVENTS": np.nan,
        "EVENTS_PER_PREDICTOR": np.nan,
        "LOW_EVENT_WARNING": False,
        "SIG_COVAR": "{}",
        "MODEL_WARNING": "",
    }

    captured_warnings: list[str] = []

    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")

            if endpoint == "CL":
                y = np.log(model_df["CL"])
                fit = sm.OLS(y, X).fit(cov_type="HC3", use_t=True)
                result["MODEL_METHOD"] = "OLS on log(CL), HC3 robust SE"
                result["EFFECT_NAME"] = "GMR"
                result["ADJ_R2"] = safe_float(fit.rsquared_adj)

                residuals = pd.Series(fit.resid).dropna()
                if 3 <= len(residuals) <= 5000:
                    result["RESIDUAL_SHAPIRO_P"] = safe_float(
                        shapiro(residuals).pvalue
                    )

            else:
                y = model_df["ADA"]
                fit = sm.Logit(y, X).fit(disp=False, maxiter=200)
                result["MODEL_METHOD"] = "Logistic regression (maximum likelihood)"
                result["EFFECT_NAME"] = "OR"
                converged = bool(fit.mle_retvals.get("converged", False))
                result["CONVERGED"] = converged

                events = int((y == 1).sum())
                predictor_count = len(x_columns)
                result["EVENTS"] = events
                result["EVENTS_PER_PREDICTOR"] = (
                    events / predictor_count if predictor_count > 0 else np.nan
                )
                result["LOW_EVENT_WARNING"] = events < 10 * predictor_count

            captured_warnings = [str(w.message) for w in caught]

        beta = safe_float(fit.params.get(predictor, np.nan))
        se = safe_float(fit.bse.get(predictor, np.nan))
        p_value = safe_float(fit.pvalues.get(predictor, np.nan))
        confidence = fit.conf_int()
        lower_beta = safe_float(confidence.loc[predictor, 0])
        upper_beta = safe_float(confidence.loc[predictor, 1])

        effect = safe_exp(beta)
        lower = safe_exp(lower_beta)
        upper = safe_exp(upper_beta)

        status = "OK" if np.isfinite(p_value) else "NONFINITE_P_VALUE"
        if endpoint == "ADA":
            warning_text = " | ".join(captured_warnings).lower()
            unstable = (
                not bool(result.get("CONVERGED", False))
                or not np.isfinite(beta)
                or not np.isfinite(se)
                or abs(beta) > 10
                or se > 10
                or "separation" in warning_text
                or "failed to converge" in warning_text
                or "singular" in warning_text
            )
            if unstable:
                status = "LOGISTIC_UNSTABLE_OR_SEPARATED"

        result.update({
            "STATUS": status,
            "BETA": beta,
            "SE": se,
            "EFFECT": safe_float(effect),
            "PERCENT_CHANGE": safe_float((effect - 1) * 100),
            "CI_LOWER": safe_float(lower),
            "CI_UPPER": safe_float(upper),
            "P_VALUE_MODEL_RAW": p_value,
            "P_VALUE": p_value if status == "OK" else np.nan,
            "SIG_COVAR": extract_covariate_effects(
                fit=fit,
                covariates=covariates,
                endpoint=endpoint,
            ),
            "MODEL_WARNING": " | ".join(captured_warnings),
        })
        return result

    except Exception as exc:  # keep a transparent row rather than silently dropping it
        result["MODEL_WARNING"] = f"{type(exc).__name__}: {exc}"
        return result


# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

def load_data() -> tuple[pd.DataFrame, pd.DataFrame, Path, Path]:
    genotype_path = find_input_file(
        exact_name="rsid_dosage_matrix_with_alleles.csv",
        glob_pattern="rsid_dosage_matrix_with_alleles*.csv",
    )
    endpoint_path = find_input_file(
        exact_name="for_genomics_df(all_drugs).csv",
        glob_pattern="for_genomics_df(all_drugs)*.csv",
    )

    genotype_df = pd.read_csv(genotype_path)
    endpoint_df = pd.read_csv(endpoint_path)

    required_genotype = {"UID", "genomics_group"}
    missing_genotype = required_genotype - set(genotype_df.columns)
    if missing_genotype:
        raise ValueError(
            f"Genotype file is missing required columns: {sorted(missing_genotype)}"
        )

    required_endpoint = {"UID", "DRUG", "PHASE", "CL", "ADA", "SEX"}
    missing_endpoint = required_endpoint - set(endpoint_df.columns)
    if missing_endpoint:
        raise ValueError(
            f"Endpoint file is missing required columns: {sorted(missing_endpoint)}"
        )

    genotype_df = genotype_df.loc[genotype_df["UID"].notna()].copy()
    genotype_df["UID"] = normalize_uid(genotype_df["UID"])

    endpoint_df = endpoint_df.loc[endpoint_df["UID"].notna()].copy()
    endpoint_df["UID"] = normalize_uid(endpoint_df["UID"])
    endpoint_df["DRUG"] = endpoint_df["DRUG"].astype(str).str.strip().str.lower()
    endpoint_df["PHASE_ORIGINAL"] = endpoint_df["PHASE"].astype(str)
    endpoint_df["PHASE"] = (
        endpoint_df["PHASE"]
        .astype(str)
        .str.strip()
        .str.upper()
        .str.split("_")
        .str[0]
    )

    if "WEIGHT" not in endpoint_df.columns:
        if "WT" not in endpoint_df.columns:
            raise ValueError("Endpoint file needs WEIGHT or WT.")
        endpoint_df["WEIGHT"] = endpoint_df["WT"]

    if "ALBUMIN" not in endpoint_df.columns:
        if "ALB" not in endpoint_df.columns:
            raise ValueError("Endpoint file needs ALBUMIN or ALB.")
        endpoint_df["ALBUMIN"] = endpoint_df["ALB"]

    for column in ["CL", "ADA", "SEX", "WEIGHT", "ALBUMIN"]:
        endpoint_df[column] = numeric(endpoint_df[column])

    endpoint_df = endpoint_df.loc[
        (endpoint_df["DRUG"] == DRUG)
        & endpoint_df["PHASE"].isin(["IND", "MAINT"])
    ].copy()

    # Validate genotype uniqueness. Duplicate identical UIDs are safely collapsed;
    # conflicting duplicate genotype rows are rejected.
    rsid_columns = [c for c in genotype_df.columns if re.match(r"^rs\d+", c)]
    if not rsid_columns:
        raise ValueError("No rsID dosage columns were found in the genotype file.")

    duplicate_uid = genotype_df["UID"].duplicated(keep=False)
    if duplicate_uid.any():
        conflict = (
            genotype_df.loc[duplicate_uid, ["UID"] + rsid_columns]
            .groupby("UID", dropna=False)
            .nunique(dropna=False)
            .gt(1)
            .any(axis=1)
        )
        if conflict.any():
            bad_uids = conflict[conflict].index.tolist()[:10]
            raise ValueError(
                "Conflicting duplicate genotype rows were found for UID(s): "
                f"{bad_uids}"
            )
        genotype_df = genotype_df.drop_duplicates("UID", keep="first")

    return genotype_df, endpoint_df, genotype_path, endpoint_path


def prepare_phase_table(endpoint_df: pd.DataFrame, phase: str) -> pd.DataFrame:
    if phase == "ALL":
        phase_df = endpoint_df.loc[
            endpoint_df["PHASE"].isin(["IND", "MAINT"])
        ].copy()
    else:
        phase_df = endpoint_df.loc[endpoint_df["PHASE"] == phase].copy()

    # Existing analysis definition: one patient-level summary per phase/ALL.
    patient_df = (
        phase_df.groupby("UID", as_index=False)
        .agg(
            CL=("CL", "mean"),
            ADA=("ADA", "max"),
            SEX=("SEX", "first"),
            WEIGHT=("WEIGHT", "mean"),
            ALBUMIN=("ALBUMIN", "mean"),
            N_SOURCE_ROWS=("UID", "size"),
        )
    )

    # CL must be positive for logarithmic analysis. Invalid CL is retained as NaN
    # so that ADA analysis can still use the same patient.
    patient_df.loc[patient_df["CL"] <= 0, "CL"] = np.nan
    return patient_df


# -----------------------------------------------------------------------------
# Main analysis
# -----------------------------------------------------------------------------

def run_analysis() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    genotype_df, endpoint_df, genotype_path, endpoint_path = load_data()
    rsid_columns = [c for c in genotype_df.columns if re.match(r"^rs\d+", c)]

    result_rows: list[dict[str, Any]] = []

    for phase in PHASES:
        patient_df = prepare_phase_table(endpoint_df, phase)
        total_phase_uids = int(patient_df["UID"].nunique())

        for genotype_column in rsid_columns:
            rsid = base_rsid(genotype_column)
            reference_allele, variant_allele = parse_alleles(genotype_column)

            one_genotype = genotype_df[["UID", genotype_column]].copy()
            one_genotype = one_genotype.rename(
                columns={genotype_column: "GENOTYPE_DOSAGE"}
            )

            for genetic_model in GENETIC_MODELS:
                modeled_genotype, predictor, group_0_label, group_1_label = (
                    build_genetic_predictor(one_genotype, genetic_model)
                )

                merged = patient_df.merge(
                    modeled_genotype[["UID", "GENOTYPE_DOSAGE", predictor]],
                    on="UID",
                    how="inner",
                    validate="one_to_one",
                )

                for endpoint in ENDPOINTS:
                    covariates = ["SEX", "WEIGHT", "ALBUMIN"]
                    if endpoint == "CL":
                        covariates.append("ADA")

                    required = [endpoint, predictor] + covariates
                    model_df = merged.dropna(subset=required).copy()

                    if endpoint == "CL":
                        model_df = model_df.loc[model_df["CL"] > 0].copy()
                    else:
                        model_df = model_df.loc[model_df["ADA"].isin([0, 1])].copy()

                    descriptive_df = merged.dropna(subset=[endpoint, predictor]).copy()
                    if endpoint == "CL":
                        descriptive_df = descriptive_df.loc[
                            descriptive_df["CL"] > 0
                        ].copy()
                    else:
                        descriptive_df = descriptive_df.loc[
                            descriptive_df["ADA"].isin([0, 1])
                        ].copy()

                    group_estimates = format_group_estimates(
                        descriptive_df,
                        endpoint=endpoint,
                        predictor=predictor,
                        genetic_model=genetic_model,
                    )

                    eligible, eligibility_status, counts = model_eligibility(
                        model_df,
                        predictor=predictor,
                        genetic_model=genetic_model,
                    )

                    outcome_unique = model_df[endpoint].nunique(dropna=True)
                    if outcome_unique < 2:
                        eligible = False
                        eligibility_status = "NO_OUTCOME_VARIATION"

                    base_row: dict[str, Any] = {
                        "DRUG": DRUG,
                        "PHASE": phase,
                        "END_POINT": endpoint,
                        "GENETIC_MODEL": genetic_model,
                        "RSID": rsid,
                        "GENE": RSID_GENE_DICT.get(rsid, ""),
                        "GENOTYPE_COLUMN": genotype_column,
                        "REFERENCE_ALLELE_DOSAGE_0": reference_allele,
                        "VARIANT_ALLELE_DOSAGE_1": variant_allele,
                        "GROUP_0_LABEL": group_0_label,
                        "GROUP_1_LABEL": group_1_label,
                        "GROUP_ESTIMATES": group_estimates,
                        "GROUP_COUNTS_MODEL_COMPLETE": json.dumps(
                            counts, ensure_ascii=False
                        ),
                        "TOTAL_PHASE_UIDS": total_phase_uids,
                        "ENDPOINT_AVAILABLE_N": int(len(descriptive_df)),
                        "MODEL_N": int(len(model_df)),
                        "DATA_AVAILABILITY": (
                            f"{len(model_df)}/{total_phase_uids} "
                            f"({100 * len(model_df) / total_phase_uids:.1f}%)"
                            if total_phase_uids > 0
                            else "0/0 (NA)"
                        ),
                        "COVARIATES": ", ".join(covariates),
                        "FDR_FAMILY": (
                            f"{DRUG}|{phase}|{endpoint}|{genetic_model}"
                        ),
                        "FDR_METHOD": "Benjamini-Hochberg",
                        "FDR_ALPHA": FDR_ALPHA,
                        "P_VALUE_FDR": np.nan,
                        "N_TESTS_IN_FDR_FAMILY": 0,
                        "SIGNIFICANT_FDR_0_05": False,
                    }

                    if not eligible:
                        base_row.update({
                            "STATUS": eligibility_status,
                            "MODEL_METHOD": "",
                            "BETA": np.nan,
                            "SE": np.nan,
                            "EFFECT_NAME": (
                                "GMR" if endpoint == "CL" else "OR"
                            ),
                            "EFFECT": np.nan,
                            "PERCENT_CHANGE": np.nan,
                            "CI_LOWER": np.nan,
                            "CI_UPPER": np.nan,
                            "P_VALUE": np.nan,
                            "P_VALUE_MODEL_RAW": np.nan,
                            "ADJ_R2": np.nan,
                            "RESIDUAL_SHAPIRO_P": np.nan,
                            "CONVERGED": np.nan,
                            "EVENTS": (
                                int((model_df["ADA"] == 1).sum())
                                if endpoint == "ADA"
                                else np.nan
                            ),
                            "EVENTS_PER_PREDICTOR": np.nan,
                            "LOW_EVENT_WARNING": (
                                bool(
                                    endpoint == "ADA"
                                    and int((model_df["ADA"] == 1).sum())
                                    < 10 * (1 + len(covariates))
                                )
                            ),
                            "SIG_COVAR": "{}",
                            "MODEL_WARNING": "",
                        })
                    else:
                        base_row.update(
                            fit_association_model(
                                model_df=model_df,
                                endpoint=endpoint,
                                predictor=predictor,
                                covariates=covariates,
                            )
                        )

                    result_rows.append(base_row)

    results = pd.DataFrame(result_rows)

    # FDR is applied independently within PHASE x END_POINT x GENETIC_MODEL,
    # across the candidate rsIDs that produced a valid p-value.
    fdr_group_cols = ["PHASE", "END_POINT", "GENETIC_MODEL"]
    for _, sub_df in results.groupby(fdr_group_cols, dropna=False):
        valid_mask = (
            sub_df["STATUS"].eq("OK")
            & pd.to_numeric(sub_df["P_VALUE"], errors="coerce").notna()
        )
        valid_index = sub_df.index[valid_mask]
        n_tests = len(valid_index)
        results.loc[sub_df.index, "N_TESTS_IN_FDR_FAMILY"] = n_tests

        if n_tests > 0:
            adjusted = multipletests(
                results.loc[valid_index, "P_VALUE"].astype(float).values,
                alpha=FDR_ALPHA,
                method="fdr_bh",
            )[1]
            results.loc[valid_index, "P_VALUE_FDR"] = adjusted

    results["SIGNIFICANT_FDR_0_05"] = (
        pd.to_numeric(results["P_VALUE_FDR"], errors="coerce") <= FDR_ALPHA
    )

    results["_PHASE_ORDER"] = results["PHASE"].map(PHASE_ORDER)
    results["_ENDPOINT_ORDER"] = results["END_POINT"].map(ENDPOINT_ORDER)
    results["_MODEL_ORDER"] = results["GENETIC_MODEL"].map(MODEL_ORDER)
    results = results.sort_values(
        [
            "_PHASE_ORDER",
            "_ENDPOINT_ORDER",
            "_MODEL_ORDER",
            "P_VALUE_FDR",
            "P_VALUE",
            "RSID",
        ],
        na_position="last",
    ).drop(columns=["_PHASE_ORDER", "_ENDPOINT_ORDER", "_MODEL_ORDER"])

    results["P_VALUE_VALID_FOR_FDR"] = np.where(
        results["STATUS"].eq("OK"),
        results["P_VALUE"],
        np.nan,
    )

    family_summary = (
        results.groupby(
            ["DRUG", "PHASE", "END_POINT", "GENETIC_MODEL", "FDR_FAMILY"],
            as_index=False,
        )
        .agg(
            TOTAL_CANDIDATE_RSIDS=("RSID", "nunique"),
            N_VALID_TESTS=("N_TESTS_IN_FDR_FAMILY", "max"),
            N_FDR_SIGNIFICANT=("SIGNIFICANT_FDR_0_05", "sum"),
            MIN_RAW_P=("P_VALUE_VALID_FOR_FDR", "min"),
            MIN_FDR_P=("P_VALUE_FDR", "min"),
        )
    )

    metadata = {
        "genotype_file": str(genotype_path),
        "endpoint_file": str(endpoint_path),
        "output_directory": str(OUTPUT_DIR),
        "drug": DRUG,
        "phases": list(PHASES),
        "endpoints": list(ENDPOINTS),
        "genetic_models": list(GENETIC_MODELS),
        "fdr_family_definition": (
            "Separate BH-FDR within each PHASE x END_POINT x GENETIC_MODEL "
            "across candidate rsIDs with valid p-values"
        ),
        "cl_model": "log(CL) ~ genotype + SEX + WEIGHT + ALBUMIN + ADA",
        "ada_model": "ADA ~ genotype + SEX + WEIGHT + ALBUMIN (maximum-likelihood logistic regression; unstable/separated fits excluded from FDR)",
        "all_phase_definition": (
            "Patient-level mean CL/continuous covariates and maximum ADA "
            "across IND and MAINT; not a repeated-measures interaction model"
        ),
        "minimum_binary_group_n": MIN_BINARY_GROUP_N,
        "minimum_additive_observed_level_n": MIN_ADDITIVE_LEVEL_N,
        "minimum_total_n": MIN_TOTAL_N,
        "fdr_alpha": FDR_ALPHA,
    }

    return results, family_summary, metadata


def save_outputs(
    results: pd.DataFrame,
    family_summary: pd.DataFrame,
    metadata: dict[str, Any],
) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    combined_path = OUTPUT_DIR / "infliximab_all_genetic_models_stratified_fdr.csv"
    significant_path = OUTPUT_DIR / "infliximab_fdr_significant_results.csv"
    family_path = OUTPUT_DIR / "infliximab_fdr_family_summary.csv"
    metadata_path = OUTPUT_DIR / "analysis_metadata.json"

    results.to_csv(combined_path, index=False, encoding="utf-8-sig")
    results.loc[results["SIGNIFICANT_FDR_0_05"]].to_csv(
        significant_path,
        index=False,
        encoding="utf-8-sig",
    )
    family_summary.to_csv(family_path, index=False, encoding="utf-8-sig")

    for genetic_model in GENETIC_MODELS:
        model_path = OUTPUT_DIR / f"infliximab_{genetic_model}_stratified_fdr.csv"
        results.loc[results["GENETIC_MODEL"] == genetic_model].to_csv(
            model_path,
            index=False,
            encoding="utf-8-sig",
        )

    with metadata_path.open("w", encoding="utf-8") as file:
        json.dump(metadata, file, ensure_ascii=False, indent=2)

    print("Analysis completed.")
    print(f"Combined results : {combined_path}")
    print(f"Significant only : {significant_path}")
    print(f"FDR family summary: {family_path}")
    print(f"Metadata         : {metadata_path}")
    print()

    display_columns = [
        "PHASE",
        "END_POINT",
        "GENETIC_MODEL",
        "N_VALID_TESTS",
        "N_FDR_SIGNIFICANT",
        "MIN_RAW_P",
        "MIN_FDR_P",
    ]
    print(family_summary[display_columns].to_string(index=False))


if __name__ == "__main__":
    analysis_results, fdr_family_summary, analysis_metadata = run_analysis()
    save_outputs(
        results=analysis_results,
        family_summary=fdr_family_summary,
        metadata=analysis_metadata,
    )
