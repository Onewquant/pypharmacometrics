"""Sensitivity analyses for FDR-significant PGx hits.

For every CL row of the main analysis (03) that is either FDR-significant
or belongs to the pre-designated lead finding (see LEAD_TARGETS), this
script re-evaluates the genotype effect with:
  - leave-one-out over the variant-group subjects (min/max p across
    re-fits; checks that no single subject drives the signal)
  - Mann-Whitney U test on the raw endpoint (covariate-free,
    rank-based check)
  - HC3 heteroskedasticity-robust standard errors (same ANCOVA)
  - exclusion of patients without any observed concentration sample
    (their CL EBEs are shrunk toward covariate-predicted typical values,
    so the association is re-tested on observation-informed CL only)

Run 03_pgx_ancova_fdr.py first.

Output:
  - paper_works_new/output/Table_pgx_sensitivity.csv
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as sps

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
data_dir = f"{prj_dir}/paper_works_new/data"
output_dir = f"{prj_dir}/paper_works_new/output"
datacheck_path = (
    f"{prj_dir}/results/modeling_df_covar/"
    "infliximab_integrated_datacheck(covar)(for pda).csv"
)


def get_no_sample_uids():
    """Patients without any observed concentration sample
    (same counting logic as Table 1 'Patients with samples')."""
    dc = pd.read_csv(datacheck_path)
    dc["UID"] = dc["UID"].astype(str)
    samp = dc[
        (dc["MDV"] != 1)
        & (~((dc["DV"].astype(str) == "0.0") & (dc["TIME"] == 0)))
    ]
    return set(dc["UID"]) - set(samp["UID"])

COMPARISONS = {
    "HOM_vs_OTHERS": ([2], [0, 1]),
    "CARRIER_vs_NONCARRIER": ([1, 2], [0]),
}

# Rows always carried through the robustness checks, regardless of whether
# they reach the FDR threshold. Under the exploratory reporting frame
# (option A, 2026-08) the lead association is reported with its sensitivity
# analyses even though q > 0.05.
LEAD_TARGETS = [("rs1061622", "HOM_vs_OTHERS"), ("rs396991", "HOM_vs_OTHERS")]


def build_phase_df(ep_df, phase):
    if phase == "OVERALL":
        cond = ep_df["PHASE"].isin(["IND", "MAINT"])
    else:
        cond = ep_df["PHASE"] == phase

    uid_df = (
        ep_df[cond]
        .groupby("UID", as_index=False)
        .agg(
            CL=("CL", "mean"),
            ADA=("ADA", "max"),
            SEX=("SEX", "first"),
            WEIGHT=("WEIGHT", "mean"),
            ALBUMIN=("ALBUMIN", "mean"),
        )
    )
    uid_df = uid_df[uid_df["CL"].isna() | (uid_df["CL"] > 0)].copy()
    uid_df["LOG_CL"] = np.log(uid_df["CL"])
    return uid_df


def fit_group_p(df, y_col, covariates, cov_type=None):
    x_cols = ["GROUP"] + covariates
    model_df = df[[y_col] + x_cols].dropna()

    if model_df["GROUP"].nunique() < 2 or len(model_df) < len(x_cols) + 2:
        return np.nan

    X = sm.add_constant(model_df[x_cols], has_constant="add")

    try:
        if cov_type:
            fit = sm.OLS(model_df[y_col], X).fit(cov_type=cov_type)
        else:
            fit = sm.OLS(model_df[y_col], X).fit()
        return fit.pvalues.get("GROUP", np.nan)
    except Exception:
        return np.nan


rsid_df = pd.read_csv(f"{data_dir}/rsid_dosage_matrix_with_alleles.csv")
rsid_df = rsid_df[~rsid_df["UID"].isna()].copy()
rsid_df["UID"] = rsid_df["UID"].map(lambda x: str(x).split(".")[0])

ep_df = pd.read_csv(f"{data_dir}/for_genomics_df(all_drugs).csv")
ep_df["UID"] = ep_df["UID"].astype(str)
ep_df["PHASE"] = ep_df["PHASE"].map(lambda x: x.split("_")[0])
ep_df = ep_df[ep_df["DRUG"] == "infliximab"].copy()
ep_df["WEIGHT"] = ep_df["WT"]
ep_df["ALBUMIN"] = ep_df["ALB"]

main_df = pd.read_csv(f"{output_dir}/Table_pgx_ancova_fdr_results.csv")
main_df["RS"] = main_df["RSID"].str.split("(").str[0]

is_cl = main_df["END_POINT"] == "CL"
is_sig = main_df["P_VALUE_FDR"] < 0.05
is_lead = pd.Series(False, index=main_df.index)
for rs, comp in LEAD_TARGETS:
    is_lead |= (main_df["RS"] == rs) & (main_df["COMPARISON"] == comp)

sig_df = main_df[is_cl & (is_sig | is_lead)].copy()
sig_df["TRIGGER"] = np.where(
    sig_df["P_VALUE_FDR"] < 0.05, "FDR-significant", "pre-designated lead"
)

if len(sig_df) == 0:
    print("No rows selected for sensitivity analysis.")

no_sample_uids = get_no_sample_uids()

rows = []

for _, hit in sig_df.iterrows():
    phase = hit["PHASE"]
    comparison = hit["COMPARISON"]
    rsid = hit["RSID"]
    y_col = "LOG_CL" if hit["MODEL_SCALE"] == "log(CL)" else "CL"

    variant_dos, reference_dos = COMPARISONS[comparison]
    covariates = ["SEX", "WEIGHT", "ALBUMIN", "ADA"]

    uid_df = build_phase_df(ep_df, phase)

    geno_df = rsid_df[["UID", rsid]].rename(columns={rsid: "DOS"})
    geno_df["GROUP"] = np.where(
        geno_df["DOS"].isin(variant_dos),
        1,
        np.where(geno_df["DOS"].isin(reference_dos), 0, np.nan),
    )

    d = uid_df.merge(geno_df[["UID", "GROUP"]], on="UID", how="inner")
    d = d.dropna(subset=["GROUP", y_col, "SEX", "WEIGHT", "ALBUMIN", "ADA"])

    variant_cl = d.loc[d["GROUP"] == 1, "CL"]
    reference_cl = d.loc[d["GROUP"] == 0, "CL"]

    # leave-one-out over variant-group subjects
    loo_ps = []
    for uid in d.loc[d["GROUP"] == 1, "UID"]:
        loo_ps.append(
            fit_group_p(d[d["UID"] != uid], y_col, covariates)
        )
    loo_ps = [p for p in loo_ps if pd.notna(p)]

    mw_p = (
        sps.mannwhitneyu(variant_cl, reference_cl).pvalue
        if len(variant_cl) > 0 and len(reference_cl) > 0
        else np.nan
    )

    # exclude patients whose CL EBE is not informed by any observed sample
    d_obs = d[~d["UID"].isin(no_sample_uids)]
    excl_n = len(d) - len(d_obs)
    excl_variant_n = int(
        (d.loc[d["UID"].isin(no_sample_uids), "GROUP"] == 1).sum()
    )
    excl_p = fit_group_p(d_obs, y_col, covariates)

    rows.append({
        "PHASE": phase,
        "COMPARISON": comparison,
        "RSID": rsid,
        "GENE": hit["GENE"],
        "MODEL_SCALE": hit["MODEL_SCALE"],
        "TRIGGER": hit["TRIGGER"],
        "MAIN_P": hit["P_VALUE"],
        "MAIN_P_FDR": hit["P_VALUE_FDR"],
        "LOO_P_MIN": round(min(loo_ps), 5) if loo_ps else np.nan,
        "LOO_P_MAX": round(max(loo_ps), 5) if loo_ps else np.nan,
        "LOO_ALL_BELOW_0.05": (
            "Y" if loo_ps and max(loo_ps) < 0.05 else "N"
        ),
        "MANN_WHITNEY_P": round(mw_p, 5) if pd.notna(mw_p) else np.nan,
        "HC3_ROBUST_P": round(
            fit_group_p(d, y_col, covariates, cov_type="HC3"), 5
        ),
        "EXCL_NO_SAMPLE_P": round(excl_p, 5) if pd.notna(excl_p) else np.nan,
        "EXCL_NO_SAMPLE_N": excl_n,
        "EXCL_NO_SAMPLE_VARIANT_N": excl_variant_n,
        "VARIANT_N": int((d["GROUP"] == 1).sum()),
        "REFERENCE_N": int((d["GROUP"] == 0).sum()),
    })

result_df = pd.DataFrame(rows)
result_df.to_csv(
    f"{output_dir}/Table_pgx_sensitivity.csv", index=False, encoding="utf-8-sig"
)

if len(result_df) > 0:
    print(result_df.to_string(index=False))
