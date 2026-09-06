"""Table 1 - Baseline characteristics (analytic cohort / infliximab PopPK
modeling cohort).

Standalone rewrite of paper_works/demotable_paper.py for the
infliximab popPK + PGx manuscript (adalimumab column dropped; the pooled
"analytic cohort" column is kept for the eligibility flow context).

2026-09-06 revision (to match the simplified Figure 1):
  - the infliximab column is built from the NONMEM *estimation* dataset
    (97 patients; the maintenance-starting patient with a single
    concentration, UID 38339532 / ID 97, is not in it), not from the
    98-patient "(for pda)" dataset
  - "Treatment phase" rows use the phase-level endpoint dataset
    (for_genomics_df(all_drugs).csv, PHASE == IND) - the same source that
    gives induction n = 83 in Figure 1 - instead of the TIME == 0
    observation heuristic, which misses maintenance-starting patients in
    the "(for pda)" datasets (their initial concentration row is absent).
    Analytic cohort: induction data for any drug; infliximab column:
    infliximab induction data.

Inputs:
  - results/modeling_df_covar/{drug}_integrated_datacheck(covar)(for pda).csv
  - C:/Users/ilma0/NONMEMProjects/IBDPGX/infliximab_integrated_modeling_df_dayscale.csv
  - paper_works_new/data/for_genomics_df(all_drugs).csv  (PHASE per UID/DRUG)

Output:
  - paper_works_new/output/Table1_demographics.csv
"""

from datetime import datetime

import numpy as np
import pandas as pd

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
results_dir = f"{prj_dir}/results"
nonmem_dir = "C:/Users/ilma0/NONMEMProjects/IBDPGX"
output_dir = f"{prj_dir}/paper_works_new/output"

TL_CUTOFF = 5

LAB_COLS = {
    "ALB": "Albumin (g/dL)",
    "AST": "AST (IU/L)",
    "ALT": "ALT (IU/L)",
    "CRP": "CRP (mg/dL)",
    "FCAL": "Fecal calprotectin (mg/kg)",
    "CREATININE": "Serum creatinine (mg/dL)",
}


def n_pct(n, denom):
    return f"{n} ({round(100 * n / denom, 2)})" if denom else f"{n} (NA)"


def mean_sd(x, digits=2):
    x = pd.Series(x).astype(float)
    return f"{round(np.mean(x), digits)} ({round(np.std(x), digits)})"


def load_analytic_cohort():
    """Pooled infliximab + adalimumab datacheck files -> analytic cohort."""
    dfs = []
    for drug in ["infliximab", "adalimumab"]:
        df = pd.read_csv(
            f"{results_dir}/modeling_df_covar/"
            f"{drug}_integrated_datacheck(covar)(for pda).csv"
        )
        dfs.append(df)

    md_df = pd.concat(dfs)[list(dfs[0].columns)].sort_values(["ID", "DATETIME"])
    md_df["IBD_TYPE"] = md_df["IBD_TYPE"].map({"CD": 0, "UC": 1})
    md_df["AGE"] = md_df.apply(
        lambda x: int(
            (
                datetime.strptime(x["DATETIME"], "%Y-%m-%d")
                - datetime.strptime(x["AGE"], "%Y-%m-%d")
            ).days
            / 365.25
        ),
        axis=1,
    )
    md_df["SEX"] = md_df["SEX"].map({"남": 0, "여": 1})
    md_df["ROUTE"] = md_df["ROUTE"].map({"IV": 1, "SC": 2, ".": "."})
    md_df["UID"] = md_df["UID"].astype(str)
    return md_df, "UID"


def load_infliximab_cohort():
    """NONMEM estimation dataset (97 patients) + UID via the datacheck map."""
    md_df = pd.read_csv(
        f"{nonmem_dir}/infliximab_integrated_modeling_df_dayscale.csv"
    )
    id_map = (
        pd.read_csv(
            f"{results_dir}/modeling_df_covar/"
            "infliximab_integrated_datacheck(covar)(for pda).csv",
            usecols=["ID", "UID"],
        )
        .drop_duplicates("ID")
        .set_index("ID")["UID"]
        .astype(str)
    )
    md_df["UID"] = md_df["ID"].map(id_map)
    assert md_df["UID"].notna().all(), "unmapped NONMEM ID in estimation dataset"
    return md_df, "ID"


def load_induction_uids():
    """UIDs with induction-phase data (Figure 1 / attrition definition)."""
    ep = pd.read_csv(f"{prj_dir}/paper_works_new/data/for_genomics_df(all_drugs).csv")
    ep["UID"] = ep["UID"].astype(str)
    ind = ep[ep["PHASE"].astype(str).str.split("_").str[0] == "IND"]
    return {
        "any": set(ind["UID"]),
        "infliximab": set(ind[ind["DRUG"] == "infliximab"]["UID"]),
    }


def summarize(md_df, id_col, totals, ind_uids):
    md_df = md_df.copy()
    md_df["Pediatric"] = (md_df["AGE"] < 19).astype(int)
    md_df["BMI"] = md_df["WT"] / ((md_df["HT"] / 100) ** 2)

    first_df = md_df.drop_duplicates(subset=[id_col])
    subtotal_n = md_df[id_col].nunique()

    row = {}
    row["Demographic categorical variables, n (%)"] = ""

    if totals["patients"] is None:
        totals["patients"] = subtotal_n
    row["Subtotal patients"] = n_pct(subtotal_n, totals["patients"])
    row["Female"] = n_pct(int((first_df["SEX"] == 1).sum()), subtotal_n)
    row["Pediatric (Age < 19)"] = n_pct(
        int((first_df["Pediatric"] == 1).sum()), subtotal_n
    )

    row["Demographic continuous variables, Mean (SD)"] = ""
    age_series = md_df[md_df["MDV"] == 1].drop_duplicates([id_col])["AGE"]
    row["Age at the 1st Dose"] = mean_sd(age_series)
    row["Height"] = mean_sd(first_df["HT"])
    row["Weight"] = mean_sd(first_df["WT"])
    row["BMI"] = mean_sd(first_df["BMI"])

    row["Laboratory test, mean (SD)"] = ""
    for col, label in LAB_COLS.items():
        row[label] = mean_sd(first_df[col].dropna())

    row["Diagnosis, n (%)"] = ""
    row["CD"] = n_pct(int((first_df["IBD_TYPE"] == 0).sum()), subtotal_n)
    row["UC"] = n_pct(int((first_df["IBD_TYPE"] == 1).sum()), subtotal_n)

    row["Treatment phase, n (%)"] = ""
    # induction-phase data available (same definition as Figure 1 / S5)
    induction_n = md_df[md_df["UID"].isin(ind_uids)][id_col].nunique()
    row["Whole phases"] = n_pct(induction_n, subtotal_n)
    row["Maintenance only"] = n_pct(subtotal_n - induction_n, subtotal_n)

    row["Blood Sampling, n (%)"] = ""
    sampling_df = md_df[
        (md_df["MDV"] != 1)
        & (~((md_df["DV"].astype(str) == "0.0") & (md_df["TIME"] == 0)))
    ].copy()
    sampling_desc = sampling_df.groupby(id_col).agg(DV_COUNT=("DV", "count"))

    if totals["samples"] is None:
        totals["samples"] = len(sampling_df)
        totals["patients_with_samples"] = len(sampling_desc)

    dv_vals = sampling_df["DV"].astype(float)
    row["Total samples"] = n_pct(len(sampling_df), totals["samples"])
    row["Patients with samples"] = n_pct(
        len(sampling_desc), totals["patients_with_samples"]
    )
    row[f"TL < {TL_CUTOFF}"] = n_pct(int((dv_vals < TL_CUTOFF).sum()), len(dv_vals))
    row[f"TL >= {TL_CUTOFF}"] = n_pct(int((dv_vals >= TL_CUTOFF).sum()), len(dv_vals))

    ada_pos_n = md_df[md_df["ADA"] != 0].drop_duplicates([id_col])[id_col].nunique()
    row["Anti-drug antibody"] = n_pct(ada_pos_n, subtotal_n)
    row["Samples/person, mean (SD)"] = mean_sd(sampling_desc["DV_COUNT"])

    row["Dosing route, n (%)"] = ""
    dosing_df = md_df[md_df["MDV"] == 1]
    route_counts = dosing_df["CMT"].value_counts()
    total_routes = int(route_counts.sum())
    row["Total routes"] = n_pct(total_routes, total_routes)
    row["Subcutaneous"] = n_pct(int(route_counts.get(1, 0)), total_routes)
    row["Intravenous"] = n_pct(int(route_counts.get(2, 0)), total_routes)

    return row


totals = {"patients": None, "samples": None, "patients_with_samples": None}

analytic_df, analytic_id = load_analytic_cohort()
ifx_df, ifx_id = load_infliximab_cohort()
ind_uids = load_induction_uids()

table = pd.DataFrame({
    "Characteristics": list(
        summarize(analytic_df, analytic_id, dict(totals), ind_uids["any"]).keys()
    ),
})

totals = {"patients": None, "samples": None, "patients_with_samples": None}
analytic_row = summarize(analytic_df, analytic_id, totals, ind_uids["any"])
ifx_row = summarize(ifx_df, ifx_id, totals, ind_uids["infliximab"])

table["Analytic cohort"] = table["Characteristics"].map(analytic_row)
table["Infliximab PopPK modeling cohort"] = table["Characteristics"].map(ifx_row)

table.to_csv(
    f"{output_dir}/Table1_demographics.csv", index=False, encoding="utf-8-sig"
)

print(table.to_string(index=False))
