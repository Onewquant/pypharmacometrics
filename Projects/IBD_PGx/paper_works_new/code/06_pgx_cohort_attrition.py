"""PGx cohort attrition numbers for the eligibility flow chart (Figure 1).

Reproduces every patient count in the PGx branch of the flow chart and
flags the special-case UIDs that need manual confirmation:
  - infliximab cohort (98)
  - genotype data unavailable after variant QC / sample-ID matching (1)
  - PGx cohort (97)
  - per-phase analysis populations (phase data -> phase-specific CL ->
    covariate-complete), for IND / MAINT / ALL

Wording note: the manuscript defines the analytic cohort as WGS-performed
patients, so the PGx-branch exclusion must NOT be labeled "without WGS".
Confirmed reasons (2026-08):
  - UID 17439372: sequencing sample could not be matched to the patient
    ID (ID-matching failure) -> excluded from the genotype matrix
  - UID 35093356: included in popPK modeling (NONMEM ID 78) but the
    phase-level derived dataset row is entirely empty - phase-specific
    estimates could not be derived (maintenance-only patient, phase
    anchor dates not determinable)

Output:
  - paper_works_new/output/Table_pgx_attrition.csv
"""

import pandas as pd

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
data_dir = f"{prj_dir}/paper_works_new/data"
output_dir = f"{prj_dir}/paper_works_new/output"
datacheck_path = (
    f"{prj_dir}/results/modeling_df_covar/"
    "infliximab_integrated_datacheck(covar)(for pda).csv"
)

rsid_df = pd.read_csv(f"{data_dir}/rsid_dosage_matrix_with_alleles.csv")
n_geno_samples_total = len(rsid_df)
n_geno_unmatched = int(rsid_df["UID"].isna().sum())
rsid_df = rsid_df[~rsid_df["UID"].isna()].copy()
rsid_df["UID"] = rsid_df["UID"].map(lambda x: str(x).split(".")[0])
geno_uids = set(rsid_df["UID"])

ep_df = pd.read_csv(f"{data_dir}/for_genomics_df(all_drugs).csv")
ep_df["UID"] = ep_df["UID"].astype(str)
ep_df["PHASE"] = ep_df["PHASE"].map(lambda x: x.split("_")[0])
ifx_df = ep_df[ep_df["DRUG"] == "infliximab"].copy()
ifx_df["WEIGHT"] = ifx_df["WT"]
ifx_df["ALBUMIN"] = ifx_df["ALB"]

ifx_uids = set(ifx_df["UID"])
no_geno_uids = sorted(ifx_uids - geno_uids)
pgx_uids = ifx_uids & geno_uids

# patients without any observed concentration sample (Table 1 logic)
dc = pd.read_csv(datacheck_path)
dc["UID"] = dc["UID"].astype(str)
samp = dc[
    (dc["MDV"] != 1)
    & (~((dc["DV"].astype(str) == "0.0") & (dc["TIME"] == 0)))
]
no_sample_uids = sorted(set(dc["UID"]) - set(samp["UID"]))

rows = [
    {"STEP": "Infliximab cohort (popPK)", "N": len(ifx_uids), "NOTE": ""},
    {
        "STEP": "Excluded: sequencing sample not matched to patient ID",
        "N": len(no_geno_uids),
        "NOTE": f"UID {', '.join(no_geno_uids)} - ID-matching failure "
                f"(confirmed 2026-08; genotype matrix: "
                f"{n_geno_samples_total} samples, "
                f"{n_geno_unmatched} unmatched to a patient UID)",
    },
    {"STEP": "PGx analysis cohort", "N": len(pgx_uids), "NOTE": ""},
]

covar_cols = ["SEX", "WEIGHT", "ALBUMIN", "ADA"]

for phase in ["IND", "MAINT", "ALL"]:
    cond = (
        ifx_df["PHASE"].isin(["IND", "MAINT"])
        if phase == "ALL" else ifx_df["PHASE"] == phase
    )
    d = (
        ifx_df[cond]
        .groupby("UID", as_index=False)
        .agg(
            CL=("CL", "mean"),
            ADA=("ADA", "max"),
            SEX=("SEX", "first"),
            WEIGHT=("WEIGHT", "mean"),
            ALBUMIN=("ALBUMIN", "mean"),
        )
    )
    d = d[d["UID"].isin(pgx_uids)]

    n_phase = d["UID"].nunique()
    d_cl = d.dropna(subset=["CL"])
    d_cl = d_cl[d_cl["CL"] > 0]
    n_cl = d_cl["UID"].nunique()
    n_full = d_cl.dropna(subset=covar_cols)["UID"].nunique()
    no_cl_uids = sorted(set(d["UID"]) - set(d_cl["UID"]))

    rows.append({
        "STEP": f"[{phase}] phase data available",
        "N": n_phase,
        "NOTE": "",
    })
    rows.append({
        "STEP": f"[{phase}] phase-specific CL available (analysis N, CL)",
        "N": n_cl,
        "NOTE": (
            f"without phase-specific CL: UID {', '.join(no_cl_uids)}"
            if no_cl_uids else ""
        ),
    })
    rows.append({
        "STEP": f"[{phase}] covariate-complete (SEX/WT/ALB/ADA)",
        "N": n_full,
        "NOTE": "",
    })

rows.append({
    "STEP": "Reference: patients without any observed concentration sample",
    "N": len(no_sample_uids),
    "NOTE": f"UID {', '.join(no_sample_uids)} - CL EBEs are "
            f"shrinkage-driven; handled as a sensitivity analysis (04)",
})

result_df = pd.DataFrame(rows)
result_df.to_csv(
    f"{output_dir}/Table_pgx_attrition.csv", index=False, encoding="utf-8-sig"
)

print(result_df.to_string(index=False))
