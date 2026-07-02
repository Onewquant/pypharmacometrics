import pandas as pd
import numpy as np

from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests


rsid_gene_dict = {
    "rs9828223": "CD96",
    "rs2097432": "HLA-DQA1",
    "rs396991": "FCGR3A",
    "rs1800629": "TNFα (TNF)",
    "rs4149570": "TNFRSF1A",
    "rs3397": "TNFRSF1B",
    "rs1061624": "TNFRSF1B",
    "rs5030728": "TLR4",
    "rs3804099": "TLR2",
    "rs10499563": "IL6",
    "rs2275913": "IL17A",
    "rs1800872": "IL10",
    "rs3024505": "IL10",
    "rs361525": "TNFα (TNF)",
    "rs767455": "TNFα (TNF)",
    "rs1061622": "TNFα (TNF)",
    "rs765249238": "SLCO2A1",
    "rs776813259": "SLCO2A1",
}


def fdr_adjust(pvals):
    pvals = np.array(pvals, dtype=float)
    mask = ~np.isnan(pvals)

    adj = np.full(len(pvals), np.nan)

    if mask.sum() > 0:
        adj[mask] = multipletests(
            pvals[mask],
            method="fdr_bh"
        )[1]

    return adj


def fmt_median_iqr(x):
    x = pd.Series(x).dropna()
    n = len(x)

    if n == 0:
        return "NA"

    median = x.median()
    q1 = x.quantile(0.25)
    q3 = x.quantile(0.75)

    return f"{median:.3f} ({q1:.3f}-{q3:.3f}), n={n}"


def fmt_binary_count(x):
    x = pd.Series(x).dropna()
    n = len(x)

    if n == 0:
        return "NA"

    count = int((x == 1).sum())
    pct = count / n * 100

    return f"{count}/{n} ({pct:.1f}%)"


def run_kruskal_three_groups(df, y_col, group_col):
    g0 = df.loc[df[group_col] == 0, y_col].dropna()
    g1 = df.loc[df[group_col] == 1, y_col].dropna()
    g2 = df.loc[df[group_col] == 2, y_col].dropna()

    if len(g0) == 0 or len(g1) == 0 or len(g2) == 0:
        return np.nan

    try:
        return kruskal(g0, g1, g2).pvalue
    except Exception:
        return np.nan


prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
resource_dir = f"{prj_dir}/gene_pd_cor"
output_dir = f"{prj_dir}/gene_pd_cor"

rsid_df = pd.read_csv(f"{resource_dir}/rsid_dosage_matrix_with_alleles.csv")
rsid_df = rsid_df[~rsid_df["UID"].isna()].copy()
rsid_df["UID"] = rsid_df["UID"].map(lambda x: str(x).split(".")[0])

rsid_list = list(rsid_df.loc[:, "genomics_group":].columns)[1:]

ep_df = pd.read_csv(f"{resource_dir}/for_genomics_df(all_drugs).csv")
ep_df["UID"] = ep_df["UID"].astype(str)
ep_df["PHASE"] = ep_df["PHASE"].map(lambda x: x.split("_")[0])

if "WT" in ep_df.columns and "WEIGHT" not in ep_df.columns:
    ep_df["WEIGHT"] = ep_df["WT"]

if "ALB" in ep_df.columns and "ALBUMIN" not in ep_df.columns:
    ep_df["ALBUMIN"] = ep_df["ALB"]

result_rows = []

group_col = "GENOTYPE_DOSAGE"

for drug in ["infliximab", "adalimumab"]:
    for phase in ["IND", "MAINT", "ALL"]:

        if phase == "ALL":
            phase_cond = ep_df["PHASE"].isin(["IND", "MAINT"])
        else:
            phase_cond = ep_df["PHASE"] == phase

        med_ep_df = ep_df[
            (ep_df["DRUG"] == drug) & phase_cond
        ].copy()

        uid_ep_df = (
            med_ep_df
            .groupby("UID", as_index=False)
            .agg(
                CL=("CL", "mean"),
                ADA=("ADA", "max"),
                SEX=("SEX", "first"),
                WEIGHT=("WEIGHT", "mean"),
                ALBUMIN=("ALBUMIN", "mean"),
            )
        )

        for rsid in rsid_list:

            geno_df = rsid_df[["UID", rsid]].copy()
            geno_df = geno_df.rename(columns={rsid: group_col})
            geno_df[group_col] = pd.to_numeric(
                geno_df[group_col],
                errors="coerce"
            )

            analysis_df = uid_ep_df.merge(
                geno_df[["UID", group_col]],
                on="UID",
                how="inner"
            )

            total_n = analysis_df["UID"].nunique()

            if drug == "adalimumab":
                endpoint_list = ["CL"]
            else:
                endpoint_list = ["ADA", "CL"]

            for ep_col in endpoint_list:

                tmp_df = analysis_df.dropna(
                    subset=["UID", group_col, ep_col]
                ).copy()

                tmp_df = tmp_df[tmp_df[group_col].isin([0, 1, 2])].copy()

                group_counts = tmp_df[group_col].value_counts()
                n0 = group_counts.get(0, 0)
                n1 = group_counts.get(1, 0)
                n2 = group_counts.get(2, 0)

                # 세 군 모두 n >= 8인 경우만 분석
                if n0 < 8 or n1 < 8 or n2 < 8:
                    continue

                available_n = tmp_df["UID"].nunique()

                data_availability = (
                    f"{available_n}/{total_n} ({available_n / total_n * 100:.1f}%)"
                    if total_n > 0 else "0/0 (NA)"
                )

                valid_values = tmp_df[ep_col].dropna().unique()

                if len(valid_values) == 0:
                    continue

                if set(valid_values).issubset({0, 1, 0.0, 1.0}):
                    endpoint_type = "binary"

                    group0_est = fmt_binary_count(
                        tmp_df.loc[tmp_df[group_col] == 0, ep_col]
                    )
                    group1_est = fmt_binary_count(
                        tmp_df.loc[tmp_df[group_col] == 1, ep_col]
                    )
                    group2_est = fmt_binary_count(
                        tmp_df.loc[tmp_df[group_col] == 2, ep_col]
                    )

                    # binary endpoint에 Kruskal-Wallis를 적용할 수도 있지만,
                    # 해석은 연속형 endpoint보다 약합니다.
                    p_value = run_kruskal_three_groups(
                        tmp_df,
                        y_col=ep_col,
                        group_col=group_col
                    )

                else:
                    endpoint_type = "continuous"

                    group0_est = fmt_median_iqr(
                        tmp_df.loc[tmp_df[group_col] == 0, ep_col]
                    )
                    group1_est = fmt_median_iqr(
                        tmp_df.loc[tmp_df[group_col] == 1, ep_col]
                    )
                    group2_est = fmt_median_iqr(
                        tmp_df.loc[tmp_df[group_col] == 2, ep_col]
                    )

                    p_value = run_kruskal_three_groups(
                        tmp_df,
                        y_col=ep_col,
                        group_col=group_col
                    )

                result_rows.append({
                    "DRUG": drug,
                    "PHASE": phase,
                    "RSID": rsid,
                    "GENE": rsid_gene_dict.get(rsid.split("(")[0], ""),
                    "END_POINT": ep_col,
                    "ENDPOINT_TYPE": endpoint_type,
                    "DATA_AVAILABILITY": data_availability,
                    "GENOTYPE_0_ESTIMATE": group0_est,
                    "GENOTYPE_1_ESTIMATE": group1_est,
                    "GENOTYPE_2_ESTIMATE": group2_est,
                    "P_VALUE": p_value,
                    "P_VALUE_FDR": np.nan,
                    "GENOTYPE_0_N": n0,
                    "GENOTYPE_1_N": n1,
                    "GENOTYPE_2_N": n2,
                })


result_df = pd.DataFrame(result_rows)

if len(result_df) > 0:
    for keys, sub_df in result_df.groupby(["DRUG", "PHASE", "END_POINT"]):
        result_df.loc[sub_df.index, "P_VALUE_FDR"] = fdr_adjust(
            sub_df["P_VALUE"].values
        )

    result_df = result_df[
        [
            "DRUG",
            "PHASE",
            "RSID",
            "GENE",
            "END_POINT",
            "ENDPOINT_TYPE",
            "DATA_AVAILABILITY",
            "GENOTYPE_0_ESTIMATE",
            "GENOTYPE_1_ESTIMATE",
            "GENOTYPE_2_ESTIMATE",
            "P_VALUE",
            "P_VALUE_FDR",
            "GENOTYPE_0_N",
            "GENOTYPE_1_N",
            "GENOTYPE_2_N",
        ]
    ]

    result_df = result_df.sort_values(
        ["P_VALUE_FDR", "P_VALUE"],
        na_position="last"
    )

result_df.to_csv(
    f"{output_dir}/pgx_kruskal_wallis_three_genotype_groups_min8_results.csv",
    index=False,
    encoding="utf-8-sig"
)