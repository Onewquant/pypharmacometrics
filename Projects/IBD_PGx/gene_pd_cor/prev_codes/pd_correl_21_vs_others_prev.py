import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, fisher_exact
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

def fmt_mean_ci(x):
    x = pd.Series(x).dropna()
    n = len(x)

    if n == 0:
        return "NA"

    mean = x.mean()
    se = x.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan

    if n > 1:
        lcl = mean - 1.96 * se
        ucl = mean + 1.96 * se
        return f"{mean:.3f} ({lcl:.3f}-{ucl:.3f}), n={n}"
    else:
        return f"{mean:.3f}, n={n}"


def fmt_binary_count(x):
    x = pd.Series(x).dropna()
    n = len(x)

    if n == 0:
        return "NA"

    count = int((x == 1).sum())
    pct = count / n * 100

    return f"{count}/{n} ({pct:.1f}%)"


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


prj_dir = "/Projects/IBD_PGx"
resource_dir = f"{prj_dir}/gene_pd_cor"
output_dir = f"{prj_dir}/gene_pd_cor"

pid_df = pd.read_csv(f"{resource_dir}/pid_df.csv")

rsid_df = pd.read_csv(f"{resource_dir}/rsid_dosage_matrix_with_alleles.csv")
rsid_df = rsid_df[~rsid_df["UID"].isna()].copy()
rsid_df["UID"] = rsid_df["UID"].map(lambda x: str(x).split(".")[0])
rsid_df = rsid_df.rename(columns={"sid": "s"})

rsid_list = list(rsid_df.loc[:, "genomics_group":].columns)[1:]

ep_df = pd.read_csv(f"{resource_dir}/for_genomics_df(all_drugs).csv")
ep_df["UID"] = ep_df["UID"].astype(str)
ep_df["PHASE"] = ep_df["PHASE"].map(lambda x: x.split("_")[0])
ep_df["IBDTYPE(UC)"] = ep_df["IBDTYPE"].map({"CD": 0, "UC": 1})

sex_dict = (
    ep_df[~ep_df["SEX"].isna()]
    .drop_duplicates(["UID"])[["UID", "SEX"]]
    .copy()
)

sex_dict = sex_dict.set_index("UID")
sex_dict.index.name = None
sex_dict["SEX"] = sex_dict["SEX"].astype(int)
sex_dict = sex_dict["SEX"].to_dict()

ep_df["SEX"] = ep_df["UID"].map(sex_dict)


result_rows = []

for drug in ["infliximab", "adalimumab"]:
    for phase in ["IND", "MAINT", "ALL"]:
        for rsid in rsid_list:

            if phase == "ALL":
                phase_cond = ep_df["PHASE"].isin(["IND", "MAINT"])
            else:
                phase_cond = ep_df["PHASE"] == phase

            med_ep_df = ep_df[
                (ep_df["DRUG"] == drug) & phase_cond
            ].copy()

            # =====================================================
            # Carrier vs Non-carrier, Dominant model
            # carrier: dosage == 1 or 2
            # non: dosage == 0
            # =====================================================
            carrier_uids = list(rsid_df[rsid_df[rsid].isin([1, 2])]["UID"])
            non_uids = list(rsid_df[rsid_df[rsid] == 0]["UID"])

            group_df = pd.DataFrame({
                "UID": carrier_uids + non_uids,
                "GROUP": (
                    ["carrier"] * len(carrier_uids) +
                    ["non"] * len(non_uids)
                )
            })

            for ep_col in ["SEX", "IBDTYPE(UC)", "ADA", "CL"]:

                actv_ep_df = med_ep_df[["UID", ep_col]].copy()
                actv_ep_df = actv_ep_df.dropna(subset=["UID"])

                tmp_df = group_df.merge(
                    actv_ep_df,
                    on="UID",
                    how="left"
                )

                valid_values = tmp_df[ep_col].dropna().unique()

                total_n = tmp_df["UID"].nunique()
                available_n = tmp_df.dropna(subset=[ep_col])["UID"].nunique()

                data_availability = (
                    f"{available_n}/{total_n} ({available_n / total_n * 100:.1f}%)"
                    if total_n > 0 else "0/0 (NA)"
                )

                if len(valid_values) == 0:
                    continue

                # =====================================================
                # Binary endpoint
                # =====================================================
                if set(valid_values).issubset({0, 1, 0.0, 1.0}):
                    endpoint_type = "binary"

                    agg_df = (
                        tmp_df
                        .groupby(["UID", "GROUP"], as_index=False)[ep_col]
                        .max()
                    )

                    table = pd.crosstab(
                        agg_df["GROUP"],
                        agg_df[ep_col]
                    ).reindex(
                        index=["carrier", "non"],
                        columns=[0, 1],
                        fill_value=0
                    )

                    try:
                        overall_p = fisher_exact(table.values)[1]
                    except Exception:
                        overall_p = np.nan

                    carrier_est = fmt_binary_count(
                        agg_df.loc[agg_df["GROUP"] == "carrier", ep_col]
                    )

                    non_est = fmt_binary_count(
                        agg_df.loc[agg_df["GROUP"] == "non", ep_col]
                    )

                # =====================================================
                # Continuous endpoint
                # =====================================================
                elif len(valid_values) > 2:
                    endpoint_type = "continuous"

                    agg_df = (
                        tmp_df
                        .groupby(["UID", "GROUP"], as_index=False)[ep_col]
                        .mean()
                    )

                    carrier_x = agg_df.loc[
                        agg_df["GROUP"] == "carrier",
                        ep_col
                    ].dropna()

                    non_x = agg_df.loc[
                        agg_df["GROUP"] == "non",
                        ep_col
                    ].dropna()

                    try:
                        overall_p = mannwhitneyu(
                            carrier_x,
                            non_x,
                            alternative="two-sided"
                        ).pvalue
                    except Exception:
                        overall_p = np.nan

                    carrier_est = fmt_mean_ci(carrier_x)
                    non_est = fmt_mean_ci(non_x)

                else:
                    continue

                result_rows.append({
                    "DRUG": drug,
                    "PHASE": phase,
                    "RSID": rsid,
                    "GENE": rsid_gene_dict.get(rsid.split("(")[0], ""),
                    "END_POINT": ep_col,
                    "ENDPOINT_TYPE": endpoint_type,
                    "DATA_AVAILABILITY": data_availability,
                    "CARRIER_ESTIMATE": carrier_est,
                    "NON_CARRIER_ESTIMATE": non_est,
                    "P_VALUE": overall_p
                })


result_df = pd.DataFrame(result_rows)

result_df["P_VALUE_FDR"] = np.nan

group_cols = ["DRUG", "PHASE", "END_POINT"]

for keys, sub_df in result_df.groupby(group_cols):
    adj_p = fdr_adjust(sub_df["P_VALUE"].values)
    result_df.loc[sub_df.index, "P_VALUE_FDR"] = adj_p

result_df = result_df[
    [
        "DRUG",
        "PHASE",
        "RSID",
        "GENE",
        "END_POINT",
        "ENDPOINT_TYPE",
        "DATA_AVAILABILITY",
        "CARRIER_ESTIMATE",
        "NON_CARRIER_ESTIMATE",
        "P_VALUE",
        "P_VALUE_FDR"
    ]
]

result_df = result_df.sort_values("P_VALUE_FDR")

result_df.to_csv(
    f"{output_dir}/pgx_endpoint_carrier_vs_noncarrier_results.csv",
    index=False,
    encoding="utf-8-sig"
)