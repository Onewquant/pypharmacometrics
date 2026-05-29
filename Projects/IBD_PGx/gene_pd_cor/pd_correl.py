import pandas as pd
import numpy as np
import glob
import os
from scipy.stats import kruskal, mannwhitneyu, fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

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
        adj[mask] = multipletests(pvals[mask], method="fdr_bh")[1]
    return adj

def pairwise_posthoc(df, ep_col, endpoint_type):
    pairs = [
        ("hom", "het"),
        ("hom", "non"),
        ("het", "non")
    ]

    posthoc = []

    for g1, g2 in pairs:
        x1 = df.loc[df["GROUP"] == g1, ep_col].dropna()
        x2 = df.loc[df["GROUP"] == g2, ep_col].dropna()

        if len(x1) == 0 or len(x2) == 0:
            p = np.nan

        elif endpoint_type == "continuous":
            try:
                p = mannwhitneyu(x1, x2, alternative="two-sided").pvalue
            except:
                p = np.nan

        elif endpoint_type == "binary":
            table = pd.crosstab(
                df.loc[df["GROUP"].isin([g1, g2]), "GROUP"],
                df.loc[df["GROUP"].isin([g1, g2]), ep_col]
            )

            table = table.reindex(index=[g1, g2], columns=[0, 1], fill_value=0)

            try:
                p = fisher_exact(table.values)[1]
            except:
                p = np.nan

        posthoc.append({
            "comparison": f"{g1} vs {g2}",
            "p": p
        })

    posthoc_df = pd.DataFrame(posthoc)
    posthoc_df["p_adj_fdr"] = fdr_adjust(posthoc_df["p"].values)

    summary = "; ".join([
        f"{r['comparison']}: p={r['p']:.4g}, FDR={r['p_adj_fdr']:.4g}"
        if not pd.isna(r["p"]) else f"{r['comparison']}: NA"
        for _, r in posthoc_df.iterrows()
    ])

    return summary


prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx'
resource_dir = f'{prj_dir}/gene_pd_cor'
output_dir = f"{prj_dir}/gene_pd_cor"

pid_df = pd.read_csv(f"{resource_dir}/pid_df.csv")

rsid_df = pd.read_csv(f"{resource_dir}/rsid_dosage_matrix_with_alleles.csv")
rsid_df = rsid_df[~rsid_df['UID'].isna()].copy()
rsid_df['UID'] = rsid_df['UID'].map(lambda x:str(x).split('.')[0])
rsid_df = rsid_df.rename(columns={'sid':'s'})
rsid_list = list(rsid_df.loc[:,'genomics_group':].columns)[1:]
rsid_count = len(rsid_list)

ep_df = pd.read_csv(f"{resource_dir}/for_genomics_df(all_drugs).csv")
ep_df['UID'] = ep_df['UID'].astype(str)
ep_df['PHASE'] = ep_df['PHASE'].map(lambda x:x.split('_')[0])
ep_df['IBDTYPE(UC)'] = ep_df['IBDTYPE'].map({'CD':0,'UC':1})
sex_dict = ep_df[~ep_df['SEX'].isna()].drop_duplicates(['UID'])[['UID','SEX']].copy()
sex_dict = sex_dict.set_index(['UID'])
sex_dict.index.name = None
sex_dict['SEX'] = sex_dict['SEX'].astype(int)
sex_dict = sex_dict['SEX'].to_dict()
ep_df['SEX'] = ep_df['UID'].map(sex_dict)

# ep_df.drop_duplicates(['UID','SEX'])
# ep_df[['NAME','SEX']]
# ep_df.columns
result_rows = []

for drug in ['infliximab', 'adalimumab']:
    for phase in ['IND', 'MAINT', 'ALL']:
        for rsid in rsid_list:

            if phase == 'ALL':
                phase_cond = ep_df['PHASE'].isin(['IND', 'MAINT'])
            else:
                phase_cond = ep_df['PHASE'] == phase

            med_ep_df = ep_df[(ep_df['DRUG'] == drug) & phase_cond].copy()

            hom_uids = list(rsid_df[rsid_df[rsid] == 2]['UID'])
            het_uids = list(rsid_df[rsid_df[rsid] == 1]['UID'])
            non_uids = list(rsid_df[rsid_df[rsid] == 0]['UID'])

            group_df = pd.DataFrame({
                "UID": hom_uids + het_uids + non_uids,
                "GROUP": (
                    ["hom"] * len(hom_uids) +
                    ["het"] * len(het_uids) +
                    ["non"] * len(non_uids)
                )
            })

            for ep_col in ['SEX', 'IBDTYPE(UC)', 'ADA', 'CL']:

                actv_ep_df = med_ep_df[['UID', ep_col]].copy()
                actv_ep_df = actv_ep_df.dropna(subset=['UID'])

                tmp_df = group_df.merge(actv_ep_df, on="UID", how="left")

                valid_values = tmp_df[ep_col].dropna().unique()

                total_n = tmp_df['UID'].nunique()
                available_n = tmp_df.dropna(subset=[ep_col])['UID'].nunique()
                data_availability = (
                    f"{available_n}/{total_n} ({available_n / total_n * 100:.1f}%)"
                    if total_n > 0 else "0/0 (NA)"
                )

                if len(valid_values) == 0:
                    continue

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
                    ).reindex(index=["hom", "het", "non"], columns=[0, 1], fill_value=0)

                    try:
                        overall_p = chi2_contingency(table.values)[1]
                    except:
                        overall_p = np.nan

                    hom_est = fmt_binary_count(agg_df.loc[agg_df["GROUP"] == "hom", ep_col])
                    het_est = fmt_binary_count(agg_df.loc[agg_df["GROUP"] == "het", ep_col])
                    non_est = fmt_binary_count(agg_df.loc[agg_df["GROUP"] == "non", ep_col])

                elif len(valid_values) > 2:
                    endpoint_type = "continuous"

                    agg_df = (
                        tmp_df
                        .groupby(["UID", "GROUP"], as_index=False)[ep_col]
                        .mean()
                    )

                    hom_x = agg_df.loc[agg_df["GROUP"] == "hom", ep_col].dropna()
                    het_x = agg_df.loc[agg_df["GROUP"] == "het", ep_col].dropna()
                    non_x = agg_df.loc[agg_df["GROUP"] == "non", ep_col].dropna()

                    try:
                        overall_p = kruskal(hom_x, het_x, non_x).pvalue
                    except:
                        overall_p = np.nan

                    hom_est = fmt_mean_ci(hom_x)
                    het_est = fmt_mean_ci(het_x)
                    non_est = fmt_mean_ci(non_x)

                else:
                    continue

                posthoc_summary = pairwise_posthoc(
                    agg_df,
                    ep_col,
                    endpoint_type
                )

                result_rows.append({
                    "DRUG": drug,
                    "PHASE": phase,
                    "RSID": rsid,
                    "GENE": rsid_gene_dict.get(rsid.split('(')[0], ""),
                    "END_POINT": ep_col,
                    "ENDPOINT_TYPE": endpoint_type,
                    "DATA_AVAILABILITY": data_availability,
                    "HOM_ESTIMATE": hom_est,
                    "HET_ESTIMATE": het_est,
                    "NON_ESTIMATE": non_est,
                    "P_VALUE": overall_p,
                    "POSTHOC_SUMMARY": posthoc_summary
                })

result_df = pd.DataFrame(result_rows)

result_df["P_VALUE_FDR"] = np.nan

group_cols = ["DRUG", "PHASE", "END_POINT"]

for keys, sub_df in result_df.groupby(group_cols):
    adj_p = fdr_adjust(sub_df["P_VALUE"].values)
    result_df.loc[sub_df.index, "P_VALUE_FDR"] = adj_p

result_df = result_df[["DRUG", "PHASE", "RSID", "GENE", "END_POINT", "ENDPOINT_TYPE", "DATA_AVAILABILITY", "HOM_ESTIMATE", "HET_ESTIMATE", "NON_ESTIMATE", "P_VALUE", "P_VALUE_FDR", "POSTHOC_SUMMARY"]]
result_df = result_df.sort_values(['P_VALUE_FDR'])

result_df.to_csv(
    f"{output_dir}/pgx_endpoint_group_comparison_results.csv",
    index=False,
    encoding="utf-8-sig"
)


