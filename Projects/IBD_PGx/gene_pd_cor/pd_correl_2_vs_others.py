import pandas as pd
import numpy as np

import statsmodels.api as sm
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
        adj[mask] = multipletests(pvals[mask], method="fdr_bh")[1]

    return adj


def fmt_mean_ci(x):
    x = pd.Series(x).dropna()
    n = len(x)

    if n == 0:
        return "NA"

    mean = x.mean()
    se = x.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan

    if n > 1:
        return f"{mean:.3f} ({mean - 1.96 * se:.3f}-{mean + 1.96 * se:.3f}), n={n}"
    return f"{mean:.3f}, n={n}"


def fmt_binary_count(x):
    x = pd.Series(x).dropna()
    n = len(x)

    if n == 0:
        return "NA"

    count = int((x == 1).sum())
    return f"{count}/{n} ({count / n * 100:.1f}%)"


def get_model_result(df, y_col, x_cols, endpoint_type, group_col):
    model_df = df[[y_col] + x_cols].dropna().copy()

    if len(model_df) < len(x_cols) + 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df)

    if model_df[y_col].nunique() < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df)

    if model_df[group_col].nunique() < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df)

    X = sm.add_constant(model_df[x_cols], has_constant="add")
    y = model_df[y_col]

    try:
        if endpoint_type == "binary":
            fit = sm.Logit(y, X).fit(disp=False)

            beta = fit.params.get(group_col, np.nan)
            se = fit.bse.get(group_col, np.nan)
            pval = fit.pvalues.get(group_col, np.nan)

            effect = np.exp(beta)
            ci_lower = np.exp(beta - 1.96 * se)
            ci_upper = np.exp(beta + 1.96 * se)

        else:
            fit = sm.OLS(y, X).fit()

            effect = fit.params.get(group_col, np.nan)
            se = fit.bse.get(group_col, np.nan)
            pval = fit.pvalues.get(group_col, np.nan)

            ci_lower = effect - 1.96 * se
            ci_upper = effect + 1.96 * se

        return effect, se, ci_lower, ci_upper, pval, len(model_df)

    except Exception:
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df)


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

group_col = "HOM_GROUP"

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
            geno_df = geno_df.rename(columns={rsid: "GENOTYPE_DOSAGE"})

            # homozygote variant: dosage == 2
            # others: dosage == 0 or 1
            geno_df[group_col] = np.where(
                geno_df["GENOTYPE_DOSAGE"] == 2,
                1,
                np.where(geno_df["GENOTYPE_DOSAGE"].isin([0, 1]), 0, np.nan)
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

                # 각 그룹 n >= 8 필터
                group_counts = tmp_df[group_col].value_counts()
                hom_n = group_counts.get(1, 0)
                other_n = group_counts.get(0, 0)

                if hom_n < 8 or other_n < 8:
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
                    effect_name = "OR"

                    hom_est = fmt_binary_count(
                        tmp_df.loc[tmp_df[group_col] == 1, ep_col]
                    )
                    other_est = fmt_binary_count(
                        tmp_df.loc[tmp_df[group_col] == 0, ep_col]
                    )

                else:
                    endpoint_type = "continuous"
                    effect_name = "BETA"

                    hom_est = fmt_mean_ci(
                        tmp_df.loc[tmp_df[group_col] == 1, ep_col]
                    )
                    other_est = fmt_mean_ci(
                        tmp_df.loc[tmp_df[group_col] == 0, ep_col]
                    )

                if ep_col == "ADA":
                    covariates = ["SEX", "WEIGHT", "ALBUMIN"]
                elif drug == "adalimumab":
                    covariates = ["SEX", "WEIGHT", "ALBUMIN"]
                else:
                    covariates = ["SEX", "WEIGHT", "ALBUMIN", "ADA"]

                x_cols = [group_col] + covariates

                effect, se, ci_lower, ci_upper, p_value, model_n = get_model_result(
                    tmp_df,
                    y_col=ep_col,
                    x_cols=x_cols,
                    endpoint_type=endpoint_type,
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
                    "HOM_ESTIMATE": hom_est,
                    "OTHER_ESTIMATE": other_est,
                    "EFFECT_NAME": effect_name,
                    "EFFECT": effect,
                    "SE": se,
                    "CI_LOWER": ci_lower,
                    "CI_UPPER": ci_upper,
                    "P_VALUE": p_value,
                    "P_VALUE_FDR": np.nan,
                    "COVARIATES": ", ".join(covariates),
                    "MODEL_N": model_n,
                    "HOM_N": hom_n,
                    "OTHER_N": other_n,
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
            "HOM_ESTIMATE",
            "OTHER_ESTIMATE",
            "EFFECT_NAME",
            "EFFECT",
            "SE",
            "CI_LOWER",
            "CI_UPPER",
            "P_VALUE",
            "P_VALUE_FDR",
            "COVARIATES",
            "MODEL_N",
            "HOM_N",
            "OTHER_N",
        ]
    ]

    result_df = result_df.sort_values(
        ["P_VALUE_FDR", "P_VALUE"],
        na_position="last"
    )

result_df.to_csv(
    f"{output_dir}/pgx_homozygote_vs_others_effect_size_min8_results.csv",
    index=False,
    encoding="utf-8-sig"
)