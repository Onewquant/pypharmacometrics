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


def round_or_nan(x, digits=4):
    if pd.isna(x):
        return np.nan
    return round(float(x), digits)


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


def fmt_dosage_distribution(df, y_col):
    result = {}

    for dosage in [0, 1, 2]:
        frag = df.loc[df["GENOTYPE_DOSAGE"] == dosage, y_col]

        if len(frag.dropna()) == 0:
            result[dosage] = "NA"
            continue

        valid_values = frag.dropna().unique()

        if set(valid_values).issubset({0, 1, 0.0, 1.0}):
            result[dosage] = fmt_binary_count(frag)
        else:
            result[dosage] = fmt_mean_ci(frag)

    return result


def get_model_result(df, y_col, x_cols, endpoint_type, genotype_col, alpha=0.05):
    model_df = df[[y_col] + x_cols].dropna().copy()

    fail_effect_name = (
        "OR_PER_ALLELE" if endpoint_type == "binary"
        else "BETA_PER_ALLELE"
    )

    if len(model_df) < len(x_cols) + 2:
        return fail_effect_name, np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df), {}

    if model_df[y_col].nunique() < 2:
        return fail_effect_name, np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df), {}

    if model_df[genotype_col].nunique() < 2:
        return fail_effect_name, np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df), {}

    X = sm.add_constant(model_df[x_cols], has_constant="add")
    y = model_df[y_col]

    try:
        if endpoint_type == "binary":
            fit = sm.Logit(y, X).fit(disp=False)

            beta = fit.params.get(genotype_col, np.nan)
            se = fit.bse.get(genotype_col, np.nan)
            pval = fit.pvalues.get(genotype_col, np.nan)

            effect_name = "OR_PER_ALLELE"
            effect = np.exp(beta)
            ci_lower = np.exp(beta - 1.96 * se)
            ci_upper = np.exp(beta + 1.96 * se)

        else:
            fit = sm.OLS(y, X).fit()

            beta = fit.params.get(genotype_col, np.nan)
            se = fit.bse.get(genotype_col, np.nan)
            pval = fit.pvalues.get(genotype_col, np.nan)

            effect_name = "BETA_PER_ALLELE"
            effect = beta
            ci_lower = beta - 1.96 * se
            ci_upper = beta + 1.96 * se

        sig_covar = {}

        for term in x_cols:
            if term == genotype_col:
                continue

            term_beta = fit.params.get(term, np.nan)
            term_se = fit.bse.get(term, np.nan)
            term_p = fit.pvalues.get(term, np.nan)

            if pd.isna(term_p) or term_p >= alpha:
                continue

            if endpoint_type == "binary":
                term_effect_type = "OR"
                term_effect = np.exp(term_beta)
                term_ci_lower = np.exp(term_beta - 1.96 * term_se)
                term_ci_upper = np.exp(term_beta + 1.96 * term_se)
            else:
                term_effect_type = "BETA"
                term_effect = term_beta
                term_ci_lower = term_beta - 1.96 * term_se
                term_ci_upper = term_beta + 1.96 * term_se

            sig_covar[term] = {
                "effect_type": term_effect_type,
                "effect": round_or_nan(term_effect, 4),
                "ci": (
                    round_or_nan(term_ci_lower, 4),
                    round_or_nan(term_ci_upper, 4),
                ),
                "p": round_or_nan(term_p, 5),
            }

        return (
            effect_name,
            effect,
            se,
            ci_lower,
            ci_upper,
            pval,
            len(model_df),
            sig_covar,
        )

    except Exception:
        return fail_effect_name, np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df), {}


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

genotype_col = "GENOTYPE_DOSAGE"

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
            geno_df = geno_df.rename(columns={rsid: genotype_col})
            geno_df[genotype_col] = pd.to_numeric(
                geno_df[genotype_col],
                errors="coerce"
            )

            analysis_df = uid_ep_df.merge(
                geno_df[["UID", genotype_col]],
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
                    subset=["UID", genotype_col, ep_col]
                ).copy()

                dosage_counts = tmp_df[genotype_col].value_counts().to_dict()

                dosage0_n = int(dosage_counts.get(0, 0))
                dosage1_n = int(dosage_counts.get(1, 0))
                dosage2_n = int(dosage_counts.get(2, 0))

                # Additive model 최소 조건:
                # 전체 model n >= 16, genotype dosage가 최소 2개 이상 존재
                if tmp_df[genotype_col].nunique() < 2:
                    continue

                if len(tmp_df) < 16:
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
                else:
                    endpoint_type = "continuous"

                dosage_estimates = fmt_dosage_distribution(tmp_df, ep_col)

                if ep_col == "ADA":
                    covariates = ["SEX", "WEIGHT", "ALBUMIN"]
                elif drug == "adalimumab":
                    covariates = ["SEX", "WEIGHT", "ALBUMIN"]
                else:
                    covariates = ["SEX", "WEIGHT", "ALBUMIN", "ADA"]

                x_cols = [genotype_col] + covariates

                (
                    effect_name,
                    effect,
                    se,
                    ci_lower,
                    ci_upper,
                    p_value,
                    model_n,
                    sig_covar,
                ) = get_model_result(
                    tmp_df,
                    y_col=ep_col,
                    x_cols=x_cols,
                    endpoint_type=endpoint_type,
                    genotype_col=genotype_col,
                    alpha=0.05,
                )

                result_rows.append({
                    "DRUG": drug,
                    "PHASE": phase,
                    "RSID": rsid,
                    "GENE": rsid_gene_dict.get(rsid.split("(")[0], ""),
                    "END_POINT": ep_col,
                    "ENDPOINT_TYPE": endpoint_type,
                    "DATA_AVAILABILITY": data_availability,
                    "DOSAGE_0_ESTIMATE": dosage_estimates.get(0, "NA"),
                    "DOSAGE_1_ESTIMATE": dosage_estimates.get(1, "NA"),
                    "DOSAGE_2_ESTIMATE": dosage_estimates.get(2, "NA"),
                    "EFFECT_NAME": effect_name,
                    "EFFECT": effect,
                    "SE": se,
                    "CI_LOWER": ci_lower,
                    "CI_UPPER": ci_upper,
                    "P_VALUE": p_value,
                    "P_VALUE_FDR": np.nan,
                    "COVARIATES": ", ".join(covariates),
                    "SIG_COVAR": sig_covar,
                    "MODEL_N": model_n,
                    "DOSAGE_0_N": dosage0_n,
                    "DOSAGE_1_N": dosage1_n,
                    "DOSAGE_2_N": dosage2_n,
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
            "DOSAGE_0_ESTIMATE",
            "DOSAGE_1_ESTIMATE",
            "DOSAGE_2_ESTIMATE",
            "EFFECT_NAME",
            "EFFECT",
            "SE",
            "CI_LOWER",
            "CI_UPPER",
            "P_VALUE",
            "P_VALUE_FDR",
            "COVARIATES",
            "SIG_COVAR",
            "MODEL_N",
            "DOSAGE_0_N",
            "DOSAGE_1_N",
            "DOSAGE_2_N",
        ]
    ]

    result_df = result_df.sort_values(
        ["P_VALUE_FDR", "P_VALUE"],
        na_position="last"
    )

result_df.to_csv(
    f"{output_dir}/pgx_additive_genetic_model_effect_size_results.csv",
    index=False,
    encoding="utf-8-sig"
)