import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats as sps
from statsmodels.stats.multitest import multipletests


# Analysis frame (per professor's comments):
#   - infliximab only (adalimumab excluded; no PK model data)
#   - data split into phase subgroups: IND / MAINT / ALL (IND+MAINT pooled)
#   - within each subgroup, genotype effect on each endpoint is tested with
#     ANCOVA (CL) or logistic regression (ADA):
#       ADA ~ GROUP + SEX + WEIGHT + ALBUMIN
#       CL  ~ GROUP + SEX + WEIGHT + ALBUMIN + ADA
#   - CL is analyzed on both raw scale (BETA) and log scale (GMR,
#     log-normal CL assumption)
#   - two genotype grouping schemes:
#       HOM_vs_OTHERS:         dosage 2  vs dosage 0/1
#       CARRIER_vs_NONCARRIER: dosage 1/2 vs dosage 0
#   - FDR (Benjamini-Hochberg) is applied within each
#     (PHASE x END_POINT x CL scale x comparison) subgroup, across the
#     genotypes tested in that subgroup (multiple testing over genotypes)

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
    "rs767455": "TNFRSF1A",
    "rs1061622": "TNFRSF1B",
    "rs765249238": "SLCO2A1",
    "rs776813259": "SLCO2A1",
}

GROUP_COL = "GENO_GROUP"
MIN_GROUP_N = 8

COMPARISONS = {
    "HOM_vs_OTHERS": {
        "variant_dosages": [2],
        "reference_dosages": [0, 1],
        "variant_label": "HOM(2)",
        "reference_label": "OTHERS(0,1)",
    },
    "CARRIER_vs_NONCARRIER": {
        "variant_dosages": [1, 2],
        "reference_dosages": [0],
        "variant_label": "CARRIER(1,2)",
        "reference_label": "NONCARRIER(0)",
    },
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
        lcl = mean - 1.96 * se
        ucl = mean + 1.96 * se
        return f"{mean:.3f} ({lcl:.3f}-{ucl:.3f}), n={n}"

    return f"{mean:.3f}, n={n}"


def fmt_geomean_ci(x):
    x = pd.Series(x).dropna()
    x = x[x > 0]
    n = len(x)

    if n == 0:
        return "NA"

    log_x = np.log(x)
    mean_log = log_x.mean()
    se_log = log_x.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    gmean = np.exp(mean_log)

    if n > 1:
        lcl = np.exp(mean_log - 1.96 * se_log)
        ucl = np.exp(mean_log + 1.96 * se_log)
        return f"{gmean:.3f} ({lcl:.3f}-{ucl:.3f}), n={n}"

    return f"{gmean:.3f}, n={n}"


def fmt_binary_count(x):
    x = pd.Series(x).dropna()
    n = len(x)

    if n == 0:
        return "NA"

    count = int((x == 1).sum())
    pct = count / n * 100
    return f"{count}/{n} ({pct:.1f}%)"


def get_model_result(df, y_col, x_cols, model_type, group_col=GROUP_COL, alpha=0.05):
    nan_result = (
        np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, len(df[[y_col] + x_cols].dropna()), {},
    )

    model_df = df[[y_col] + x_cols].dropna().copy()

    if len(model_df) < len(x_cols) + 2:
        return nan_result

    if model_df[y_col].nunique() < 2 or model_df[group_col].nunique() < 2:
        return nan_result

    X = sm.add_constant(model_df[x_cols], has_constant="add")
    y = model_df[y_col]

    try:
        if model_type == "logistic":
            fit = sm.Logit(y, X).fit(disp=False)
        else:
            # standard ANCOVA (plain OLS) for both raw and log CL, so the
            # two scales differ only in the distributional assumption
            fit = sm.OLS(y, X).fit()

        beta = fit.params.get(group_col, np.nan)
        se = fit.bse.get(group_col, np.nan)
        pval = fit.pvalues.get(group_col, np.nan)

        # residual diagnostics (OLS only): Shapiro-Wilk normality p and skewness
        if model_type == "logistic":
            resid_shapiro_p = np.nan
            resid_skew = np.nan
        else:
            resid = np.asarray(fit.resid, dtype=float)
            resid_shapiro_p = (
                sps.shapiro(resid)[1] if 3 <= len(resid) <= 5000 else np.nan
            )
            resid_skew = sps.skew(resid) if len(resid) > 2 else np.nan

        # logistic -> OR, log-scale OLS -> GMR, raw OLS -> BETA
        if model_type in ("logistic", "log_ols"):
            effect = np.exp(beta)
            ci_lower = np.exp(beta - 1.96 * se)
            ci_upper = np.exp(beta + 1.96 * se)
        else:
            effect = beta
            ci_lower = beta - 1.96 * se
            ci_upper = beta + 1.96 * se

        sig_covar = {}

        for term in x_cols:
            if term == group_col:
                continue

            term_beta = fit.params.get(term, np.nan)
            term_se = fit.bse.get(term, np.nan)
            term_p = fit.pvalues.get(term, np.nan)

            if pd.isna(term_p) or term_p >= alpha:
                continue

            if model_type == "logistic":
                effect_type = "OR"
                term_effect = np.exp(term_beta)
                term_ci = (
                    np.exp(term_beta - 1.96 * term_se),
                    np.exp(term_beta + 1.96 * term_se),
                )
            elif model_type == "log_ols":
                effect_type = "GMR"
                term_effect = np.exp(term_beta)
                term_ci = (
                    np.exp(term_beta - 1.96 * term_se),
                    np.exp(term_beta + 1.96 * term_se),
                )
            else:
                effect_type = "BETA"
                term_effect = term_beta
                term_ci = (
                    term_beta - 1.96 * term_se,
                    term_beta + 1.96 * term_se,
                )

            sig_covar[term] = {
                "effect_type": effect_type,
                "effect": round_or_nan(term_effect, 4),
                "ci": (round_or_nan(term_ci[0], 4), round_or_nan(term_ci[1], 4)),
                "p": round_or_nan(term_p, 5),
            }

        return (
            effect, se, ci_lower, ci_upper, pval,
            resid_shapiro_p, resid_skew, len(model_df), sig_covar,
        )

    except Exception:
        return (
            np.nan, np.nan, np.nan, np.nan, np.nan,
            np.nan, np.nan, len(model_df), {},
        )


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
ep_df = ep_df[ep_df["DRUG"] == "infliximab"].copy()

if "WT" in ep_df.columns and "WEIGHT" not in ep_df.columns:
    ep_df["WEIGHT"] = ep_df["WT"]

if "ALB" in ep_df.columns and "ALBUMIN" not in ep_df.columns:
    ep_df["ALBUMIN"] = ep_df["ALB"]


result_rows = []

for phase in ["IND", "MAINT", "ALL"]:

    if phase == "ALL":
        phase_cond = ep_df["PHASE"].isin(["IND", "MAINT"])
    else:
        phase_cond = ep_df["PHASE"] == phase

    med_ep_df = ep_df[phase_cond].copy()

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

    uid_ep_df = uid_ep_df[uid_ep_df["CL"].isna() | (uid_ep_df["CL"] > 0)].copy()
    uid_ep_df["LOG_CL"] = np.log(uid_ep_df["CL"])

    for comparison, comp_def in COMPARISONS.items():

        for rsid in rsid_list:

            geno_df = rsid_df[["UID", rsid]].copy()
            geno_df = geno_df.rename(columns={rsid: "GENOTYPE_DOSAGE"})

            geno_df[GROUP_COL] = np.where(
                geno_df["GENOTYPE_DOSAGE"].isin(comp_def["variant_dosages"]),
                1,
                np.where(
                    geno_df["GENOTYPE_DOSAGE"].isin(comp_def["reference_dosages"]),
                    0,
                    np.nan,
                ),
            )

            analysis_df = uid_ep_df.merge(
                geno_df[["UID", GROUP_COL]],
                on="UID",
                how="inner",
            )

            total_n = analysis_df["UID"].nunique()

            # (END_POINT, model_type, y_col)
            model_specs = [
                ("ADA", "logistic", "ADA"),
                ("CL", "raw_ols", "CL"),
                ("CL", "log_ols", "LOG_CL"),
            ]

            for ep_col, model_type, y_col in model_specs:

                tmp_df = analysis_df.dropna(
                    subset=["UID", GROUP_COL, y_col]
                ).copy()

                group_counts = tmp_df[GROUP_COL].value_counts()
                variant_n = group_counts.get(1, 0)
                reference_n = group_counts.get(0, 0)

                if variant_n < MIN_GROUP_N or reference_n < MIN_GROUP_N:
                    continue

                available_n = tmp_df["UID"].nunique()
                data_availability = (
                    f"{available_n}/{total_n} ({available_n / total_n * 100:.1f}%)"
                    if total_n > 0 else "0/0 (NA)"
                )

                variant_vals = tmp_df.loc[tmp_df[GROUP_COL] == 1, ep_col]
                reference_vals = tmp_df.loc[tmp_df[GROUP_COL] == 0, ep_col]

                if model_type == "logistic":
                    effect_name = "OR"
                    scale_label = "binary"
                    covariates = ["SEX", "WEIGHT", "ALBUMIN"]
                    variant_est = fmt_binary_count(variant_vals)
                    reference_est = fmt_binary_count(reference_vals)
                elif model_type == "log_ols":
                    effect_name = "GMR"
                    scale_label = "log(CL)"
                    covariates = ["SEX", "WEIGHT", "ALBUMIN", "ADA"]
                    variant_est = fmt_geomean_ci(variant_vals)
                    reference_est = fmt_geomean_ci(reference_vals)
                else:
                    effect_name = "BETA"
                    scale_label = "raw CL"
                    covariates = ["SEX", "WEIGHT", "ALBUMIN", "ADA"]
                    variant_est = fmt_mean_ci(variant_vals)
                    reference_est = fmt_mean_ci(reference_vals)

                x_cols = [GROUP_COL] + covariates

                (
                    effect,
                    se,
                    ci_lower,
                    ci_upper,
                    p_value,
                    resid_shapiro_p,
                    resid_skew,
                    model_n,
                    sig_covar,
                ) = get_model_result(
                    tmp_df,
                    y_col=y_col,
                    x_cols=x_cols,
                    model_type=model_type,
                    alpha=0.05,
                )

                percent_change = (
                    (effect - 1) * 100
                    if pd.notna(effect) and effect_name == "GMR"
                    else np.nan
                )

                result_rows.append({
                    "DRUG": "infliximab",
                    "PHASE": phase,
                    "COMPARISON": comparison,
                    "RSID": rsid,
                    "GENE": rsid_gene_dict.get(rsid.split("(")[0], ""),
                    "END_POINT": ep_col,
                    "MODEL_SCALE": scale_label,
                    "DATA_AVAILABILITY": data_availability,
                    "VARIANT_GROUP": comp_def["variant_label"],
                    "REFERENCE_GROUP": comp_def["reference_label"],
                    "VARIANT_ESTIMATE": variant_est,
                    "REFERENCE_ESTIMATE": reference_est,
                    "EFFECT_NAME": effect_name,
                    "EFFECT": effect,
                    "PERCENT_CHANGE": percent_change,
                    "SE": se,
                    "CI_LOWER": ci_lower,
                    "CI_UPPER": ci_upper,
                    "P_VALUE": p_value,
                    "P_VALUE_FDR": np.nan,
                    "RESID_SHAPIRO_P": resid_shapiro_p,
                    "RESID_SKEW": resid_skew,
                    "RESID_NORMALITY_OK": (
                        "" if pd.isna(resid_shapiro_p)
                        else ("Y" if resid_shapiro_p >= 0.05 else "N")
                    ),
                    "COVARIATES": ", ".join(covariates),
                    "SIG_COVAR": sig_covar,
                    "MODEL_N": model_n,
                    "VARIANT_N": variant_n,
                    "REFERENCE_N": reference_n,
                })


result_df = pd.DataFrame(result_rows)

if len(result_df) > 0:
    # FDR within each phase-subgroup analysis, across the genotypes tested:
    # one correction family per (PHASE x END_POINT x MODEL_SCALE x COMPARISON)
    fdr_keys = ["PHASE", "END_POINT", "MODEL_SCALE", "COMPARISON"]

    for keys, sub_df in result_df.groupby(fdr_keys):
        result_df.loc[sub_df.index, "P_VALUE_FDR"] = fdr_adjust(
            sub_df["P_VALUE"].values
        )

    result_df = result_df.sort_values(
        ["PHASE", "END_POINT", "MODEL_SCALE", "COMPARISON",
         "P_VALUE_FDR", "P_VALUE"],
        na_position="last",
    )

result_df.to_csv(
    f"{output_dir}/pgx_infliximab_phase_subgroup_fdr_results.csv",
    index=False,
    encoding="utf-8-sig",
)

if len(result_df) > 0:
    sig_df = result_df[result_df["P_VALUE_FDR"] < 0.05]

    print(f"Total rows: {len(result_df)}")
    print(f"FDR-significant rows (q < 0.05): {len(sig_df)}")

    if len(sig_df) > 0:
        print(sig_df[
            ["PHASE", "COMPARISON", "RSID", "GENE", "END_POINT",
             "MODEL_SCALE", "EFFECT_NAME", "EFFECT",
             "P_VALUE", "P_VALUE_FDR"]
        ].to_string(index=False))
else:
    print("No results produced.")

sig_df.to_csv(
    f"{output_dir}/pgx_infliximab_phase_subgroup_fdr_significant_results.csv",
    index=False,
    encoding="utf-8-sig",)