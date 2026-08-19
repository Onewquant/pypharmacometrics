import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


# Reflects professor's feedback (2026-07):
#   - adalimumab has no PK model data -> infliximab only
#   - FDR correction must be done separately within each
#     (PHASE x ENDPOINT) stratum, across the genotypes tested in that stratum
#     (not pooled across drug/phase/endpoint)
#   - "whole phase" should not be a naive pooled model; induction and
#     maintenance are analyzed as separate strata (method 1), and a
#     supplementary GENOTYPE x PHASE interaction model (method 2) is
#     provided to test whether the effect differs across phases
#   - CL is assumed log-normal, so CL is modeled on the log scale
#     (residuals ~ normal), reported as GMR; ADA is modeled as OR
#   - ADA is included as a covariate in the CL model

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

GROUP_COL = "HOM_GROUP"
PHASE_DUMMY_COL = "PHASE_MAINT"
INTERACTION_COL = "INTERACTION"


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


def fit_endpoint_model(y, X, endpoint_type, cluster_groups=None):
    if endpoint_type == "binary":
        if cluster_groups is not None:
            return sm.Logit(y, X).fit(
                disp=False, cov_type="cluster", cov_kwds={"groups": cluster_groups}
            )
        return sm.Logit(y, X).fit(disp=False)

    if cluster_groups is not None:
        return sm.OLS(y, X).fit(
            cov_type="cluster", cov_kwds={"groups": cluster_groups}
        )
    return sm.OLS(y, X).fit(cov_type="HC3")


def extract_term(fit, term):
    beta = fit.params.get(term, np.nan)
    se = fit.bse.get(term, np.nan)
    p = fit.pvalues.get(term, np.nan)

    if pd.isna(beta) or pd.isna(se):
        return np.nan, np.nan, np.nan, np.nan, np.nan

    effect = np.exp(beta)
    ci_lower = np.exp(beta - 1.96 * se)
    ci_upper = np.exp(beta + 1.96 * se)

    return effect, se, ci_lower, ci_upper, p


def run_stratified_model(df, y_col, x_cols, endpoint_type, group_col=GROUP_COL):
    model_df = df[[y_col] + x_cols].dropna().copy()

    if len(model_df) < len(x_cols) + 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df)

    if model_df[y_col].nunique() < 2 or model_df[group_col].nunique() < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df)

    X = sm.add_constant(model_df[x_cols], has_constant="add")
    y = model_df[y_col]

    try:
        fit = fit_endpoint_model(y, X, endpoint_type)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(model_df)

    effect, se, ci_lower, ci_upper, p = extract_term(fit, group_col)
    return effect, se, ci_lower, ci_upper, p, len(model_df)


def run_interaction_model(df, y_col, covariates, endpoint_type, group_col=GROUP_COL):
    work_df = df.copy()
    work_df[INTERACTION_COL] = work_df[group_col] * work_df[PHASE_DUMMY_COL]

    x_cols = [group_col, PHASE_DUMMY_COL, INTERACTION_COL] + covariates
    model_df = work_df[["UID", y_col] + x_cols].dropna().copy()

    nan_main = (np.nan,) * 5
    nan_result = nan_main, nan_main, len(model_df)

    if len(model_df) < len(x_cols) + 3:
        return nan_result

    if model_df[y_col].nunique() < 2:
        return nan_result

    if model_df[group_col].nunique() < 2 or model_df[PHASE_DUMMY_COL].nunique() < 2:
        return nan_result

    X = sm.add_constant(model_df[x_cols], has_constant="add")
    y = model_df[y_col]

    try:
        fit = fit_endpoint_model(
            y, X, endpoint_type, cluster_groups=model_df["UID"]
        )
    except Exception:
        return nan_result

    # main effect = HOM_GROUP effect at induction (PHASE_MAINT == 0, reference)
    main_effect = extract_term(fit, group_col)
    # interaction effect = how much the HOM_GROUP effect changes at maintenance vs induction
    interaction_effect = extract_term(fit, INTERACTION_COL)

    return main_effect, interaction_effect, len(model_df)


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


def build_uid_phase_df(ep_df, phase):
    med_ep_df = ep_df[ep_df["PHASE"] == phase].copy()

    uid_df = (
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
    uid_df["PHASE"] = phase
    return uid_df


def get_genotype_group(rsid):
    geno_df = rsid_df[["UID", rsid]].copy()
    geno_df = geno_df.rename(columns={rsid: "GENOTYPE_DOSAGE"})

    # homozygote variant: dosage == 2, others: dosage == 0 or 1
    geno_df[GROUP_COL] = np.where(
        geno_df["GENOTYPE_DOSAGE"] == 2,
        1,
        np.where(geno_df["GENOTYPE_DOSAGE"].isin([0, 1]), 0, np.nan),
    )
    return geno_df[["UID", GROUP_COL]]


ENDPOINT_LIST = ["ADA", "CL"]

# ---------------------------------------------------------------------------
# 1) Stratified analysis (method 1): induction and maintenance analyzed
#    separately; FDR correction performed within each PHASE x ENDPOINT stratum
# ---------------------------------------------------------------------------

stratified_rows = []

for phase in ["IND", "MAINT"]:

    uid_df = build_uid_phase_df(ep_df, phase)
    uid_df = uid_df[uid_df["CL"].isna() | (uid_df["CL"] > 0)].copy()
    uid_df["LOG_CL"] = np.log(uid_df["CL"])

    for rsid in rsid_list:

        geno_df = get_genotype_group(rsid)
        analysis_df = uid_df.merge(geno_df, on="UID", how="inner")
        total_n = analysis_df["UID"].nunique()

        for ep_col in ENDPOINT_LIST:

            y_col = "LOG_CL" if ep_col == "CL" else "ADA"
            endpoint_type = "continuous_log" if ep_col == "CL" else "binary"

            tmp_df = analysis_df.dropna(subset=["UID", GROUP_COL, y_col]).copy()

            group_counts = tmp_df[GROUP_COL].value_counts()
            hom_n = group_counts.get(1, 0)
            other_n = group_counts.get(0, 0)

            if hom_n < 8 or other_n < 8:
                continue

            available_n = tmp_df["UID"].nunique()
            data_availability = (
                f"{available_n}/{total_n} ({available_n / total_n * 100:.1f}%)"
                if total_n > 0 else "0/0 (NA)"
            )

            if ep_col == "ADA":
                effect_name = "OR"
                covariates = ["SEX", "WEIGHT", "ALBUMIN"]
                hom_est = fmt_binary_count(tmp_df.loc[tmp_df[GROUP_COL] == 1, "ADA"])
                other_est = fmt_binary_count(tmp_df.loc[tmp_df[GROUP_COL] == 0, "ADA"])
            else:
                effect_name = "GMR"
                covariates = ["SEX", "WEIGHT", "ALBUMIN", "ADA"]
                hom_est = fmt_geomean_ci(tmp_df.loc[tmp_df[GROUP_COL] == 1, "CL"])
                other_est = fmt_geomean_ci(tmp_df.loc[tmp_df[GROUP_COL] == 0, "CL"])

            x_cols = [GROUP_COL] + covariates

            effect, se, ci_lower, ci_upper, p_value, model_n = run_stratified_model(
                tmp_df, y_col=y_col, x_cols=x_cols, endpoint_type=endpoint_type
            )

            percent_change = (
                (effect - 1) * 100 if pd.notna(effect) and ep_col == "CL" else np.nan
            )

            stratified_rows.append({
                "DRUG": "infliximab",
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
                "PERCENT_CHANGE": percent_change,
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


stratified_df = pd.DataFrame(stratified_rows)

if len(stratified_df) > 0:
    for keys, sub_df in stratified_df.groupby(["PHASE", "END_POINT"]):
        stratified_df.loc[sub_df.index, "P_VALUE_FDR"] = fdr_adjust(
            sub_df["P_VALUE"].values
        )

    stratified_df = stratified_df.sort_values(
        ["PHASE", "END_POINT", "P_VALUE_FDR", "P_VALUE"], na_position="last"
    )

stratified_df.to_csv(
    f"{output_dir}/pgx_hom_vs_others_logCL_infliximab_stratified_results.csv",
    index=False,
    encoding="utf-8-sig",
)


# ---------------------------------------------------------------------------
# 2) Supplementary GENOTYPE x PHASE interaction model (method 2): pools
#    induction + maintenance, subject-clustered SE, tests whether the
#    genotype effect differs by phase. FDR correction performed separately
#    for the main-effect p-values and the interaction p-values, each within
#    its own ENDPOINT stratum.
# ---------------------------------------------------------------------------

all_df = pd.concat(
    [build_uid_phase_df(ep_df, phase) for phase in ["IND", "MAINT"]],
    ignore_index=True,
)
all_df = all_df[all_df["CL"].isna() | (all_df["CL"] > 0)].copy()
all_df["LOG_CL"] = np.log(all_df["CL"])
all_df[PHASE_DUMMY_COL] = (all_df["PHASE"] == "MAINT").astype(int)

interaction_rows = []

for rsid in rsid_list:

    geno_df = get_genotype_group(rsid)
    analysis_df = all_df.merge(geno_df, on="UID", how="inner")
    total_n = analysis_df["UID"].nunique()

    for ep_col in ENDPOINT_LIST:

        y_col = "LOG_CL" if ep_col == "CL" else "ADA"
        endpoint_type = "continuous_log" if ep_col == "CL" else "binary"

        tmp_df = analysis_df.dropna(subset=["UID", GROUP_COL, y_col]).copy()

        group_counts = tmp_df[GROUP_COL].value_counts()
        hom_n = group_counts.get(1, 0)
        other_n = group_counts.get(0, 0)

        if hom_n < 8 or other_n < 8:
            continue

        available_n = tmp_df["UID"].nunique()
        data_availability = (
            f"{available_n}/{total_n} ({available_n / total_n * 100:.1f}%)"
            if total_n > 0 else "0/0 (NA)"
        )

        ind_df = tmp_df[tmp_df["PHASE"] == "IND"]
        maint_df = tmp_df[tmp_df["PHASE"] == "MAINT"]

        if ep_col == "ADA":
            effect_name = "OR"
            covariates = ["SEX", "WEIGHT", "ALBUMIN"]
            ind_hom_est = fmt_binary_count(ind_df.loc[ind_df[GROUP_COL] == 1, "ADA"])
            ind_other_est = fmt_binary_count(ind_df.loc[ind_df[GROUP_COL] == 0, "ADA"])
            maint_hom_est = fmt_binary_count(maint_df.loc[maint_df[GROUP_COL] == 1, "ADA"])
            maint_other_est = fmt_binary_count(maint_df.loc[maint_df[GROUP_COL] == 0, "ADA"])
        else:
            effect_name = "GMR"
            covariates = ["SEX", "WEIGHT", "ALBUMIN", "ADA"]
            ind_hom_est = fmt_geomean_ci(ind_df.loc[ind_df[GROUP_COL] == 1, "CL"])
            ind_other_est = fmt_geomean_ci(ind_df.loc[ind_df[GROUP_COL] == 0, "CL"])
            maint_hom_est = fmt_geomean_ci(maint_df.loc[maint_df[GROUP_COL] == 1, "CL"])
            maint_other_est = fmt_geomean_ci(maint_df.loc[maint_df[GROUP_COL] == 0, "CL"])

        main_effect, interaction_effect, model_n = run_interaction_model(
            tmp_df, y_col=y_col, covariates=covariates, endpoint_type=endpoint_type
        )

        main_eff, main_se, main_ci_lower, main_ci_upper, main_p = main_effect
        (
            inter_eff, inter_se, inter_ci_lower, inter_ci_upper, inter_p,
        ) = interaction_effect

        interaction_rows.append({
            "DRUG": "infliximab",
            "PHASE": "ALL_interaction",
            "RSID": rsid,
            "GENE": rsid_gene_dict.get(rsid.split("(")[0], ""),
            "END_POINT": ep_col,
            "ENDPOINT_TYPE": endpoint_type,
            "DATA_AVAILABILITY": data_availability,
            "IND_HOM_ESTIMATE": ind_hom_est,
            "IND_OTHER_ESTIMATE": ind_other_est,
            "MAINT_HOM_ESTIMATE": maint_hom_est,
            "MAINT_OTHER_ESTIMATE": maint_other_est,
            "EFFECT_NAME": effect_name,
            # main effect: HOM_GROUP effect at induction (reference phase)
            "MAIN_EFFECT_IND": main_eff,
            "MAIN_SE": main_se,
            "MAIN_CI_LOWER": main_ci_lower,
            "MAIN_CI_UPPER": main_ci_upper,
            "MAIN_P_VALUE": main_p,
            "MAIN_P_VALUE_FDR": np.nan,
            # interaction effect: ratio of (HOM_GROUP effect at maintenance) to
            # (HOM_GROUP effect at induction) -- tests if the effect differs by phase
            "INTERACTION_EFFECT": inter_eff,
            "INTERACTION_SE": inter_se,
            "INTERACTION_CI_LOWER": inter_ci_lower,
            "INTERACTION_CI_UPPER": inter_ci_upper,
            "INTERACTION_P_VALUE": inter_p,
            "INTERACTION_P_VALUE_FDR": np.nan,
            "COVARIATES": ", ".join(covariates + [PHASE_DUMMY_COL]),
            "MODEL_N": model_n,
            "HOM_N": hom_n,
            "OTHER_N": other_n,
        })


interaction_df = pd.DataFrame(interaction_rows)

if len(interaction_df) > 0:
    for ep_col, sub_df in interaction_df.groupby("END_POINT"):
        interaction_df.loc[sub_df.index, "MAIN_P_VALUE_FDR"] = fdr_adjust(
            sub_df["MAIN_P_VALUE"].values
        )
        interaction_df.loc[sub_df.index, "INTERACTION_P_VALUE_FDR"] = fdr_adjust(
            sub_df["INTERACTION_P_VALUE"].values
        )

    interaction_df = interaction_df.sort_values(
        ["END_POINT", "INTERACTION_P_VALUE_FDR", "INTERACTION_P_VALUE"],
        na_position="last",
    )

interaction_df.to_csv(
    f"{output_dir}/pgx_hom_vs_others_logCL_infliximab_interaction_results.csv",
    index=False,
    encoding="utf-8-sig",
)

print(f"Stratified results: {len(stratified_df)} rows")
print(f"Interaction results: {len(interaction_df)} rows")
