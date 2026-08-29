"""Paper-ready core tables and supplementary tables.

Consolidates the raw analysis outputs (01-06) into the numbered tables
used in the manuscript. Run 01-06 first.

Reporting frame (decision 2026-08, option A): the rs1061622 association is
reported as an exploratory finding - the effect estimate is consistent
across analysis periods but does not reach the FDR-adjusted significance
threshold.

Outputs -> paper_works_new/core_fig_tab/
  Table1_baseline_characteristics.csv
  Table2_popPK_parameters.csv                (provisional - see note)
  Table3_variant_characteristics.csv
  Table4_CL_association_overall.csv
  Table5_rs1061622_across_periods.csv
  SupplTableS1_CL_association_by_period.csv
  SupplTableS2_ADA_association.csv
  SupplTableS3_sensitivity_analyses.csv
  SupplTableS4_CL_association_original_scale.csv
  SupplTableS5_cohort_attrition.csv
"""

import os
import shutil

import pandas as pd

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
out_dir = f"{prj_dir}/paper_works_new/output"
cft_dir = f"{prj_dir}/paper_works_new/core_fig_tab"

LEAD_RSID = "rs1061622"
PERIOD_LABEL = {
    "OVERALL": "Overall treatment",
    "IND": "Induction phase",
    "MAINT": "Maintenance phase",
}
PERIOD_ORDER = ["OVERALL", "IND", "MAINT"]
MODEL_LABEL = {
    "HOM_vs_OTHERS": "Recessive",
    "CARRIER_vs_NONCARRIER": "Dominant",
}


def fmt_est(effect, lo, hi, digits=2):
    if pd.isna(effect):
        return "NA"
    return f"{effect:.{digits}f} ({lo:.{digits}f}-{hi:.{digits}f})"


def fmt_p(p):
    if pd.isna(p):
        return "NA"
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def genotype_labels(rsid_full):
    """rs1061622(0=T,1=G) -> (rs, ref, alt, [TT, TG, GG])."""
    rs = rsid_full.split("(")[0]
    ref = alt = "?"
    if "(" in rsid_full:
        body = rsid_full[rsid_full.find("(") + 1: rsid_full.rfind(")")]
        for part in body.split(","):
            if part.startswith("0="):
                ref = part[2:]
            elif part.startswith("1="):
                alt = part[2:]
    return rs, ref, alt, [ref + ref, ref + alt, alt + alt]


res = pd.read_csv(f"{out_dir}/Table_pgx_ancova_fdr_results.csv")
geno = pd.read_csv(f"{out_dir}/Table_genotype_summary.csv")
vqc = pd.read_csv(f"{out_dir}/Table_variant_qc.csv")
attr = pd.read_csv(f"{out_dir}/Table_pgx_attrition.csv")

res["RS"] = res["RSID"].str.split("(").str[0]

# ---------------------------------------------------------------- Table 1
shutil.copy(f"{out_dir}/Table1_demographics.csv",
            f"{cft_dir}/Table1_baseline_characteristics.csv")

# ---------------------------------------------------------------- Table 2
# Transcribed from the final NONMEM run and the bootstrap summary.
# NOTE: these estimates predate the body-weight data correction;
# re-estimation is pending (see README and the manuscript note).
popPK_rows = [
    ("Typical values", "", "", ""),
    ("CL (L/day)", "0.310 (9)", "0.295 (6)", "0.291 (0.267-0.330)"),
    ("  ADA on CL", "-", "0.712 (22)", "0.763 (0.421-1.054)"),
    ("  Sex on CL", "-", "-0.113 (68)", "-0.098 (-0.225-0.049)"),
    ("  Weight on CL", "-", "0.484 (47)", "0.514 (0.224-0.945)"),
    ("Vc (L)", "3.990 (12)", "4.490 (7)", "4.418 (3.944-4.856)"),
    ("  Albumin on Vc", "-", "1.410 (27)", "1.443 (0.607-2.083)"),
    ("  Weight on Vc", "-", "0.800 (35)", "0.919 (0.191-1.460)"),
    ("Vp (L)", "0.898 (37)", "0.407 (26)", "0.450 (0.306-0.689)"),
    ("Q (L/day)", "0.0646 FIX", "0.0646 FIX", "0.065 (0.065-0.065)"),
    ("F1", "0.667 FIX", "0.667 FIX", "0.667 (0.667-0.667)"),
    ("Ka (day-1)", "0.058 (11)", "0.055 (7)", "0.055 (0.045-0.095)"),
    ("Interindividual variability", "", "", ""),
    ("IIV on CL (CV%)", "43.5 (21) [6]", "27.3 (13) [8]", "27.5 (20.7-37.4)"),
    ("Residual variability", "", "", ""),
    ("Proportional error", "0.381 (8)", "0.380 (8)", "0.386 (0.337-0.440)"),
    ("Additive error", "0 FIX", "0 FIX", "0.000 (0.000-0.000)"),
]
pd.DataFrame(popPK_rows, columns=[
    "Parameter",
    "Base model (OFV = -855.775), Estimate (%RSE) [%shrinkage]",
    "Final model (OFV = -935.665), Estimate (%RSE) [%shrinkage]",
    "Bootstrap, Median (5th-95th percentile)",
]).to_csv(f"{cft_dir}/Table2_popPK_parameters.csv",
          index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- Table 3
tested = set(res["RS"])
t3 = []
for _, g in geno.iterrows():
    rs = g["RSID"]
    full = rs + str(g["ALLELE_CODING"])
    _, ref, alt, labels = genotype_labels(full)
    counts = str(g["IFX_COHORT_GENO_0/1/2"]).split("/")
    qcr = vqc[vqc["RSID"] == rs]
    passed = qcr["PASS"].iloc[0] if len(qcr) else "N"
    reason = qcr["REASON"].iloc[0] if len(qcr) else ""
    t3.append({
        "rsID": rs,
        "Gene": g["GENE"],
        "Alleles (ref>alt)": f"{ref}>{alt}",
        "Genotype counts": ", ".join(
            f"{labels[i]} {counts[i]}" for i in range(3)
        ),
        "MAF": round(float(g["IFX_COHORT_MAF"]), 3),
        "HWE P": ("NA" if pd.isna(g["IFX_COHORT_HWE_P"])
                  else round(float(g["IFX_COHORT_HWE_P"]), 3)),
        "Passed variant QC": "Yes" if passed == "Y" else "No",
        "Reason for exclusion": "" if passed == "Y" else reason,
        "Tested in association analysis": "Yes" if rs in tested else "No",
    })
pd.DataFrame(t3).to_csv(f"{cft_dir}/Table3_variant_characteristics.csv",
                        index=False, encoding="utf-8-sig")


def association_table(df, periods, effect_label):
    rows = []
    for period in periods:
        for comp in ["HOM_vs_OTHERS", "CARRIER_vs_NONCARRIER"]:
            sub = df[(df["PHASE"] == period) & (df["COMPARISON"] == comp)]
            for _, r in sub.sort_values("P_VALUE").iterrows():
                _, ref, alt, labels = genotype_labels(r["RSID"])
                if comp == "HOM_vs_OTHERS":
                    grp = f"{labels[2]} vs {labels[0]}+{labels[1]}"
                else:
                    grp = f"{labels[1]}+{labels[2]} vs {labels[0]}"
                rows.append({
                    "Analysis period": PERIOD_LABEL[period],
                    "Genetic model": MODEL_LABEL[comp],
                    "rsID": r["RS"],
                    "Gene": r["GENE"],
                    "Comparison": grp,
                    "n (variant group)": int(r["VARIANT_N"]),
                    "n (reference group)": int(r["REFERENCE_N"]),
                    effect_label: fmt_est(
                        r["EFFECT"], r["CI_LOWER"], r["CI_UPPER"]
                    ),
                    "P value": fmt_p(r["P_VALUE"]),
                    "FDR q value": fmt_p(r["P_VALUE_FDR"]),
                })
    return pd.DataFrame(rows)


cl_log = res[(res["END_POINT"] == "CL") & (res["MODEL_SCALE"] == "log(CL)")]
cl_raw = res[(res["END_POINT"] == "CL") & (res["MODEL_SCALE"] == "raw CL")]
ada = res[res["END_POINT"] == "ADA"]

association_table(cl_log, ["OVERALL"], "GMR (95% CI)").to_csv(
    f"{cft_dir}/Table4_CL_association_overall.csv",
    index=False, encoding="utf-8-sig")

association_table(cl_log, ["IND", "MAINT"], "GMR (95% CI)").to_csv(
    f"{cft_dir}/SupplTableS1_CL_association_by_period.csv",
    index=False, encoding="utf-8-sig")

association_table(ada, PERIOD_ORDER, "OR (95% CI)").to_csv(
    f"{cft_dir}/SupplTableS2_ADA_association.csv",
    index=False, encoding="utf-8-sig")

association_table(
    cl_raw, PERIOD_ORDER, "Adjusted difference (95% CI), L/day"
).to_csv(f"{cft_dir}/SupplTableS4_CL_association_original_scale.csv",
         index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- Table 5
lead = res[(res["RS"] == LEAD_RSID)
           & (res["COMPARISON"] == "HOM_vs_OTHERS")
           & (res["END_POINT"] == "CL")]

t5 = []
for period in PERIOD_ORDER:
    lg = lead[(lead["PHASE"] == period) & (lead["MODEL_SCALE"] == "log(CL)")]
    rw = lead[(lead["PHASE"] == period) & (lead["MODEL_SCALE"] == "raw CL")]
    if len(lg) == 0:
        continue
    lg = lg.iloc[0]
    row = {
        "Analysis period": PERIOD_LABEL[period],
        "n (GG)": int(lg["VARIANT_N"]),
        "n (TT+TG)": int(lg["REFERENCE_N"]),
        "Geometric mean CL, GG (95% CI), L/day": lg["VARIANT_ESTIMATE"],
        "Geometric mean CL, TT+TG (95% CI), L/day": lg["REFERENCE_ESTIMATE"],
        "GMR (95% CI)": fmt_est(lg["EFFECT"], lg["CI_LOWER"], lg["CI_UPPER"]),
        "P value": fmt_p(lg["P_VALUE"]),
        "FDR q value": fmt_p(lg["P_VALUE_FDR"]),
        "Shapiro-Wilk P (log-scale residuals)": fmt_p(lg["RESID_SHAPIRO_P"]),
    }
    if len(rw):
        rw = rw.iloc[0]
        row["Original-scale adjusted difference (95% CI), L/day"] = fmt_est(
            rw["EFFECT"], rw["CI_LOWER"], rw["CI_UPPER"], digits=3)
        row["Original-scale FDR q value"] = fmt_p(rw["P_VALUE_FDR"])
    t5.append(row)

pd.DataFrame(t5).to_csv(f"{cft_dir}/Table5_rs1061622_across_periods.csv",
                        index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- Suppl S3
sens_path = f"{out_dir}/Table_pgx_sensitivity.csv"
s3 = pd.DataFrame([{
    "Note": "Script 04 runs sensitivity analyses only for FDR-significant "
            "CL associations; none were significant under the final "
            "analysis frame, so no rows were produced. Use "
            "04_pgx_sensitivity.py with a relaxed trigger to reproduce the "
            "robustness checks reported in the manuscript."
}])
if os.path.exists(sens_path):
    try:
        sens = pd.read_csv(sens_path)
    except pd.errors.EmptyDataError:
        sens = pd.DataFrame()
    if len(sens):
        sens["RS"] = sens["RSID"].str.split("(").str[0]
        sub = sens[sens["RS"] == LEAD_RSID]
        if len(sub):
            s3 = pd.DataFrame({
                "Analysis period": sub["PHASE"].map(PERIOD_LABEL),
                "CL scale": sub["MODEL_SCALE"],
                "Primary P value": sub["MAIN_P"].map(fmt_p),
                "Leave-one-out P range": [
                    f"{a:.3f}-{b:.3f}"
                    for a, b in zip(sub["LOO_P_MIN"], sub["LOO_P_MAX"])
                ],
                "All leave-one-out P < 0.05": sub["LOO_ALL_BELOW_0.05"],
                "Mann-Whitney P": sub["MANN_WHITNEY_P"].map(fmt_p),
                "HC3 robust P": sub["HC3_ROBUST_P"].map(fmt_p),
                "P excluding patients without observed concentrations":
                    sub["EXCL_NO_SAMPLE_P"].map(fmt_p),
            })
s3.to_csv(f"{cft_dir}/SupplTableS3_sensitivity_analyses.csv",
          index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- Suppl S5
attr.to_csv(f"{cft_dir}/SupplTableS5_cohort_attrition.csv",
            index=False, encoding="utf-8-sig")

print("core_fig_tab tables written:")
for f in sorted(os.listdir(cft_dir)):
    if f.endswith(".csv"):
        print(f"  {f}  ({len(pd.read_csv(f'{cft_dir}/{f}'))} rows)")
