"""Build the two Word deliverables for the manuscript.

  manuscript/Methods_and_Results_FINAL.docx
      Methods + Results + a Discussion draft for the pharmacogenomic part.
      The popPK sections are carried over verbatim from the user-reviewed
      version (output/Methods_and_results_reviewed.docx); every
      pharmacogenomic paragraph is replaced with the EBE-based text, and
      TLR4 is removed from the candidate-gene list (rs5030728 is absent
      from the genotype matrix).
      Format (2026-09-07, per user request): plain manuscript layout like
      the reviewed version - no document title, no editorial notes, no
      to-do section; "[Methods]" / "[Results]" / "[Discussion]" as bold
      section headers and bold subsection titles. Editorial notes and open
      items live in manuscript/Methods_and_Results_FINAL.md and README.md.

  core_fig_tab/[IFX_POPPK]_core_fig_tab_FINAL.docx
      Table 1-5, Figure 1-3, Supplementary Table S1-S5 and Supplementary
      Figure S1-S3 with captions, built from the CSV/PNG files in
      core_fig_tab/ so the numbers cannot drift from the analysis output.

Run 01-09 first. Requires python-docx.
"""

import os

import pandas as pd
from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
cft_dir = f"{prj_dir}/paper_works_new/core_fig_tab"
ms_dir = f"{prj_dir}/paper_works_new/manuscript"

BODY_PT = 11
NOTE_COLOR = RGBColor(0x8A, 0x60, 0x00)


def new_doc():
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(BODY_PT)
    style.paragraph_format.space_after = Pt(6)
    return doc


def h(doc, text, level=1):
    """Bold in-line heading, as in the reviewed manuscript draft.

    level 1 -> "[Section]" (Methods / Results / Discussion)
    level 2 -> bold subsection title
    """
    p = doc.add_paragraph()
    label = f"[{text}]" if level == 1 else text
    run = p.add_run(label)
    run.bold = True
    run.font.size = Pt(BODY_PT)
    p.paragraph_format.space_before = Pt(12 if level == 1 else 8)
    p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.keep_with_next = True
    return p


def para(doc, text, italic=False, size=None):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.italic = italic
    if size:
        run.font.size = Pt(size)
    return p


def bullets(doc, items):
    for it in items:
        doc.add_paragraph(it, style="List Bullet")


def note(doc, text):
    """Editorial note to the co-authors, not manuscript body text."""
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.italic = True
    run.font.size = Pt(9)
    run.font.color.rgb = NOTE_COLOR
    return p


def add_csv_table(doc, csv_path, font_size=8, first_col_width=None):
    df = pd.read_csv(csv_path, dtype=str).fillna("")
    table = doc.add_table(rows=1, cols=len(df.columns))
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    hdr = table.rows[0].cells
    for j, col in enumerate(df.columns):
        hdr[j].text = ""
        run = hdr[j].paragraphs[0].add_run(str(col))
        run.bold = True
        run.font.size = Pt(font_size)
    for _, row in df.iterrows():
        cells = table.add_row().cells
        for j, col in enumerate(df.columns):
            cells[j].text = ""
            run = cells[j].paragraphs[0].add_run(str(row[col]))
            run.font.size = Pt(font_size)
    if first_col_width:
        for r in table.rows:
            r.cells[0].width = Inches(first_col_width)
    return table


def add_figure(doc, png_path, width_in=6.2):
    doc.add_picture(png_path, width=Inches(width_in))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER


def caption(doc, label, text, size=9):
    p = doc.add_paragraph()
    r1 = p.add_run(label + " ")
    r1.bold = True
    r1.font.size = Pt(size)
    r2 = p.add_run(text)
    r2.font.size = Pt(size)
    return p


# ===========================================================================
# 1) Methods and Results
# ===========================================================================
doc = new_doc()

h(doc, "Methods", 1)

h(doc, "Study Population and Data Collection", 2)
para(doc,
     "This retrospective population pharmacokinetic (PopPK) study included "
     "patients with inflammatory bowel disease (IBD) who received anti-tumor "
     "necrosis factor (anti-TNF) inhibitor therapy for at least one day at a "
     "tertiary referral hospital. Patients without available anti-TNF "
     "inhibitor concentration measurements were excluded from the study "
     "population, and those without whole-genome sequencing (WGS) data were "
     "further excluded to define the overall analytic cohort (n = 139). From "
     "this cohort, patients without infliximab administration records (n = 41) "
     "and one patient whose data started in the maintenance phase with only a "
     "single concentration measurement, which was used to initialize the model "
     "and could not inform individual parameter estimation, were excluded "
     "(n = 42 in total); the infliximab PopPK modeling cohort therefore "
     "comprised 97 patients.")
para(doc,
     "For the pharmacogenomic analyses, one additional patient was excluded "
     "because the corresponding sequencing sample was removed during genotype "
     "quality control, leaving 96 patients in the pharmacogenomic analysis "
     "cohort. Pharmacogenomic analyses were conducted separately for the "
     "overall treatment period (n = 96), the induction phase (n = 83), and the "
     "maintenance phase (n = 96), based on the availability of the "
     "corresponding clearance estimates and complete covariate data (Figure 1, "
     "Supplementary Table S5).")
para(doc,
     "Demographic characteristics, laboratory test results, therapeutic drug "
     "monitoring (TDM) measurements, anti-drug antibody (ADA) status, dosing "
     "records, and treatment-related variables were extracted from the "
     "electronic medical record system.")

h(doc, "Genotyping, Quality Control, and Variant Processing", 2)
para(doc,
     "Peripheral blood-derived genomic DNA samples were processed for "
     "whole-genome sequencing (WGS) by Macrogen Inc. (Seoul, Republic of "
     "Korea). Sequencing libraries were prepared using the TruSeq DNA "
     "PCR-Free (350 bp insert) library preparation kit and subjected to "
     "library quality control (QC) according to the manufacturer\u2019s and "
     "sequencing facility\u2019s protocols. Library concentration and fragment "
     "size distributions were evaluated using qPCR and Agilent 2100 "
     "Bioanalyzer analysis, respectively. For sequencing on the NovaSeq X "
     "platform, libraries with concentrations >5 nM satisfied the QC criteria "
     "for downstream sequencing.")
para(doc,
     "Initial variant calling and preprocessing were performed by the "
     "Precision Medicine Center, and variant call format (VCF) files were "
     "generated after application of variant quality score recalibration "
     "(VQSR) filtering and Hardy\u2013Weinberg equilibrium (HWE) filtering "
     "(P \u2265 1\u00d710\u207b\u2076). Additional sample- and variant-level QC "
     "procedures were subsequently performed using Hail (version 0.2) in "
     "Python.")
para(doc,
     "For sample-level QC, call rate, genotype quality (GQ), sequencing depth "
     "(DP), and transition/transversion (Ti/Tv) ratios were evaluated. Samples "
     "with low call rates (<0.95) were considered for exclusion. Ti/Tv ratios "
     "were calculated using biallelic single nucleotide variants (SNVs) with "
     "non-reference genotypes. Additional QC metrics, including heterozygosity "
     "assessment, relatedness estimation, and principal component analysis "
     "(PCA)-based outlier detection, were reviewed according to "
     "recommendations from the Precision Medicine Center.")
para(doc,
     "Variant-level QC was performed after sample QC. Variants with call rates "
     "<0.95 were excluded. Minor allele frequency (MAF) filtering was applied "
     "to exclude rare variants with MAF <0.01, and analyses were restricted to "
     "biallelic variants.")
para(doc,
     "For the candidate-variant association analyses, an additional "
     "variant-level quality control step was applied within the "
     "pharmacogenomic analysis cohort using genotype data only, without "
     "reference to any endpoint. Variants with a MAF below 0.05 or an HWE "
     "exact test P value below 0.05 were excluded. Of the 17 candidate "
     "variants, 14 passed this filter; three were excluded (IL10 rs3024505, "
     "MAF 0.036; SLCO2A1 rs765249238, MAF 0.016 and HWE P = 0.016; SLCO2A1 "
     "rs776813259, monomorphic in the cohort) (Table 3).")

h(doc, "Population Pharmacokinetic Modeling", 2)
para(doc,
     "Population pharmacokinetic analysis was performed using nonlinear "
     "mixed-effects modeling with NONMEM. Parameter estimation was conducted "
     "using the first-order conditional estimation method with interaction "
     "(FOCE-I). Perl-speaks-NONMEM and R software were used for model "
     "automation, diagnostics, and visualization.")
para(doc,
     "Infliximab pharmacokinetics were described using an integrated "
     "intravenous and subcutaneous two-compartment model with first-order "
     "absorption and first-order elimination, implemented using ADVAN4 "
     "TRANS4.")
para(doc, "The structural model included:")
bullets(doc, ["clearance (CL)",
              "central volume of distribution (V2)",
              "peripheral volume of distribution (V3)",
              "intercompartmental clearance (Q)",
              "absorption rate constant (Ka)",
              "relative bioavailability (F1)"])
para(doc,
     "Intercompartmental clearance (Q) was fixed at 0.0646 L/day, and relative "
     "bioavailability (F1) was fixed at 0.667.")

h(doc, "Interindividual Variability and Residual Error Model", 2)
para(doc,
     "Interindividual variability (IIV) was modeled using an exponential error "
     "model assuming log-normal parameter distributions, where Pi represents "
     "the individual parameter estimate, TVP represents the typical population "
     "value, and \u03b7i represents a random effect with mean 0 and variance "
     "\u03c9\u00b2. In the final model, interindividual variability was "
     "estimated only for clearance, whereas variability terms for V2, V3, Ka, "
     "and F1 were fixed to zero.")
para(doc,
     "Residual unexplained variability was evaluated using additive, "
     "proportional, and combined residual error models. The final residual "
     "variability model consisted of a proportional error model, while the "
     "additive residual error component was fixed near zero.")

h(doc, "Covariate Analysis", 2)
para(doc,
     "Potential covariates evaluated included ADA status, sex, body weight, "
     "albumin, age, inflammatory markers, renal function markers, and "
     "disease-related variables. Continuous covariates were incorporated using "
     "normalized power models, with body weight and albumin normalized to the "
     "median values of the study population (67.95 kg and 4.5 g/dL, "
     "respectively). Categorical covariates were modeled using proportional "
     "relationships.")

h(doc, "Model Evaluation", 2)
para(doc, "Model performance was evaluated using:")
bullets(doc, ["goodness-of-fit plots",
              "visual predictive checks (VPCs)",
              "shrinkage assessment",
              "parameter precision estimates",
              "nonparametric bootstrap analysis"])
para(doc,
     "Visual predictive checks were based on 200 simulated datasets, with "
     "observations grouped into eight time bins of equal count. Bootstrap "
     "resampling was performed to estimate median parameter values and "
     "5th\u201395th percentile confidence intervals.")

h(doc, "Individual Clearance Estimation for Pharmacogenomic Analysis", 2)
para(doc,
     "Individual clearance (CL) values were derived from the final PopPK model "
     "as empirical Bayes estimates (EBEs) obtained during model estimation. "
     "Because body weight, serum albumin, and ADA status varied over time, "
     "record-level individual CL values were generated by combining each "
     "patient\u2019s EBE of the random effect on clearance with the "
     "time-varying covariate values (CL(t) = TVCL(covariates at t) \u00d7 "
     "exp(\u03b7_EBE)), and phase-specific individual CL was summarized as the "
     "median of the record-level values within each analysis period. The "
     "induction-phase window extended from the first infliximab dose to the "
     "third dose, and the maintenance-phase window from the start of "
     "maintenance dosing to the last dose administered on or before the "
     "one-year assessment (10\u201312 months after the start of maintenance "
     "dosing), or to the end of follow-up when no one-year assessment was "
     "available; the overall treatment period combined both windows. "
     "Time-varying covariates, including ADA status, were taken from the "
     "records within each window, so ADA positivity in the pharmacogenomic "
     "analyses refers to a positive result within the corresponding analysis "
     "period. For one patient without any observed concentration in the "
     "estimation dataset, no EBE could be derived and the patient was excluded "
     "from the CL analyses.")

h(doc, "Pharmacogenomic Association Analysis", 2)
para(doc,
     "Candidate variants in genes previously implicated in anti-TNF "
     "pharmacokinetics or treatment response (TNF, TNFRSF1A, TNFRSF1B, FCGR3A, "
     "TLR2, IL6, IL10, IL17A, CD96, HLA-DQA1, and SLCO2A1) were extracted from "
     "the quality-controlled WGS data. Genotypes were coded as the dosage of "
     "the coded allele (0, 1, or 2).")
para(doc,
     "Because infliximab CL was assumed to follow a log-normal distribution, "
     "the primary CL analysis was performed on the logarithmic scale using "
     "analysis of covariance (ANCOVA), with effect estimates reported as "
     "geometric mean ratios (GMRs) with 95% confidence intervals (CIs). "
     "Consistent with the covariates retained in the final PopPK model, the CL "
     "model was adjusted for sex, body weight, serum albumin, and ADA status "
     "(log CL ~ genotype group + sex + weight + albumin + ADA). ADA positivity "
     "was analyzed using logistic regression adjusted for sex, body weight, "
     "and serum albumin, with associations reported as odds ratios (ORs) with "
     "95% CIs.")
para(doc,
     "For each analysis period, two genetic models were prespecified: a "
     "recessive model comparing homozygotes for the coded allele with all "
     "other genotypes, and a dominant model comparing carriers of the coded "
     "allele with non-carriers. All variants for which both comparison groups "
     "were non-empty were tested; no minimum group-size threshold was applied. "
     "Accordingly, 11 variants were evaluable under the recessive model (three "
     "of the 14 QC-passing variants had no homozygous carriers) and 14 under "
     "the dominant model. Within each analysis stratum defined by analysis "
     "period, endpoint, and genetic model, P values across the tested variants "
     "were adjusted for multiple comparisons using the Benjamini\u2013Hochberg "
     "false discovery rate (FDR) procedure. An FDR-adjusted q value < 0.05 was "
     "considered statistically significant.")
para(doc,
     "Robustness of the association with the smallest P value was examined "
     "using prespecified sensitivity analyses: ANCOVA on the original CL scale "
     "with the same covariate adjustment and FDR procedure, leave-one-out "
     "re-estimation in which each subject of the variant genotype group was "
     "sequentially excluded, a covariate-free Mann\u2013Whitney U test, ANCOVA "
     "with heteroskedasticity-robust (HC3) standard errors, and exclusion of "
     "patients whose individual CL estimates were not informed by observed "
     "infliximab concentration measurements. Residual normality of the ANCOVA "
     "models was assessed using the Shapiro\u2013Wilk test. Analyses were "
     "performed in Python (pandas, statsmodels, SciPy).")

h(doc, "Results", 1)

h(doc, "Patient Characteristics", 2)
para(doc,
     "The final analytic cohort included 139 patients, including 95 patients "
     "(68.4%) with Crohn\u2019s disease and 44 (31.6%) with ulcerative "
     "colitis; 29 patients (20.9%) were female and 17 (12.2%) were pediatric "
     "patients (<19 years old). The mean age was 32.6 years, and the mean body "
     "weight was 66.9 kg. The infliximab PopPK modeling cohort included 97 "
     "patients, whereas the adalimumab cohort included 52 patients.")
para(doc,
     "A total of 567 serum concentration samples were available for PopPK "
     "analysis, including 393 infliximab samples from the PopPK modeling "
     "cohort (69.3%) and 173 adalimumab samples (30.5%). The mean number of "
     "samples per patient was 4.08. Anti-drug antibodies were detected at any "
     "time during follow-up in 10 patients (7.2%) in the analytic cohort and in "
     "9 patients (9.3%) in the infliximab PopPK modeling cohort. Baseline "
     "demographic and clinical characteristics are summarized in Table 1.")

h(doc, "Population Pharmacokinetic Model", 2)
para(doc,
     "An integrated intravenous and subcutaneous two-compartment model "
     "adequately described infliximab concentration-time profiles (Figure 2). "
     "The final model incorporated ADA status, sex, and body weight on "
     "clearance, and albumin and body weight on central volume of "
     "distribution. The final population parameter estimates are summarized in "
     "Table 2.")
para(doc, "Typical population estimates in the final model were as follows:")
bullets(doc, ["CL: 0.295 L/day", "V2: 4.49 L", "V3: 0.407 L",
               "Ka: 0.055 day\u207b\u00b9"])
para(doc,
     "Intercompartmental clearance (Q) and relative bioavailability (F1) were "
     "fixed at 0.0646 L/day and 0.667, respectively. Interindividual "
     "variability was retained only for clearance in the final model.")
h(doc, "Covariate Effects", 2)
bullets(doc, [
    "Body weight was positively associated with both clearance and central "
    "volume of distribution.",
    "Lower albumin levels were associated with increased central volume of "
    "distribution.",
    "ADA positivity was associated with increased infliximab clearance.",
    "Sex showed a modest effect on clearance, although the associated "
    "parameter estimate demonstrated relatively high uncertainty.",
])
para(doc,
     "After covariate inclusion, the estimated interindividual variability for "
     "clearance decreased from 43.5% in the base model to 27.3% in the final "
     "model.")

h(doc, "Model Evaluation and Bootstrap Analysis", 2)
para(doc,
     "Goodness-of-fit plots and visual predictive checks demonstrated "
     "acceptable agreement between observed and model-predicted "
     "concentrations (Supplementary Figures S1 and S2). Bootstrap analysis "
     "demonstrated acceptable parameter "
     "stability and robustness. Median bootstrap parameter estimates were "
     "comparable to the final model estimates, supporting the adequacy of the "
     "final PopPK model.")

h(doc, "Pharmacogenomic Association Analysis", 2)
para(doc,
     "Of the 17 candidate variants, 14 passed variant-level quality control "
     "and were carried into the association analyses (Table 3). A total of 225 "
     "prespecified tests were performed (150 for clearance and 75 for ADA "
     "positivity); 26 ADA tests were not estimable because no ADA-positive "
     "patient was observed within the induction-phase window, leaving 199 "
     "estimable tests. No variant was significantly associated with infliximab "
     "clearance or with ADA positivity after FDR correction in any analysis "
     "period (Table 4, Figure 3, Supplementary Tables S1, S2 and S4).")
para(doc,
     "The smallest adjusted P value for clearance was observed for the FCGR3A "
     "rs396991 variant under the recessive model, in which CC homozygotes "
     "(n = 5) showed higher adjusted clearance than A-allele carriers, with a "
     "consistent direction across the three analysis periods (maintenance "
     "phase: GMR 1.36, 95% CI 1.09\u20131.71, P = 0.009, q = 0.098; overall "
     "treatment period: GMR 1.32, 95% CI 1.06\u20131.64, P = 0.014, q = 0.157; "
     "induction phase: GMR 1.35, 95% CI 1.05\u20131.72, P = 0.020, q = 0.218) "
     "(Table 5, Supplementary Figure S3). This was the only variant with a "
     "nominally significant clearance association in the entire analysis. In "
     "leave-one-out analyses, however, nominal significance was not retained "
     "when any single CC homozygote was excluded (maximum P = 0.112 in the "
     "maintenance phase and 0.141 over the overall treatment period), "
     "indicating dependence on a small number of observations (Supplementary "
     "Table S3). The corresponding dominant-model comparison was not "
     "significant (overall treatment period, AC + CC versus AA: GMR 1.08, "
     "95% CI 0.98\u20131.20, P = 0.134). Given the small number of homozygotes "
     "and the absence of FDR-adjusted significance, this signal should be "
     "regarded as hypothesis-generating only.")
para(doc,
     "No variant was associated with ADA positivity in the maintenance phase "
     "or over the overall treatment period (smallest P = 0.076, smallest "
     "q > 0.99; Supplementary Table S2). Six patients in the pharmacogenomic "
     "analysis cohort were ADA-positive within the analysis periods; three "
     "additional patients first became ADA-positive after the one-year "
     "maintenance window and were therefore classified as ADA-negative in "
     "these analyses (Table 1 reports ADA positivity at any time during "
     "follow-up).")
para(doc,
     "Original-scale ANCOVA yielded consistent conclusions (Supplementary "
     "Table S4). Residual normality was imperfect on both scales but "
     "substantially less deviated on the logarithmic scale "
     "(Shapiro\u2013Wilk P = 0.0006\u20130.006 for log-scale models versus "
     "P < 10\u207b\u2076 for original-scale models), supporting the log-scale "
     "analysis as primary.")

h(doc, "Discussion", 1)
para(doc,
     "In this candidate-gene analysis of an infliximab PopPK cohort, none of "
     "the 14 variants that passed quality control was associated with "
     "individual infliximab clearance or with ADA positivity after correction "
     "for multiple comparisons. This included TNFRSF1B rs1061622, a variant "
     "previously reported in relation to anti-TNF treatment response, which "
     "showed no association with clearance in any analysis period (Table 4, "
     "Figure 3). The variant with the smallest P value, FCGR3A "
     "rs396991, is mechanistically plausible: it encodes the FcγRIIIa "
     "p.Phe158Val substitution, and FcγRIIIa mediates IgG binding on effector "
     "cells, so an effect on the disposition of an IgG1 monoclonal antibody "
     "cannot be excluded a priori. The observed direction was consistent "
     "across the induction phase, the maintenance phase, and the overall "
     "treatment period, and the association persisted in covariate-free and "
     "robust-standard-error analyses. Nevertheless, only five patients were CC "
     "homozygotes, and nominal significance was lost whenever any one of them "
     "was excluded, so the finding is reported as hypothesis-generating rather "
     "than as evidence of an effect.")
para(doc,
     "Two features of the analysis are worth noting. First, individual "
     "clearance was taken from empirical Bayes estimates of the final PopPK "
     "model and combined with time-varying body weight, serum albumin, and ADA "
     "status, so the pharmacokinetic endpoint reflects each patient\u2019s "
     "observed concentrations and the covariate values prevailing in each "
     "treatment period rather than a single static value. Second, "
     "variant-level quality control (MAF and HWE) was applied using genotype "
     "data alone, before any endpoint was examined, and multiplicity was "
     "controlled within each analysis stratum; both choices were prespecified "
     "and reduce the chance that a reported association reflects analytic "
     "flexibility.")
para(doc,
     "The main limitation is statistical power. The pharmacogenomic cohort "
     "comprised 96 patients, homozygote groups under the recessive model "
     "ranged from 1 to 32 patients, and only six patients were ADA-positive "
     "within the analysis periods, with none in the induction-phase window; "
     "the analysis was therefore able "
     "to detect only comparatively large effects, and the ADA endpoint was "
     "essentially uninformative. Interindividual variability in clearance had "
     "already been reduced from 43.5% to 27.3% by the covariates retained in "
     "the PopPK model, leaving a modest amount of unexplained variability for "
     "genetic factors to explain in a cohort of this size. The findings should "
     "therefore be interpreted as an absence of detectable candidate-variant "
     "effects in this cohort rather than as evidence that such effects do not "
     "exist, and the FCGR3A rs396991 signal warrants evaluation in a larger "
     "population.")

ms_docx = f"{ms_dir}/Methods_and_Results_FINAL.docx"
doc.save(ms_docx)
print(f"saved: {os.path.relpath(ms_docx, prj_dir)}")

# ===========================================================================
# 2) Core figures and tables
# ===========================================================================
doc = new_doc()
doc.add_heading("Infliximab PopPK + pharmacogenomics \u2014 core figures and "
                "tables", level=0)
note(doc, "Revised 2026-09-03. Every table is generated from the analysis "
          "output in core_fig_tab/ (code 01\u201309); individual clearance "
          "values are empirical Bayes estimates from the final PopPK model.")

CORE = [
    ("figure", "Figure 1.", "Patient eligibility flow chart.",
     "Anti-TNF inhibitors included infliximab, adalimumab, and ustekinumab. "
     "Patients without available anti-TNF concentration measurements or "
     "without whole-genome sequencing data were excluded to define the "
     "analytic cohort (n = 139). Patients without infliximab administration "
     "records or without an evaluable infliximab concentration for population "
     "pharmacokinetic (PopPK) model estimation were excluded (n = 42), leaving "
     "the infliximab PopPK modeling cohort (n = 97); one further patient was "
     "excluded from the pharmacogenomic analyses because the corresponding "
     "sequencing sample was removed during genotype quality control (n = 96). Analyses were performed separately for the "
     "overall treatment period, the induction phase, and the maintenance "
     "phase. Screening-stage counts (n = X,XXX) to be completed from the EMR "
     "extraction.",
     # 2026-09-07: user-drawn version replaces the matplotlib flow chart
     "Figure1_eligibility_flowchart_revised.png", 6.2),
    ("figure", "Figure 2.", "Structural model diagram.",
     "Integrated intravenous and subcutaneous two-compartment model with "
     "first-order absorption and first-order elimination (ADVAN4 TRANS4) used "
     "to describe infliximab pharmacokinetics.",
     "Figure2_structural_model.png", 5.6),
    ("figure", "Figure 3.",
     "Candidate variants and infliximab clearance over the overall treatment "
     "period.",
     "Adjusted geometric mean ratios (GMRs) with 95% confidence intervals "
     "from analysis of covariance on log-transformed clearance, adjusted for "
     "sex, body weight, serum albumin, and anti-drug antibody status, under "
     "the recessive model (upper panel; homozygotes for the coded allele "
     "versus all other genotypes, 11 variants) and the dominant model (lower "
     "panel; carriers versus non-carriers, 14 variants). Variants are ordered "
     "by P value; group sizes are given in parentheses as variant/reference. "
     "The dashed line marks a GMR of 1. P values were adjusted for multiple "
     "comparisons within each genetic model using the Benjamini\u2013Hochberg "
     "false discovery rate procedure; no association reached the prespecified "
     "threshold of q < 0.05.",
     "Figure3_forest_CL_overall.png", 5.2),
    ("table", "Table 1.",
     "Baseline characteristics of the analytic cohort and the infliximab PopPK modeling cohort.",
     "Categorical variables are presented as n (%) and continuous variables "
     "as mean (standard deviation). Laboratory values are the first available "
     "measurement per patient. Whole phases, patients with induction-phase "
     "data; maintenance only, patients whose data began in the maintenance "
     "phase. ADA, anti-drug antibody; ADA positivity denotes a positive result "
     "at any time during follow-up. TL, trough level.",
     "Table1_baseline_characteristics.csv", 8),
    ("table", "Table 2.",
     "Base and final population pharmacokinetic model parameter estimates.",
     "Values are estimates with relative standard error (%RSE) in parentheses "
     "and shrinkage in square brackets. Bootstrap results are medians with "
     "5th\u201395th percentiles. CL, clearance; Vc, central volume of "
     "distribution; Vp, peripheral volume of distribution; Q, "
     "intercompartmental clearance; F1, relative bioavailability; Ka, "
     "absorption rate constant; IIV, interindividual variability; ADA, "
     "anti-drug antibody.",
     "Table2_popPK_parameters.csv", 8),
    ("table", "Table 3.",
     "Characteristics of the candidate variants in the pharmacogenomic "
     "analysis cohort (n = 96).",
     "Genotype counts are given for each genotype in the coded-allele order. "
     "Variants with a minor allele frequency (MAF) below 0.05 or a "
     "Hardy\u2013Weinberg equilibrium (HWE) exact test P value below 0.05 "
     "were excluded from the association analyses. Variants without "
     "homozygous carriers of the coded allele could not be evaluated under "
     "the recessive model.",
     "Table3_variant_characteristics.csv", 7.5),
    ("table", "Table 4.",
     "Association between candidate variants and infliximab clearance over "
     "the overall treatment period.",
     "Analysis of covariance on log-transformed clearance adjusted for sex, "
     "body weight, serum albumin, and anti-drug antibody status. Effect "
     "estimates are geometric mean ratios (GMRs) with 95% confidence "
     "intervals. P values were adjusted for multiple comparisons within each "
     "genetic model using the Benjamini\u2013Hochberg false discovery rate "
     "procedure (11 variants under the recessive model, 14 under the dominant "
     "model). No association reached the prespecified threshold of q < 0.05.",
     "Table4_CL_association_overall.csv", 7.5),
    ("table", "Table 5.",
     "FCGR3A rs396991 and infliximab clearance across analysis periods "
     "(recessive model).",
     "The variant with the smallest P value in the analysis. Geometric mean "
     "clearance values are shown for CC homozygotes and for patients carrying "
     "the A allele, with adjusted geometric mean ratios from analysis of "
     "covariance on log-transformed clearance. The leave-one-out P range "
     "gives the smallest and largest P values obtained when each CC "
     "homozygote was excluded in turn. Results of the corresponding "
     "original-scale analysis and the Shapiro\u2013Wilk test of residual "
     "normality are shown for comparison.",
     "Table5_rs396991_across_periods.csv", 7),
]

SUPPL = [
    ("figure", "Supplementary Figure S1.",
     "Goodness-of-fit plots for the final infliximab population "
     "pharmacokinetic model.",
     "Observed infliximab concentrations versus individual (top left) and population (top right) predicted concentrations, and conditional weighted residuals (CWRES) versus population predicted concentrations (bottom left) and time after the first dose (bottom right). Black lines are lines of identity or zero, red lines are locally weighted regression trends, and dashed lines mark CWRES of ±2. Concentrations are in µg/mL and time in days.",
     "SupplFigureS1_GOF_model89.png", 5.8),
    ("figure", "Supplementary Figure S2.",
     "Visual predictive check for the final infliximab population "
     "pharmacokinetic model.",
     "Circles are observed infliximab concentrations (µg/mL) versus time after the first dose (days). Lines are the 5th, 50th (red), and 95th percentiles of the observed data; shaded areas are the 95% confidence intervals of the corresponding percentiles obtained from 200 simulated datasets, with observations grouped into eight time bins of equal count.",
     "SupplFigureS2_VPC_model89.png", 6.2),
    ("figure", "Supplementary Figure S3.",
     "Infliximab clearance by FCGR3A rs396991 genotype.",
     "Individual model-estimated clearance by genotype in the maintenance "
     "phase (left) and over the overall treatment period (right). Circles "
     "represent individual patients; horizontal bars and whiskers denote the "
     "geometric mean and its 95% confidence interval. The y-axis is on a "
     "logarithmic scale. Brackets show the adjusted geometric mean ratio "
     "(GMR) with 95% confidence interval, the P value, and the FDR-adjusted q "
     "value for the recessive comparison (CC versus AA + AC), from analysis "
     "of covariance on log-transformed clearance adjusted for sex, body "
     "weight, serum albumin, and anti-drug antibody status.",
     "SupplFigureS3_CL_by_rs396991.png", 6.2),
    ("table", "Supplementary Table S1.",
     "Association between candidate variants and infliximab clearance in the "
     "induction and maintenance phases.",
     "Presented as in Table 4.",
     "SupplTableS1_CL_association_by_period.csv", 7),
    ("table", "Supplementary Table S2.",
     "Association between candidate variants and anti-drug antibody "
     "positivity.",
     "Logistic regression adjusted for sex, body weight, and serum albumin; "
     "odds ratios (ORs) with 95% confidence intervals. Induction-phase "
     "estimates are not available because no ADA-positive patient was "
     "observed within the induction-phase window.",
     "SupplTableS2_ADA_association.csv", 7),
    ("table", "Supplementary Table S3.",
     "Sensitivity analyses for the recessive-model association of FCGR3A "
     "rs396991 with infliximab clearance.",
     "Leave-one-out re-estimation, covariate-free Mann\u2013Whitney U test, "
     "HC3 robust standard errors, and exclusion of patients whose clearance "
     "estimates were not informed by observed concentration measurements.",
     "SupplTableS3_sensitivity_analyses.csv", 7),
    ("table", "Supplementary Table S4.",
     "Association between candidate variants and infliximab clearance on the "
     "original (untransformed) scale.",
     "Presented as in Table 4, with adjusted differences in L/day instead of "
     "geometric mean ratios.",
     "SupplTableS4_CL_association_original_scale.csv", 7),
    ("table", "Supplementary Table S5.",
     "Cohort attrition for the pharmacogenomic analyses.",
     "Per-period analysis populations and the reason for each exclusion.",
     "SupplTableS5_cohort_attrition.csv", 7),
]


def emit(doc, items):
    for kind, label, title, cap, fname, size in items:
        path = f"{cft_dir}/{fname}"
        if not os.path.exists(path):
            note(doc, f"[missing file: {fname}]")
            continue
        caption(doc, label, title + " " + cap)
        if kind == "figure":
            add_figure(doc, path, width_in=size)
        else:
            add_csv_table(doc, path, font_size=size)
        doc.add_paragraph()


h(doc, "Main figures and tables", 1)
emit(doc, CORE)
doc.add_page_break()
h(doc, "Supplementary material", 1)
emit(doc, SUPPL)

cft_docx = f"{cft_dir}/[IFX_POPPK]_core_fig_tab_FINAL.docx"
doc.save(cft_docx)
print(f"saved: {os.path.relpath(cft_docx, prj_dir)}")
