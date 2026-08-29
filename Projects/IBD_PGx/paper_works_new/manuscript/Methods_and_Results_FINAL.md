# Methods and Results — FINAL (option A: exploratory reporting)

> ## ⚠️ 중요: 아래 PGx 수치는 재산출 예정 (2026-08-30 확인)
>
> 이 문서의 PGx 수치는 `sim90`에서 나온 개인 CL을 기반으로 합니다.
> 그런데 `sim90`(90.mod)은 `$SIM ... ONLYSIM` + `$OMEGA 0.0711`로
> **개인 ETA를 난수 생성**한 시뮬레이션이어서, 그 CL은 환자의 실제
> 개인 추정치(EBE)가 아닙니다 (sim ETA vs 89.phi EBE 상관 −0.08).
>
> 실제 EBE로 확인하면 rs1061622 GG vs TT+TG는 **GMR ≈ 1.00, p = 0.99**로
> 연관성이 없습니다. 즉 아래의 GMR 1.29–1.30은 난수 배정에서 비롯된
> 값일 가능성이 높습니다.
>
> **다음 세션 예정 작업**: `95.mod`(run 89 EBE 주입)로 시뮬레이션 재실행
> → 파생 데이터·표·그림 전체 재생성 → 이 문서 수치 갱신.
> 코호트도 97 → **96명**으로 바뀝니다 (UID 38339532는 유지기 시작
> + 관찰농도 1건이라 추정모델에서 제외되어 EBE 없음).
> 절차: `../RERUN_STEPS.md`
>
> **따라서 이 문서는 아직 원고에 반영하지 마십시오.**
> 구조·서술 방식(옵션 A 프레이밍)만 참고용으로 유효합니다.

Reporting frame confirmed 2026-08: the TNFRSF1B rs1061622 association is
reported as an **exploratory finding**. The effect estimate is consistent
in direction and magnitude across all three analysis periods but does not
reach the FDR-adjusted significance threshold.

All numbers below are reproduced by `code/01`–`08`; sources are noted in
square brackets. Replaces the PGx sections of
`output/Methods_and_results_reviewed.docx`.

---

## [Methods]

### Study Population and Data Collection

This retrospective population pharmacokinetic (PopPK) study included
patients with inflammatory bowel disease (IBD) who received anti-tumor
necrosis factor (anti-TNF) inhibitor therapy for at least one day at a
tertiary referral hospital. Patients without available anti-TNF inhibitor
concentration measurements were excluded from the study population, and
those without whole-genome sequencing (WGS) data were further excluded to
define the overall analytic cohort (n = 139). From this cohort, patients
with infliximab administration records were selected to form the
infliximab PopPK modeling cohort (n = 98), which was used for development
of the infliximab PopPK model. For the subsequent pharmacogenomic
analyses, one patient was excluded because the corresponding sequencing
sample was removed during genotype quality control and was therefore
absent from the final quality-controlled variant dataset, leaving 97
patients in the pharmacogenomic analysis cohort. Pharmacogenomic analyses
were conducted separately for the overall treatment period (n = 97), the
induction phase (n = 83), and the maintenance phase (n = 97), based on
the availability of the corresponding clearance estimates and complete
covariate data (Figure 1). [Table_pgx_attrition.csv]

Demographic characteristics, laboratory test results, therapeutic drug
monitoring (TDM) measurements, anti-drug antibody (ADA) status, dosing
records, and treatment-related variables were extracted from the
electronic medical record system.

### Genotyping, Quality Control, and Variant Processing

*(unchanged from the reviewed draft up to the end of variant-level QC;
append the paragraph below)*

For the candidate-variant association analyses, an additional
variant-level quality control step was applied within the pharmacogenomic
analysis cohort using genotype data only, without reference to any
endpoint. Variants with a minor allele frequency (MAF) below 0.05 or a
Hardy–Weinberg equilibrium (HWE) exact test P value below 0.05 were
excluded. Of the 17 candidate variants, 14 passed this filter; three were
excluded (IL10 rs3024505, MAF 0.036; SLCO2A1 rs765249238, MAF 0.016 and
HWE P = 0.016; SLCO2A1 rs776813259, monomorphic in the cohort)
(Table 3). [Table_variant_qc.csv]

### Population Pharmacokinetic Modeling

*(unchanged from the reviewed draft)*

### Pharmacogenomic Association Analysis

Candidate variants in genes previously implicated in anti-TNF
pharmacokinetics or treatment response (TNF, TNFRSF1A, TNFRSF1B, FCGR3A,
TLR2, IL6, IL10, IL17A, CD96, HLA-DQA1, and SLCO2A1) were extracted
from the quality-controlled WGS data. Genotypes were coded as the dosage
of the coded allele (0, 1, or 2).

Pharmacogenomic association analyses were performed in the 97 patients
included in the pharmacogenomic analysis cohort. Individual empirical
Bayes estimates of infliximab clearance (CL) derived from the final PopPK
model were used as the pharmacokinetic endpoint, and ADA positivity was
evaluated as an additional endpoint. Separate analyses were conducted for
the overall treatment period (n = 97), the induction phase (n = 83), and
the maintenance phase (n = 97), according to the availability of the
corresponding clearance estimates and complete covariate data.

Because infliximab CL was assumed to follow a log-normal distribution,
the primary CL analysis was performed on the logarithmic scale using
analysis of covariance (ANCOVA), with effect estimates reported as
geometric mean ratios (GMRs) with 95% confidence intervals (CIs).
Consistent with the covariates retained in the final PopPK model, the CL
model was adjusted for sex, body weight, serum albumin, and ADA status
(log CL ~ genotype group + sex + weight + albumin + ADA). ADA positivity
was analyzed using logistic regression adjusted for sex, body weight, and
serum albumin, with associations reported as odds ratios (ORs) with
95% CIs.

For each analysis period, two genetic models were prespecified: a
recessive model comparing homozygotes for the coded allele with all other
genotypes, and a dominant model comparing carriers of the coded allele
with non-carriers. All variants for which both comparison groups were
non-empty were tested; no minimum group-size threshold was applied.
Accordingly, 11 variants were evaluable under the recessive model (three
of the 14 QC-passing variants had no homozygous carriers) and 14 under
the dominant model. Within each analysis stratum defined by analysis
period, endpoint, and genetic model, P values across the tested variants
were adjusted for multiple comparisons using the Benjamini–Hochberg false
discovery rate (FDR) procedure. An FDR-adjusted q value < 0.05 was
considered statistically significant.

Robustness of the leading CL association was evaluated using prespecified
sensitivity analyses: ANCOVA on the original CL scale using the same
covariate adjustment and FDR procedure; leave-one-out re-estimation with
sequential exclusion of each subject in the variant genotype group; a
covariate-free Mann–Whitney U test; ANCOVA with heteroskedasticity-robust
(HC3) standard errors; and exclusion of patients whose individual CL
estimates were not informed by observed infliximab concentration
measurements. Residual normality of the ANCOVA models was assessed using
the Shapiro–Wilk test.

---

## [Results]

### Patient Characteristics

*(unchanged from the reviewed draft — Table 1)*

### Population Pharmacokinetic Model / Covariate Effects / Model Evaluation

*(unchanged from the reviewed draft — Table 2, Figure 2)*

### Pharmacogenomic Association Analysis

Of the 17 candidate variants, 14 passed variant-level quality control and
were carried into the association analyses (Table 3). No variant reached
the FDR-adjusted significance threshold for either endpoint in any
analysis period (Table 4, Supplementary Tables S1–S2).

The strongest signal was observed for the TNFRSF1B rs1061622 variant
under the recessive model, in which GG homozygotes had higher adjusted
infliximab CL than patients with TT or TG genotypes. The magnitude and
direction of this association were consistent across all three analysis
periods: overall treatment (n = 8 GG vs 89 TT/TG; GMR 1.30, 95% CI
1.08–1.56; P = 0.008; q = 0.084), maintenance phase (n = 8 vs 89; GMR
1.29, 95% CI 1.07–1.55; P = 0.009; q = 0.098), and induction phase
(n = 5 vs 78; GMR 1.30, 95% CI 1.05–1.61; P = 0.017; q = 0.190)
(Table 5, Figure 3). Corresponding geometric mean CL values in the
overall treatment period were 0.342 L/day (95% CI 0.285–0.410) in GG
homozygotes and 0.262 L/day (95% CI 0.248–0.277) in TT/TG patients.
After FDR adjustment across the 11 variants tested under the recessive
model, none of these associations met the prespecified significance
threshold.

ANCOVA performed on the original CL scale yielded consistent estimates
(overall treatment, adjusted difference +0.079 L/day, 95% CI
0.024–0.135, q = 0.069; maintenance phase, +0.078 L/day, 95% CI
0.022–0.134, q = 0.085) (Supplementary Table S4). Residuals from the
log-scale ANCOVA models were consistent with normality (Shapiro–Wilk
P = 0.38–0.69), whereas residuals from the original-scale models showed
significant deviations from normality (P < 0.05 in all periods),
supporting the use of log-transformed CL as the primary endpoint.

In sensitivity analyses of the rs1061622 association, the nominal
association persisted under several alternative specifications
(Supplementary Table S3). In the overall treatment and maintenance
analyses, nominal significance was retained after sequential exclusion of
each GG homozygote (leave-one-out P range 0.002–0.039 and 0.002–0.046,
respectively), in covariate-free Mann–Whitney U tests (P = 0.008 and
0.015), using HC3 robust standard errors (P = 0.022 and 0.026), and
after excluding the five patients whose individual CL estimates were not
informed by observed infliximab concentration measurements (P = 0.010 and
0.012); all eight GG homozygotes had observed concentration measurements.
In the induction-phase analysis, which included only five GG homozygotes,
the leave-one-out P range extended above 0.05 (0.002–0.084), indicating
greater sensitivity to individual observations.

No variant was associated with ADA positivity after FDR correction in any
analysis period (Supplementary Table S2). Under the dominant model, the
smallest adjusted P value for CL was observed for IL6 rs10499563 in the
overall treatment period (GMR 0.86, 95% CI 0.77–0.97; P = 0.012;
q = 0.164), which likewise did not meet the significance threshold
(Table 4).

---

## Notes for the Discussion

- Frame the rs1061622 signal as **exploratory and hypothesis-generating**:
  consistent effect size (GMR 1.29–1.30) across three analysis periods
  and robust to several sensitivity analyses, but not significant after
  multiplicity correction.
- **Power is the principal limitation**: only eight GG homozygotes were
  available. With this group size the analysis had approximately 80%
  power to detect a GMR of about 1.29, i.e. the study was powered only
  for an effect at least as large as the one observed.
- The maintenance-phase and overall-treatment analyses **share the same
  eight GG homozygotes** and are therefore not independent replications;
  the induction-phase estimate (n = 5) is the only partially independent
  comparison.
- The signal was **specific to the recessive model**; no association was
  seen under the dominant model, and ADA positivity showed no association
  with any variant.
- TNFRSF1B rs1061622 (TNFR2 M196R) has prior literature in anti-TNF
  response; a TNFR2-mediated target-mediated drug disposition mechanism
  is a plausible interpretation but cannot be established from these data.
- **External validation in an independent cohort with a larger number of
  GG homozygotes is required.**

## Open items to resolve before submission

1. **PopPK model re-estimation.** The final model parameters in Table 2
   were estimated before the body-weight data correction (one record,
   834.0 kg, a decimal-point error for 83.4 kg). The simulation used for
   the individual CL estimates has been re-run with corrected data, but
   the model itself has not been re-estimated. Table 2 should be
   regenerated after re-estimation, and the PGx analyses re-run on the
   resulting CL estimates.
2. **Inclusion rule for sparsely sampled patients.** One patient with
   only two observed concentrations has the highest CL in the cohort and
   contributes to the reference group. A prespecified rule for the
   minimum number of observed concentrations required for a patient's CL
   estimate to enter the association analysis should be agreed on
   methodological grounds.
3. **Screening-stage counts in Figure 1** (total treated, exclusions 1–2)
   are placeholders and need the EMR extraction numbers.
4. Confirm with the sequencing center whether sample `23-B02281_EB-01`
   was sequenced, to refine the Exclusion 4 wording if needed.
5. **TLR4 removed from the candidate gene list.** The reviewed draft listed
   TLR4, but no TLR4 variant (rs5030728) is present in the genotype matrix;
   the 17 candidate variants span 11 genes without TLR4. Verified 2026-08.
