# 원고 교체용 섹션 — 최종 확정본 (옵션 A, 2026-08 교수님 결정 반영)

`output/Methods_and_results_reviewed.docx`(가장 최신 원고)를 기준으로,
**수치·서술이 바뀌는 섹션만** 교체용으로 작성한 문서입니다. 여기 없는
섹션(Genotyping/QC, PopPK Modeling, IIV/Residual, Covariate, Model
Evaluation, PopPK Results 등)은 기존 docx 그대로 사용하면 됩니다.

근거 수치: `core_fig_tab/` 및 `output/` (2026-08-23 실행, 03~08 스크립트로 재현).

**[TODO — 원고 확정 전 처리]**
1. Figure 1의 스크리닝 단계 수치(n = X,XXX 3곳)를 EMR 추출 기록으로 채울 것
2. **Table 2 (popPK 파라미터)는 체중 오류 교정 전 추정치** — 재추정 여부
   교수님 결정 대기. 재추정 시 Table 2와 popPK Results 문단의 수치 갱신 필요
3. rs5030728(TLR4)은 Methods의 후보 유전자 나열에 있으나 실제 매트릭스에
   없음 — 후보 목록 문장에서 TLR4 제외 여부 확인

---

## [Methods] Study Population and Data Collection — 1문단 교체

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
analyses, one patient was excluded because the corresponding genotype
data had been removed during quality control, leaving 97 patients in the
pharmacogenomic analysis cohort. Pharmacogenomic analyses were conducted
separately for the overall treatment period (n = 97), the induction phase
(n = 83), and the maintenance phase (n = 97), based on the availability
of phase-specific clearance estimates and complete covariate data
(Figure 1).

## [Methods] Pharmacogenomic Association Analysis — 전체 교체

Candidate variants in genes previously implicated in anti-TNF
pharmacokinetics or treatment response (TNF, TNFRSF1A, TNFRSF1B, FCGR3A,
TLR2, IL6, IL10, IL17A, CD96, HLA-DQA1, and SLCO2A1) were extracted from
the quality-controlled WGS data. Genotypes were coded according to the
dosage of the minor (coded) allele (0, 1, or 2). Variant-level quality
control was applied on genotype data alone, independently of any
endpoint: variants with a minor allele frequency (MAF) below 0.05 or a
Hardy–Weinberg equilibrium (HWE) exact-test P value below 0.05 in the
pharmacogenomic analysis cohort were excluded. Of the 17 candidate
variants, three were excluded (rs3024505, MAF 0.036; rs765249238, MAF
0.016 with HWE P = 0.016; rs776813259, monomorphic), leaving 14 variants
for analysis (Table 3).

Pharmacogenomic association analyses were performed in the 97 patients of
the pharmacogenomic analysis cohort. Individual empirical Bayes estimates
of infliximab clearance (CL) derived from the final PopPK model were used
as the pharmacokinetic endpoint, and anti-drug antibody (ADA) positivity
was evaluated as an additional endpoint. Analyses were conducted
separately for the overall treatment period (n = 97), the induction phase
(n = 83), and the maintenance phase (n = 97), according to the
availability of phase-specific clearance estimates and complete covariate
data.

Because infliximab CL was assumed to follow a log-normal distribution,
the primary CL analysis was performed on the logarithmic scale using
analysis of covariance (ANCOVA), with effect estimates reported as
geometric mean ratios (GMRs) with 95% confidence intervals (CIs).
Consistent with the covariates retained in the final PopPK model, the CL
model was adjusted for sex, body weight, serum albumin, and ADA status
(log CL ~ genotype group + sex + weight + albumin + ADA). ADA positivity
was analyzed using logistic regression adjusted for sex, body weight, and
serum albumin, with associations reported as odds ratios (ORs) with 95%
CIs.

For each analysis period, two genetic models were prespecified: a
recessive model comparing homozygotes for the coded allele with all other
genotypes, and a dominant model comparing carriers of the coded allele
with non-carriers. No minimum genotype-group size was imposed; every
variant with both genotype groups represented was tested (11 variants
under the recessive model and 14 under the dominant model). Within each
analysis stratum defined by analysis period, endpoint, and genetic model,
P values across the tested variants were adjusted for multiple
comparisons using the Benjamini–Hochberg false discovery rate (FDR)
procedure, and an FDR-adjusted q value <0.05 was considered statistically
significant.

The robustness of the lead association was further examined using
sensitivity analyses, including ANCOVA on the original CL scale with the
same covariate adjustment and FDR procedure, leave-one-out re-estimation
across subjects in the variant genotype group, a covariate-free
Mann–Whitney U test, ANCOVA with heteroskedasticity-robust (HC3) standard
errors, and exclusion of patients whose individual CL estimates were not
informed by observed infliximab concentration measurements. Residual
normality of the ANCOVA models was assessed using the Shapiro–Wilk test.

## [Results] Pharmacogenomic Association Analysis — 전체 교체

No candidate variant reached the FDR-adjusted significance threshold for
infliximab CL or ADA positivity in any analysis period (Table 4;
Supplementary Tables S1 and S2).

The strongest and most consistent signal was observed for the TNFRSF1B
rs1061622 GG genotype under the recessive model. GG homozygotes showed
approximately 29–30% higher covariate-adjusted CL than patients with TT
or TG genotypes, with effect estimates that were consistent in direction
and magnitude across all three analysis periods: GMR 1.30 (95% CI
1.08–1.56; P = 0.008, q = 0.084) for the overall treatment period (8 GG
vs 89 TT/TG), GMR 1.29 (95% CI 1.07–1.55; P = 0.009, q = 0.098) for the
maintenance phase (8 vs 89), and GMR 1.30 (95% CI 1.05–1.61; P = 0.017,
q = 0.190) for the induction phase (5 vs 78) (Table 5, Figure 3). In the
overall treatment period, the geometric mean CL was 0.342 L/day (95% CI
0.285–0.410) in GG homozygotes and 0.262 L/day (95% CI 0.248–0.277) in
TT/TG patients.

Analysis on the original CL scale yielded consistent estimates (adjusted
difference +0.079 L/day, q = 0.069 for the overall treatment period;
+0.078 L/day, q = 0.085 for the maintenance phase; Supplementary Table
S4). Residuals from the log-scale ANCOVA models were consistent with
normality (Shapiro–Wilk P = 0.38–0.69), whereas residuals from the
original-scale models deviated significantly from normality (P ≤ 0.016),
supporting the use of log-transformed CL as the primary endpoint.

Although the rs1061622 association did not remain significant after FDR
adjustment, the nominal association was robust in sensitivity analyses of
the overall-treatment and maintenance-phase estimates (Supplementary
Table S3). The association was retained after sequential exclusion of
each GG homozygote in leave-one-out analyses (maximum P = 0.039 and
0.046, respectively), in covariate-free Mann–Whitney U tests (P = 0.008
and 0.015), with HC3 robust standard errors (P = 0.022 and 0.026), and
after excluding the five patients whose individual CL estimates were not
informed by observed concentration measurements (P = 0.010 and 0.012).
All eight GG homozygotes had observed infliximab concentration
measurements.

No variant was associated with ADA positivity after FDR correction
(Supplementary Table S2), and no association with infliximab CL was
observed under the dominant model (smallest q = 0.16, for IL6 rs10499563
in the overall treatment period; Table 4), indicating that the rs1061622
signal was restricted to the recessive model. Given the small number of
GG homozygotes and the multiplicity burden of the candidate-variant
screen, the rs1061622 finding should be regarded as exploratory and
hypothesis-generating, requiring confirmation in an independent cohort.

## [Discussion에 반영할 포인트] (초안 아님 — 집필 가이드)

- **주 메시지 프레임**: "TNFR2 M196R(rs1061622) 열성 유전형에서 infliximab
  청