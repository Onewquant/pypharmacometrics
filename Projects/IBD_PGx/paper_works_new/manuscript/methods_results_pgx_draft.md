# Draft — Pharmacogenomic Analysis (Methods & Results sections)

영어 초안. `materials/Methods and results.docx`의 popPK Methods 뒤에 이어붙일
PGx 통계분석 섹션과, Results에 추가할 PGx 결과 단락입니다.
숫자는 `output/Table_pgx_ancova_fdr_results.csv` / `Table_pgx_sensitivity.csv`
기준 (2026-08-20 실행).

---

## [Methods] Pharmacogenomic Association Analysis

Candidate variants in genes previously implicated in anti-TNF
pharmacokinetics or treatment response (TNF, TNFRSF1A, TNFRSF1B, FCGR3A,
TLR2, TLR4, IL6, IL10, IL17A, CD96, HLA-DQA1, SLCO2A1) were extracted from
the quality-controlled WGS data. Genotypes were coded as the dosage of the
minor (coded) allele (0, 1, 2).

Of the 98 patients in the infliximab cohort, one patient was excluded
because the sequencing sample could not be matched to the patient
identifier, leaving 97 patients in the pharmacogenomic analysis cohort.
Phase-specific analysis populations comprised patients with available
phase-specific clearance estimates and complete covariate data
(induction, n = 83; maintenance, n = 92; whole treatment period, n = 96);
for one patient, phase-specific estimates could not be derived.

The association between genotype and infliximab pharmacokinetics was
evaluated using the individual empirical Bayes estimates of clearance (CL)
from the final population pharmacokinetic model and ADA positivity as
endpoints. Because infliximab CL was assumed to follow a log-normal
distribution, CL was analyzed on the logarithmic scale using analysis of
covariance (ANCOVA); results are reported as geometric mean ratios (GMR)
with 95% confidence intervals. ADA positivity was analyzed with logistic
regression and reported as odds ratios. Consistent with the covariates
retained in the final population pharmacokinetic model, the CL model was
adjusted for sex, body weight, serum albumin, and ADA status
(log CL ~ genotype group + sex + weight + albumin + ADA), and the ADA model
was adjusted for sex, body weight, and serum albumin.

Analyses were stratified by treatment phase (induction, maintenance, and
the whole treatment period). Within each phase, two genetic contrasts were
prespecified: a recessive contrast (homozygous variant carriers vs all
others) and a dominant contrast (variant carriers vs non-carriers).
Contrasts were tested only when both groups included at least eight
patients. Within each stratum (phase × endpoint × contrast), p-values
across the tested variants were adjusted for multiple testing using the
Benjamini–Hochberg false discovery rate (FDR) procedure; an FDR-adjusted
q < 0.05 was considered significant.

For significant associations, robustness was assessed by (i)
leave-one-out re-estimation over the variant-group subjects, (ii) the
covariate-free Mann–Whitney U test, (iii) heteroskedasticity-robust
(HC3) standard errors, and (iv) exclusion of patients without any
observed concentration sample, whose empirical Bayes clearance estimates
are shrunk toward covariate-predicted typical values. Residual normality
of the ANCOVA models was assessed with the Shapiro–Wilk test.

## [Results] Pharmacogenomic Association Analysis (draft)

Among the candidate variants tested, the TNFRSF1B rs1061622 GG genotype
(recessive contrast; n = 8 GG vs 84–88 others) was significantly
associated with higher infliximab clearance. On the log scale, GG
homozygotes showed approximately 29–31% higher CL than TT/TG patients
(maintenance phase: GMR 1.29, 95% CI 1.08–1.55, q = 0.048; whole period:
GMR 1.31, 95% CI 1.09–1.57, q = 0.031). The raw-scale ANCOVA gave
consistent results (maintenance: +0.080 L/day, q = 0.020; whole period:
+0.084 L/day, q = 0.010). The induction-phase contrast was not testable
(GG n < 8). Residuals of the log-scale models satisfied normality
(Shapiro–Wilk p = 0.62–0.90) whereas raw-scale residuals did not
(p < 0.01), supporting the log-normal CL assumption.

The association was robust in sensitivity analyses: significance was
retained after excluding any single GG subject (leave-one-out maximum
p = 0.040 on the log scale), in the covariate-free Mann–Whitney test
(p = 0.006–0.013), and with HC3 robust standard errors (p = 0.012–0.026).
All eight GG patients had observed concentration samples; excluding the
five patients whose clearance estimates were not informed by any observed
concentration did not attenuate the association (log-scale p = 0.005–0.009).
No variant was significantly associated with ADA positivity, and no
association reached significance under the dominant (carrier) contrast,
indicating a recessive effect pattern.

### Figure 1 (flow chart) PGx 분기 문구

analytic cohort(139) 정의가 "WGS 시행 환자"이므로, PGx 분기의 제외 사유를
"without WGS"로 쓰면 본문과 모순됨. 아래 문구/수치 사용 (근거:
`output/Table_pgx_attrition.csv`, 06 스크립트로 재현):

> Infliximab cohort (n=98)
> → Excluded: sequencing sample not matched to patient ID (n=1)
> → PGx analysis cohort (n=97)
> → Phase-specific analysis: induction n=83 / maintenance n=92 / whole n=96
>   (phase-specific CL unavailable: maintenance 5, whole 1)

확정된 사유 (2026-08 확인):
- UID 17439372: 시퀀싱 샘플 ↔ 환자 ID 매칭 실패 → 유전형 매트릭스 미포함
- UID 35093356: popPK 모델링에는 포함(NONMEM ID 78, 샘플 2건)되었으나
  maintenance-only 환자로 phase 기준일 산출이 안 되어 phase별 파생
  데이터셋 행이 전부 결측 → phase-specific estimates not derivable

### 해석/한계 (Discussion에 반영할 포인트)

- GG군이 8명으로 작아 exploratory finding으로 프레이밍하고 외부 검증
  필요성을 명시할 것.
- Maintenance와 whole-period 결과는 동일한 GG 8명을 공유하므로 독립적
  재현이 아님을 명시 (whole-period는 보조 결과로 제시).
- TNFR2(TNFRSF1B) M196R(rs1061622)은 anti-TNF 반응성 문헌이 축적된
  변이로, TNFR2-매개 약물 처리/타겟 매개 소실(TMDD) 관점의 기전 논의
  가능.
- Carrier/additive 모델 및 Kruskal-Wallis에서는 비유의 → 열성 패턴
  특이적 신호임을 투명하게 기술.
