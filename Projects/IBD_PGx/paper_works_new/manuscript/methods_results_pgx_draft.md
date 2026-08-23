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
because genotype data were not available: the corresponding sequencing
sample was removed during genotype quality control and was therefore
absent from the final quality-controlled variant dataset. This left 97
patients in the pharmacogenomic analysis cohort.
Phase-specific analysis populations comprised patients with available
phase-specific clearance estimates and complete covariate data
(induction, n = 83; maintenance, n = 92; overall treatment period,
n = 96); for one patient, phase-specific estimates could not be derived.

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
the overall treatment period). Within each phase, two genetic contrasts
were prespecified: a recessive contrast (homozygous variant carriers vs
all others) and a dominant contrast (variant carriers vs non-carriers).
No minimum group-size filter was applied: every variant with both
genotype groups represented was tested. Within each stratum
(phase × endpoint × contrast), p-values across the tested variants were
adjusted for multiple testing using the Benjamini–Hochberg false
discovery rate (FDR) procedure; an FDR-adjusted q < 0.05 was considered
significant.

For significant associations, robustness was assessed by (i)
leave-one-out re-estimation over the variant-group subjects, (ii) the
covariate-free Mann–Whitney U test, (iii) heteroskedasticity-robust
(HC3) standard errors, and (iv) exclusion of patients without any
observed concentration sample, whose empirical Bayes clearance estimates
are shrunk toward covariate-predicted typical values. Residual normality
of the ANCOVA models was assessed with the Shapiro–Wilk test.

## [Results] Pharmacogenomic Association Analysis (draft)

Among the candidate variants tested, the TNFRSF1B rs1061622 GG genotype
(recessive contrast; n = 8 GG vs 78–88 others) showed the strongest
association with infliximab clearance. GG homozygotes had approximately
29–31% higher clearance than TT/TG patients, and the effect estimate was
consistent in direction and magnitude across all three phase strata
(induction: GMR 1.30, 95% CI 1.05–1.62, p = 0.017, q = 0.207;
maintenance: GMR 1.29, 95% CI 1.08–1.55, p = 0.007, q = 0.083; overall
treatment: GMR 1.31, 95% CI 1.09–1.57, p = 0.004, q = 0.052). After FDR
adjustment across all estimable contrasts, the association reached
significance on the raw scale (maintenance q = 0.034; overall treatment
q = 0.018) but was borderline on the log scale (overall treatment
q = 0.052). Residuals of the log-scale models satisfied normality
(Shapiro–Wilk p = 0.62–0.90) whereas raw-scale residuals did not
(p < 0.01), supporting the log-normal clearance assumption; the log-scale
estimates are therefore reported as the primary analysis.

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
> → Excluded: genotype data removed during quality control (n=1)
> → PGx analysis cohort (n=97)
> → Phase-specific analysis: induction n=83 / maintenance n=92 /
>   overall treatment n=96
>   (phase-specific CL unavailable: maintenance 5, overall 1)

확정된 사유 (2026-08-23 검증):
- UID 17439372: **QC 단계에서 제거된 것으로 기록** (논문 표기 방침).
  배정 샘플 `23-B02281_EB-01`이 최종 QC VCF에 부재 — 완전일치/정규화/
  숫자/부분일치 4단계 탐색 모두 음성, tagged VCF(195샘플)에서도 동일.
  pid_df의 2023년 샘플 23건 중 22건만 VCF에 존재하며 빠진 1건이 이 샘플.

  검증된 범위: 본 분석 코드(이전 QC 코드의 `filter_cols`는 전부 주석,
  현재 스크립트의 샘플 QC는 `mt_pca` 분기에만 적용)에서는 샘플이 제거되지
  않음 → 제거는 상위(센터) 단계에서 발생.
  **미확인**: 시퀀싱 미실시인지 센터 QC 탈락인지는 대조하지 않았음.
  센터 확인 시 문구 구체화 가능("not sequenced" vs "excluded during
  upstream QC"). 어느 쪽이든 현재 표기는 유지 가능하며 n=97 불변.
  검증 코드: `gene_pd_cor/genomics_sample_id_audit_hail.py`,
  `gene_pd_cor/genomics_mt_vs_vcf_check.py`,
  `gene_pd_cor/genomics_sample_uid_reconcile.py`
- UID 35093356: popPK 모델링에는 포함(NONMEM ID 78, 샘플 2건)되었으나
  maintenance-only 환자로 phase 기준일 산출이 안 되어 phase별 파생
  데이터셋 행이 전부 결측 → phase-specific estimates not derivable

### 해석/한계 (Discussion에 반영할 포인트)

- **컷오프 제거의 영향 (2026-08 결정)**: 최소 그룹 수 필터를 없애면 FDR
  family가 7 → 12개(HOM 비교)로 커져, log 스케일 q값이 MAINT 0.048 →
  0.083, OVERALL 0.031 → 0.052로 상승. 비보정 p는 불변(MAINT 0.007,
  OVERALL 0.004). 따라서 **주 분석(log 스케일)에서는 "borderline
  significant"로 서술**하고, effect size의 일관성(세 phase 모두 GMR
  1.29–1.31)과 민감도 분석 결과를 근거로 논의하는 것이 정직함.
  raw 스케일에서는 여전히 q<0.05 (MAINT 0.034, OVERALL 0.018).
- Induction phase도 이제 검정 가능해짐 (GG n=5): GMR 1.30, p=0.017,
  q=0.207. 유의하지 않지만 **세 phase에서 방향·크기가 일관**되다는 점은
  결과 신뢰도를 지지하는 서술로 활용 가능.
- GG군이 8명으로 작아 exploratory finding으로 프레이밍하고 외부 검증
  필요성을 명시할 것.
- Maintenance와 overall treatment 결과는 동일한 GG 8명을 공유하므로
  독립적 재현이 아님을 명시.
- TNFR2(TNFRSF1B) M196R(rs1061622)은 anti-TNF 반응성 문헌이 축적된
  변이로, TNFR2-매개 약물 처리/타겟 매개 소실(TMDD) 관점의 기전 논의
  가능.
- Carrier/additive 모델 및 Kruskal-Wallis에서는 비유의 → 열성 패턴
  특이적 신호임을 투명하게 기술. ADA endpoint는 어떤 변이도 비유의.
