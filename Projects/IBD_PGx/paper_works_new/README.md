# paper_works_new — Infliximab popPK + Pharmacogenomic 논문 작업 폴더

논문 골자: **IBD 환자에서 infliximab 집단약동학 모델링 + 후보 유전변이의
PK(CL)/ADA 연관성 분석** (adalimumab 제외 확정, 2026-07 교수님 코멘트 반영).

기존 `paper_works/` 자료를 재정리하고, PGx 분석을 최종 프레임
(infliximab only, phase 층화, phase×endpoint×contrast별 FDR 보정)으로
재작성한 폴더입니다.

## 폴더 구조

```
paper_works_new/
├── code/         분석·표·그림 생성 스크립트 (번호 순서대로 실행)
├── data/         분석 입력 스냅샷 (gene_pd_cor에서 복사)
├── materials/    기존 paper_works에서 가져온 원고·그림·참고자료
├── output/       스크립트가 생성하는 표·그림 (재실행으로 재현 가능)
└── manuscript/   원고 초안 텍스트
```

## 실행 방법

프로젝트 루트 venv 사용. 실행 순서 의존성: 04, 05는 03 이후에 실행.

```
C:/Users/ilma0/PycharmProjects/pypharmacometrics/venv/Scripts/python.exe -X utf8 <script>
```

## code/ ↔ 논문 요소 매핑

| 스크립트 | 산출물 (output/) | 논문 요소 |
|---|---|---|
| `01_table1_demographics.py` | `Table1_demographics.csv` | Table 1 — Baseline characteristics (analytic / infliximab cohort). 기존 `demotable_paper.py`의 독립 실행형 재작성, adalimumab 컬럼 제거. 기존 원고 수치(139/98명)와 일치 확인됨. |
| `02_table_genotype_summary.py` | `Table_genotype_summary.csv` | 후보 변이 요약표 — 유전형 분포(0/1/2), coded allele freq, MAF, HWE exact p (전체 genotyped 138명 / infliximab cohort 97명). |
| `03_pgx_ancova_fdr.py` | `Table_pgx_ancova_fdr_results.csv` | **PGx 본분석.** infliximab only, phase(IND/MAINT/OVERALL) 층화, endpoint(CL/ADA), CL은 raw(BETA)+log(GMR) 두 스케일, 비교는 HOM_vs_OTHERS / CARRIER_vs_NONCARRIER. FDR은 (phase×endpoint×스케일×비교)별로 변이들에 대해 BH 보정. 최소 그룹 수 컷오프 없음. 잔차 Shapiro/skew 컬럼 포함. |
| `04_pgx_sensitivity.py` | `Table_pgx_sensitivity.csv` | 유의 결과의 견고성 — leave-one-out(min/max p), Mann-Whitney, HC3 robust p, 관찰 농도 샘플 없는 환자(5명) 제외 재분석. |
| `05_figure_cl_by_genotype.py` | `Figure_CL_by_rs1061622.png/.pdf` | PGx Figure — 유전형(TT/TG/GG)별 개인 CL 산점도 + 기하평균(95% CI) + GMR/q 주석, 2패널(MAINT/Whole). |
| `06_pgx_cohort_attrition.py` | `Table_pgx_attrition.csv` | Figure 1 PGx 분기 수치 산출 — 98 → 97(유전형 매트릭스 탈락 1명) → phase별 83/92/96, 특이 UID 목록 포함. |

## 핵심 결과 (2026-08-22 실행 기준, 최소 그룹 수 컷오프 제거 후)

**rs1061622 (TNFRSF1B, TNFR2 M196R) GG homozygote → infliximab CL 상승**

| Phase | log(CL) GMR (95% CI) | p | q | raw q | GG n |
|---|---|---|---|---|---|
| IND | 1.30 (1.05–1.61) | 0.017 | 0.207 | 0.144 | 5 |
| MAINT | 1.29 (1.08–1.55) | 0.007 | 0.083 | **0.034** | 8 |
| OVERALL | 1.31 (1.09–1.57) | 0.004 | 0.052 | **0.018** | 8 |

- 세 phase 모두 GMR 1.29–1.31로 방향·크기 일관
- log 모델 잔차 정규성 만족(Shapiro p 0.62–0.90), raw는 기각 → **log 모델이 주 분석**
- 주 분석(log) q는 borderline (OVERALL 0.052), raw 스케일은 q<0.05
- 민감도: LOO 전 케이스 p<0.05, Mann-Whitney p 0.006–0.013, HC3 p 0.012–0.026,
  무샘플 5명 제외 후에도 유지 (GG 8명 전원 관찰 농도 보유)
- ADA endpoint 및 carrier 비교는 전부 비유의 (열성 패턴 특이적)
- 주의: MAINT와 OVERALL은 같은 GG 8명 공유 → 독립 재현 아님

### 컷오프 관련 이력 (중요)

기존 코드는 그룹당 최소 8명 필터를 썼으나, 근거가 임의적이고 핵심 결과의
GG군이 정확히 8명이라 post-hoc 비판 소지가 있어 **2026-08 컷오프를
제거**함 (`MIN_GROUP_N = 1`). 영향:

- FDR family: HOM 비교 7 → 12개, CARRIER 비교 14 → 16개로 증가
- 비보정 p는 불변, **q값만 상승** (log MAINT 0.048→0.083, OVERALL 0.031→0.052)
- IND phase가 검정 가능해짐 (GG n=5)

컷오프를 되돌리려면 `03_pgx_ancova_fdr.py`의 `MIN_GROUP_N`만 수정.

## PGx cohort attrition (Figure 1 분기, 06 스크립트로 재현)

```
Infliximab cohort 98
 └─ 제외 1명 (UID 17439372): 유전형 QC 단계에서 제거 (2026-08-23 확정)
    → Figure 1 표기: "genotype data removed during quality control (n=1)"
      ("without WGS"로 쓰면 본문의 analytic cohort 정의와 모순됨)
PGx cohort 97
 ├─ IND     83 (15명은 induction phase 데이터 없음)
 ├─ MAINT   92 (5명 phase별 CL 없음: 18898880, 35093356, 37291334, 37366865, 38241008)
 └─ OVERALL 96 (35093356만 CL 전무)
공변량(SEX/WT/ALB/ADA) 결측 탈락: 0명
```

### 특이 UID 사유 (2026-08 확인 완료)

1. **UID 17439372**: **유전형 QC 단계에서 제거된 것으로 기록** (방침 확정).
   배정 샘플 `23-B02281_EB-01`이 최종 QC VCF에 없음. ID 표기 차이 문제가
   **아님** — tagged VCF(195샘플) 대상 완전일치/정규화/숫자/부분일치 4단계
   탐색 모두 음성. 미매칭 57샘플은 pid_df와 교집합 0건인 별개 코호트 샘플.
   - 검증된 범위: 본 분석 코드에서는 샘플이 제거되지 않음 (이전 QC 코드의
     `filter_cols`는 전부 주석 처리, 현재 스크립트의 샘플 QC는 `mt_pca`
     분기에만 적용되어 dosage matrix 경로인 `mt`은 무영향)
   - **미확인 사항**: 시퀀싱 미실시인지 센터 QC 탈락인지는 대조하지 않음.
     센터 확인 시 문구를 구체화할 수 있으나 n=97은 불변.
   - 검증 코드: `genomics_sample_id_audit_hail.py`(Hail),
     `genomics_mt_vs_vcf_check.py`(mt↔VCF 헤더 대조),
     `genomics_sample_uid_reconcile.py`(UID 병합/감사),
     `check_vcf_header_samples.py`(파일 단독 판정)
2. **UID 35093356**: popPK 모델링에는 포함(NONMEM ID 78, 투여 23건·샘플
   2건)되었으나 maintenance-only 환자로 phase 기준일 산출 불가 →
   phase별 파생 데이터셋(for_genomics_df) 행 전부 결측. 논문 표기:
   "phase-specific estimates could not be derived (n=1)".

## materials/ (기존 paper_works에서 복사)

- `Methods and results.docx` — popPK Methods/Results 기존 초안
- `[IFX_POPPK]_core_fig_tab.docx` — Table 1, Table 2(파라미터 추정치) 등 핵심 표
- `Figure1_eligibility flow chart.png`, `Figure2_structural model diagram.png`
- `bootstrap_summary_corrected.xlsx`, `bootstrap_results - 복사본.csv` — 부트스트랩
- `rsid_genotype_summary.xlsx` — 기존 유전형 요약 (02 스크립트가 CSV로 재현)
- `QC python코드.txt`, `기본적용된 QC 확인...txt` — WGS QC (Hail) 기록

## manuscript/

- `methods_results_pgx_draft.md` — PGx 통계분석 Methods + Results 영어 초안,
  Discussion 반영 포인트 포함. `Methods and results.docx`에 병합해서 사용.

## 원본 데이터/코드 위치 (이 폴더 밖)

- popPK 모델링 코드: `../modeling_codes_infliximab/`
- NONMEM 데이터셋: `C:/Users/ilma0/NONMEMProjects/IBDPGX/`
- PGx 분석 원본(개발 이력 포함): `../gene_pd_cor/`
- 주의: `../gene_pd_cor/` 및 기존 스크립트들의 `rsid_gene_dict`에는
  rs1061622가 TNF로 잘못 주석되어 있음 → 실제는 **TNFRSF1B** (본 폴더
  코드에는 수정 반영됨). rs767455도 TNFRSF1A로 수정함 (dbSNP 확인 권장).
