# paper_works_new — Infliximab popPK + Pharmacogenomic 논문 작업 폴더

논문 골자: **IBD 환자에서 infliximab 집단약동학 모델링 + 후보 유전변이의
PK(CL)/ADA 연관성 분석**

**보고 방향: (A) exploratory finding 확정 (2026-08)**
rs1061622(TNFRSF1B) 연관성은 세 분석기간에서 효과크기가 일관되나
다중검정 보정 후 유의수준에 도달하지 못함 → 탐색적 발견으로 서술.

## 바로 쓸 수 있는 것 (논문 작성용)

| 항목 | 위치 |
|---|---|
| **원고 본문 (Methods/Results 최종본)** | `manuscript/Methods_and_Results_FINAL.md` |
| **Figure/Table 캡션 초안** | `core_fig_tab/CAPTIONS.md` |
| **본문 표·그림 파일** | `core_fig_tab/Table1~5, Figure1~3` |
| **보충자료 표** | `core_fig_tab/SupplTableS1~S5` |

`manuscript/Methods_and_Results_FINAL.md`가 기준 문서입니다. 기존
`output/Methods_and_results_reviewed.docx`는 **이전 수치(PGx 96명,
MAINT 92명, q=0.048 유의)** 기준이므로 PGx 관련 문단을 위 파일 내용으로
교체해야 합니다. popPK 부분은 기존 docx 그대로 사용 가능합니다.

## 폴더 구조

```
paper_works_new/
├── code/          분석·표·그림 생성 스크립트 (01~08, 번호순 실행)
├── core_fig_tab/  ★ 논문용 최종 표·그림 + 캡션
├── manuscript/    ★ 원고 본문 최종본, 교수님 보고 메일 기록
├── data/          분석 입력 스냅샷
├── materials/     기존 paper_works에서 가져온 원본 자료
└── output/        중간 산출물 (raw 결과 CSV, 사용자 개정 docx 등)
```

## 실행 방법

프로젝트 루트 venv 사용. 01~06 실행 후 07(표), 08(Figure 1) 실행.

```
C:/Users/ilma0/PycharmProjects/pypharmacometrics/venv/Scripts/python.exe -X utf8 code/<script>
```

## code/ ↔ 논문 요소 매핑

| 스크립트 | 산출물 | 논문 요소 |
|---|---|---|
| `01_table1_demographics.py` | `Table1_demographics.csv` | Table 1 (원고 수치와 일치 확인됨) |
| `02_table_genotype_summary.py` | `Table_genotype_summary.csv` | Table 3 재료 (유전형 분포/MAF/HWE) |
| `03_pgx_ancova_fdr.py` | `Table_pgx_ancova_fdr_results.csv`, `Table_variant_qc.csv` | **PGx 본분석** + 변이 QC |
| `04_pgx_sensitivity.py` | `Table_pgx_sensitivity.csv` | Suppl S3 (견고성) |
| `05_figure_cl_by_genotype.py` | `Figure_CL_by_rs1061622.png/pdf` | **Figure 3** |
| `06_pgx_cohort_attrition.py` | `Table_pgx_attrition.csv` | Suppl S5, Figure 1 수치 |
| `07_core_tables.py` | `core_fig_tab/Table1~5, SupplS1~S5` | **논문용 표 전체** |
| `08_figure1_flowchart.py` | `core_fig_tab/Figure1_*.png/pdf` | **Figure 1** |

## 확정된 분석 프레임

- infliximab만 (adalimumab 제외), phase 층화: **IND / MAINT / OVERALL**
- 모델: `ADA ~ GROUP+SEX+WEIGHT+ALBUMIN` (logistic),
  `CL ~ GROUP+SEX+WEIGHT+ALBUMIN+ADA` (일반 OLS ANCOVA)
- **log(CL)이 주 분석** (잔차 정규성 log만 통과), raw는 보조
- 최소 그룹 수 컷오프 **없음** (`MIN_GROUP_N=1`)
- 변이 QC 사전 필터: **MAF ≥ 0.05, HWE exact p ≥ 0.05**
  → 3개 제외 (rs3024505, rs765249238, rs776813259), 14개 통과
  → 열성 모델 검정 가능 11개(동형접합 0명인 3개 제외), 우성 모델 14개
- FDR(BH)은 (phase × endpoint × CL스케일 × 유전모델) 층 안에서 변이들에 대해 보정

## 핵심 결과

**rs1061622 (TNFRSF1B, TNFR2 M196R) GG vs TT+TG, log(CL), family m=11**

| Phase | GMR (95% CI) | p | q | GG n |
|---|---|---|---|---|
| Overall | 1.30 (1.08–1.56) | 0.008 | 0.084 | 8 |
| Maintenance | 1.29 (1.07–1.55) | 0.009 | 0.098 | 8 |
| Induction | 1.30 (1.05–1.61) | 0.017 | 0.190 | 5 |

- **FDR 유의(q<0.05) 0건** → exploratory finding으로 보고
- 세 기간 모두 GMR 1.29–1.30으로 방향·크기 일관
- 민감도(Overall/Maintenance): LOO 전 케이스 p<0.05, Mann-Whitney
  p=0.008/0.015, HC3 p=0.022/0.026, 무샘플 5명 제외 p=0.010/0.012
  (Induction은 GG 5명이라 LOO 최대 p=0.084로 취약)
- ADA endpoint, 우성 모델 전부 비유의 (열성 패턴 특이적)
- BH rank-1이라 q = p × m. q<0.05엔 m ≤ 6 필요(현재 11)

## Cohort attrition (Figure 1)

```
Analytic cohort 139
 └─ Exclusion 3: No infliximab records 41
Infliximab PopPK modeling cohort 98
 └─ Exclusion 4: 유전형 QC 단계에서 제거 1 (UID 17439372)
PGx analysis cohort 97
 ├─ Overall     97
 ├─ Induction   83  (15명은 induction phase 데이터 없음)
 └─ Maintenance 97
```

## 제출 전 해결할 항목

1. **popPK 모델 재추정** — Table 2 추정치는 체중 오류(834kg) 교정 **이전**
   값. 시뮬레이션은 교정 데이터로 재실행했으나 파라미터는 미재추정.
   재추정 후 Table 2 갱신 + PGx 재분석 필요.
2. **희소 샘플 환자 포함 규칙** — 관찰 농도 2건뿐인 환자가 코호트 최대 CL로
   참조군에 포함됨. 최소 관찰 농도 수 기준을 **결과 보기 전에** 확정할 것.
3. **Figure 1 스크리닝 단계 수치** — 상단 n=X,XXX 3곳은 EMR 추출 수치 필요.
4. 센터에 `23-B02281_EB-01` 시퀀싱 여부 확인 (Exclusion 4 문구 구체화용).
5. VPC/GOF 플롯을 NONMEM(`C:/Users/ilma0/NONMEMProjects/IBDPGX/`)에서 가져오기.

## 데이터 수정 이력 (모두 결과와 무관하게 타당한 수정)

1. **체중 834kg 오류** (UID 25269024, NONMEM ID 43): 83.4→834.0 소수점 오류.
   NONMEM 모델링 데이터셋 원본에도 존재. 교정 및 시뮬레이션 재실행 완료.
2. **PK가 PD 창에 연동된 버그**: `pd_endpoint_data_for_genomics2.py`에서 CL
   추출이 DE_DATE(1년차 PD 측정일 종속)에 묶여, PD 추적 없는 환자 5명의 CL이
   통째로 결측. PK를 PD와 분리 → MAINT 92→97, OVERALL 96→97.
   감사파일 `results/PKPD_EDA/GENOMICS/pk_window_audit.csv`
3. **파이프라인 `(for pda)` 브랜치 미갱신**: `added_filename_str` 플래그로
   두 갈래 파일셋이 생성되는데 PGx가 읽는 쪽만 낡아 있었음 (해결).
4. **최소 그룹 수 컷오프 8명 제거**: 기준이 임의적이고 핵심 결과의 GG군이
   정확히 8명이라 post-hoc 비판 소지 → 제거 (family 7→12→11).

### q값 변화 이력

| 단계 | MAINT q | OVERALL q |
|---|---|---|
| ① 최초 보고 (컷오프8, 체중오류有) | 0.048 | 0.031 |
| ② 컷오프 제거 (family 7→12) | 0.083 | 0.052 |
| ③ 체중 교정 + PD창 분리 (n 증가) | 0.107 | 0.091 |
| ④ 변이 QC 필터 (m 12→11) | **0.098** | **0.084** |

②는 비보정 p 불변, family 크기만 12/7=1.71배 → q도 정확히 1.71배 증가.

## 유전자 주석 수정

rs1061622 = **TNFRSF1B** (기존 자료의 TNF 표기는 오류), rs767455 = TNFRSF1A.
`gene_pd_cor/`의 옛 스크립트와 `rsid_genotype_summary.xlsx`에는 아직 잔존.

## 유전형 샘플 추적 코드 (`../gene_pd_cor/`)

| 파일 | 용도 |
|---|---|
| `genomics_sample_id_audit_hail.py` | Hail: VCF 샘플 목록 추출 + 타깃 4단계 탐색 |
| `genomics_mt_vs_vcf_check.py` | Hail: `mt` ↔ VCF 헤더 대조 (필터 영향 판정) |
| `genomics_sample_uid_reconcile.py` | 샘플ID↔UID 정규화 병합 + 감사 |
| `check_vcf_header_samples.py` | VCF 파일 단독 판정 (Hail 불필요) |

## 원본 데이터/코드 위치

- popPK 모델링 코드: `../modeling_codes_infliximab/`
- NONMEM: `C:/Users/ilma0/NONMEMProjects/IBDPGX/`
- PGx 분석 개발 이력: `../gene_pd_cor/`
- PD/PK 파생 데이터 생성: `../pd_analysis_codes(new2)/`
