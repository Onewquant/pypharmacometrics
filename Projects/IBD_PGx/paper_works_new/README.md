# paper_works_new — Infliximab popPK + Pharmacogenomic 논문 작업 폴더

논문 골자: **IBD 환자에서 infliximab 집단약동학 모델링 + 후보 유전변이의
PK(CL)/ADA 연관성 분석**

**상태: 원고 검토 반영 갱신본 준비 완료 (2026-09-07)** — 결론 변화 없음.
2026-09-07 변경: Figure 1 단순화(Excl 4 → Excl 3 병합, 98 박스 삭제),
Table 1을 97명·Figure 1 기준 phase로 재산출, ADA 양성 정의(분석 창 vs 추적
전체) 명시, Suppl Fig S1 GOF·S2 VPC 추가, rs1061622 전용 자료(Results 문단·
Suppl S6·Suppl Fig S4·S3 행) 삭제. 보충자료 최종 Table S1–S5, Figure S1–S3.
메일: `manuscript/email_to_professor_20260907.md`, 첨부 `for_professor_20260907/`

**이전 상태: 최종 수치 확정 + 원고/표·그림 전면 갱신 완료 (2026-09-03)**
개인 CL 산출을 sim90(난수 ETA) → sim95(run 89 EBE 주입)로 교정 완료.
최종 결론: **FDR 보정 후 유의한 변이 없음** (225개 검정, q<0.05 = 0건).
이전에 보고했던 rs1061622(TNFRSF1B) 신호는 난수 ETA의 산물로 확정 —
EBE 기반에서는 GMR 1.05, p=0.579로 완전 소멸. 최소 p는 rs396991(FCGR3A)
q=0.098이며 LOO에서 1명만 빼도 유의성 소실.

## 교수님 발송 자료

| 항목 | 위치 |
|---|---|
| **메일 본문 (2026-09-07, 갱신본)** | `manuscript/email_to_professor_20260907.md` |
| **첨부 4개 모음 (2026-09-07)** | `for_professor_20260907/` |
| 메일 본문 (2026-09-03, 발송용 초안 — 실제 발송 여부는 사용자 확인) | `manuscript/email_to_professor_20260903_send.md` |
| 첨부 4개 모음 (2026-09-03) | `for_professor_20260903/` |

첨부: `Methods_and_Results_FINAL.docx`,
`[IFX_POPPK]_core_fig_tab_FINAL.docx`,
`Figure3_forest_CL_overall.pdf`, `Figure1_eligibility_flowchart.pdf`

## 바로 쓸 수 있는 것 (논문 작성용)

| 항목 | 위치 |
|---|---|
| **원고 Word (Methods/Results/Discussion)** | `manuscript/Methods_and_Results_FINAL.docx` |
| 원고 마크다운 기준본 | `manuscript/Methods_and_Results_FINAL.md` |
| **표·그림 Word (캡션 포함)** | `core_fig_tab/[IFX_POPPK]_core_fig_tab_FINAL.docx` |
| **Figure/Table 캡션** | `core_fig_tab/CAPTIONS.md` |
| **본문 표·그림 파일** | `core_fig_tab/Table1~5, Figure1~3` |
| **보충자료** | `core_fig_tab/SupplTableS1~S5, SupplFigureS1~S3` |

기존 `output/Methods_and_results_reviewed.docx`와
`manuscript_optionA_final_sections.md`는 **폐기된 수치(sim90 난수 CL
기반, q=0.048 유의 등)** 이므로 참조하지 말 것. popPK 문단은
`Methods_and_Results_FINAL.docx`에 그대로 옮겨져 있습니다.

## 표·그림 구성 (2026-09-03 확정)

| 항목 | 내용 |
|---|---|
| Table 1–2 | 환자 특성 / popPK 파라미터 (재추정으로 확정) |
| Table 3 | 후보 변이 특성 (유전형 분포·MAF·HWE·QC) |
| Table 4 | 전체 기간 CL 연관성 (열성 11 + 우성 14) |
| **Table 5** | rs396991 시기별 결과 + LOO 범위 |
| Figure 1 | Eligibility flow chart (98→97→96) |
| Figure 2 | Structural model diagram |
| **Figure 3** | 전체 변이 forest plot (GMR·95% CI·p·q) |
| Suppl S1–S5 | 시기별 CL / ADA / 민감도(rs396991) / 원스케일 / 코호트 |
| Suppl Fig S1–S2 | popPK GOF / VPC (run 89) |
| **Suppl Fig S3** | rs396991 유전형별 CL 산점도 |
| ~~Suppl S6, Suppl Fig S4~~ | rs1061622 전용 자료 — 2026-09-07 삭제 (비유의, 별도 보고 근거 없음) |

## 폴더 구조

```
paper_works_new/
├── code/          분석·표·그림·docx 생성 스크립트 (01~10, 번호순 실행)
├── core_fig_tab/  ★ 논문용 최종 표·그림 + 캡션
├── manuscript/    ★ 원고 본문 최종본, 교수님 보고 메일 기록
├── data/          분석 입력 스냅샷
├── materials/     기존 paper_works에서 가져온 원본 자료
├── output/        중간 산출물 (raw 결과 CSV, 폐기된 구 docx 등)
└── for_professor_20260903/   교수님 발송 첨부 4개
```

## 실행 방법

프로젝트 루트 venv 사용. 01~06 실행 후 07(표), 08(Figure 1),
09(Figure 3 forest), 10(docx 생성) 순서로 실행.

```
C:/Users/ilma0/PycharmProjects/pypharmacometrics/venv/Scripts/python.exe -X utf8 code/<script>
```

## code/ ↔ 논문 요소 매핑

| 스크립트 | 산출물 | 논문 요소 |
|---|---|---|
| `01_table1_demographics.py` | `Table1_demographics.csv` | Table 1 (2026-09-06: IFX 열을 추정용 97명 데이터셋으로, Treatment phase는 `for_genomics_df` PHASE 기준으로 Figure 1과 일치시킴) |
| `02_table_genotype_summary.py` | `Table_genotype_summary.csv` | Table 3 재료 (유전형 분포/MAF/HWE) |
| `03_pgx_ancova_fdr.py` | `Table_pgx_ancova_fdr_results.csv`, `Table_variant_qc.csv` | **PGx 본분석** + 변이 QC |
| `04_pgx_sensitivity.py` | `Table_pgx_sensitivity.csv` | Suppl S3 (견고성) |
| `05_figure_cl_by_genotype.py` | `SupplFigureS3_*.png/pdf` | **Suppl Figure S3** (S1 GOF·S2 VPC는 NONMEM run 89 PDF를 수동 반입, PNG는 PyMuPDF 300 dpi 래스터) |
| `06_pgx_cohort_attrition.py` | `Table_pgx_attrition.csv` | Suppl S5, Figure 1 수치 |
| `07_core_tables.py` | `core_fig_tab/Table1~5, SupplS1~S5` | **논문용 표 전체** |
| `08_figure1_flowchart.py` | `core_fig_tab/Figure1_*.png/pdf` | **Figure 1** |
| `09_figure3_forest.py` | `core_fig_tab/Figure3_forest_*.png/pdf` | **Figure 3** |
| `10_build_manuscript_docx.py` | 원고 docx + 표·그림 docx | **Word 산출물** (python-docx 필요) |

## 확정된 분석 프레임

- infliximab만 (adalimumab 제외), phase 층화: **IND / MAINT / OVERALL**
- 모델: `ADA ~ GROUP+SEX+WEIGHT+ALBUMIN` (logistic),
  `CL ~ GROUP+SEX+WEIGHT+ALBUMIN+ADA` (일반 OLS ANCOVA)
- **log(CL)이 주 분석** (log-normal CL 가정; 잔차 정규성은 불완전하나 raw보다 우수), raw는 보조
- 최소 그룹 수 컷오프 **없음** (`MIN_GROUP_N=1`)
- 변이 QC 사전 필터: **MAF ≥ 0.05, HWE exact p ≥ 0.05**
  → 3개 제외 (rs3024505, rs765249238, rs776813259), 14개 통과
  → 열성 모델 검정 가능 11개(동형접합 0명인 3개 제외), 우성 모델 14개
- FDR(BH)은 (phase × endpoint × CL스케일 × 유전모델) 층 안에서 변이들에 대해 보정

## 핵심 결과 (EBE 기반 최종)

**FDR 보정 후 유의한 변이 없음** (225개 검정 중 q<0.05 = 0건, ADA 최소 q>0.99)

**비보정 p<0.05는 rs396991(FCGR3A) 열성비교 하나뿐** (전체 6행, 세 phase ×
log/raw). CC n=5, 가설생성 수준으로만 보고.

| Phase | GMR (95% CI) | p | q | LOO 최대 p |
|---|---|---|---|---|
| Maintenance | 1.36 (1.09–1.71) | 0.009 | 0.098 | 0.112 |
| Overall | 1.32 (1.06–1.64) | 0.014 | 0.157 | 0.141 |
| Induction | 1.35 (1.05–1.72) | 0.020 | 0.218 | 0.184 |

- **LOO에서 1명만 빼도 유의성 소실** → 소수 관측 의존. 과대 해석 금지
- CC 5명 개인 CL: 0.293/0.329/0.336/0.424/0.547 (최고값 1명 영향 큼)
- 우성 모델은 비유의 (OVERALL GMR 1.08, p=0.134)
- ADA endpoint 전부 비유의 (최소 p=0.076, q>0.99). PGx 코호트 ADA 양성
  6명, 도입기 창에는 양성 0명이라 ADA 검정 75개 중 26개는 추정 불가
  (실제 산출 검정 199개 = CL 150 + ADA 49)

**rs1061622 (TNFRSF1B)**: 이전 신호는 sim90 난수 ETA의 산물로 확정.
EBE 기반 OVERALL GMR 1.05 (0.88–1.27), p=0.579 (GG 기하평균 0.295 vs
TT+TG 0.284 L/day).

- 검증: sim95 implied ETA vs 89.phi EBE 상관 **r=1.000**
- 잔차 정규성: log 스케일도 불완전(Shapiro p 0.0006–0.006)하나
  raw(p<1e-6)보다 훨씬 나음 → log 주 분석 유지

## Cohort attrition (Figure 1)

```
Analytic cohort 139
 └─ Exclusion 3 (n=42): No infliximab records 41
     + 유지기 시작 + 관찰농도 1건 → 추정모델 제외, EBE 없음 (UID 38339532) 1
     (2026-09-06: 기존 Exclusion 4를 Exclusion 3에 통합, 98 중간 박스 삭제)
Infliximab PopPK modeling cohort 97
 └─ Exclusion 4: 유전형 QC 단계에서 제거 (UID 17439372) (n=1)
PGx 분석 코호트 96
 ├─ Overall     96
 ├─ Induction   83  (15명은 induction phase 데이터 없음)
 └─ Maintenance 96
```

## 제출 전 해결할 항목

1. ~~popPK 모델 재추정~~ → **완료/불필요 확정**: run 89 재추정 결과가 기존과
   동일 (체중 오류 레코드는 다음 관찰까지 1,200일 간격으로 우도 기여 0).
   Table 2 그대로 유효.
2. ~~희소 샘플 환자 포함 규칙~~ → **해소**: EBE 전환으로 관찰 1건 환자
   (38339532)는 추정모델에서 자연 제외. 관찰 2건 환자(35093356)의 CL은
   이제 실측 기반 EBE라 문제였던 극단값(0.617)도 사라짐.
3. **Figure 1 스크리닝 단계 수치** — 상단 n=X,XXX 3곳은 EMR 추출 수치 필요
   (`08_figure1_flowchart.py`의 최상단 박스 문자열).
4. 센터에 `23-B02281_EB-01` 시퀀싱 여부 확인 (Exclusion 문구 구체화용).
5. ~~VPC/GOF 플롯 추가~~ → **완료 (2026-09-07)**: `SupplFigureS1_GOF_model89`,
   `SupplFigureS2_VPC_model89` (run/vpc_89: 200 samples, 8 bins by count).
   VPC PDF 축 제목이 xpose 기본값이라 재출력 권장.
6. ~~Figure 3 대상 결정~~ → **전체 변이 forest plot으로 확정**(09 스크립트).
   유전형별 산점도는 Suppl Figure S3(rs396991)만 유지 (rs1061622 산점도는 2026-09-07 삭제).
7. ~~교수님 보고~~ → **메일 본문·첨부 준비 완료**
   (`manuscript/email_to_professor_20260903_send.md`). 발송 및 회신 대기.
   질의 5건: negative finding 방향, rs396991 본문 비중, forest plot 구성,
   Figure 1 스크리닝 수치 출처, 센터 문의 진행 여부.

## 데이터 수정 이력 (모두 결과와 무관하게 타당한 수정)

0. **개인 CL 소스 교정 (2026-09, 최종 결론을 바꾼 수정)**: 기존 sim90은
   `$SIM ONLYSIM` + Ω=0.0711로 개인 ETA를 난수 생성 → 개인 CL이 환자
   EBE와 무관했음(상관 −0.08). run 89 재추정의 EBE를 데이터 컬럼(EBE1)으로
   주입한 sim95로 교체. 관련: `../modeling_codes_infliximab/prepare_sim95_ebe.py`,
   `run/95.mod`, sim95 검증 r=1.000.

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
**⑤ (최종) EBE 전환 후 위 수치 전부 폐기** — rs1061622는 p=0.58~0.69로
신호 자체가 소멸. 위 표는 "난수 CL 기반에서의 이력"으로만 의미 있음.

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
