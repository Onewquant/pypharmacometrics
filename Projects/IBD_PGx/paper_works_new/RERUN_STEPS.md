# 재실행 절차 (체중 교정 반영 + EBE 기반 CL로 전환)

> **✅ 완료 (2026-09-02)**: 89 재추정 → 95.mod 검증 → sim95 생성 →
> 파이프라인 재실행 → 논문 자료 갱신까지 전부 완료됨.
> 최종 결과: FDR 유의 변이 없음. 이 문서는 이력 보존용.

작성 2026-08. 두 가지를 함께 해결하는 절차입니다.

1. **체중 교정 반영** — 최종모델 파라미터가 834 kg 오류 이전 추정치였음
2. **CL 산출 방식 수정** — 기존 `sim90`은 `$SIM ... ONLYSIM` + `$OMEGA 0.0711`로
   개인 ETA를 **난수로 새로 뽑았음**. 따라서 sim90의 개인 CL은 환자 EBE와
   무관했음 (sim ETA vs `93.phi` EBE 상관 **−0.08**, 평균 절대차 0.30).
   PGx endpoint로는 부적합.

시변 공변량(체중·알부민·ADA)에 따른 phase별 CL 변화를 반영한다는 원래
설계 의도는 그대로 유지합니다. 바꾸는 것은 **개인차 항만** 난수 →
재추정 EBE입니다:

```
CL_record = TVCL(공변량_t) x exp(EBE_i)      <- run 95
CL_record = TVCL(공변량_t) x exp(난수_i)      <- 기존 run 90
```

---

## 준비된 파일

| 파일 | 내용 |
|---|---|
| `run/94.mod` | 최종모델(93.mod) 재추정용. 데이터·구조 동일, TABLE만 94로 변경 |
| `modeling_codes_infliximab/prepare_sim95_ebe.py` | 94 결과로 `95.mod` + EBE 컬럼 데이터셋 자동 생성 |
| `pd_analysis_codes(new2)/pd_endpoint_data_for_genomics2.py` | `sim90` → `sim95` 읽도록 수정 완료 |

입력 데이터는 이미 교정본입니다 (확인됨):
- `infliximab_integrated_modeling_df_dayscale.csv` — 체중 max 119.5 kg
- `infliximab_integrated_simulation_df.csv` — 체중 max 119.5 kg

---

## 실행 순서

### 1) NONMEM: 최종모델 재추정

```
run/94.mod  실행   (Pirana/PsN 등 평소 방식대로)
```

산출: `94.lst`, `94.ext`, `94.phi`, `patab94` 등

확인할 것:
- OFV가 기존 −935.665와 크게 다르지 않은지 (체중 1건 교정이므로 소폭 변화 예상)
- 수렴 및 `$COV` 성공 여부
- THETA 추정치 (특히 CL_WT — 기존 0.4843)

### 2) EBE 기반 시뮬레이션 구성

```
venv python  Projects/IBD_PGx/modeling_codes_infliximab/prepare_sim95_ebe.py
```

하는 일:
- `94.ext`에서 최종 THETA를 읽어 `95.mod`의 `$THETA`에 FIX로 기록
  (기존 `90.mod`는 2025-11 값이 박혀 있었고 CL_WT가 0.515로 최종치와 6.3% 차이)
- `94.phi`의 posthoc ETA(1)을 시뮬 데이터셋에 `EBE1` 컬럼으로 추가
  → `infliximab_integrated_simulation_df_ebe.csv`
- `95.mod` 생성: `CL = TVCL * EXP(EBE1 + ETA(1))`, `$OMEGA` CL은 `0 FIX`

실행 시 THETA 변화량과 EBE 범위를 출력하므로 눈으로 확인 가능합니다.

### 3) NONMEM: 시뮬레이션

```
run/95.mod  실행
```

산출: `sim95`

### 4) 파생 데이터 및 논문 자료 재생성

```
# PD/PK 파생 데이터 (sim95 사용)
cd <pypharmacometrics 루트>
set PYTHONPATH=C:/Users/ilma0/PycharmProjects/pypharmacometrics
venv python "Projects/IBD_PGx/pd_analysis_codes(new2)/pd_endpoint_data_for_genomics2.py"

# for_genomics_df 사본 동기화
copy results/PKPD_EDA/GENOMICS/for_genomics_df(all_drugs).csv  gene_pd_cor/
copy results/PKPD_EDA/GENOMICS/for_genomics_df(all_drugs).csv  paper_works_new/data/

# 논문 표·그림 전체 재생성
venv python paper_works_new/code/01_table1_demographics.py
venv python paper_works_new/code/02_table_genotype_summary.py
venv python paper_works_new/code/03_pgx_ancova_fdr.py
venv python paper_works_new/code/04_pgx_sensitivity.py
venv python paper_works_new/code/05_figure_cl_by_genotype.py
venv python paper_works_new/code/06_pgx_cohort_attrition.py
venv python paper_works_new/code/07_core_tables.py
venv python paper_works_new/code/08_figure1_flowchart.py
```

4단계는 제가 대신 실행해 드릴 수 있습니다.

---

## 예상되는 결과 변화 — 미리 알아두실 것

재추정 후 EBE로 전환하면 **rs1061622 연관성이 사라질 가능성이 높습니다.**

구 데이터의 최종모델 EBE(`93.phi`)로 직접 확인한 결과:

| 기준 | GG(n=8) vs TT+TG(n=88) | p |
|---|---|---|
| 진짜 EBE (`93.phi`) | GMR ≈ **1.00** | **0.99** |
| sim90 난수 ETA (현재 논문 수치의 근원) | GMR ≈ 1.25 | 0.017 |

GG 8명의 실제 개인 ETA는 −0.124 ~ +0.114로 0 근처에 고르게 분포하며
참조군과 차이가 없었습니다. 즉 현재 보고 중인 GMR 1.29~1.30은 환자
약동학이 아니라 난수 배정에서 비롯된 것으로 보입니다.

이 점은 그동안의 관찰과도 일관됩니다:
- 환자 1~2명 편입/제외에 결과가 크게 흔들림
- 체중 오류 교정이 결과를 눈에 띄게 움직임 (난수 ETA에 공변량 예측값만
  곱해지는 구조라 공변량 변화가 CL에 직접 반영됨)
- 관찰 농도 2건인 환자가 코호트 최대 CL을 가짐

따라서 재실행 후에는 **PGx 결과 부분을 다시 정리해야 하고, 교수님께도
CL 산출 방식 오류와 그 영향을 보고**하실 필요가 있습니다. 검증 근거는
이 문서와 `paper_works_new/manuscript/`에 남겨두었습니다.

---

## 참고: run 88은 최종모델이 아님

2026-08-22에 실행된 `88.mod`는 CL 공변량이 BMI
(`CLBMI = (BMI/24)**THETA(10)`)이고 sex·weight가 없는 다른 구조입니다
(OFV −926.102 vs 최종 −935.665). 탐색용 실행으로 보이며, 체중 교정 후
최종모델 재추정은 아직 수행되지 않았습니다.
