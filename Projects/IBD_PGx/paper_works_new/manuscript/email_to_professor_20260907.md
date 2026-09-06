# 교수님 보고 메일 (2026-09-07, 원고 검토 반영본 송부)

> 2026-09-03 발송 메일(`email_to_professor_20260903_send.md`)의 후속입니다.
> 분석 결과와 결론은 그대로이고, 원고 검토 과정에서 정리한 변경 사항을
> 요약해 갱신본을 보내는 상황입니다.

**첨부파일 (4개)** — `for_professor_20260907/`

1. `Methods_and_Results_FINAL.docx` — Methods / Results / Discussion(참고용) 초안
2. `[IFX_POPPK]_core_fig_tab_FINAL.docx` — Table 1–5, Figure 1–3,
   Supplementary Table S1–S5, Supplementary Figure S1–S3 (캡션 포함)
3. `Figure1_eligibility_flowchart.pdf` — 본문 Figure 1 (단순화 반영)
4. `Figure3_forest_CL_overall.pdf` — 본문 Figure 3

**메일 제목(안)**: [IBD infliximab PGx] 원고 검토 반영 갱신본 송부
(Figure 1·Table 1 정리, 보충 GOF/VPC 추가, rs1061622 전용 자료 삭제)

---

교수님 안녕하십니까? 백승환입니다.

지난 3일 보내드린 원고와 표·그림을 다시 검토하면서 정리한 내용을 반영한
갱신본을 보내드립니다. **분석 결과와 결론(FDR 보정 후 유의 변이 없음,
rs396991은 가설 생성 수준)은 변한 것이 없고**, 아래는 모두 코호트
기술·표 구성·보충자료에 관한 정리입니다.

---

## 1. Figure 1 단순화와 코호트 기술 통일

지난 버전은 infliximab 코호트 98명 → popPK 모델링 97명(유지기 시작 +
농도 1건 환자 제외) → PGx 96명의 3단계였습니다. 98명 단계에서 1명만
빠지는 박스를 따로 두는 것이 번거로워, 이 1명을 "infliximab 기록 없음"
41명과 합쳐 **분석 코호트 139명 → 제외 42명 → infliximab PopPK 모델링
코호트 97명 → PGx 96명**으로 정리했습니다. Methods의 대상자 문단도
같은 흐름으로 고쳤고, 제외 사유(첫 농도가 모델 초기화에만 쓰여 개인
추정에 기여하지 못함)는 한 문장으로 남겼습니다.

## 2. Table 1 재산출

Figure 1을 바꾸면서 Table 1을 점검했더니 두 가지가 어긋나 있어 함께
고쳤습니다.

- **infliximab 열이 98명 기준**이었습니다. 추정에 쓴 97명 데이터셋으로
  다시 만들어 열 제목을 "Infliximab PopPK modeling cohort"로 했습니다.
  검체 394 → 393건, ADA 양성 9.2 → 9.3% 등 소폭 변동입니다.
- **Treatment phase 행이 Figure 1과 맞지 않았습니다.** 기존 표는 유지기
  시작 환자를 4명(infliximab) / 5명(분석 코호트)으로 집계했는데, Figure 1의
  유도기 분석 대상은 83명이라 유지기 시작이 14명이어야 합니다. 원인은
  표 산출용 데이터셋에 유지기 시작 환자의 첫 농도 행이 없어 판정
  규칙이 이들을 놓친 것이었습니다. Figure 1과 같은 기준(유도기 자료
  유무)으로 다시 산출해 **infliximab 83 / 14, 분석 코호트 117 / 22**가
  되었습니다. 분석 코호트 열은 "어느 약물이든 유도기 자료가 있는 환자"를
  Whole phases로 잡았습니다.

## 3. ADA 양성 수의 정의 차이 명시

Table 1의 ADA 양성 9명과 PGx 결과의 6명이 어긋나 보여 산출 경로를
확인했습니다. Table 1은 **추적기간 전체에서 한 번이라도 양성**, PGx 분석은
**분석 창(유도기: 첫 투여~3회차 투여, 유지기: 유지 투여 시작~1년차 평가
전 마지막 투여) 안에서의 양성**이었습니다. 3명은 1년차 창이 끝난 뒤
(1,161일·1,456일·2,165일 이후) 양성으로 전환한 환자여서 PGx 분석에서는
음성으로 들어간 것이 맞습니다. 두 수치가 모두 옳으므로 값은 두고,
Methods에 분석 창 정의를 한 문장 추가하고 Results·Discussion·Table 1
캡션에 "분석 기간 내 / 추적기간 전체" 한정어를 넣어 구분했습니다.

## 4. 보충자료: GOF·VPC 추가, rs1061622 전용 자료 삭제

- **Supplementary Figure S1(GOF), S2(VPC)** 를 새로 넣었습니다. 최종
  모델(run 89) 기준이며 VPC는 PsN으로 200회 시뮬레이션, 8개 구간입니다.
  Methods 모델 평가 항목에 VPC 설정 한 문장을, Results 모델 평가 문단에
  참조를 추가했습니다. 기존 유전형별 산점도는 S3으로 밀렸습니다.
- **rs1061622(TNFRSF1B) 전용 자료를 뺐습니다.** 14개 변이 모두 비유의인데
  이 변이만 Results 문단, Suppl Table S6, Suppl Figure S4, 민감도 분석
  행을 따로 갖는 것은 이전 분석 이력의 잔재일 뿐 근거가 없다고 판단했습니다.
  Table 4·Figure 3·Suppl Table S1에는 다른 변이와 같이 들어 있습니다.
  Discussion에는 "항-TNF 반응 관련 문헌이 있는 TNFRSF1B rs1061622도
  연관성이 없었다" 한 문장만 참고용으로 넣어 두었습니다(Discussion은
  제 담당이 아니어서 편집자 메모로 표시했습니다).
- 보충자료 최종 구성: **Supplementary Table S1–S5, Supplementary Figure
  S1–S3**.

## 5. 그 밖의 소소한 정리

- Results 환자 특성 문단의 불릿 2개(여성·소아 수)를 문단에 합쳤습니다.
- Methods 변이 QC 문단에서 "HWE P < 1×10⁻⁶ 필터링은 상위 전처리에서 이미
  적용" 문장은 앞 문단과 중복이라 삭제했습니다.
- 국문 원고(Methods_and_Results_FINAL_Korean.docx)도 같은 내용으로
  맞춰 두었습니다. 필요하시면 함께 보내드리겠습니다.

---

## 6. 확인 부탁드릴 점

1. **Table 1 분석 코호트 열의 Treatment phase 정의**: "어느 약물이든
   유도기 자료가 있음"을 Whole phases로 잡았는데(117 / 22), 이 정의가
   괜찮으실지요. 마땅하지 않다면 이 행을 빼는 것도 방법입니다.
2. **VPC 그림의 축 제목**이 현재 xpose 기본값("Independent variable /
   Dependent variable")입니다. 시간(일)·농도(µg/mL)로 다시 출력해 교체할
   예정입니다.
3. 지난 메일에서 여쭌 다섯 가지(negative finding 방향, rs396991 본문 비중,
   Figure 3 패널 구성, Figure 1 스크리닝 수치 출처, 정밀의료센터 문의)는
   회신 주시는 대로 반영하겠습니다.

검토 부탁드립니다. 감사합니다.

백승환 드림
