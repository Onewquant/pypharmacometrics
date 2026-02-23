import pandas as pd
import numpy as np
from datetime import datetime, timedelta

prj_name = 'VPA'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource/ysy_req'
output_dir = f"{prj_dir}/results"

# sim_df = pd.read_csv(f'{output_dir}/sim114.csv')
dose_df = pd.read_csv(f'{resource_dir}/drug/VPA_dose_df.csv')
# dose_df = pd.read_csv('./prepdata/amikacin_dose.csv')
alt_lab_df = pd.read_csv(f"{resource_dir}/LAB_exam/VPA/VPA_LAB_ALT.csv")
ast_lab_df = pd.read_csv(f"{resource_dir}/LAB_exam/VPA/VPA_LAB_AST.csv")
tbil_lab_df = pd.read_csv(f"{resource_dir}/LAB_exam/VPA/VPA_LAB_T.B.csv")


dose_df = dose_df.rename(columns={'UID':'ID'})
alt_lab_df = alt_lab_df.rename(columns={'VALUE':'ALT'})
ast_lab_df = ast_lab_df.rename(columns={'VALUE':'AST'})
tbil_lab_df = tbil_lab_df.rename(columns={'VALUE':'T.B'})

"""
# (1) ALT > 120 and T.B >2.4
# (2) ALT/AST > 200
"""
# ==== 설정 (너 데이터 컬럼명에 맞춰 바꿔) ====
ID_COL = 'ID'
DOSE_DATE_COL = 'DATE'   # 투약일 컬럼
GAP_DAYS = 30            # 이보다 큰 공백이면 episode 끊기
WASHOUT_DAYS = 30        # risk window 연장

# ==== 1) 투약일 episode 생성 ====
dose = dose_df[[ID_COL, DOSE_DATE_COL]].copy()
dose[DOSE_DATE_COL] = pd.to_datetime(dose[DOSE_DATE_COL])
dose = dose.dropna(subset=[ID_COL, DOSE_DATE_COL]).sort_values([ID_COL, DOSE_DATE_COL])

dose['prev_dt'] = dose.groupby(ID_COL)[DOSE_DATE_COL].shift(1)
dose['gap'] = (dose[DOSE_DATE_COL] - dose['prev_dt']).dt.days

# 새 episode 시작: 첫 기록이거나 gap > GAP_DAYS
dose['new_ep'] = ((dose['prev_dt'].isna()) | (dose['gap'] > GAP_DAYS)).astype(int)
dose['ep_no'] = dose.groupby(ID_COL)['new_ep'].cumsum()

episode_df = (dose.groupby([ID_COL, 'ep_no'])
              .agg(EP_START=(DOSE_DATE_COL, 'min'),
                   EP_END=(DOSE_DATE_COL, 'max'),
                   N_DOSE=(DOSE_DATE_COL, 'count'))
              .reset_index())

episode_df['RW_END'] = episode_df['EP_END'] + pd.Timedelta(days=WASHOUT_DAYS)

episode_df.head()

LAB_ID_COL = 'ID'
LAB_DATE_COL = 'DATE'
ALT_COL = 'ALT'
AST_COL = 'AST'
TBIL_COL = 'T.B'   # 네 파일에서는 T.B.를 어떤 컬럼명으로 저장했는지에 맞춰 수정

alt = alt_lab_df[[LAB_ID_COL, LAB_DATE_COL, ALT_COL]].copy()
ast = ast_lab_df[[LAB_ID_COL, LAB_DATE_COL, AST_COL]].copy()
tbi = tbil_lab_df[[LAB_ID_COL, LAB_DATE_COL, TBIL_COL]].copy()

def normalize_id(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip()
    # 10000456.0 -> 10000456
    s = s.str.replace(r'\.0$', '', regex=True)
    return s

dose_df[ID_COL] = normalize_id(dose_df[ID_COL])
alt_lab_df[LAB_ID_COL] = normalize_id(alt_lab_df[LAB_ID_COL])
ast_lab_df[LAB_ID_COL] = normalize_id(ast_lab_df[LAB_ID_COL])
tbil_lab_df[LAB_ID_COL] = normalize_id(tbil_lab_df[LAB_ID_COL])

for df in [alt, ast, tbi]:
    df[LAB_DATE_COL] = pd.to_datetime(df[LAB_DATE_COL])

# 같은 ID/DATE에 여러 값이 있으면 (최대값 등)으로 하나로 정리 (보수적으로 MAX 권장)
alt = alt.groupby([LAB_ID_COL, LAB_DATE_COL], as_index=False)[ALT_COL].max()
ast = ast.groupby([LAB_ID_COL, LAB_DATE_COL], as_index=False)[AST_COL].max()
tbi = tbi.groupby([LAB_ID_COL, LAB_DATE_COL], as_index=False)[TBIL_COL].max()

lab_wide = alt.merge(ast, on=[LAB_ID_COL, LAB_DATE_COL], how='outer') \
              .merge(tbi, on=[LAB_ID_COL, LAB_DATE_COL], how='outer')

lab_wide = lab_wide.sort_values([LAB_ID_COL, LAB_DATE_COL])
lab_wide.head()
lab_wide['ID'] = normalize_id(lab_wide['ID'])

# 3) DILI flag 계산 (너 기준 반영)
# =========================
# Rule 1: ALT > 120 and TBil > 2.4
rule1 = (lab_wide[ALT_COL] > 120) & (lab_wide[TBIL_COL] > 2.4)

# Rule 2: ALT > 200 or AST > 200
rule2 = (lab_wide[ALT_COL] > 200) | (lab_wide[AST_COL] > 200)

lab_wide['DILI_FLAG'] = rule1 | rule2

# (선택) 어떤 규칙으로 걸렸는지 확인용
lab_wide['DILI_RULE'] = np.select(
    [rule1 & rule2, rule1, rule2],
    ['RULE1+RULE2', 'RULE1', 'RULE2'],
    default=''
)

# =========================
# 4) episode risk window 안에서 DILI 발생일(최초) 찾기
# =========================
episode_df[ID_COL] = normalize_id(episode_df[ID_COL])
lab_wide[LAB_ID_COL] = normalize_id(lab_wide[LAB_ID_COL])
tmp = episode_df.merge(lab_wide, left_on=ID_COL, right_on=LAB_ID_COL, how='left')

tmp = tmp[
    (tmp[LAB_DATE_COL] >= tmp['EP_START']) &
    (tmp[LAB_DATE_COL] <= tmp['RW_END'])
].copy()

dili_by_ep = (tmp[tmp['DILI_FLAG']]
              .sort_values([ID_COL, 'ep_no', LAB_DATE_COL])
              .groupby([ID_COL, 'ep_no'], as_index=False)
              .first()[[ID_COL, 'ep_no', LAB_DATE_COL, 'DILI_RULE']])

dili_by_ep = dili_by_ep.rename(columns={LAB_DATE_COL: 'DILI_DATE'})

episode_df = episode_df.merge(dili_by_ep, on=[ID_COL, 'ep_no'], how='left')
episode_df['DILI'] = episode_df['DILI_DATE'].notna().astype(int)

# =========================
# 5) 결과 확인
# =========================
print(episode_df[['ID','ep_no','EP_START','EP_END','RW_END','N_DOSE','DILI','DILI_DATE','DILI_RULE']].head(30))

# 환자 기준으로 "어떤 episode든 DILI 있으면 1" 요약도 필요하면:
patient_dili = (episode_df.groupby(ID_COL, as_index=False)
                .agg(DILI=('DILI','max'),
                     FIRST_DILI_DATE=('DILI_DATE','min')))
print(patient_dili.head(20))

patient_dili['DILI'].mean(), patient_dili['DILI'].sum(), patient_dili.shape[0]
# tmp = episode_df[episode_df['DILI']==1]
# tmp['DILI_RULE'].value_counts()
# tmp.to_csv(f'{output_dir}/DILI_event.csv', index=False, encoding='utf-8-sig')
episode_df.to_csv(f'{output_dir}/DILI_event(0_1).csv', index=False, encoding='utf-8-sig')
# patient_dili['DILI'].mean(), patient_dili['DILI'].sum(), patient_dili.shape[0]