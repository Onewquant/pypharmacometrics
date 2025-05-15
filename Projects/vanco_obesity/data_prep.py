from tools import *
from pynca.tools import *

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/vanco_obesity"

res_cols = ['ID','DATETIME','DV','MDV','CMT','AMT','RATE']

id_df = pd.read_excel(f"{prj_dir}/resource/ID(2020-2024).xlsx")
dose_df = pd.read_excel(f"{prj_dir}/resource/Dose(2020-2024).xlsx")
conc_df = pd.read_excel(f"{prj_dir}/resource/VCMconc(2020-2024).xlsx")
# lab_df1 = pd.read_excel(f"{prj_dir}/resource/Lab(2020-2024).xlsx", sheet_name="약제")
# lab_df2 = pd.read_excel(f"{prj_dir}/resource/Lab(2020-2024).xlsx", sheet_name="진검약리")
# lab_df = pd.concat([lab_df1, lab_df2], ignore_index=True)
# lab_df.to_csv(f"{prj_dir}/resource/Lab(2020-2024)_Total.csv", index=False, encoding='utf-8-sig')
lab_df = pd.read_csv(f"{prj_dir}/resource/Lab(2020-2024)_Total.csv")

# dose_df.columns
dose_df = dose_df[['등록번호','약품 오더일자','[함량단위환산] 1회 처방량','[실처방] 경로','수행시간']].copy()
dose_df = dose_df[dose_df['[실처방] 경로'].isin(['IV','MIV','IVS'])].copy()
dose_df = dose_df[~dose_df['수행시간'].isna()].copy()
dose_df['수행시간'] = dose_df['수행시간'].map(lambda x:x.split(', '))
dose_df = dose_df.explode('수행시간')
dose_df['약품 오더일자'] = dose_df['약품 오더일자'].astype(str)
dose_df['수행시간'] = dose_df.apply(lambda x:x['약품 오더일자']+' '+x['수행시간'] if len(x['수행시간'].split(' '))<2 else x['수행시간'], axis=1)

dose_df = dose_df.rename(columns={'등록번호':'ID','수행시간':'DATETIME','[함량단위환산] 1회 처방량':'AMT'})
dose_df = dose_df[['ID','DATETIME','AMT']].copy()
# dose_df.to_csv(f"{prj_dir}/resource/dose_prep.csv", index=False)
dose_df = dose_df[dose_df['DATETIME'].map(lambda x:x.split('/')[-1]).isin(['Y','DR','Z','O'])].copy()
dose_df['DATETIME'] = dose_df['DATETIME'].map(lambda x:x.split('/')[0])
dose_df = dose_df.sort_values(['ID','DATETIME'])
dose_df['DV'] = '.'
dose_df['MDV'] = 1
dose_df['AMT'] = dose_df['AMT'].map(lambda x: x.replace('mg',''))
dose_df['RATE'] = dose_df['AMT']
dose_df['CMT'] = 1
dose_df = dose_df[res_cols].copy()

conc_df = conc_df[['등록번호','채혈일시','검사결과']].copy()
conc_df = conc_df.rename(columns={'등록번호':'ID', '채혈일시':'DATETIME','검사결과':'DV'})
conc_df = conc_df.sort_values(['ID','DATETIME'])
conc_df['DATETIME'] = conc_df['DATETIME'].astype(str).map(lambda x:x[:-3])
conc_df['MDV'] = '.'
conc_df['AMT'] = '.'
conc_df['RATE'] = '.'
conc_df['CMT'] = 1
conc_df = conc_df[res_cols].copy()

md_df = pd.concat([dose_df, conc_df], ignore_index=True)
md_df = md_df.sort_values(['ID','DATETIME'], ignore_index=True)

###########################

## 필요함수들

def calculate_time(group):
    base_time_row = group[group['MDV'] == 1]
    if base_time_row.empty:
        group['TIME'] = np.nan
        return group
    base_time = base_time_row.iloc[0]['DATETIME']
    group['TIME'] = (group['DATETIME'] - base_time).dt.total_seconds() / 3600
    return group

def remove_late_conc(group):
    # MDV가 0인 데이터만 필터링
    conc_rows = group[group['MDV'] == 0]
    if conc_rows.empty:
        return group

    # 최초 MDV=0 값의 TIME을 찾기
    first_dv_time = conc_rows['TIME'].min()

    # 60시간 후의 시간 계산 (TIME 기준으로)
    limit_time = first_dv_time + 60  # 60시간 이후

    # limit_time을 초과한 모든 데이터 삭제
    return group[group['TIME'] <= limit_time]

pd.to_datetime(md_df['DATETIME'], errors='coerce')
md_df = md_df.sort_values(['ID', 'DATETIME'], ignore_index=True)


# MDV가 0인 'ID' 제거: 투약 없는 데이터 제거
md_df['MDV'] = md_df['MDV'].replace('.',0).astype(int)
id_all_0 = md_df.groupby('ID')['MDV'].apply(lambda x: (x == 0).all())
id_to_remove = id_all_0[id_all_0].index
md_df = md_df[~md_df['ID'].isin(id_to_remove)].copy()

# --- TIME 계산 (위의 코드에서 TIME 계산) ---
# 먼저 TIME을 계산
md_df['DATETIME'] = pd.to_datetime(md_df['DATETIME'], errors='coerce')
md_df = md_df.groupby('ID', as_index=False).apply(calculate_time).copy()

# --- 농도 60시간 초과 제거 (TIME 기준으로): 첫 TDM만 cycle만 남기려고 ---

# 'ID' 열이 확실히 사라지지 않도록 'groupby'에서 as_index=False 사용
md_df = md_df.groupby('ID', as_index=False).apply(remove_late_conc).copy()

# 각 ID에 대해 `DV=0`, `MDV=0`, `TIME=0` 추가
zero_md_df = md_df.groupby(['ID'], as_index=False).agg({'DATETIME':'min'})
zero_dict = {'DV':0.0,'MDV':'.','CMT':1,'AMT':'.','RATE':'.','TIME':0.0}
for k, v in zero_dict.items():
    zero_md_df[k] = v
md_df = pd.concat([zero_md_df, md_df],ignore_index=True)

# Data sorting
md_df = md_df.sort_values(by=['ID', 'TIME', 'MDV'], ascending=[True, True, False]).drop_duplicates(['ID','TIME','DV','MDV','AMT'])

# 투약 간격이 일정치 않은 사람들 확인위해 따로 저장
irreg_interval_df = list()
for inx, df_frag in (md_df[md_df['MDV']==1]).groupby(['ID']):
    eval_ds = df_frag['TIME'].diff().dropna()
    cb_ratio = eval_ds/(eval_ds.shift(1).bfill())
    ac_ratio = (eval_ds.shift(-1).ffill())/eval_ds
    time_interval_peak_exist = (cb_ratio/ac_ratio > 3).sum()
    if time_interval_peak_exist > 0:
        irreg_interval_df.append(inx)

irreg_interval_df = md_df[md_df['ID'].isin(irreg_interval_df)].reset_index(drop=True)
irreg_interval_df.to_csv(f'{prj_dir}/resource/irreg_df.csv', index=False)

"""
# (1) Dose 없이 DV가 0이 아닌 측정 데이터가 먼저 나오는 사람들이 있다 (TIME < 0)
# (2) 쌩뚱맞은 첫 Dose 한참 이후 갑자기 Dose 나오는 사람들 있다 (TIME > 500 등등)
"""

final_cols = ['ID', 'TIME', 'DV', 'MDV', 'CMT', 'AMT', 'RATE']
md_df = md_df[final_cols]


# df_id 데이터에 성별을 M/F로 맵핑, AGE, SEX, HT, WT 계산
id_df.rename(columns={'성별': 'SEX', '당시나이': 'AGE', '몸무게': 'WT', '키': 'HT','등록번호':'ID','검사일':'DATETIME'}, inplace=True)
id_df["ID"] = id_df["ID"].astype(int)
id_df["DATETIME"] = pd.to_datetime(id_df["DATETIME"], errors="coerce")
id_df['SEX'] = id_df['SEX'].map({'M': 1, 'F': 0})

# 2. BMI 계산
id_df['BMI'] = round(id_df['WT'] / (id_df['HT'] / 100) ** 2, 2)
id_df['OB'] = np.where(id_df['BMI'] >= 30, 1, 0)  # BMI ≥ 30 → OB=1
id_df['LBM'] = np.where(
    id_df['SEX'] == 1,
    round(0.407 * id_df['WT'] + 0.267 * id_df['HT'] - 19.2, 2),
    round(0.252 * id_df['WT'] + 0.473 * id_df['HT'] - 48.3, 2)
)

# 날짜
lab_df.rename(columns={'등록번호':'ID','검사일자':'DATETIME','검사코드':'LAB_CODE','검사명':'LAB','검사결과':"VALUE"},inplace=True)
lab_df["ID"] = lab_df["ID"].astype(int)
lab_df["DATETIME"] = pd.to_datetime(lab_df["DATETIME"], errors="coerce")
lab_df["LAB"] = lab_df["LAB"].map(lambda x:x.split('(')[0].strip())
lab_df["VALUE"] = lab_df["VALUE"].map(lambda x:str(x).replace('>','').replace('<','').split(' (')[0].strip()).replace('-',np.nan).replace('+/-',np.nan)
val_list = list()
for v in lab_df["VALUE"]:
    try: val_list.append(float(v))
    except: val_list.append(np.nan)
lab_df["VALUE"] = val_list

# lab code <-> lab name 변환 dict 생성
lab_rename_dict = dict([(row['LAB_CODE'],row['LAB']) for inx, row in lab_df[['LAB_CODE','LAB']].drop_duplicates(['LAB_CODE']).iterrows()])

"""
# {'L2002': 'RBC', 'L2001': 'WBC', 'L201218': 'ANC', 'L2009': 'PLT', 'L2004': 'Hct', 'L2003': 'Hb', 'L6108': 'Bilirubin', 'L6109': 'Blood', 'L6110': 'Urobilinogen', 'L6111': 'Nitrite', 'L6112': 'WBC stick', 'L61301': 'RBC', 'L61302': 'WBC', 'L61306': 'Bacteria', 'L61021': 'Color', 'L61022': 'Turbidity', 'L61309': 'Others', 'L6103': 'Specific Gravity', 'L6104': 'pH', 'L6105': 'Protein', 'L6106': 'Glucose', 'L6107': 'Ketone', 'L8186': 'CRP', 'L8136': 'Bilirubin, total', 'L81354': 'eGFR-Cockcroft-Gault', 'L81353': 'eGFR-Schwartz', 'L81352': 'eGFR-CKD-EPI', 'L81351': 'eGFR', 'L8135': 'Creatinine', 'L8134': 'BUN', 'L8164': 'Uric Acid', 'L8166': 'Alkaline Phosphotase', 'L8163': 'Albumin', 'L8162': 'Protein, total', 'L8153': 'ALT', 'L8152': 'AST', 'L8165': 'Cholesterol'}
소변검사 여부 확인
"""


# 가장 가까운 검사일자 찾는 함수
mlab_df = lab_df.merge(id_df[['ID','DATETIME']].rename(columns={'DATETIME':'TDM_DT'}), on=['ID'], how='left')
nearTDM_lab_df = mlab_df[mlab_df['DATETIME'] <= mlab_df['TDM_DT']].sort_values(['ID','DATETIME']).groupby(['ID','LAB'], as_index=False).agg({'DATETIME':'last','VALUE':'last'})
nearTDM_lab_df = nearTDM_lab_df.pivot(index=['ID'],columns=['LAB'],values='VALUE')
nearTDM_lab_df = nearTDM_lab_df.reset_index(drop=False)
nearTDM_lab_df.index.name = None
nearTDM_lab_df.columns.name = None


# lab_results = id_df.apply(lambda row:mlab_df get_latest_lab_result(row, lab_df), axis=1)
# idlab_df = pd.concat([id_df, lab_results], axis=1)
# idlab_df.rename(columns=lab_rename_dict, inplace=True)
# id_df.columns
# md_df와 df_idlab을 'ID' 기준으로 병합
final_df = md_df.merge(id_df[['ID','SEX','AGE','HT','WT','BMI','OB','LBM']], how='left', on=['ID']).merge(nearTDM_lab_df, how='left', on=['ID'])
# nearTDM_lab_df['ID'].iloc[0]
# md_df['ID'].iloc[0]
# final_df.columns
final_df = final_df.rename(columns={'Bilirubin':'TBil','eGFR-CKD-EPI':'eGFR_CE','eGFR-Schwartz':'eGFR_Shw','eGFR':'eGFR_MDRD','Alkaline Phosphotase':'ALP','Albumin':'ALB','Creatinine':'Cr', 'eGFR-Cockcroft-Gault':'eGFR_CG'})

final_df.columns
# 중복 컬럼 제거
final_df = final_df[['ID', 'TIME', 'DV', 'MDV', 'CMT', 'AMT', 'RATE', 'SEX', 'AGE', 'HT', 'WT', 'BMI',
                     'OB', 'LBM', 'eGFR_CE', 'eGFR_Shw', 'AST', 'WBC', 'PLT', 'eGFR_MDRD', 'CRP', 'TBil',
                     'ALP', 'ALB', 'Hct', 'pH', 'RBC', 'ALT', 'ANC', 'Cr', 'BUN', 'eGFR_CG', 'Hb']]

final_df.to_csv(f'{prj_dir}/resource/final_df.csv', index=False)

# 결과 출력
# print(final_df.head(50))
