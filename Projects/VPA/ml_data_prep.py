from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
import glob
from scipy.stats import spearmanr

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'VPA'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource/ysy_req'
output_dir = f"{prj_dir}/results"


dili_df = pd.read_csv(f"{output_dir}/DILI_event(0_1).csv")
dili_df = dili_df.rename(columns={'ID':'UID'})
dili_df['UID'] = dili_df['UID'].astype(str)
# dili_df = dili_df.drop_duplicates(['UID'])
# dili_df = dili_df[dili_df['DILI'] > 0].copy()
# dili_df[['EP_START','EP_END']]
# dili_df.columns
print(f"Total DILI cases: {len(dili_df['UID'].drop_duplicates())}")
# except_dili_df = dili_df[dili_df['DILI_DT']<24].copy()
# excepted_dili_ids = set(except_dili_df['UID'])
# print(f"-(Patients with DILI < 24hr): {excepted_dili_ids}")
# dili_df = dili_df[dili_df['DILI_DT']>=24].copy()
# print(f"DILI cases (Filtered): {len(dili_df['UID'].drop_duplicates())}")


demo_df = pd.read_csv(f"{resource_dir}/demo/VPA_DEMO_SEX,AGE.csv")
demo_df = demo_df.rename(columns={'ID':'UID'})
demo_df['UID'] = demo_df['UID'].astype(str)
demo_df['SEX'] = demo_df['SEX'].map({'남':0, '여':1})
#
# uid_df_dict = {}
#
# for i in range(1,5):
#     res_filename = f'251202_재식별 파일_DILI_{i}.csv'
#     uid_df_dict[i] = pd.read_csv(f'{resource_dir}/re_identification/{res_filename}')
#     uid_df_dict[i] = uid_df_dict[i].rename(columns={'Deidentification_ID':'UID','환자번호':'ID'})
#     uid_df_dict[i]['ID'] = uid_df_dict[i]['ID'].map(lambda x:x.split('-')[0])
#     uid_df_dict[i] = uid_df_dict[i][['ID','UID']].copy()
#
# fdf.columns
# fdf = fdf.rename(columns={'항목명':'LAB','항목값':'VALUE','환자번호':'ID','간호기록 작성일자':'DATE'})[['ID','DATE', 'LAB','VALUE']]
# # fdf = fdf.rename(columns={'BS':'LAB',})[['ID','DATE', 'LAB','VALUE']]
# # lab_df = pd.read_csv(lab_file)
# fdf['ID'] = fdf['ID'].map(lambda x:x.split('-')[0])
#
# for i in range(1,5):
#     # i = 2
#     fdf_test = fdf.merge(uid_df_dict[i], on=['ID'], how='left')
#     uid_pct = ((~(fdf_test['UID'].isna())).sum()) / len(fdf_test)
#     print(uid_pct)
#     if uid_pct > 0.7:
#         # raise ValueError
#         # lab_df_test
#         fdf_test['ID'] = fdf_test['UID']
#         fdf_test = fdf_test.drop(['UID'], axis=1)
#         fdf_test.to_csv(f"{resource_dir}/BS_exam/VPA_BS_HT.csv", index=False, encoding='utf-8-sig')
#         print(f"파일 치환 완료 | {resource_dir}/BS_exam/VPA_BS_HT.csv")
#         break
#     else:
#         continue
#
# fdf = pd.read_csv(f"{resource_dir}/BS_exam/VPA_BS_HT.csv")


bs_files = glob.glob(f"{resource_dir}/BS_exam/VPA_BS_*.csv")
for finx, fpath in enumerate(bs_files): #break
    print(f"({finx}) Append / {fpath} ")
    fdf = pd.read_csv(fpath)
    if type(fdf['ID'].iloc[0])==float:
        fdf['ID']=fdf['ID'].map(lambda x:str(x).split('.')[0])
    if finx==0:
        bsize_df = fdf.copy()
    else:
        bsize_df = pd.concat([bsize_df, fdf])

# bsize_df[bsize_df['LAB']=='WT']
bsize_df = bsize_df.rename(columns={'ID':'UID'})
bsize_df['LAB'] = bsize_df['LAB'].map({'BMI':'BMI','몸무게':'WT','키':'HT'})
bsize_df['UID'] = bsize_df['UID'].astype(str)
bsize_df = bsize_df.reset_index(drop=True)
for inx, row in bsize_df.iterrows():
    wt_value = str(row['VALUE']).split('(')[0]
    try:
        wt_value = float(wt_value)
    except:
        wt_value = np.nan
    bsize_df.at[inx,'VALUE'] = wt_value
bsize_df['VALUE'] = bsize_df['VALUE'].astype(float)
bsize_df = bsize_df.dropna(subset=['VALUE'])
# bsize_df['VALUE'].map(lambda x:str(x).split('(')[0])

# lab_list_df = pd.read_csv(f"{output_dir}/vpa_lablist_df.csv")
# total_lab_cols = lab_list_df['LAB'].to_list()
total_lab_cols = {'WBC', 'ANC', 'Lymphocyte', 'RBC', 'Hb', 'PLT', 'Platelet', 'Na', 'K', 'Ca', 'Ca, total', 'P', 'Cl',
                'Mg', 'Protein', 'Protein (mg/dL)', 'T. Protein', 'Glucose', 'Albumin', 'CRP', 'hsCRP', 'ESR', 'BUN',
                'Creatinine', 'Cr', 'Cr (S)', 'SCr', 'eGFR-CKD-EPI', 'eGFR-Cockcroft-Gault', 'eGFR-MDRD', 'AST',
                'AST(GOT)', 'ALT', 'ALT(GPT)', 'Alk. phos', 'Alk. phos.', 'γ-GT', 'Bilirubin, total', 'T. Bil.',
                'PT INR', 'PT sec', 'aPTT', 'B2-MG(S)', 'Cystatin C', 'HCO3-'}

res_lab_cols = ['DATE','WBC', 'ANC', 'Lymphocyte', 'RBC', 'Hb', 'PLT', 'Na', 'K', 'Ca', 'P', 'Cl', 'Mg', 'Protein',
                    'Glucose', 'Albumin', 'CRP', 'ESR', 'BUN', 'CREATININE', 'eGFR-CKD-EPI', 'eGFR-Cockcroft-Gault',
                    'eGFR-MDRD', 'AST', 'ALT', 'ALP', 'γ-GT', 'TBIL', 'PT INR', 'PT sec', 'aPTT', 'B2-MG(S)',
                    'Cystatin C', 'HCO3-']
# lab_result_df = list()

# print(f"# 각 환자별 날짜별 랩수치 파악 시작\n")
#
lab_files = glob.glob(f"{resource_dir}/LAB_exam/VPA/*.csv")
for finx, fpath in enumerate(lab_files): #break
    print(f"({finx}) Append / {fpath} ")
    fdf = pd.read_csv(fpath)
    if finx==0:
        lab_df = fdf.copy()
    else:
        lab_df = pd.concat([lab_df, fdf])

lab_df = lab_df.rename(columns={'ID':'UID'})
lab_df['UID'] = lab_df['UID'].astype(str)
#
#     pid = fpath.split('(')[-1].split('_')[0]
#     pname = fpath.split('_')[-1].split(')')[0]
#
#     print(f"({finx}) {pname} / {pid}")
#
#     # fdf.columns
#     fdf = pd.read_excel(fpath)
#     fdf['DATETIME'] = fdf[['보고일', '오더일']].max(axis=1)
#     fdf['LAB'] = fdf['검사명']
#     fdf['VALUE'] = fdf['검사결과']
#     fdf['VALUE'] = pd.to_numeric(fdf['VALUE'], errors='coerce') # 숫자 아닌 랩은 NAN으로 변환
#     fdf = fdf.drop_duplicates(['DATETIME','LAB'], keep='last', ignore_index=True)
#     # fpv_df = fdf.pivot_table(index='DATETIME', columns='LAB', values='VALUE', aggfunc='min').reset_index(drop=False).fillna(method='ffill')
#     fpv_df = fdf.pivot_table(index='DATETIME', columns='LAB', values='VALUE', aggfunc='min').reset_index(drop=False)
#     # fpv_df = fdf.pivot(index='DATETIME', columns='LAB', values='VALUE').reset_index(drop=False).fillna(method='ffill')
#
#     # fpv_df = fdf.pivot_table(index='DATETIME', columns='LAB', values='VALUE', aggfunc='min').reset_index(drop=False)
#     fpv_df.columns.name = None
#     fpv_df['UID'] = [pid]*len(fpv_df)
#     ind_lab_cols = list(fpv_df.columns)

# totlab_df = pd.read_csv(f"{output_dir}/totlab_df.csv")
# totlab_df = totlab_df.rename(columns={"DATETIME":'DATE'})

# comed_df = pd.read_csv(f"{resource_dir}/cm/VPA_CM_DIAG.csv")
# # comed_df['DATE_LIST'].iloc[0]
# comed_df['DATE_LIST'] = comed_df['DATE_LIST'].map(lambda x:x.replace("('",'').replace("',)",'').replace("')",'').split("', '"))
# comed_df = comed_df.explode('DATE_LIST')
# comed_df = comed_df.rename(columns={"DATE_LIST":'DATE'})
# comed_dict = pd.read_excel(f"{resource_dir}/[AMK_AKI_ML_DATA]/CCM_index_table.xlsx")
# comed_dict = comed_dict.set_index('CAT_NUM')['CCM_CAT']
# comed_dict.index.name = None
# comed_dict = dict(comed_dict)
# comed_df['DATE_LIST'].columns
cm_df = pd.read_csv(f"{resource_dir}/cm/VPA_CM_DIAG.csv")
cm_df = cm_df.rename(columns={'ID':'UID'})
cm_df['UID'] = cm_df['UID'].astype(str)
# cm_dict = pd.read_excel(f"{resource_dir}/cm/CM_index_table.xlsx")
# cm_dict = cm_dict.set_index('CAT_NUM')['CM_CAT']
# cm_dict.index.name = None
# cm_dict = dict(cm_dict)

vs_df = pd.read_csv(f"{resource_dir}/VS_exam/VPA_VS.csv")
vs_df = vs_df.rename(columns={'ID':'UID'})
vs_df['UID'] = vs_df['UID'].astype(str)
# vs_df = vs_df.drop_duplicates(['UID','DATE','VS_TYPE'],keep='last')

# proc_df = pd.read_csv(f"{output_dir}/final_procedure_df.csv")
# proc_df['PROC_CATNUM'].unique()
# proc_dict = {1:'CARDIAC SURGERY',2:'SURGERY',3:'ICU',4:'MECHANICAL VENT'}
#
# micbio_df = pd.read_csv(f"{output_dir}/final_microbiology_df.csv")
# micbio_dict = {1:"Pseudomonas aeruginosa",
#              2:"Acinetobacter baumannii",
#              3:"Klebsiella pneumoniae",
#              4:"Escherichia coli",
#              5:"Enterobacter species",
#              6:"Serratia marcescens",
#              7:"Mycobacterium tuberculosis",
#              8:"NonTuberculous Mycobacterium",
#              9:"Others",
#              10:"Empirical"}

# losmot_df = pd.read_csv(f"{output_dir}/final_locmotality_df.csv")
dose_df = pd.read_csv(f'{resource_dir}/drug/VPA_dose_df.csv')
dose_df = dose_df.rename(columns={'DOSE':'AMT'})
dose_df['UID'] = dose_df['UID'].astype(str)

ml_df = dose_df.drop_duplicates(['UID'])[['UID']].reset_index(drop=True)
ml_res_df = list()
no_lab_pids = list()
no_dose_pids = list()
no_demo_pids = list()
for inx, row in ml_df.iterrows(): #break
    # raise ValueError
    ml_row_dict = dict()
    uid = row['UID']
    print(f"({inx}) {uid}")

    # 개별 dataframe 로드
    uid_dili_rows = dili_df[dili_df['UID'] == uid].copy()
    uid_demo_df = demo_df[demo_df['UID'] == str(uid)].copy()
    uid_bsize_df = bsize_df[(bsize_df['UID']==uid)&(~bsize_df['LAB'].isna())].copy()

    uid_dose_df = dose_df[dose_df['UID']==uid].copy()
    uid_dose_df['AMT'] = uid_dose_df['AMT'].astype(float)

    if len(uid_demo_df)==0:
        no_demo_pids.append(uid)
        print(f"/ No demo data")
        continue
    uid_demo_last_row = uid_demo_df.iloc[-1]
    # uid_bsize_row = uid_bsize_df
    uid_lab_df = lab_df[lab_df['UID']==uid].copy()
    # uid_lab_df['LAB'] = uid_lab_df['검사명']Z
    # uid_lab_df['VALUE'] = uid_lab_df['검사결과']
    # uid_lab_df['VALUE'] = pd.to_numeric(uid_lab_df['VALUE'], errors='coerce')  # 숫자 아닌 랩은 NAN으로 변환
    uid_lab_df = uid_lab_df.drop_duplicates(['DATE', 'LAB'], keep='last', ignore_index=True)
    # uid_lab_df = uid_lab_df.drop(['보고일', '오더일','검사명','검사결과', '참고치', '직전결과'],axis=1)


    uid_lab_df = uid_lab_df.pivot_table(index='DATE', columns='LAB', values='VALUE', aggfunc='min').reset_index(drop=False)
    uid_lab_df.columns.name=None

    raw_lab_cols = set(list(uid_lab_df.columns))
    # uid_lab_df[uid_lab_df['LAB']=='CRP']
    existing_cols = total_lab_cols.intersection(raw_lab_cols)
    not_exist_cols = total_lab_cols.difference(raw_lab_cols)
    # pivoted_uid_lab_df[['DATE'] + list(existing_cols)]
    for lab_col in not_exist_cols:
        uid_lab_df[lab_col] = np.nan
    uid_lab_df = uid_lab_df[['DATE']+list(total_lab_cols)].copy()

    # total_lab_cols = {'WBC', 'ANC', 'Lymphocyte', 'RBC', 'Hb', 'PLT', 'Platelet', 'Na', 'K', 'Ca', 'Ca, total', 'P',
    #                   'Cl',
    #                   'Mg', 'Protein', 'Protein (mg/dL)', 'T. Protein', 'Glucose', 'Albumin', 'CRP', 'hsCRP', 'ESR',
    #                   'BUN',
    #                   'Creatinine', 'Cr', 'Cr (S)', 'SCr', 'eGFR-CKD-EPI', 'eGFR-Cockcroft-Gault', 'eGFR-MDRD', 'AST',
    #                   'AST(GOT)', 'ALT', 'ALT(GPT)', 'Alk. phos', 'Alk. phos.', 'γ-GT', 'Bilirubin, total', 'T. Bil.',
    #                   'PT INR', 'PT sec', 'aPTT', 'B2-MG(S)', 'Cystatin C', 'HCO3-'}
    uid_lab_df['PLT'] = uid_lab_df[['PLT', 'Platelet']].max(axis=1)
    uid_lab_df['Ca'] = uid_lab_df[['Ca', 'Ca, total']].max(axis=1)
    uid_lab_df['Protein'] = uid_lab_df[['Protein','Protein (mg/dL)','T. Protein']].max(axis=1)
    uid_lab_df['CRP'] = uid_lab_df[['CRP','hsCRP']].max(axis=1)
    uid_lab_df['CREATININE'] = uid_lab_df[['Creatinine', 'Cr', 'Cr (S)', 'SCr']].max(axis=1)
    uid_lab_df['AST'] = uid_lab_df[['AST', 'AST(GOT)']].max(axis=1)
    uid_lab_df['ALT'] = uid_lab_df[['ALT', 'ALT(GPT)']].max(axis=1)
    uid_lab_df['ALP'] = uid_lab_df[['Alk. phos', 'Alk. phos.']].max(axis=1)
    uid_lab_df['TBIL'] = uid_lab_df[['Bilirubin, total', 'T. Bil.']].max(axis=1)
    uid_lab_df = uid_lab_df.drop(['Platelet','Ca, total','Protein (mg/dL)','T. Protein','hsCRP','Creatinine', 'Cr', 'Cr (S)', 'SCr','AST(GOT)','ALT(GPT)','Alk. phos', 'Alk. phos.','Bilirubin, total', 'T. Bil.'],axis=1)
    uid_lab_df = uid_lab_df[res_lab_cols].copy()

    if len(uid_lab_df)==0:
        no_lab_pids.append(uid)
        print(f"/ No lab data")
        continue

    # raise ValueError
    # uid_lab_df = totlab_df[totlab_df['UID']==uid].copy()

    # uid_comed_df = comed_df[comed_df['UID']==uid].copy()
    uid_cm_df = cm_df[cm_df['UID']==uid].copy()
    uid_vs_df = vs_df[vs_df['UID']==uid].copy()
    # uid_proc_df = proc_df[proc_df['UID']==uid].copy()
    # uid_micbio_df = micbio_df[micbio_df['UID']==uid].copy()
    # uid_losmot_df = losmot_df[losmot_df['UID']==uid].copy()

    # DILI 발생여부 기록

    # DILI 발생 한 사람 (DILI 발생 시점이 첫 dosing 24시간 이후 일때만 포함되어 있음)
    if len(uid_dili_rows[uid_dili_rows['DILI']>0]):
        dili_occurrence = 1
        uid_cycle_row = uid_dili_rows[uid_dili_rows['DILI']>0].iloc[0]
        dili_occurrence_date = uid_cycle_row['DILI_DATE']
    else:
        # uid_dili_rows['N_DOSE']
        dili_occurrence = 0
        uid_cycle_row = uid_dili_rows.sort_values(['N_DOSE'], ascending=False).iloc[0]
        dili_occurrence_date = uid_cycle_row['EP_END'] # DOSING 가장 오래한 CYCLE로 선정

    # DILI 발생일 (or DOSING을 가장 오래한 CYCLE의 마지막 날) 기준으로 1일 전에 DOSING 기록이 존재하는 환자만 남김

    ml_dose_df = uid_dose_df[uid_dose_df['DATE'] < dili_occurrence_date].copy()
    if len(ml_dose_df) == 0:
        no_dose_pids.append(uid)
        print(f"/ No dose data")
        continue

    # BASELINE Features 정리
    min30_bl_date = (datetime.strptime(dili_occurrence_date, '%Y-%m-%d') - timedelta(days=30)).strftime('%Y-%m-%d')
    min14_bl_date = (datetime.strptime(dili_occurrence_date, '%Y-%m-%d') - timedelta(days=14)).strftime('%Y-%m-%d')
    min7_bl_date = (datetime.strptime(dili_occurrence_date, '%Y-%m-%d') - timedelta(days=7)).strftime('%Y-%m-%d')
    min3_bl_date = (datetime.strptime(dili_occurrence_date, '%Y-%m-%d') - timedelta(days=3)).strftime('%Y-%m-%d')
    min_bl_date = (datetime.strptime(dili_occurrence_date, '%Y-%m-%d') - timedelta(days=3)).strftime('%Y-%m-%d')
    max_bl_date = (datetime.strptime(dili_occurrence_date, '%Y-%m-%d') - timedelta(days=1)).strftime('%Y-%m-%d')
    # max_bl_dt = (datetime.strptime(dosing_start_date, '%Y-%m-%dT%H:%M') - timedelta(hours=24)).strftime('%Y-%m-%dT%H:%M')

    # 데이터 중간 처리
    # uid_cm_df['CM_PERIOD_FROM_DATE'] = uid_cm_df['CM_PERIOD_FROM'].map(lambda x: (datetime.strptime(max_bl_date, '%Y-%m-%d') - timedelta(days=x)).strftime('%Y-%m-%d'))
    # uid_cm_df.columns
    # uid_proc_df['CM_PERIOD_FROM_DATE'] = uid_proc_df['PROC_PERIOD_FROM'].map(lambda x: (datetime.strptime(max_bl_date, '%Y-%m-%d') - timedelta(days=x)).strftime('%Y-%m-%d'))

    # ml feature 기록용 데이터 필터링
    ml_bsize_df = uid_bsize_df[(uid_bsize_df['DATE'] < dili_occurrence_date)&(uid_bsize_df['DATE'] >= min7_bl_date)].copy()

    uid_demo_df['AGE'] = int((datetime.strptime(max_bl_date, '%Y-%m-%d') - datetime.strptime(uid_demo_df['BIRTHDATE'].iloc[0], '%Y-%m-%d')).total_seconds()/86400/365.25)


    # raise ValueError
    # bsize 데이터 정상화 되면 uid_demo_df에다가 WT, HT 붙여서 사용해야함
    ml_demo_df_last_row = uid_demo_df.iloc[-1]
    if len(ml_bsize_df[ml_bsize_df['LAB']=='WT'])>0:
        ml_demo_df_last_row['WT'] = ml_bsize_df[ml_bsize_df['LAB']=='WT']['VALUE'].iloc[-1]
    elif (len(ml_bsize_df[ml_bsize_df['LAB']=='HT'])>0) and (len(ml_bsize_df[ml_bsize_df['LAB']=='BMI'])>0):
        bmi = ml_bsize_df[ml_bsize_df['LAB'] == 'BMI']['VALUE'].iloc[-1]
        ht = ml_bsize_df[ml_bsize_df['LAB'] == 'HT']['VALUE'].iloc[-1]
        ml_demo_df_last_row['WT'] = bmi * (ht/100)**2
    else:
        ml_demo_df_last_row['WT'] = np.nan

        # ml_lab_df = uid_lab_df[(uid_lab_df['DATE'] <= max_bl_date)&(uid_lab_df['DATE'] >= min_bl_date)].copy()
    ml_lab_df = uid_lab_df[(uid_lab_df['DATE'] <= max_bl_date)&(uid_lab_df['DATE'] >= min30_bl_date)].copy()
    # raise ValueError
    # ml_comed_df = uid_comed_df[(uid_comed_df['DATE'] <= max_bl_date)&(uid_comed_df['DATE'] >= min_bl_date)].copy()
    ml_cm_df = uid_cm_df[(uid_cm_df['DATE'] >= min30_bl_date)&(uid_cm_df['DATE'] <= max_bl_date)].copy()
    ml_vs_df = uid_vs_df[(uid_vs_df['DATE'] <= max_bl_date)&(uid_vs_df['DATE'] >= min30_bl_date)].copy()
    # ml_proc_df = uid_proc_df[(uid_proc_df['DATE'] >= uid_proc_df['CM_PERIOD_FROM_DATE'])&(uid_proc_df['DATE'] <= max_bl_date)].copy()
    # ml_micbio_df = uid_micbio_df[(uid_micbio_df['DATE'] <= max_bl_date)&(uid_micbio_df['DATE'] >= min14_bl_date)].copy()
    # uid_losmot_df = uid_micbio_df[(uid_micbio_df['DATE'] <= max_bl_date)&(uid_micbio_df['DATE'] >= min14_bl_date)].copy()
    # uid_dose_df['AMT']


    # first_dose_dt = ml_dose_df.iloc[0]['DATETIME_ORI']
    # last_dose_dt = ml_dose_df.iloc[-1]['DATETIME_ORI']
    first_dose_dt = uid_dose_df.iloc[0]['DATE']
    last_dose_dt = ml_dose_df.iloc[-1]['DATE']

    ## ml feature 기록
    ml_row_dict['UID'] = uid

    # Demographics
    # for c in ['ID','AGE','SEX','WT','HT']:
    for c in ['UID','AGE','SEX','WT']:
        ml_row_dict[c] = ml_demo_df_last_row[c]
    # ml_row_dict['BMI'] = round(ml_row_dict['WT']/((ml_row_dict['HT']/100)**2),3)


    # DOSE (현재 WT 정보가 정확치 않아서 WT 로 계산된 것은 주석처리되어 있음)
    ml_row_dict['LAST_DOSE'] = ml_dose_df.iloc[-1]['AMT']
    ml_row_dict['LAST_DAILY_DOSE'] = ml_dose_df.groupby('DATE').agg({'AMT':'sum'}).iloc[-1]['AMT']
    ml_row_dict['LAST_DOSE_perWT'] = round(ml_row_dict['LAST_DOSE']/ml_row_dict['WT'],3)
    ml_row_dict['LAST_DAILY_DOSE_perWT'] = round(ml_row_dict['LAST_DAILY_DOSE']/ml_row_dict['WT'],3)
    # ml_row_dict['DURATION(DAYS)'] = (ml_dose_df.iloc[-1]['TIME'] - ml_dose_df.iloc[0]['TIME'] + 1)/24
    ml_row_dict['DURATION(DAYS)'] = round((datetime.strptime(ml_dose_df.iloc[-1]['DT_DOSE'].split(':')[0], '%Y-%m-%dT%H') - datetime.strptime(ml_dose_df.iloc[0]['DT_DOSE'].split(':')[0], '%Y-%m-%dT%H')).total_seconds()/86400,3)
    ml_row_dict['CUM_DOSE'] = ml_dose_df['AMT'].sum()
    ml_row_dict['MEAN_DAILY_DOSE'] = round(ml_row_dict['CUM_DOSE']/ml_row_dict['DURATION(DAYS)'],3)
    ml_row_dict['MEAN_DAILY_DOSE_perWT'] = round(ml_row_dict['CUM_DOSE']/ml_row_dict['WT'], 3)


    ## COMORBIDITY
    ml_row_dict['CM_COUNT'] = len(ml_cm_df.drop_duplicates(['DIAGNOSIS']))

    ## CONCOMITANT MEDICATION
    # uid_uniq_comed_inx_list = list(ml_comed_df['CCM_CATNUM'].unique())
    # for comed_inx, comed_name in comed_dict.items():
    #     if comed_inx in uid_uniq_comed_inx_list:
    #         ml_row_dict[comed_name] = 1
    #     else:
    #         ml_row_dict[comed_name] = 0
    #
    # if len(ml_comed_df)==0:
    #     ml_row_dict['AT LEAST ONE NEPHROTOXIC AGENT'] = 0
    # else:
    #     ml_row_dict['AT LEAST ONE NEPHROTOXIC AGENT'] = int(ml_comed_df['NEPHTOX_DRUG_YN'].max())

    ## LAB
    try:
        uid_last_lab_row = ml_lab_df.fillna(method='ffill').fillna(method='bfill').iloc[-1]
    except:
        uid_last_lab_row = {c:np.nan for c in res_lab_cols}
    for lab_col in res_lab_cols[1:]:
        ml_row_dict[lab_col] = uid_last_lab_row[lab_col]

    # VITAL
    for vs_type in ['SBP','DBP','HR','BT']:
        for minmax_type in ['MIN','MAX']:
            vs_str = f"{vs_type}_{minmax_type}"
            try:
                ml_row_dict[vs_str] = ml_vs_df[ml_vs_df['VS_TYPE']==vs_str].iat[-1,-1]
            except:
                ml_row_dict[vs_str] = np.nan

    # # PROCEDURE
    # uid_uniq_proc_inx_list = list(ml_proc_df['PROC_CATNUM'].unique())
    # for proc_inx, proc_name in proc_dict.items():
    #     if proc_inx in uid_uniq_proc_inx_list:
    #         ml_row_dict[proc_name] = 1
    #     else:
    #         ml_row_dict[proc_name] = 0
    #
    # # MICROBIOLOGY
    # uid_uniq_micbio_inx_list = list(ml_micbio_df['MICBIO_CATNUM'].unique())
    # for micbio_inx, micbio_name in micbio_dict.items():
    #     if micbio_inx in uid_uniq_micbio_inx_list:
    #         ml_row_dict[micbio_name] = 1
    #     else:
    #         ml_row_dict[micbio_name] = 0

    # OUTCOME
    ml_row_dict['AKI_OCCURRENCE'] = dili_occurrence
    ml_row_dict['AKI_DATE'] = dili_occurrence_date
    ml_row_dict['BL_LAB_DATE'] = uid_last_lab_row['DATE']
    ml_row_dict['FIRST_DOSE_DATE'] = first_dose_dt
    ml_row_dict['LAST_DOSE_DATE'] = last_dose_dt

    # raise ValueError
    ml_res_df.append(ml_row_dict)
    # raise ValueError

ml_res_df = pd.DataFrame(ml_res_df)
# ml_res_df['UID']
# ml_res_df = pd.read_csv(f"{output_dir}/final_mlres_data.csv")
ml_res_df2 = ml_res_df[(ml_res_df['DURATION(DAYS)']>=1)].copy()
ml_res_df2.to_csv(f"{output_dir}/final_mlres_data.csv", index=False, encoding='utf-8-sig')
# ml_res_df['AT LEAST ONE NEPHROTOXIC AGENT'].unique()
# ml_res_df['AT LEAST ONE NEPHROTOXIC AGENT'] = ml_res_df['AT LEAST ONE NEPHROTOXIC AGENT'].map(lambda x: int(x) if not np.isnan(x) else np.nan)
print(f"Patients without lab data: {len(no_lab_pids)} / {set(no_lab_pids)}")
print(f"Patients without dose data: {len(no_dose_pids)} / {set(no_dose_pids)}")


# str(set(modeling_df['UID']) - set(ml_res_df2['UID'])).replace(', ',"', '").replace('{',"[").replace('}',"]")




# ml_res_df['AKI_OCCURRENCE'].sum()
# uid_demo_df.columns
# aki_df['UID'].drop_duplicates()
# aki_df[aki_df['AKI_DT'] > 24]['UID'].drop_duplicates()

# set(ml_res_df[ml_res_df[['WBC', 'ANC', 'Lymphocyte', 'RBC', 'Hb', 'PLT', 'Na', 'K', 'Ca', 'P', 'Cl']].sum(axis=1)==0]['UID'])