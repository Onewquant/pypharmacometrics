from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

aki_df = pd.read_csv(f"{output_dir}/amk_aki.csv")
aki_df = aki_df.rename(columns={'ID':'UID'})
print(f"Total AKI cases: {len(aki_df['UID'].drop_duplicates())}")
except_aki_df = aki_df[aki_df['AKI_DT']<24].copy()
excepted_aki_ids = set(except_aki_df['UID'])
print(f"-(Patients with AKI < 24hr): {excepted_aki_ids}")
aki_df = aki_df[aki_df['AKI_DT']>=24].copy()
print(f"AKI cases (Filtered): {len(aki_df['UID'].drop_duplicates())}")

# aki_df.columns
# modeling_df.columns
modeling_df = pd.read_csv(f"{nonmem_dir}/amk_modeling_df_covar.csv")
modeling_df['DATE'] = modeling_df['DATETIME_ORI'].map(lambda x:x.split('T')[0])



lab_list_df = pd.read_csv(f"{output_dir}/amk_lablist_df.csv")
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

# for finx, fpath in enumerate(lab_files): #break
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

comed_df = pd.read_csv(f"{output_dir}/final_comed_df.csv")
# comed_df['DATE_LIST'].iloc[0]
comed_df['DATE_LIST'] = comed_df['DATE_LIST'].map(lambda x:x.replace("('",'').replace("',)",'').replace("')",'').split("', '"))
comed_df = comed_df.explode('DATE_LIST')
comed_df = comed_df.rename(columns={"DATE_LIST":'DATE'})
comed_dict = pd.read_excel(f"{resource_dir}/[AMK_AKI_ML_DATA]/CCM_index_table.xlsx")
comed_dict = comed_dict.set_index('CAT_NUM')['CCM_CAT']
comed_dict.index.name = None
comed_dict = dict(comed_dict)
# comed_df['DATE_LIST'].columns
cm_df = pd.read_csv(f"{output_dir}/final_comorbidity_df.csv")
cm_dict = pd.read_excel(f"{resource_dir}/[AMK_AKI_ML_DATA]/CM_index_table.xlsx")
cm_dict = cm_dict.set_index('CAT_NUM')['CM_CAT']
cm_dict.index.name = None
cm_dict = dict(cm_dict)

vs_df = pd.read_csv(f"{output_dir}/final_vs_data.csv")
# vs_df = vs_df.drop_duplicates(['UID','DATE','VS_TYPE'],keep='last')

proc_df = pd.read_csv(f"{output_dir}/final_procedure_df.csv")
proc_df['PROC_CATNUM'].unique()
proc_dict = {1:'CARDIAC SURGERY',2:'SURGERY',3:'ICU',4:'MECHANICAL VENT'}

micbio_df = pd.read_csv(f"{output_dir}/final_microbiology_df.csv")
micbio_dict = {1:"Pseudomonas aeruginosa",
             2:"Acinetobacter baumannii",
             3:"Klebsiella pneumoniae",
             4:"Escherichia coli",
             5:"Enterobacter species",
             6:"Serratia marcescens",
             7:"Mycobacterium tuberculosis",
             8:"NonTuberculous Mycobacterium",
             9:"Others",
             10:"Empirical"}

losmot_df = pd.read_csv(f"{output_dir}/final_locmotality_df.csv")

ml_df = modeling_df.drop_duplicates(['UID'])[['UID']].reset_index(drop=True)
ml_res_df = list()
no_lab_pids = list()
no_dose_pids = list()
for inx, row in ml_df.iterrows():
    # raise ValueError
    ml_row_dict = dict()
    uid = row['UID']
    print(f"({inx}) {uid}")


    # 개별 dataframe 로드
    uid_aki_rows = aki_df[aki_df['UID'] == uid].copy()

    uid_modeling_df = modeling_df[modeling_df['UID']==uid].copy()
    uid_dose_df = uid_modeling_df[(uid_modeling_df['MDV'] > 0)&(uid_modeling_df['DV']=='.')].copy()
    uid_dose_df['AMT'] = uid_dose_df['AMT'].astype(float)
    # uid_demo_last_row = uid_demo_df.iloc[-1]
    lab_fpath = glob.glob(f'{resource_dir}/lab\\AMK_lab({uid}_*).xlsx')[0]
    uid_lab_df = pd.read_excel(lab_fpath)
    uid_lab_df['DATE'] = uid_lab_df[['보고일', '오더일']].max(axis=1)
    uid_lab_df['LAB'] = uid_lab_df['검사명']
    uid_lab_df['VALUE'] = uid_lab_df['검사결과']
    uid_lab_df['VALUE'] = pd.to_numeric(uid_lab_df['VALUE'], errors='coerce')  # 숫자 아닌 랩은 NAN으로 변환
    uid_lab_df = uid_lab_df.drop_duplicates(['DATE', 'LAB'], keep='last', ignore_index=True)
    uid_lab_df = uid_lab_df.drop(['보고일', '오더일','검사명','검사결과', '참고치', '직전결과'],axis=1)


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

    uid_comed_df = comed_df[comed_df['UID']==uid].copy()
    uid_cm_df = cm_df[cm_df['UID']==uid].copy()
    uid_vs_df = vs_df[vs_df['UID']==uid].copy()
    uid_proc_df = proc_df[proc_df['UID']==uid].copy()
    uid_micbio_df = micbio_df[micbio_df['UID']==uid].copy()
    uid_losmot_df = losmot_df[losmot_df['UID']==uid].copy()

    # AKI 발생여부 기록

    # AKI 발생 한 사람 (AKI 발생 시점이 첫 dosing 24시간 이후 일때만 포함되어 있음)
    if len(uid_aki_rows):
        aki_occurrence = 1
        aki_occurrence_time = uid_aki_rows.iloc[0]['AKI_DT']
        aki_occurrence_dt = uid_aki_rows.iloc[0]['AKI_DATETIME']
        aki_occurrence_date = aki_occurrence_dt.split('T')[0]
    else:
        aki_occurrence = 0
        aki_occurrence_time = np.nan
        aki_occurrence_dt = uid_modeling_df.iloc[-1]['DATETIME_ORI']
        aki_occurrence_date = aki_occurrence_dt.split('T')[0]

    # BASELINE Features 정리
    min30_bl_date = (datetime.strptime(aki_occurrence_date, '%Y-%m-%d') - timedelta(days=30)).strftime('%Y-%m-%d')
    min14_bl_date = (datetime.strptime(aki_occurrence_date, '%Y-%m-%d') - timedelta(days=14)).strftime('%Y-%m-%d')
    min7_bl_date = (datetime.strptime(aki_occurrence_date, '%Y-%m-%d') - timedelta(days=7)).strftime('%Y-%m-%d')
    min3_bl_date = (datetime.strptime(aki_occurrence_date, '%Y-%m-%d') - timedelta(days=3)).strftime('%Y-%m-%d')
    min_bl_date = (datetime.strptime(aki_occurrence_date, '%Y-%m-%d') - timedelta(days=3)).strftime('%Y-%m-%d')
    max_bl_date = (datetime.strptime(aki_occurrence_date, '%Y-%m-%d') - timedelta(days=1)).strftime('%Y-%m-%d')

    # 데이터 중간 처리
    uid_cm_df['CM_PERIOD_FROM_DATE'] = uid_cm_df['CM_PERIOD_FROM'].map(lambda x: (datetime.strptime(max_bl_date, '%Y-%m-%d') - timedelta(days=x)).strftime('%Y-%m-%d'))
    uid_proc_df['CM_PERIOD_FROM_DATE'] = uid_proc_df['PROC_PERIOD_FROM'].map(lambda x: (datetime.strptime(max_bl_date, '%Y-%m-%d') - timedelta(days=x)).strftime('%Y-%m-%d'))

    # ml feature 기록용 데이터 필터링
    ml_demo_df = uid_modeling_df[(uid_modeling_df['DATE'] <= aki_occurrence_date)&(uid_modeling_df['DATE'] >= min7_bl_date)].copy()
    ml_demo_df_last_row = ml_demo_df.iloc[-1]

    # ml_lab_df = uid_lab_df[(uid_lab_df['DATE'] <= max_bl_date)&(uid_lab_df['DATE'] >= min_bl_date)].copy()
    ml_lab_df = uid_lab_df[(uid_lab_df['DATE'] <= max_bl_date)&(uid_lab_df['DATE'] >= min14_bl_date)].copy()
    ml_comed_df = uid_comed_df[(uid_comed_df['DATE'] <= max_bl_date)&(uid_comed_df['DATE'] >= min_bl_date)].copy()
    ml_cm_df = uid_cm_df[(uid_cm_df['DATE'] >= uid_cm_df['CM_PERIOD_FROM_DATE'])&(uid_cm_df['DATE'] <= max_bl_date)].copy()
    ml_vs_df = uid_vs_df[(uid_vs_df['DATE'] <= max_bl_date)&(uid_vs_df['DATE'] >= min30_bl_date)].copy()
    ml_proc_df = uid_proc_df[(uid_proc_df['DATE'] >= uid_proc_df['CM_PERIOD_FROM_DATE'])&(uid_proc_df['DATE'] <= max_bl_date)].copy()
    ml_micbio_df = uid_micbio_df[(uid_micbio_df['DATE'] <= max_bl_date)&(uid_micbio_df['DATE'] >= min14_bl_date)].copy()
    uid_losmot_df = uid_micbio_df[(uid_micbio_df['DATE'] <= max_bl_date)&(uid_micbio_df['DATE'] >= min14_bl_date)].copy()

    ml_dose_df = uid_dose_df[(uid_dose_df['DATE'] <= max_bl_date)].copy()
    if len(ml_dose_df)==0:
        no_dose_pids.append(uid)
        print(f"/ No dose data")
        continue

    # first_dose_dt = ml_dose_df.iloc[0]['DATETIME_ORI']
    # last_dose_dt = ml_dose_df.iloc[-1]['DATETIME_ORI']
    first_dose_dt = ml_dose_df.iloc[0]['DATE']
    last_dose_dt = ml_dose_df.iloc[-1]['DATE']

    ## ml feature 기록
    ml_row_dict['UID'] = uid

    # Demographics
    for c in ['ID','AGE','SEX','WT','HT']:
        ml_row_dict[c] = ml_demo_df_last_row[c]
    ml_row_dict['BMI'] = round(ml_row_dict['WT']/((ml_row_dict['HT']/100)**2),3)

    # DOSE
    ml_row_dict['LAST_DOSE'] = ml_dose_df.iloc[-1]['AMT']
    ml_row_dict['LAST_DAILY_DOSE'] = ml_dose_df.groupby('DATE').agg({'AMT':'sum'}).iloc[-1]['AMT']
    ml_row_dict['LAST_DOSE_perWT'] = ml_row_dict['LAST_DOSE']/ml_row_dict['WT']
    ml_row_dict['LAST_DAILY_DOSE_perWT'] = ml_row_dict['LAST_DAILY_DOSE']/ml_row_dict['WT']
    ml_row_dict['DURATION(DAYS)'] = (ml_dose_df.iloc[-1]['TIME'] - ml_dose_df.iloc[0]['TIME'] + 1)/24
    ml_row_dict['CUM_DOSE'] = ml_dose_df['AMT'].sum()
    ml_row_dict['MEAN_DAILY_DOSE'] = ml_row_dict['CUM_DOSE']/ml_row_dict['DURATION(DAYS)']
    ml_row_dict['MEAN_DAILY_DOSE_perWT'] = ml_row_dict['CUM_DOSE']/ml_row_dict['WT']

    # COMORBIDITY
    uid_uniq_cm_inx_list = list(ml_cm_df['CM_CATNUM'].unique())
    for cm_inx, cm_name in cm_dict.items():
        if cm_inx in uid_uniq_cm_inx_list:
            ml_row_dict[cm_name] = 1
        else:
            ml_row_dict[cm_name] = 0

    # CONCOMITANT MEDICATION
    # raise ValueError
    # ml_comed_df.columns
    uid_uniq_comed_inx_list = list(ml_comed_df['CCM_CATNUM'].unique())
    for comed_inx, comed_name in comed_dict.items():
        if comed_inx in uid_uniq_comed_inx_list:
            ml_row_dict[comed_name] = 1
        else:
            ml_row_dict[comed_name] = 0

    if len(ml_comed_df)==0:
        ml_row_dict['AT LEAST ONE NEPHROTOXIC AGENT'] = 0
    else:
        ml_row_dict['AT LEAST ONE NEPHROTOXIC AGENT'] = int(ml_comed_df['NEPHTOX_DRUG_YN'].max())

    # LAB
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

    # PROCEDURE
    uid_uniq_proc_inx_list = list(ml_proc_df['PROC_CATNUM'].unique())
    for proc_inx, proc_name in proc_dict.items():
        if proc_inx in uid_uniq_proc_inx_list:
            ml_row_dict[proc_name] = 1
        else:
            ml_row_dict[proc_name] = 0

    # MICROBIOLOGY
    uid_uniq_micbio_inx_list = list(ml_micbio_df['MICBIO_CATNUM'].unique())
    for micbio_inx, micbio_name in micbio_dict.items():
        if micbio_inx in uid_uniq_micbio_inx_list:
            ml_row_dict[micbio_name] = 1
        else:
            ml_row_dict[micbio_name] = 0

    # OUTCOME
    ml_row_dict['AKI_OCCURRENCE'] = aki_occurrence
    ml_row_dict['AKI_DATE'] = aki_occurrence_date
    ml_row_dict['BL_LAB_DATE'] = uid_last_lab_row['DATE']
    ml_row_dict['FIRST_DOSE_DATE'] = first_dose_dt
    ml_row_dict['LAST_DOSE_DATE'] = last_dose_dt

    # raise ValueError
    ml_res_df.append(ml_row_dict)
    # raise ValueError

ml_res_df = pd.DataFrame(ml_res_df)
ml_res_df.to_csv(f"{output_dir}/final_mlres_data.csv", index=False, encoding='utf-8-sig')
# ml_res_df['AT LEAST ONE NEPHROTOXIC AGENT'].unique()
# ml_res_df['AT LEAST ONE NEPHROTOXIC AGENT'] = ml_res_df['AT LEAST ONE NEPHROTOXIC AGENT'].map(lambda x: int(x) if not np.isnan(x) else np.nan)
print(f"Patients without lab data: {len(no_lab_pids)} / {set(no_lab_pids)}")
print(f"Patients without dose data: {len(no_dose_pids)} / {set(no_dose_pids)}")

# ml_res_df['AKI_OCCURRENCE'].sum()
# uid_demo_df.columns
# aki_df['UID'].drop_duplicates()
# aki_df[aki_df['AKI_DT'] > 24]['UID'].drop_duplicates()