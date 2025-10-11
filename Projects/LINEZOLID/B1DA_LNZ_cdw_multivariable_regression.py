from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

prj_name = 'LNZ'
prj_dir = './Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"



# dose_df['DAILY_DOSE'] = dose_df['DOSE'] * 24/dose_df['INTERVAL'].map(lambda x:float(x.replace('q','').replace('h','')))
# dose_df['DOSAGE_REGIMEN'] = (dose_df['DOSE'].astype(str)+'_'+dose_df['INTERVAL'])
#
# str(list(dose_df['INTERVAL'].drop_duplicates().sort_values().reset_index(drop=True))).replace("', '",", ")
# str(list(dose_df['DOSE'].astype(int).drop_duplicates().sort_values().reset_index(drop=True))).replace("', '",", ")
# str(list(dose_df['DAILY_DOSE'].astype(int).drop_duplicates().sort_values().reset_index(drop=True))).replace("', '",", ")
# dose_df['DOSAGE_REGIMEN'].drop_duplicates().reset_index(drop=True)

ptinfo_df = pd.read_csv(f"{output_dir}/lnz_final_ptinfo_df.csv")
lab_df = pd.read_csv(f"{output_dir}/lnz_final_lab_df.csv")
# lab_df = pd.read_csv(f"{output_dir}/lnz_final_lab_df2.csv")
bodysize_df = pd.read_csv(f"{output_dir}/lnz_final_bodysize_df.csv")
dose_df = pd.read_csv(f"{output_dir}/lnz_final_dose_df.csv")
surv_res_df = pd.read_csv(f"{output_dir}/b1da_lnz_surv_res_df.csv")
dose_df['UID'].drop_duplicates()
# surv_res_df['UID'].drop_duplicates()
# surv_res_df['']
# dose_df[dose_df['UID']==155505674746153]
# lab_df['ANC']

# lab_df 처리
raw_pd_list = ['ANC (em)', 'Hct', 'Hct (em)', 'Hct(i-STAT)', 'Hematocrit (ICU용)', 'Hematocrit (마취과용)', 'Hematocrit (응급실용)', 'PCT', 'PDW', 'PLT (em)', 'PT (%)', 'PT (INR)', 'PT (INR) (em)', 'PT (sec)', 'Plasma cell', 'Lymphocyte', 'MCH', 'MCHC', 'MCV', 'MPV', 'Metamyelocyte', 'Mixing test (PT, aPTT 제외)', 'Monocyte', 'Myelocyte', 'Normoblast', 'Other', 'Band neutrophil', 'Basophil',  'CBC (em) (differential count) RDW제외', 'Blast', 'Eosinophil', 'Eosinophil count',  '절대단구수', '절대림프구수', 'Promyelocyte', 'Prothrombin time (%) : MIX', 'Prothrombin time (INR) : MIX', 'Prothrombin time (sec) : MIX', 'RBC', 'RBC (em)', 'RDW(CV)', 'RDW(SD)', 'Reticulocyte', 'Segmented neutrophil', 'WBC (em)', 'WBC stick', 'WBC stick (em)','Hb (em)', 'Hb(i-STAT)', 'Hemoglobin (ICU용)', 'Hemoglobin (마취과용)', 'Hemoglobin (응급실용)', 'Immature cell', 'Joint WBC stick (em)', 'Activated PTT : MIX', 'Atypical lymphocyte',  'BE', 'BE(i-STAT)','aPTT', 'aPTT (em)',]
pd_list = ['ANC','WBC',  'Hb', 'PLT', 'Lactate', 'pH']
etc_list1 = ['Creatinine (random urine)', 'Creatinine (urine, em)', 'Creatinine Clerarance (24hrs urine, Ccr)','Creatinine (24hrs urine) (g/day)', 'Creatinine (24hrs urine) (mg/dL)','Glucose (간이혈당기용)', 'Urine creatinine (24hrs urine)' ]
etc_list2 = ['FBS (serum)', 'Fibrinogen', 'Microalbumin (random urine)', 'Microalbumin/Creatinine ratio', 'O₂CT', 'O₂SAT', 'O₂SAT(i-STAT)', 'Potasium(i-STAT)', 'Protein (random urine)', 'Protein/Creatinine ratio',  'Sodium(i-STAT)', 'TotalCO₂(i-STAT)',  'i- Calcium (i-STAT)', 'pCO₂', 'pCO₂(i-STAT)','pO₂', 'pO₂(i-STAT)',]

raw_pd_list = ['ANC (em)', 'Hct', 'Hct (em)', 'Hct(i-STAT)', 'Hematocrit (ICU용)', 'Hematocrit (마취과용)', 'Hematocrit (응급실용)', 'PCT', 'PDW', 'PLT (em)', 'PT (%)', 'PT (INR)', 'PT (sec)', 'Plasma cell', 'Lymphocyte', 'MCH', 'MCHC', 'MCV', 'MPV', 'Metamyelocyte', 'Mixing test (PT, aPTT 제외)', 'Monocyte', 'Myelocyte', 'Normoblast', 'Other', 'Band neutrophil', 'Basophil',  'CBC (em) (differential count) RDW제외', 'Blast', 'Eosinophil', 'Eosinophil count',  '절대단구수', '절대림프구수', 'Promyelocyte', 'Prothrombin time (%) : MIX', 'Prothrombin time (INR) : MIX', 'Prothrombin time (sec) : MIX', 'RBC', 'RBC (em)', 'RDW(CV)', 'RDW(SD)', 'Reticulocyte', 'Segmented neutrophil', 'WBC (em)', 'Hb (em)', 'Hb(i-STAT)', 'Hemoglobin (ICU용)', 'Hemoglobin (마취과용)', 'Hemoglobin (응급실용)', 'Immature cell', 'Activated PTT : MIX', 'Atypical lymphocyte',  'BE', 'BE(i-STAT)','aPTT', 'aPTT (em)',]
etc_list1 = []
etc_list2 = ['Fibrinogen', 'O₂CT', 'O₂SAT', 'O₂SAT(i-STAT)', 'Potasium(i-STAT)',  'Sodium(i-STAT)', 'TotalCO₂(i-STAT)', 'i- Calcium (i-STAT)', 'pCO₂', 'pCO₂(i-STAT)','pO₂', 'pO₂(i-STAT)',]

lab_df['ANC'] = lab_df[['ANC', 'ANC (em)']].min(axis=1)
lab_df['PLT'] = lab_df[['PLT', 'PLT (em)']].min(axis=1)
lab_df['WBC'] = lab_df[['WBC', 'WBC (em)', ]].min(axis=1)
lab_df['Hb'] = lab_df[['Hb', 'Hb (em)', ]].min(axis=1)

lab_df.drop(raw_pd_list,axis=1)
lab_df = lab_df.drop(raw_pd_list+etc_list1+etc_list2,axis=1)

lab_df['ALT'] = lab_df[['ALT(GPT) (em)','GPT (ALT)']].max(axis=1)
lab_df['AST'] = lab_df[['AST(GOT) (em)','GOT (AST)']].max(axis=1)
lab_df['GGT'] = lab_df[['GGT', 'GGT (em)']].max(axis=1)
lab_df['ALB'] = lab_df[['Albumin', 'Albumin (em)']].max(axis=1)
lab_df['SCR'] = lab_df[['Creatinine', 'Creatinine  (응급실용)', 'Creatinine (ICU용)','Creatinine (serum)', 'Creatinine (serum, em)', ]].max(axis=1)
lab_df['GLU'] = lab_df[['Glucose', 'Glucose (ICU용)', 'Glucose (em)', 'Glucose (마취과용)', 'Glucose (응급실용)']].max(axis=1)
lab_df['CRP'] = lab_df[['CRP (em)', 'CRP (hsCRP)', 'hsCRP (em)',]].max(axis=1)
lab_df['TPRO'] = lab_df[['Protein, total', 'Protein, total (em)',]].max(axis=1)
lab_df['TBIL'] = lab_df[['Bilirubin, total', 'Bilirubin, total  (em)', 'Bilirubin, total (NICU)',]].max(axis=1)
lab_df['Lactate'] = lab_df[['Lactate (ICU용)', 'Lactate (마취과용)', 'Lactate (응급실용)', 'Lactate, Lactic acid (em)', ]].max(axis=1)
lab_df['PROCAL'] = lab_df[['Procalcitonin', 'Procalcitonin 정량', ]].max(axis=1)
lab_df['eGFR'] = lab_df[['eGFR', ]].min(axis=1)
lab_df['eGFR-CKD-EPI'] = lab_df[['eGFR-CKD-EPI']].min(axis=1)
lab_df['eGFR-Schwartz'] = lab_df[['eGFR-Schwartz(소아)']].min(axis=1)
lab_df['eGFR-Cockcroft-Gault'] = lab_df[['eGFR-Cockcroft-Gault']].min(axis=1)
lab_df['HCO3-'] = lab_df[['HCO3-', 'HCO3-(i-STAT)', ]].min(axis=1)
lab_df['DBIL'] = lab_df[['Bilirubin, direct', ]].min(axis=1)
lab_df['ESR'] = lab_df['ESR']
lab_df['pH'] = lab_df[['pH', 'pH (i-STAT)']].min(axis=1)




# list(surv_res_df['ENDPOINT'].unique())

# uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATE'], how='left').fillna(method='ffill')

# lab_covar_list = ['ALT','AST','GGT','ALB','SCR','GLU','CRP','TPRO','TBIL','LACT','PROCAL','pH','DBIL','ESR']
lab_covar_list = ['ALT','AST','GGT','ALB','SCR','GLU','CRP','TPRO','TBIL','PROCAL','DBIL','ESR','eGFR','eGFR-CKD-EPI','HCO3-']
lab_df = lab_df[['UID', 'DATE']+lab_covar_list+pd_list].copy()
# lab_df.columns
# 정보량 반영 (상대적인 평가로, Max 정보량 가진 column 기준으로 실질 정보의 row 수 %가 50 이상인 컬럼만 남김)

lab_col_inclusion = (100*(~lab_df.isna()).sum()/((~lab_df.isna()).sum().iloc[2:].max())).sort_values(ascending=False)
# lab_col_inclusion = (100*(~lab_df.isna()).sum()/(~lab_df.isna())['SCR'].sum()).sort_values(ascending=False)


lab_covar_list = list(lab_col_inclusion[lab_col_inclusion >= 50].index)[2:] + ['Lactate','pH']
lab_df = lab_df[['UID', 'DATE']+lab_covar_list].copy()
# str(lab_covar_list).replace("', '",', ')

full_result_df = list()
count = 0
for uid, uid_df in lab_df.groupby('UID',as_index=False): #break
    uid_df = uid_df.fillna(method='ffill')
    full_result_df.append(uid_df)
lab_df = pd.concat(full_result_df)
# list(lab_df.columns)[2:]

# ep_surv_df.columns
# surv_res_df['ENDPOINT'].unique()
ep_res_dict = dict()
for endpoint, ep_surv_df in surv_res_df.groupby('ENDPOINT'):
    ep_res_df = list()

    if endpoint=='PLT':
        ep_lab_covar_list = lab_covar_list
    elif endpoint=='WBC':
        ep_lab_covar_list = lab_covar_list
    elif endpoint=='Hb':
        ep_lab_covar_list = lab_covar_list
    elif endpoint=='Lactate':
        ep_lab_covar_list = lab_covar_list
    else:
        ep_lab_covar_list = lab_covar_list

    for inx, row in ep_surv_df.iterrows():
        # raise ValueError
        uid = row['UID']
        bl_date = row['BL_DATE']
        bl_date_bf1mo = (datetime.strptime(bl_date,'%Y-%m-%d')-timedelta(days=30)).strftime('%Y-%m-%d')
        ev_date = row['EV_DATE']
        ev = row['EV']

        res_dict = {'UID':uid}

        # patient info
        ptinfo_row = ptinfo_df[ptinfo_df['UID']==uid].iloc[0]
        for ptcol in ['SEX','BIRTH_DATE']:
            if ptcol == 'BIRTH_DATE':
                res_dict['AGE'] = ptinfo_row[ptcol]
            else:
                res_dict[ptcol] = ptinfo_row[ptcol]

        # lab
        uid_lab_df = lab_df[lab_df['UID']==uid].copy()
        uid_bl_lab_df = uid_lab_df[(uid_lab_df['DATE'] >= bl_date_bf1mo) & (uid_lab_df['DATE'] <= bl_date)].fillna(method='ffill')
        bl_lab_row = uid_bl_lab_df.iloc[-1]
        # col='DATE'
        for col in ['DATE']+ep_lab_covar_list:
            if col=='DATE':
                res_dict[f'BL_{col}'] = bl_lab_row[col]
                res_dict['AGE'] = int((datetime.strptime(bl_lab_row[col],'%Y-%m-%d') - datetime.strptime(res_dict['AGE'],'%Y-%m-%d')).days/365.25)
                res_dict['ELD'] = 1 if res_dict['AGE']>=65 else 0
            else:
                res_dict[col] = bl_lab_row[col]

        # bodysize
        uid_bs_df = bodysize_df[bodysize_df['UID']==uid].copy()
        uid_bl_bs_df = uid_bs_df[(uid_bs_df['DATE'] >= bl_date_bf1mo) & (uid_bs_df['DATE'] <= bl_date)].copy()
        if len(uid_bl_bs_df)==0:
            bl_bs_row = {'DATE':bl_date,'HT':np.nan,'WT':np.nan,'BMI':np.nan}
        else:
            uid_bl_bs_df = uid_bl_bs_df.pivot_table(index="DATE", columns="BODYSIZE", values="VALUE", aggfunc="median").reset_index(drop=False).sort_values(['DATE'],ignore_index=True).fillna(method='ffill')
            uid_bl_bs_df.columns.name = None
            if len(uid_bl_bs_df.drop(['DATE'],axis=1).columns)<=1:
                for bs_col in ['WT','HT','BMI']:
                    if bs_col in uid_bl_bs_df.columns:
                        continue
                    else:
                        uid_bl_bs_df[bs_col]=np.nan
            if 'HT' not in uid_bl_bs_df.columns:
                uid_bl_bs_df['HT'] = (np.sqrt(uid_bl_bs_df['WT']/uid_bl_bs_df['BMI'])*100).map(lambda x:round(x,1))
            if 'WT' not in uid_bl_bs_df.columns:
                uid_bl_bs_df['WT'] = (np.square(uid_bl_bs_df['HT']/100) * uid_bl_bs_df['BMI']).map(lambda x:round(x,1))
            if 'BMI' not in uid_bl_bs_df.columns:
                uid_bl_bs_df['BMI'] = (uid_bl_bs_df['WT'] / np.square(uid_bl_bs_df['HT']/100)).map(lambda x:round(x,1))
            uid_bl_bs_df = uid_bl_bs_df[['DATE','HT','WT','BMI']].copy()
            bl_bs_row = uid_bl_bs_df.iloc[-1]

        for col in ['HT','WT','BMI']:
            res_dict[col] = bl_bs_row[col]

        # dose
        uid_dose_df = dose_df[dose_df['UID']==uid].copy()
        uid_cum_dose_df = uid_dose_df[(uid_dose_df['DATE'] >= bl_date) & (uid_dose_df['DATE'] <= ev_date)].copy()
        # if len(uid_cum_dose_df)>1:
        #     raise ValueError

        dosing_period = ((datetime.strptime(uid_cum_dose_df['DATE'].iloc[-1], '%Y-%m-%d') - datetime.strptime(uid_cum_dose_df['DATE'].iloc[0], '%Y-%m-%d')).days + 1)
        daily_dose = (uid_cum_dose_df['DOSE'].sum()) / dosing_period
        active_dosing_days = len(uid_cum_dose_df['DATE'].unique())
        active_daily_dose = (uid_cum_dose_df['DOSE'].sum()) / active_dosing_days


        total_dosing_period = ((datetime.strptime(uid_dose_df['DATE'].iloc[-1],'%Y-%m-%d') - datetime.strptime(uid_dose_df['DATE'].iloc[0],'%Y-%m-%d')).days + 1)
        total_daily_dose = (uid_dose_df['DOSE'].sum())/total_dosing_period

        total_active_dosing_days = len(uid_dose_df['DATE'].unique())
        total_active_daily_dose = (uid_dose_df['DOSE'].sum()) / total_active_dosing_days


        res_dict['CUM_DOSE'] = (uid_cum_dose_df['DOSE'].sum())
        res_dict['CUM_DOSE(TOTAL)'] = (uid_dose_df['DOSE'].sum())
        res_dict['DOSE_PERIOD'] = dosing_period
        res_dict['DOSE_PERIOD(ACTIVE)'] = active_dosing_days
        res_dict['DOSE_PERIOD(TOTAL)'] = total_dosing_period
        res_dict['DOSE24'] = daily_dose
        res_dict['DOSE24(ACTIVE)'] = active_daily_dose
        res_dict['DOSE24(TOTAL_ACTIVE)'] = total_active_daily_dose
        res_dict['DOSE24(TOTAL_PERIOD)'] = total_daily_dose
        res_dict['DOSE24PERWT'] = (daily_dose / res_dict['WT'])
        res_dict['DOSE24PERWT(ACTIVE)'] = (active_daily_dose / res_dict['WT'])
        res_dict['DOSE24PERWT(TOTAL_ACTIVE)'] = (total_active_daily_dose / res_dict['WT'])
        res_dict['DOSE24PERWT(TOTAL_PERIOD)'] = (total_daily_dose / res_dict['WT'])
        res_dict['DOSE_INTERVAL'] = uid_cum_dose_df['INTERVAL'].map(lambda x:float(x.replace('q','').replace('h',''))).mean()
        res_dict['DOSE_INTERVAL(TOTAL)'] = uid_dose_df['INTERVAL'].map(lambda x:float(x.replace('q','').replace('h',''))).mean()
        # raise ValueError
        res_dict['EV'] = ev
        ep_res_df.append(res_dict)
    # raise ValueError
    ep_res_df = pd.DataFrame(ep_res_df)
    ep_res_df['ENDPOINT'] = endpoint

    ep_res_dict[endpoint] = ep_res_df.copy()

for endpoint in surv_res_df['ENDPOINT'].unique():
    # endpoint = 'PLT'
    # endpoint = 'Hb'
    # endpoint = 'WBC'
    # endpoint = 'ANC'

    epreg_df = ep_res_dict[endpoint].copy()
    if not os.path.exists(f'{output_dir}/b1da'):
        os.mkdir(f'{output_dir}/b1da')
    if not os.path.exists(f'{output_dir}/b1da/datacheck'):
        os.mkdir(f'{output_dir}/b1da/datacheck')
    epreg_df.to_csv(f"{output_dir}/b1da/datacheck/b1da_lnz_datacheck_{endpoint}_df.csv", encoding='utf-8-sig', index=False)
    epreg_df = epreg_df.drop(['BL_DATE','UID','ENDPOINT'],axis=1)
    epreg_df['SEX']=epreg_df['SEX'].map({'남':0,'여':1})
    epreg_df = epreg_df.fillna(epreg_df.median(numeric_only=True))
    if not os.path.exists(f'{output_dir}/b1da/mvlreg'):
        os.mkdir(f'{output_dir}/b1da/mvlreg')
    epreg_df.to_csv(f"{output_dir}/b1da/mvlreg/b1da_lnz_mvlreg_{endpoint}_df.csv", encoding='utf-8-sig', index=False)
    print(f'lnz_mvlreg_{endpoint}_df.csv / generated')
    # ep_res_dict.keys()