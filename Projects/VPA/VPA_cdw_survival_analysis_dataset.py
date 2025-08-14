from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

prj_name = 'VPA'
prj_dir = './Projects/VPA'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

lab_df = pd.read_csv(f"{output_dir}/vpa_final_lab_df.csv")
lab_df = lab_df.rename(columns={'eGFR-Schwartz(소아)':'eGFR-Schwartz','Valproic Acid 농도':'valproate_conc','Lacosamide 농도':'lacosamide_conc'})
lab_df['Cr'] = lab_df[['Creatinine','Creatinine  (응급실용)','Creatinine (serum)','Creatinine (serum, em)']].max(axis=1)
lab_df = lab_df[['Cr','eGFR','eGFR-Cockcroft-Gault','eGFR-Schwartz','eGFR-CKD-EPI','eGFR (CKD-EPI Cys)','eGFR (CKD-EPI Cr-Cys)','valproate_conc','lacosamide_conc']].copy()
# lab_df.columns
# lab_df = lab_df[['Cr','eGFR','eGFR-Cockcroft-Gault','eGFR-Schwartz','eGFR-CKD-EPI','eGFR (CKD-EPI Cys)','eGFR (CKD-EPI Cr-Cys)','valproate_conc','lacosamide_conc']].isna().sum()
# lab_df[lab_df['eGFR'].isna()]

lab_df = pd.read_csv(f"{output_dir}/vpa_final_lab_df.csv")
lab_df = lab_df.rename(columns={'eGFR-Schwartz(소아)':'eGFR-Schwartz','Valproic Acid 농도':'valproate_conc','Lacosamide 농도':'lacosamide_conc'})
lab_df['Cr'] = lab_df[['Creatinine','Creatinine  (응급실용)','Creatinine (serum)','Creatinine (serum, em)']].max(axis=1)
lab_df = lab_df[['UID','DATE','Cr','eGFR','eGFR-Cockcroft-Gault','eGFR-Schwartz','eGFR-CKD-EPI','eGFR (CKD-EPI Cys)','eGFR (CKD-EPI Cr-Cys)','valproate_conc','lacosamide_conc']].copy()

adm_df = pd.read_excel(f"{output_dir}/merged_adm_dates.xlsx")

# adm_df.columns
no_sbase_lab_uids = list()
no_ibase_lab_uids = list()
yes_sbase_lab_uids = list()
yes_ibase_lab_uids = list()
# len(yes_sbase_lab_uids)
surv_res_df = list()
for inx, row in adm_df.iterrows(): #break
    uid = row['UID']
    uid_lab_df = lab_df[lab_df['UID']==uid].copy()

    uid_sbase_lab_df = uid_lab_df[uid_lab_df['DATE'] <= row['MIN_SINGLE_DATE']].copy()
    uid_sadm_lab_df = uid_lab_df[(uid_lab_df['DATE'] >= row['MIN_SINGLE_DATE'])&(uid_lab_df['DATE'] <= row['MAX_SINGLE_DATE'])].copy()
    if len(uid_sbase_lab_df)==0:
        print(f"({inx}) {uid} / No sbase lab value")
        no_sbase_lab_uids.append(uid)
        continue
    else:
        yes_sbase_lab_uids.append(uid)

    {'UID':uid, ''}


    iadm_min_lab_date = (datetime.strptime(row['MIN_INTSEC_DATE'],'%Y-%m-%d')-timedelta(days=30)).strftime('%Y-%m-%d')
    uid_ibase_lab_df = uid_lab_df[(uid_lab_df['DATE'] <= row['MIN_INTSEC_DATE'])&(uid_lab_df['DATE'] > iadm_min_lab_date)].copy()
    uid_iadm_lab_df = uid_lab_df[(uid_lab_df['DATE'] >= row['MIN_INTSEC_DATE'])&(uid_lab_df['DATE'] <= row['MAX_INTSEC_DATE'])].copy()
    if len(uid_ibase_lab_df)==0:
        print(f"({inx}) {uid} / No ibase lab value")
        no_ibase_lab_uids.append(uid)
        continue
    else:
        yes_ibase_lab_uids.append(uid)


    surv_res_df.append()


# len(yes_sbase_lab_uids)
# len(yes_ibase_lab_uids)
    # uid_lab_df[['UID','DATE','eGFR','Cr']]
    # row['UID']