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
# lab_df[['Cr','eGFR','eGFR-Cockcroft-Gault','eGFR-Schwartz','eGFR-CKD-EPI','eGFR (CKD-EPI Cys)','eGFR (CKD-EPI Cr-Cys)','valproate_conc','lacosamide_conc']].isna().sum().sort_values()

adm_df = pd.read_excel(f"{output_dir}/merged_adm_dates.xlsx")

# adm_df.columns
no_sbase_lab_uids = list()
no_ibase_lab_uids = list()
yes_sbase_lab_uids = list()
yes_ibase_lab_uids = list()

no_sadm_lab_uids = list()
yes_sadm_lab_uids = list()
no_iadm_lab_uids = list()
yes_iadm_lab_uids = list()
# len(yes_sbase_lab_uids)
surv_res_df = list()
for inx, row in adm_df.iterrows(): #break
    uid = row['UID']
    uid_lab_df = lab_df[lab_df['UID']==uid].copy()

    uid_sbase_lab_df = uid_lab_df[uid_lab_df['DATE'] <= row['MIN_SINGLE_DATE']].copy()
    # uid_sbase_lab_df = uid_lab_df[uid_lab_df['DATE'] < row['MIN_SINGLE_DATE']].copy()
    uid_sadm_lab_df = uid_lab_df[(uid_lab_df['DATE'] >= row['MIN_SINGLE_DATE'])&(uid_lab_df['DATE'] <= row['MAX_SINGLE_DATE'])].copy()
    if len(uid_sbase_lab_df)==0:
        print(f"({inx}) {uid} / No sbase lab value")
        no_sbase_lab_uids.append(uid)
        continue
    else:
        yes_sbase_lab_uids.append(uid)

    if len(uid_sadm_lab_df)==0:
        print(f"({inx}) {uid} / No sadm lab value")
        no_sadm_lab_uids.append(uid)
        continue
    else:
        yes_sadm_lab_uids.append(uid)

    # raise ValueError
    bl_row = uid_sbase_lab_df[~uid_sbase_lab_df['Cr'].isna()].iloc[-1]
    tar_rows = uid_sadm_lab_df[uid_sadm_lab_df['Cr'] > 1.5 * bl_row['Cr']].copy()
    single_res_dict = {'ADM_TYPE':'SINGLE',
                       'DRUG':row['SINGLE_DRUG'],
                       'UID':uid,
                       'BL_DATE':bl_row['DATE'],
                       'BL_Cr':bl_row['Cr'],
                       'BL_eGFR':bl_row['eGFR'],
                       'BL_eGFR-Schwartz':bl_row['eGFR-Schwartz'],
                       'BL_eGFR-CKD-EPI':bl_row['eGFR-CKD-EPI'],
                       'FIRST_ADM_DATE':uid_sadm_lab_df['DATE'].iloc[0],
                       'TRT':0,
                       'EV': 0 if len(tar_rows)==0 else 1,
                       'EV_DATE': uid_sadm_lab_df['DATE'].iloc[-1] if len(tar_rows)==0 else tar_rows['DATE'].iloc[0],
                       'EV_Cr': np.nan if len(tar_rows)==0 else tar_rows['Cr'].iloc[0],
                       'EV_eGFR': np.nan if len(tar_rows)==0 else tar_rows['eGFR'].iloc[0],
                       'EV_eGFR-Schwartz': np.nan if len(tar_rows)==0 else tar_rows['eGFR-Schwartz'].iloc[0],
                       'EV_eGFR-CKD-EPI': np.nan if len(tar_rows)==0 else tar_rows['eGFR-CKD-EPI'].iloc[0],
                       }

    single_res_dict['time'] = (datetime.strptime(single_res_dict['EV_DATE'],'%Y-%m-%d')-datetime.strptime(uid_sadm_lab_df['DATE'].iloc[0],'%Y-%m-%d')).total_seconds()/86400
    single_res_dict['event'] = single_res_dict['EV']
    single_res_dict['group'] = single_res_dict['TRT']

    surv_res_df.append(single_res_dict)


    iadm_min_lab_date = (datetime.strptime(row['MIN_INTSEC_DATE'],'%Y-%m-%d')-timedelta(days=30)).strftime('%Y-%m-%d')
    uid_ibase_lab_df = uid_lab_df[(uid_lab_df['DATE'] <= row['MIN_INTSEC_DATE']) & (uid_lab_df['DATE'] > iadm_min_lab_date)].copy()
    # uid_ibase_lab_df = uid_lab_df[(uid_lab_df['DATE'] < row['MIN_INTSEC_DATE'])&(uid_lab_df['DATE'] > iadm_min_lab_date)].copy()
    uid_iadm_lab_df = uid_lab_df[(uid_lab_df['DATE'] >= row['MIN_INTSEC_DATE'])&(uid_lab_df['DATE'] <= row['MAX_INTSEC_DATE'])].copy()
    if len(uid_ibase_lab_df)==0:
        print(f"({inx}) {uid} / No ibase lab value")
        no_ibase_lab_uids.append(uid)
        continue
    else:
        yes_ibase_lab_uids.append(uid)

    if len(uid_iadm_lab_df)==0:
        print(f"({inx}) {uid} / No iadm lab value")
        no_iadm_lab_uids.append(uid)
        continue
    else:
        yes_iadm_lab_uids.append(uid)


    bl_row = uid_ibase_lab_df[~uid_ibase_lab_df['Cr'].isna()].iloc[-1]
    tar_rows = uid_ibase_lab_df[uid_ibase_lab_df['Cr'] >= 1.5 * bl_row['Cr']].copy()
    intersect_res_dict = {'ADM_TYPE':'INTSEC',
                          'DRUG':'BOTH',
                          'UID':uid,
                          'BL_DATE':bl_row['DATE'],
                          'BL_Cr':bl_row['Cr'],
                          'BL_eGFR':bl_row['eGFR'],
                          'BL_eGFR-Schwartz':bl_row['eGFR-Schwartz'],
                          'BL_eGFR-CKD-EPI':bl_row['eGFR-CKD-EPI'],
                          'FIRST_ADM_DATE':uid_iadm_lab_df['DATE'].iloc[0],
                          'TRT':1,
                          'EV': 0 if len(tar_rows)==0 else 1,
                          'EV_DATE': uid_iadm_lab_df['DATE'].iloc[-1] if len(tar_rows)==0 else tar_rows['DATE'].iloc[0],
                          'EV_Cr': np.nan if len(tar_rows)==0 else tar_rows['Cr'].iloc[0],
                          'EV_eGFR': np.nan if len(tar_rows)==0 else tar_rows['eGFR'].iloc[0],
                          'EV_eGFR-Schwartz': np.nan if len(tar_rows)==0 else tar_rows['eGFR-Schwartz'].iloc[0],
                          'EV_eGFR-CKD-EPI': np.nan if len(tar_rows)==0 else tar_rows['eGFR-CKD-EPI'].iloc[0],
                       }

    intersect_res_dict['time'] = (datetime.strptime(intersect_res_dict['EV_DATE'], '%Y-%m-%d') - datetime.strptime(uid_iadm_lab_df['DATE'].iloc[0], '%Y-%m-%d')).total_seconds() / 86400
    intersect_res_dict['event'] = intersect_res_dict['EV']
    intersect_res_dict['group'] = intersect_res_dict['TRT']

    surv_res_df.append(intersect_res_dict)

surv_res_df = pd.DataFrame(surv_res_df)


## Fisher's Exact Test

single_survres_df = surv_res_df[surv_res_df['ADM_TYPE']=='SINGLE'].copy()
intsec_survres_df = surv_res_df[surv_res_df['ADM_TYPE']=='INTSEC'].copy()


single_ev_count = single_survres_df['EV'].sum()
intsec_ev_count = intsec_survres_df['EV'].sum()

from scipy.stats import fisher_exact

# 예시 2×2 분할표 데이터
#            사건O   사건X
# 그룹A       8      2
# 그룹B       1      5
table = np.array([[intsec_ev_count, len(intsec_survres_df)-intsec_ev_count],
                  [single_ev_count, len(single_survres_df)-single_ev_count]])

# fisher_exact(table, alternative='two-sided') 가능
oddsratio, p_value = fisher_exact(table, alternative='two-sided')

print(table)
print(f"Odds Ratio: {oddsratio:.3f}")
print(f"P-value: {p_value:.5f}")

# len(yes_sbase_lab_uids)
# len(yes_ibase_lab_uids)
    # uid_lab_df[['UID','DATE','eGFR','Cr']]
    # row['UID']

## Survival analysis

# surv_res_df.columns['TIME']

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
from lifelines import CoxTimeVaryingFitter

kmf = KaplanMeierFitter()

plt.figure(figsize=(10,8))
for g, gdf in surv_res_df.groupby("group"):
    # kmf = KaplanMeierFitter()
    # kmf.fit(gdf["time"], event_observed=gdf["event"], label=str(g))
    # sf = kmf.survival_function_  # S(t)
    # ax.step(sf.index, 1 - sf.values.ravel(), where="post", label=f"{g}")  # 1 - S(t)

    kmf.fit(durations=gdf["time"], event_observed=gdf["event"], label=f"Group {g}")
    kmf.plot(ci_show=True)

plt.title("Kaplan–Meier Survival Curves by Group")
plt.xlabel("Time")
plt.ylabel("Survival probability")
plt.grid(True, linestyle="--")
plt.show()

# 로그랭크 검정 (A vs B)
A = surv_res_df[surv_res_df.group==1]; B = surv_res_df[surv_res_df.group==0]
results = logrank_test(A["time"], B["time"], event_observed_A=A["event"], event_observed_B=B["event"])
print(results.summary)

# Cox 비례위험모형 (정적 공변량) / 예시: 나이(연속), 성별(범주), 군(범주)

# cph = CoxPHFitter()
# cph.fit(df_model, duration_col="time", event_col="event")
# cph.print_summary()     # HR, CI, p-value 등
# cph.plot_partial_effects_on_outcome(covariates="age", values=[40, 60, 80])
# plt.show()
# cph.check_assumptions(df_model, p_value_threshold=0.05, show_plots=False)