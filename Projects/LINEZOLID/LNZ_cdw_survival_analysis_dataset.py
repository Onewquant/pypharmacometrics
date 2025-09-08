from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

prj_name = 'LNZ'
prj_dir = './Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

lab_df = pd.read_csv(f"{output_dir}/lnz_final_lab_df.csv")
# lab_df.columns
# lab_df[['UID', 'DATE', 'ANC', 'Atypical lymphocyte', 'Band neutrophil', 'Basophil', 'Blast', 'Eosinophil', 'Hb', 'Hct', 'Immature cell', 'Lymphocyte', 'MCH', 'MCHC', 'MCV', 'MPV', 'Metamyelocyte', 'Monocyte', 'Myelocyte', 'Normoblast', 'Other', 'PCT', 'PDW', 'PLT', 'Plasma cell', 'Promyelocyte', 'RBC', 'RDW(CV)', 'RDW(SD)', 'Segmented neutrophil', 'WBC', '절대단구수', '절대림프구수']].isna().sum().sort_values()
# lab_df['PLT'].max()
# lab_df['PLT'].min()
adm_df = pd.read_excel(f"{output_dir}/merged_adm_dates.xlsx")

# adm_df.columns
no_lab_uids = list()

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

# lab_df.columns
# lab_df['WBC']


for endpoint_lab in ['PLT', 'ANC', 'Hb','WBC']:
    # endpoint_lab = 'PLT'
    # max_time_at_risk = 365 * 99999999999
    # max_time_at_risk = 365
    # adm_df.columns
    max_time_at_risk = 90
    for inx, row in adm_df.iterrows(): #break
        uid = row['UID']
        sparse_uid_lab_df = lab_df[lab_df['UID']==uid].copy()

        if len(sparse_uid_lab_df)==0:
            print(f"({inx}) {uid} / No lab value")
            no_lab_uids.append(uid)
            continue

        # uid_lab_df
        min_lab_date = sparse_uid_lab_df['DATE'].min()
        max_lab_date = sparse_uid_lab_df['DATE'].max()

        uid_lab_df = pd.DataFrame(columns=['UID', 'DATE'])
        uid_lab_df['DATE'] = pd.date_range(start=min_lab_date, end=max_lab_date).astype(str)
        uid_lab_df['UID'] = uid

        uid_lab_df = uid_lab_df.merge(sparse_uid_lab_df, on=['UID', 'DATE'], how='left').fillna(method='ffill')

        # raise ValueError

        sbase_min_lab_date = (datetime.strptime(row['MIN_SINGLE_DATE'], '%Y-%m-%d') - timedelta(days=30)).strftime('%Y-%m-%d')
        uid_sbase_lab_df = uid_lab_df[(uid_lab_df['DATE'] <= row['MIN_SINGLE_DATE'])&(uid_lab_df['DATE'] >= sbase_min_lab_date)].sort_values(['DATE'])
        # uid_sbase_lab_df = uid_lab_df[uid_lab_df['DATE'] < row['MIN_SINGLE_DATE']].sort_values(['DATE'])
        uid_sadm_lab_df = uid_lab_df[(uid_lab_df['DATE'] >= row['MIN_SINGLE_DATE'])&(uid_lab_df['DATE'] <= row['MAX_SINGLE_DATE'])].sort_values(['DATE'])
        if len(uid_sbase_lab_df)==0:
            print(f"({inx}) {uid} / {endpoint_lab} / No sbase lab value")
            no_sbase_lab_uids.append(uid)
            continue
        else:
            yes_sbase_lab_uids.append(uid)

        if len(uid_sadm_lab_df)==0:
            print(f"({inx}) {uid} / {endpoint_lab} / No sadm lab value")
            no_sadm_lab_uids.append(uid)
            continue
        else:
            yes_sadm_lab_uids.append(uid)

        # raise ValueError
        if endpoint_lab == 'PLT':
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] >= 150)]
            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                if (uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']][endpoint_lab] < 150).sum() > 0:
                    continue
            except:
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] < 150)].copy()
        elif endpoint_lab == 'ANC':
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] >= 1500)]
            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                if (uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']][endpoint_lab] < 1500).sum() > 0:
                    continue
            except:
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] < 1500)].copy()
        elif endpoint_lab == 'Hb':
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] >= 8)]
            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                if (uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']][endpoint_lab] < 8).sum() > 0:
                    continue
            except:
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] < 8)].copy()
        elif endpoint_lab == 'WBC':
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] >= 4)]
            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                if (uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']][endpoint_lab] < 4).sum() > 0:
                    continue
            except:
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] < 4)].copy()
        elif endpoint_lab == 'Lactate':
            # lactate_cond =
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] < 4)]

            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                if (uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']][endpoint_lab] >= 4).sum() > 0:
                    continue
            except:
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] >= 4)].copy()


        single_res_dict = {'UID':uid,
                           'ENDPOINT':endpoint_lab,
                           'DRUG':row['SINGLE_DRUG'],
                           'FIRST_ADM_DATE':uid_sadm_lab_df['DATE'].iloc[0],
                           'LAST_ADM_DATE':uid_sadm_lab_df['DATE'].iloc[-1],
                           'EV': 0 if len(tar_rows)==0 else 1,
                           'BL_DATE':sbl_row['DATE'],
                           f'BL_{endpoint_lab}': sbl_row[endpoint_lab],
                           'EV_DATE': uid_sadm_lab_df['DATE'].iloc[-1] if len(tar_rows)==0 else tar_rows['DATE'].iloc[0],
                           f'EV_{endpoint_lab}': np.nan if len(tar_rows)==0 else tar_rows[endpoint_lab].iloc[0],
                           }

        single_res_dict['time'] = (datetime.strptime(single_res_dict['EV_DATE'],'%Y-%m-%d')-datetime.strptime(uid_sadm_lab_df['DATE'].iloc[0],'%Y-%m-%d')).total_seconds()/86400
        if single_res_dict['time'] > max_time_at_risk:
            single_res_dict['time'] = max_time_at_risk
            single_res_dict['event'] = 0
        else:
            single_res_dict['event'] = single_res_dict['EV']
        single_res_dict['group'] = 0

        surv_res_df.append(single_res_dict)


        # ibase_min_lab_date = (datetime.strptime(row['MIN_INTSEC_DATE'],'%Y-%m-%d')-timedelta(days=30)).strftime('%Y-%m-%d')
        # uid_ibase_lab_df = uid_lab_df[(uid_lab_df['DATE'] <= row['MIN_INTSEC_DATE'])&(uid_lab_df['DATE'] > ibase_min_lab_date)].sort_values(['DATE'])
        # # uid_ibase_lab_df = uid_lab_df[(uid_lab_df['DATE'] < row['MIN_INTSEC_DATE'])&(uid_lab_df['DATE'] > iadm_min_lab_date)].sort_values(['DATE'])
        # uid_iadm_lab_df = uid_lab_df[(uid_lab_df['DATE'] >= row['MIN_INTSEC_DATE'])&(uid_lab_df['DATE'] <= row['MAX_INTSEC_DATE'])].sort_values(['DATE'])
        # if len(uid_ibase_lab_df)==0:
        #     print(f"({inx}) {uid} / No ibase lab value")
        #     no_ibase_lab_uids.append(uid)
        #     continue
        # else:
        #     yes_ibase_lab_uids.append(uid)
        #
        # if len(uid_iadm_lab_df)==0:
        #     print(f"({inx}) {uid} / No iadm lab value")
        #     no_iadm_lab_uids.append(uid)
        #     continue
        # else:
        #     yes_iadm_lab_uids.append(uid)
        #
        # if endpoint_lab=='RBC_count':
        #     bl_rows = uid_ibase_lab_df[(~uid_ibase_lab_df[endpoint_lab].isna())&(uid_ibase_lab_df[endpoint_lab] <= 5)]
        #     try:
        #         bl_row = bl_rows.iloc[-1]
        #         # Baseline에서 Cr이 이미 높아지는 경우 제외
        #         if (uid_ibase_lab_df[uid_ibase_lab_df['DATE'] > bl_row['DATE']][endpoint_lab] > 5).sum() > 0:
        #             continue
        #     except:
        #         continue
        #     tar_rows = uid_iadm_lab_df[(uid_iadm_lab_df[endpoint_lab] >= 5)].copy()
        #
        # intersect_res_dict = {'ADM_TYPE':'INTSEC',
        #                       'DRUG':'BOTH',
        #                       'UID':uid,
        #                       'BL_DATE':bl_row['DATE'],
        #                       f'BL_{endpoint_lab}':bl_row[endpoint_lab],
        #                       'FIRST_ADM_DATE':uid_iadm_lab_df['DATE'].iloc[0],
        #                       'TRT':1,
        #                       'EV': 0 if len(tar_rows)==0 else 1,
        #                       'EV_DATE': uid_iadm_lab_df['DATE'].iloc[-1] if len(tar_rows)==0 else tar_rows['DATE'].iloc[0],
        #                       f'EV_{endpoint_lab}': np.nan if len(tar_rows)==0 else tar_rows[endpoint_lab].iloc[0],
        #                       }
        #
        # intersect_res_dict['time'] = (datetime.strptime(intersect_res_dict['EV_DATE'], '%Y-%m-%d') - datetime.strptime(uid_iadm_lab_df['DATE'].iloc[0], '%Y-%m-%d')).total_seconds() / 86400
        # if intersect_res_dict['time'] > max_time_at_risk:
        #     intersect_res_dict['time'] = max_time_at_risk
        #     intersect_res_dict['event'] = 0
        # else:
        #     intersect_res_dict['event'] = intersect_res_dict['EV']
        # intersect_res_dict['group'] = intersect_res_dict['TRT']
        #
        # surv_res_df.append(intersect_res_dict)

surv_res_df = pd.DataFrame(surv_res_df).sort_values(['time','event'])
surv_res_df.to_csv(f"{output_dir}/lnz_surv_res_df.csv", encoding='utf-8-sig', index=False)
# surv_res_df


## Fisher's Exact Test

single_survres_df = surv_res_df.copy()
# single_survres_df = surv_res_df[surv_res_df['ADM_TYPE']=='SINGLE'].copy()
# intsec_survres_df = surv_res_df[surv_res_df['ADM_TYPE']=='INTSEC'].copy()

# single_survres_df['UID']


# single_ev_count = single_survres_df['EV'].sum()
# intsec_ev_count = intsec_survres_df['EV'].sum()

from scipy.stats import fisher_exact

# 예시 2×2 분할표 데이터
#            사건O   사건X
# 그룹A       8      2
# 그룹B       1      5
# table = np.array([[intsec_ev_count, len(intsec_survres_df)-intsec_ev_count],
#                   [single_ev_count, len(single_survres_df)-single_ev_count]])
# max_incidence = max(intsec_ev_count/len(intsec_survres_df), single_ev_count/len(single_survres_df))
# max_odds = max(intsec_ev_count/(len(intsec_survres_df)-intsec_ev_count), single_ev_count/(len(single_survres_df)-single_ev_count))

# fisher_exact(table, alternative='two-sided') 가능
# oddsratio, p_value = fisher_exact(table, alternative='two-sided')

# print(table)
# print(f"Odds Ratio: {oddsratio:.3f}")
# print(f"P-value: {p_value:.5f}")

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

# 로그랭크 검정 (A vs B)
# A = surv_res_df[surv_res_df.group==1]; B = surv_res_df[surv_res_df.group==0]
# results = logrank_test(A["time"], B["time"], event_observed_A=A["event"], event_observed_B=B["event"])
# print(results.summary)



kmf = KaplanMeierFitter()

fig, ax = plt.subplots(figsize=(10, 8))

for g, gdf in surv_res_df.groupby("ENDPOINT"):
    kmf.fit(durations=gdf["time"], event_observed=gdf["event"], label=f"{g}")

    # 생존 확률 S(t)
    sf = kmf.survival_function_[kmf.survival_function_.columns[0]]
    # 누적발생률 CI(t) = 1 - S(t)
    ci_curve = 1 - sf

    # 신뢰구간도 변환
    ci_bounds = 1 - kmf.confidence_interval_
    lower = ci_bounds.iloc[:, 1]  # 하한
    upper = ci_bounds.iloc[:, 0]  # 상한

    # 곡선
    ax.step(ci_curve.index, ci_curve.values, where="post", label=f"{g}")
    # 신뢰구간 음영
    ax.fill_between(ci_curve.index, lower, upper, step="post", alpha=0.2)

# 서식
ax.set_title(f"Cumulative Incidence (myelosuppression) by Group (with 95% CI)")

# ax.set_title(f"Cumulative Incidence ({endpoint_lab}) by Group (with 95% CI)\nLog-rank test (p={round(results.summary['p'].iloc[0],3)})\nContingency Table: [{table[0][0]},{table[0][1]}]/[{table[1][0]},{table[1][1]}]")
ax.set_xlabel("Time")
ax.set_ylabel(f"Cumulative Incidence (myelosuppression)")
ax.set_ylim(0, 1.2)
ax.grid(True, linestyle="--")
ax.legend(title="Group")

plt.show()


# kmf = KaplanMeierFitter()
# plt.close('all')
# fig, ax = plt.subplots(figsize=(10, 8))
#
# for g, gdf in surv_res_df.groupby("group"):
#     kmf.fit(durations=gdf["time"], event_observed=gdf["event"], label=f"{g}")
#
#     # S(t)
#     sf = kmf.survival_function_[kmf.survival_function_.columns[0]]
#     # CI(t) = 1 - S(t)
#     ci = 1 - sf
#
#     # 계단그래프(우측계단)로 플로팅
#     ax.step(ci.index, ci.values, where="post", label=f"{g}")
#
# # 레이아웃/서식
# ax.set_title("Cumulative Incidence by Group (1 - KM Survival)")
# ax.set_xlabel("Time")
# ax.set_ylabel("Cumulative Incidence")
# ax.set_ylim(0, max_odds*3)
# ax.grid(True, linestyle="--")
# ax.legend(title="Group")
# plt.show()



# kmf = KaplanMeierFitter()
#
# plt.figure(figsize=(10,8))
# for g, gdf in surv_res_df.groupby("group"):
#     # kmf = KaplanMeierFitter()
#     # kmf.fit(gdf["time"], event_observed=gdf["event"], label=str(g))
#     # sf = kmf.survival_function_  # S(t)
#     # ax.step(sf.index, 1 - sf.values.ravel(), where="post", label=f"{g}")  # 1 - S(t)
#
#     kmf.fit(durations=gdf["time"], event_observed=gdf["event"], label=f"Group {g}")
#     kmf.plot()
#     # cum_inc = 1 - kmf.survival_function_
#     # ax = cum_inc.plot()
#
# # ax.set_title("Kaplan–Meier Survival Curves (Cumulative Incidence) by Group")
# # ax.set_xlabel("Time")
# # ax.set_ylabel("Cumulative Incidence")
# # plt.grid(True, linestyle="--")
#
#
# plt.title("Kaplan–Meier Survival Curves (Cumulative Incidence) by Group")
# plt.xlabel("Time")
# plt.ylabel("Cumulative Incidence")
# plt.grid(True, linestyle="--")
# plt.show()


# Cox 비례위험모형 (정적 공변량) / 예시: 나이(연속), 성별(범주), 군(범주)

# cph = CoxPHFitter()
# cph.fit(df_model, duration_col="time", event_col="event")
# cph.print_summary()     # HR, CI, p-value 등
# cph.plot_partial_effects_on_outcome(covariates="age", values=[40, 60, 80])
# plt.show()
# cph.check_assumptions(df_model, p_value_threshold=0.05, show_plots=False)