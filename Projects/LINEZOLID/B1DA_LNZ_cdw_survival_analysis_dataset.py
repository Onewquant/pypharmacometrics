from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

prj_name = 'LNZ'
prj_dir = './Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'
# lab_df.columns
lactate_cols = ['Lactate (ICU용)', 'Lactate (마취과용)', 'Lactate (응급실용)', 'Lactate, Lactic acid (em)', ]
pH_cols = ['pH', 'pH (i-STAT)']
anc_cols = ['ANC', 'ANC (em)']
plt_cols = ['PLT', 'PLT (em)']
wbc_cols = ['WBC', 'WBC (em)', ]
hb_cols = ['Hb', 'Hb (em)', ]
total_cols = lactate_cols+pH_cols+anc_cols+plt_cols+wbc_cols+hb_cols


lab_df = pd.read_csv(f"{output_dir}/lnz_final_lab_df.csv")
# lab_df.columns
lab_df['Lactate'] = lab_df[lactate_cols].max(axis=1)
lab_df['pH'] = lab_df[pH_cols].min(axis=1)
lab_df['ANC'] = lab_df[anc_cols].min(axis=1)
lab_df['PLT'] = lab_df[plt_cols].min(axis=1)
lab_df['WBC'] = lab_df[wbc_cols].min(axis=1)
lab_df['Hb'] = lab_df[hb_cols].min(axis=1)
# lab_df['LACT']

# lab_df.columns
# lab_df[['UID', 'DATE', 'ANC', 'Atypical lymphocyte', 'Band neutrophil', 'Basophil', 'Blast', 'Eosinophil', 'Hb', 'Hct', 'Immature cell', 'Lymphocyte', 'MCH', 'MCHC', 'MCV', 'MPV', 'Metamyelocyte', 'Monocyte', 'Myelocyte', 'Normoblast', 'Other', 'PCT', 'PDW', 'PLT', 'Plasma cell', 'Promyelocyte', 'RBC', 'RDW(CV)', 'RDW(SD)', 'Segmented neutrophil', 'WBC', '절대단구수', '절대림프구수']].isna().sum().sort_values()
# lab_df['PLT'].max()
# lab_df['PLT'].min()
adm_df = pd.read_excel(f"{output_dir}/merged_adm_dates.xlsx")
# adm_df['UID'].drop_duplicates()

# adm_df.columns
no_lab_uids = list()

no_sbase_lab_uids = dict()
no_sadm_lab_uids = dict()
# yes_sbase_lab_uids = list()
no_ibase_lab_uids = list()


# yes_sadm_lab_uids = list()
# no_iadm_lab_uids = list()

not_normal_base_lab_uids = dict()
# len(yes_sbase_lab_uids)
surv_res_df = list()

# lab_df.columns
# lab_df['WBC']


for endpoint_lab in ['PLT', 'ANC', 'Hb','WBC','Lactate']:
# for endpoint_lab in ['ANC', 'Lactate']:
    # endpoint_lab = 'PLT'
    # max_time_at_risk = 365 * 99999999999
    # max_time_at_risk = 365
    # adm_df.columns
    # max_time_at_risk = 90
    max_time_at_risk = 365
    not_normal_base_lab_uids[endpoint_lab] = list()
    no_sbase_lab_uids[endpoint_lab] = list()
    no_sadm_lab_uids[endpoint_lab] = list()
    for inx, row in adm_df.iterrows(): #break
        uid = row['UID']
        # if uid==155505674746153:
        #     raise ValueError
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
        uid_sbase_lab_df = uid_lab_df[(uid_lab_df['DATE'] < row['MIN_SINGLE_DATE'])&(uid_lab_df['DATE'] >= sbase_min_lab_date)].sort_values(['DATE'])
        # uid_sbase_lab_df = uid_lab_df[uid_lab_df['DATE'] < row['MIN_SINGLE_DATE']].sort_values(['DATE'])
        uid_sadm_lab_df = uid_lab_df[(uid_lab_df['DATE'] >= row['MIN_SINGLE_DATE'])&(uid_lab_df['DATE'] <= row['MAX_SINGLE_DATE'])].sort_values(['DATE'])
        if len(uid_sbase_lab_df)==0:
            print(f"({inx}) {uid} / {endpoint_lab} / No sbase lab value")
            no_sbase_lab_uids[endpoint_lab].append(uid)
            continue
        else:
            # yes_sbase_lab_uids.append(uid)
            pass

        if len(uid_sadm_lab_df)==0:
            print(f"({inx}) {uid} / {endpoint_lab} / No sadm lab value")
            no_sadm_lab_uids[endpoint_lab].append(uid)
            continue
        else:
            # yes_sadm_lab_uids.append(uid)
            pass

        # raise ValueError
        if endpoint_lab == 'PLT':
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] >= 150)]
            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                if (uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']][endpoint_lab] < 150).sum() > 0:
                    not_normal_base_lab_uids[endpoint_lab].append(uid)
                    continue
            except:
                not_normal_base_lab_uids[endpoint_lab].append(uid)
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] < 150)].copy()
        elif endpoint_lab == 'ANC':
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] >= 1500)]
            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                if (uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']][endpoint_lab] < 1500).sum() > 0:
                    not_normal_base_lab_uids[endpoint_lab].append(uid)
                    continue
            except:
                not_normal_base_lab_uids[endpoint_lab].append(uid)
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] < 1500)].copy()
        elif endpoint_lab == 'Hb':
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] >= 8)]
            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                if (uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']][endpoint_lab] < 8).sum() > 0:
                    not_normal_base_lab_uids[endpoint_lab].append(uid)
                    continue
            except:
                not_normal_base_lab_uids[endpoint_lab].append(uid)
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] < 8)].copy()
        elif endpoint_lab == 'WBC':
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] >= 4)]
            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                if (uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']][endpoint_lab] < 4).sum() > 0:
                    not_normal_base_lab_uids[endpoint_lab].append(uid)
                    continue
            except:
                not_normal_base_lab_uids[endpoint_lab].append(uid)
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] < 4)].copy()
        elif endpoint_lab == 'Lactate':
            # lactate_cond =
            sbl_rows = uid_sbase_lab_df[(~uid_sbase_lab_df[endpoint_lab].isna()) & (uid_sbase_lab_df[endpoint_lab] < 4) & (~uid_sbase_lab_df['pH'].isna()) & (uid_sbase_lab_df['pH'] >= 7.35)]

            # sbl_row = sbl_rows.iloc[-1]
            try:
                sbl_row = sbl_rows.iloc[-1]
                # Baseline에서 Cr이 이미 높아지는 경우 제외
                uid_sbase_check_df = uid_sbase_lab_df[uid_sbase_lab_df['DATE'] > sbl_row['DATE']].copy()
                if ((uid_sbase_check_df[endpoint_lab] >= 4)&(uid_sbase_check_df['pH'] < 7.35)).sum() > 0:
                    not_normal_base_lab_uids[endpoint_lab].append(uid)
                    continue
            except:
                not_normal_base_lab_uids[endpoint_lab].append(uid)
                continue
            tar_rows = uid_sadm_lab_df[(uid_sadm_lab_df[endpoint_lab] >= 4)&(uid_sadm_lab_df['pH'] < 7.35)].copy()
        else:
            raise ValueError

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


surv_res_df = pd.DataFrame(surv_res_df).sort_values(['time','event'])
surv_res_df.to_csv(f"{output_dir}/b1da_lnz_surv_res_df.csv", encoding='utf-8-sig', index=False)

# surv_res_df['UID'].drop_duplicates()
###################################

print(f"# [Total Linezolid-Administrated]: {len(adm_df['UID'].drop_duplicates())}")

for inx, endpoint_lab in enumerate(['PLT', 'ANC', 'Hb','WBC','Lactate']):
    if inx == 0:
        no_sbase_lab_value_set = set(no_sbase_lab_uids[endpoint_lab])
        no_sadm_lab_value_set = set(no_sadm_lab_uids[endpoint_lab])
    else:
        no_sbase_lab_value_set = no_sbase_lab_value_set.intersection(set(no_sbase_lab_uids[endpoint_lab]))
        no_sadm_lab_value_set = no_sadm_lab_value_set.intersection(set(no_sadm_lab_uids[endpoint_lab]))
no_either_lab_value_set = no_sbase_lab_value_set.union(no_sadm_lab_value_set)

print(f"  (No lab value for baseline or adm period: -{len(no_either_lab_value_set)})")
# print(f"(No lab value for baseline or adm period: -{len((set(no_ibase_lab_uids)))})")
# not_normal_integ_set = set()
for inx, endpoint_lab in enumerate(['PLT', 'ANC', 'Hb','WBC','Lactate']):
    if inx==0:
        not_normal_integ_set = set(not_normal_base_lab_uids[endpoint_lab])
    else:
        not_normal_integ_set = not_normal_integ_set.intersection(set(not_normal_base_lab_uids[endpoint_lab]))
print(f"  (Not normal value in baseline period: -{len(not_normal_integ_set)})")

# for endpoint_lab in ['PLT', 'ANC', 'Hb','WBC','Lactate']:
#     print(f"Not normal baseline value ({endpoint_lab}): {len(not_normal_base_lab_uids[endpoint_lab])}")

print(f"# [Survival Analysis]: {len(surv_res_df['UID'].drop_duplicates())}")
for endpoint_lab in ['PLT', 'ANC', 'Hb','WBC','Lactate']:
    print(f"  SA_Subgroup ({endpoint_lab}): {len(surv_res_df[surv_res_df['ENDPOINT']==endpoint_lab]['UID'].drop_duplicates())}")


###################################


single_survres_df = surv_res_df.copy()

anc_df = single_survres_df[single_survres_df['ENDPOINT']=='ANC'].copy()
100 * len(anc_df[anc_df['EV']==1])/len(anc_df)

lact_df = single_survres_df[single_survres_df['ENDPOINT']=='Lactate'].copy()
100 * len(lact_df[lact_df['EV']==1])/len(lact_df)

# from lifelines import KaplanMeierFitter
# import matplotlib.pyplot as plt
#
# kmf = KaplanMeierFitter()
# fig, ax = plt.subplots(figsize=(10, 8))
#
# final_rates = []  # 최종 발생률을 저장할 리스트
#
# for g, gdf in surv_res_df.groupby("ENDPOINT"):
#     kmf.fit(durations=gdf["time"], event_observed=gdf["event"], label=f"{g}")
#
#     # 생존 확률 S(t)
#     sf = kmf.survival_function_[kmf.survival_function_.columns[0]]
#     # 누적발생률 CI(t) = 1 - S(t)
#     ci_curve = 1 - sf
#
#     # 신뢰구간도 변환
#     ci_bounds = 1 - kmf.confidence_interval_
#     lower = ci_bounds.iloc[:, 1]  # 하한
#     upper = ci_bounds.iloc[:, 0]  # 상한
#
#     # 곡선 그리기
#     ax.step(ci_curve.index, ci_curve.values, where="post", label=f"{g}")
#     ax.fill_between(ci_curve.index, lower, upper, step="post", alpha=0.2)
#
#     # 👉 최종 발생률(마지막 시점의 값)
#     final_time = ci_curve.index[-1]
#     final_value = ci_curve.iloc[-1]
#     final_rates.append({"ENDPOINT": g, "final_cumulative_incidence": float(final_value)})
#
#     # 👉 그래프 끝점에 값 표시
#     ax.text(final_time, final_value,
#             f"{final_value*100:.3f}%",  # % 단위로 표시
#             ha="left", va="center", fontsize=9)


from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import pandas as pd

kmf = KaplanMeierFitter()
fig, ax = plt.subplots(figsize=(15, 12))

final_rates = []  # 최종 발생률과 95% CI를 저장할 리스트

for g, gdf in surv_res_df.groupby("ENDPOINT"):
    kmf.fit(durations=gdf["time"], event_observed=gdf["event"], label=f"{g}")

    # 1. 생존 확률 S(t) 및 누적 발생률 CI(t)
    sf = kmf.survival_function_[kmf.survival_function_.columns[0]]
    ci_curve = 1 - sf

    # 2. 신뢰구간도 1 - S(t)로 변환
    ci_bounds = 1 - kmf.confidence_interval_
    lower_curve = ci_bounds.iloc[:, 1]  # 하한
    upper_curve = ci_bounds.iloc[:, 0]  # 상한

    # 3. 곡선과 신뢰구간 그림
    ax.step(ci_curve.index, ci_curve.values, where="post", label=f"{g}")
    ax.fill_between(ci_curve.index, lower_curve, upper_curve, step="post", alpha=0.2)

    # 4. 최종(마지막 시점) 누적 발생률과 95% CI 추출
    final_time = ci_curve.index[-1]
    final_val = float(ci_curve.iloc[-1])
    final_low = float(lower_curve.iloc[-1])
    final_high = float(upper_curve.iloc[-1])

    # 그래프 끝 지점에 값과 95% CI 표시
    ax.text(final_time, final_val,
            f"  {final_val*100:.1f}%\n  [{final_low*100:.1f}–{final_high*100:.1f}%]",
            ha="left", va="center", fontsize=12)

    # DataFrame용으로 저장
    final_rates.append({
        "ENDPOINT": g,
        "final_cumulative_incidence": final_val,
        "lower_95CI": final_low,
        "upper_95CI": final_high
    })


ax.set_title("Cumulative Incidence (ADRs) by Group (with 95% CI)",fontsize=14)
ax.set_xlabel("Time (Days)",fontsize=14)
ax.set_ylabel(f"Cumulative Incidence (ADRs)",fontsize=14)
ax.set_ylim(0, 1.1)
ax.set_xlim(0, 400)
ax.grid(True, linestyle="--")
ax.legend(title="Group", title_fontsize=12, fontsize=12)

plt.tight_layout()
plt.savefig(f"{output_dir}/B1DA_KM_plot(ADRs).png")  # PNG 파일로 저장

# 👉 최종 발생률을 DataFrame으로 정리
final_df = pd.DataFrame(final_rates)
print(final_df)