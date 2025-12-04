## 시간에 따른 모든 사람의 PD 경향성 그려보기

from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr
import statsmodels.api as sm

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## 모델링 데이터셋 로딩

df = pd.read_csv(f'{output_dir}/modeling_df_datacheck(for pda)/infliximab_integrated_datacheck(for pda).csv')
# df['DATETIME'].min()
# df['DATETIME'].max()
# df[df['ID']==1]
id_uid_dict = {row['ID']:row['UID'] for inx, row in df[['ID','UID']].drop_duplicates(['ID']).iterrows()}
df['TIME'] = df['TIME(DAY)'].copy()
uniq_df = df.drop_duplicates(['UID'])[['UID','NAME','START_INDMAINT']].copy()
ind_uids = list(uniq_df[(uniq_df['START_INDMAINT']==0)]['UID'].reset_index(drop=True))
maint_uids = list(uniq_df[(uniq_df['START_INDMAINT']==1)]['UID'].reset_index(drop=True))

## Simulation 데이터셋 로딩

sim_df = pd.read_csv(f"{nonmem_dir}/run/sim90",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
# sim_df[sim_df['ID']==1]
sim_df['ID'] = sim_df['ID'].astype(int)
sim_df['UID'] = sim_df['ID'].map(id_uid_dict)
sim_df['DT_YEAR'] = sim_df['DT_YEAR'].map(lambda x:str(int(x)).zfill(4))
sim_df['DT_MONTH'] = sim_df['DT_MONTH'].map(lambda x:str(int(x)).zfill(2))
sim_df['DT_DAY'] = sim_df['DT_DAY'].map(lambda x:str(int(x)).zfill(2))
sim_df['DATETIME'] = sim_df['DT_YEAR']+'-'+sim_df['DT_MONTH']+'-'+sim_df['DT_DAY']

# 시뮬레이션 데이터에 날짜 매칭
# sim_df[sim_df['ID']==0]
# sim_df[sim_df['ID']==1].head(50)
for uid, uid_sim_df in sim_df.groupby('ID'):
    # if uid==2:
    #     raise ValueError
    uid_base_inx = uid_sim_df[uid_sim_df['TIME']==0].index[0]
    uid_base_dt_rows = uid_sim_df.loc[:uid_base_inx,:].copy()
    uid_base_reltime = uid_base_dt_rows['TIME'].iloc[0]
    if len(uid_base_dt_rows)==0:
        raise ValueError
    cur_base_dt = uid_base_dt_rows['DATETIME'].iloc[0]
    for uid_inx, uid_row in uid_sim_df.iterrows(): #break

        # uid_row['DATETIME']
        time_delt_int_days = int(uid_row['TIME'] - uid_base_reltime)
        if uid_row['DATETIME']==cur_base_dt:
            new_dt = (datetime.strptime(uid_row['DATETIME'], '%Y-%m-%d') + timedelta(days=time_delt_int_days)).strftime('%Y-%m-%d')
        else:
            new_dt = uid_row['DATETIME']
            cur_base_dt = uid_row['DATETIME']
            uid_base_reltime = uid_row['TIME']
        sim_df.at[uid_inx,'DATETIME'] = new_dt

        # sim_df[sim_df['ID']==uid].iloc[0:50]
        # sim_df[sim_df['ID']==uid].iloc[50:100]
        # if time_delt < 1:
        #     sim_df
        # else:
        # uid_base_dt_rows = uid_sim_df[uid_sim_df['TIME']==0].copy()
        # uid_base_date = uid_base_dt_rows['DATETIME'].iloc[0]
        # uid_base_reltime = uid_base_dt_rows['TIME'].iloc[0]


        # raise ValueError
#
# sim_df['DT_MONTH']
# sim_df['DT_DAY']
# pd_df = pd.read_csv(f'{output_dir}/pd_bsize_df_ori.csv')

## PD 결과 정리한 데이터셋 로딩
# pd_res_df.columns
pd_res_df = pd.read_csv(f"{output_dir}/PKPD_EDA/PKPD_Corr0/pd_eda_df.csv")
# sim_df

## PD 결과 자체만 정리
# pd_res_df = pd_res_df.rename(columns={'BL_PD_CR':'BL_PD_PRO2(CR)','TG_PD_CR':'TG_PD_PRO2(CR)'})
pd_res_df = pd_res_df.rename(columns={'BL_PD_CR':'BL_PD_CR(PRO2)','TG_PD_CR':'TG_PD_CR(PRO2)'})
pd_res_df['BL_PD_CRP(BR)'] = (pd_res_df['BL_PD_CRP']<0.5)*1
pd_res_df['TG_PD_CRP(BR)'] = (pd_res_df['TG_PD_CRP']<0.5)*1

pd_res_df['BL_PD_FCAL250(BR)'] = (pd_res_df['BL_PD_FCAL']<250)*1
pd_res_df['TG_PD_FCAL250(BR)'] = (pd_res_df['TG_PD_FCAL']<250)*1

pd_res_df['BL_PD_FCAL150(BR)'] = (pd_res_df['BL_PD_FCAL']<150)*1
pd_res_df['TG_PD_FCAL150(BR)'] = (pd_res_df['TG_PD_FCAL']<150)*1
# pd_res_df.columns

# pd_sumres_df = pd.DataFrame(columns=['PD_EP','IND','MAINT'])
# pd_sumres_df = list()
# # for pd_ep in ['PRO2(CR)','CRP(BR)','FCAL250(BR)','FCAL150(BR)']:
# for pd_ep in ['CR(PRO2)','CRP(BR)','FCAL250(BR)','FCAL150(BR)']:
#     for trt_phase in ['IND','MAINT']:
#         if trt_phase=='IND':
#             ind_res_df = pd_res_df[pd_res_df['PHASE']==trt_phase].copy()
#         else:
#             ind_res_df = pd_res_df[pd_res_df['PHASE']!='IND'].copy()
#
#         ind_no_bl_val_df = ind_res_df[~ind_res_df[f'BL_PD_{pd_ep}'].isna()].copy()
#         ind_no_tg_val_df = ind_res_df[~ind_res_df[f'TG_PD_{pd_ep}'].isna()].copy()
#         ind_val_df = ind_res_df[(~ind_res_df[f'BL_PD_{pd_ep}'].isna())&(~ind_res_df[f'TG_PD_{pd_ep}'].isna())].copy()
#         rem_count = ind_val_df[f'TG_PD_{pd_ep}'].sum()
#
#         # ind_total_n = len(ind_no_tg_val_df)
#         ind_total_n = len(ind_val_df)
#
#         pd_sumres_df.append({'PD_EP':pd_ep, 'PHASE':trt_phase, 'STATE':'(total)', 'VAL': f"{ind_total_n}(100%)"})
#         pd_sumres_df.append({'PD_EP':pd_ep, 'PHASE':trt_phase, 'STATE':'(rem)',  'VAL': f"{rem_count}({round(100*rem_count/ind_total_n,1)}%)"})
#         pd_sumres_df.append({'PD_EP':pd_ep, 'PHASE':trt_phase, 'STATE':'(not)', 'VAL': f"{(ind_total_n-rem_count)}({round(100*(ind_total_n-rem_count)/ind_total_n,1)}%)"})
#
# pd_sumres_df = pd.DataFrame(pd_sumres_df)
# pd_sumres_df = pd_sumres_df.pivot(index=['PD_EP','STATE'],columns='PHASE',values='VAL')
# pd_sumres_df.columns.name = None
# # pd_sumres_df.index.name
# pd_sumres_df = pd_sumres_df.reset_index(drop=False).sort_values(by=['PD_EP','STATE'],ascending=[True,False])
# pd_sumres_df.to_csv(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_SIG/pd_only_analysis_df.csv", index=False, encoding='utf-8-sig')


# Exposure 비교를 위한 PK params 대표값 입력
# pd_res_df.columns
bl_demo_list = ['ALB', 'AGE', 'SEX', 'WT', 'HT', 'BMI', 'ADA']
pk_param_list = ['CL','V2','V3']
secondary_param_list = ['CUM_AUC','AUC24','CUM_DOSE','DOSE24','CMAX','ADJ_DOSE24','CUM_DOSEpWT','DOSE24pWT','CUM_ADJ_DOSEpWT','ADJ_DOSE24pWT','OPT_PCT']
for pk_param in bl_demo_list+secondary_param_list+pk_param_list:
    pd_res_df[pk_param]=np.nan
    # uid_sim_conc_df.columns
    # pd_res_df.columns
for inx, row in pd_res_df.iterrows():
    # break
    uid = row['UID']
    pname = row['NAME']
    phase = row['PHASE']
    # ibd_type = row['IBD_TYPE']

    # if uid==10875838:
    #     raise ValueError

    uid_sim_df = sim_df[sim_df['UID']==uid].copy()
    uid_sim_dose_df = uid_sim_df[uid_sim_df['MDV']==1].copy()
    uid_sim_conc_df = uid_sim_df[uid_sim_df['MDV']==0].copy()

    # sim_df.columns
    start_date = row['DS_DATE']
    end_date = row['DE_DATE']
    if (type(start_date)==float) or (type(end_date)==float):
        for pk_param in secondary_param_list:
            pd_res_df.at[inx, pk_param] = np.nan
        continue
    else:
        delta_days = (datetime.strptime(end_date, '%Y-%m-%d') - datetime.strptime(start_date, '%Y-%m-%d')).days

    # 투약 기간이 이틀밖에 안 되는 경우는 제외
    if (delta_days<=1):
        for pk_param in secondary_param_list:
            pd_res_df.at[inx, pk_param] = np.nan
        continue

    print(f"({inx}) {uid} / {pname} / {phase} / {str(start_date)} / {str(end_date)}")
    # uid_sim_conc_df.columns
    partial_conc_df = uid_sim_conc_df[(uid_sim_conc_df['DATETIME'] >= start_date)&(uid_sim_conc_df['DATETIME'] <= end_date)].copy()
    partial_dose_df = uid_sim_dose_df[(uid_sim_dose_df['DATETIME'] >= start_date)&(uid_sim_dose_df['DATETIME'] <= end_date)].copy()

    partial_conc_df = add_iAUC(partial_conc_df, time_col="TIME", conc_col="DV", id_col="ID")  # ID가 없으면 자동으로 전체 처리

    # total_days =
    # partial_dose_df.columns
    pd_res_df.at[inx, 'CUM_DOSE'] = partial_dose_df['AMT'].sum()
    pd_res_df.at[inx, 'DOSE24'] = (partial_dose_df['AMT'].sum())/(delta_days)

    pd_res_df.at[inx, 'CUM_DOSEpWT'] = pd_res_df.at[inx, 'CUM_DOSE']/(partial_dose_df['WT'].mean())
    pd_res_df.at[inx, 'DOSE24pWT'] = pd_res_df.at[inx, 'DOSE24']/(partial_dose_df['WT'].mean())
    # partial_dose_df['ROUTE']
    bioavailability = 0.667
    pd_res_df.at[inx, 'CUM_ADJ_DOSE'] = partial_dose_df.apply(lambda x:x['AMT']*bioavailability if x['ROUTE']==2 else x['AMT'], axis=1).sum()
    pd_res_df.at[inx, 'ADJ_DOSE24'] = (pd_res_df.at[inx, 'CUM_ADJ_DOSE'])/(delta_days)

    pd_res_df.at[inx, 'CUM_ADJ_DOSEpWT'] = pd_res_df.at[inx, 'CUM_ADJ_DOSE']/(partial_dose_df['WT'].mean())
    pd_res_df.at[inx, 'ADJ_DOSE24pWT'] = pd_res_df.at[inx, 'ADJ_DOSE24']/(partial_dose_df['WT'].mean())

    # pd_res_df.at[inx, 'DOSE1ADM'] = (partial_dose_df['AMT'].sum())/(len(partial_dose_df))

    pd_res_df.at[inx, 'CUM_AUC'] = partial_conc_df.iloc[-1]['cumAUC']
    pd_res_df.at[inx, 'AUC24'] = (partial_conc_df.iloc[-1]['cumAUC'])/(delta_days)


    int_day_conc_df = partial_conc_df[partial_conc_df['TIME']==(partial_conc_df['TIME'].map(int))].copy()
    pd_res_df.at[inx, 'OPT_PCT'] = ((int_day_conc_df['DV']>=5).sum())/len(int_day_conc_df)

    # partial_conc_df['DV']

    # if pd_res_df.at[inx, 'AUC24'] > 120:
    #     raise ValueError
    pd_res_df.at[inx, 'CMAX'] = partial_conc_df['DV'].max()
    pd_res_df.at[inx, 'CMAX'] = partial_conc_df['DV'].max()

    # partial_conc_df['ADA']
    # demo_c = 'ALB'
    for demo_c in bl_demo_list:
        if demo_c!='ADA':
            pd_res_df.at[inx, demo_c] = partial_conc_df[demo_c].iloc[0]
        else:
            pd_res_df.at[inx, demo_c] = partial_conc_df[demo_c].max()

    # pk_p = 'CL'
    for pk_p in pk_param_list:
        pd_res_df.at[inx, pk_p] = partial_conc_df[pk_p].median()

# pd_res_df.to_csv(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_SIG/pkpd_analysis_df.csv", index=False, encoding='utf-8-sig')

# if not os.path.exists(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups2"):
#     os.mkdir(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups2")


# phase_list = pd_res_df['PHASE'].sort_values().drop_duplicates()
# visualize_graph = False

# phase_list = ['IND','MAINT']

# pk_param_list
pd_ep_list = list({c[3:] for c in pd_res_df.columns if '_PD_' in c})

for pd_ep in pd_ep_list:         # ['PD_CR', 'PD_ABDPAIN', 'PD_PRO2', 'PD_FCAL', 'PD_RECTBLD', 'PD_STLFREQ', 'PD_TOTALSCORE', 'PD_CRP']
    # for pk_pr in pk_param_list:  # ['CUM_AUC', 'AUC24', 'CUM_DOSE', 'DOSE24', 'CMAX']
    pd_res_df[f'DELT_{pd_ep}'] = pd_res_df[f'TG_{pd_ep}'] - pd_res_df[f'BL_{pd_ep}']
    # pd_res_df

# pd_res_df.columns

id_list = ['UID','NAME', 'IBD_TYPE', 'PHASE', 'DAY_DELT']
bl_demo_res_list = bl_demo_list
pk_param_res_list = pk_param_list + secondary_param_list
pd_ep_res_list = [f"DELT_{c}" for c in pd_ep_list] + [f"TG_{c}" for c in pd_ep_list]

total_res_cols = id_list + bl_demo_res_list + pk_param_res_list + pd_ep_res_list
gene_res_df = pd_res_df[total_res_cols].copy()

# pd_res_df.to_csv(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_SIG/pkpd_analysis_df.csv", index=False, encoding='utf-8-sig')

if not os.path.exists(f"{output_dir}/PKPD_EDA/GENOMICS"):
    os.mkdir(f"{output_dir}/PKPD_EDA/GENOMICS")

gene_res_df.to_csv(f"{output_dir}/PKPD_EDA/GENOMICS/for_genomics_df.csv", index=False, encoding='utf-8-sig')

# gene_res_df['ADA'].sum()
# gene_res_df['UID'].drop_duplicates()