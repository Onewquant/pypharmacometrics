from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

# df = pd.read_csv(f"{nonmem_dir}/amk_simulation_df.csv")
# df.drop_duplicates(['ID'])
# df[df['ID']==449]
# 156, 159, 255, 263, 449, 816, 1081, 1851, 2093, 2377, 840

## DEMO Covariates Loading
# simulation_df.drop_duplicates(['ID'])
simulation_df = pd.read_csv(f"{nonmem_dir}/amk_modeling_df_covar.csv")
# simulation_df['TIME_DIFF'] = simulation_df['TIME'].diff()
# simulation_df[simulation_df['TIME_DIFF']>0]['TIME_DIFF'].sort_values()
interval = 10/60
# final_sim_df = get_model_population_sim_df(df=simulation_df, interval=interval, add_on_period=4*24)
final_sim_df = get_model_population_sim_df(df=simulation_df, interval=interval, add_on_period=0)
final_sim_df['UID'] = final_sim_df['UID'].astype(int)
# final_sim_df.to_csv(f"{output_dir}/amk_simulation_df_ori.csv",index=False, encoding='utf-8-sig')

# uniq_ids = list(final_sim_df['ID'].unique())[:30]
# final_sim_df[final_sim_df['ID'].isin(uniq_ids)].to_csv(f"{output_dir}/amk_simulation_partial_df.csv",index=False, encoding='utf-8-sig')

# final_sim_df['ID'].drop_duplicates()
# final_sim_df['TAD'] = add_time_after_dosing_column(df=final_sim_df)

## Filtered simulation dataset (aki 반영)

# final_sim_df = pd.read_csv(f"{output_dir}/amk_simulation_df.csv")
# final_sim_df = simulation_df.copy()
# final_sim_df['UID'] = final_sim_df['UID'].astype(int)

aki_df = pd.read_csv(f'{output_dir}/amk_aki.csv')
# aki_df = pd.read_csv(f'{output_dir}/amk_aki(filtered).csv')
# aki_df = aki_df.drop_duplicates(subset=['ID','COND_TYPE','BLCr_DT','AKI_RSLT'])
aki_df = aki_df.sort_values(['UID','AKI_DT']).drop_duplicates(subset=['UID'])
aki_ids = list(aki_df['UID'])

nonaki_pt_df = final_sim_df[~final_sim_df['UID'].isin(aki_ids)].copy()
nonaki_pt_df['AKI_DT'] = -1
aki_pt_df = final_sim_df[final_sim_df['UID'].isin(aki_ids)].copy()
aki_pt_df = aki_pt_df.merge(aki_df[['UID','AKI_DT']], on=['UID'],how='left')
aki_pt_df = aki_pt_df[(aki_pt_df['TIME'] <= (aki_pt_df['AKI_DT']-24))].copy()
# aki_pt_df = aki_pt_df[(aki_pt_df['TIME'] <= (aki_pt_df['AKI_DT']))].copy()
# aki_pt_df[aki_pt_df['AKI_DT'] <= 24].drop_duplicates(['UID'])
# aki_pt_df[(aki_pt_df['TIME'] <= (aki_pt_df['AKI_DT']-24))].copy()
# aki_pt_df.drop_duplicates(['UID'])
# nonaki_pt_df.drop_duplicates(['UID'])
# aki_pt_df.drop_duplicates(['UID'])
# filt_final_sim_df.drop_duplicates(['UID'])

filt_final_sim_df = pd.concat([aki_pt_df, nonaki_pt_df]).sort_values(['ID','TIME','MDV'])
filt_final_sim_df['YY'] = filt_final_sim_df['DATETIME_ORI'].map(lambda x:int(x[:4]))
filt_final_sim_df['MM'] = filt_final_sim_df['DATETIME_ORI'].map(lambda x:int(x[5:7]))
filt_final_sim_df['DD'] = filt_final_sim_df['DATETIME_ORI'].map(lambda x:int(x[8:10]))
filt_final_sim_df['HOUR'] = filt_final_sim_df['DATETIME_ORI'].map(lambda x:int(x[11:13]))
filt_final_sim_df['MINUTE'] = filt_final_sim_df['DATETIME_ORI'].map(lambda x:int(x[14:]))
filt_final_sim_df = filt_final_sim_df[['ID', 'TIME', 'TAD', 'DV', 'MDV', 'CMT', 'AMT', 'RATE', 'UID','SEX','AGE','ALB','WT','CREATININE','AKI_DT','YY','MM','DD','HOUR','MINUTE']].copy()
# filt_final_sim_df.to_csv(f"{output_dir}/amk_simulation_df.csv",index=False, encoding='utf-8-sig')
filt_final_sim_df.to_csv(f"{nonmem_dir}/amk_simulation_df.csv",index=False, encoding='utf-8-sig')

# filt_final_sim_df[filt_final_sim_df['AKI_DT']>0][['ID','UID']].drop_duplicates()
# aki_pt_df[(aki_pt_df['ID']==21)][['ID','UID','TIME','CREATININE','AKI_DT']]

# max_time_df = nonaki_pt_df.groupby('ID',as_index=False)['TIME'].agg('max')
# sns.displot(max_time_df['TIME']/24)
# aki_pt_df['AKI_DT'] = aki_pt_df['AKI_DT'].map()
# aki_pt_df[(aki_pt_df['ID']==1)][['TIME','CREATININE']]
# aki_pt_df[aki_pt_df['AKI_DT'] < 24].drop_duplicates(['ID'])

# filt_final_sim_df[filt_final_sim_df['ID']==3]

# import seaborn as sns
# sns.displot(aki_pt_df.drop_duplicates(['ID'])['AKI_DT']/24)

## AKI 발생시점보다 현재 존재하는 MAX_TIME 시점이 짧은 경우
# max_time_df = final_sim_df.groupby('ID',as_index=False)['TIME'].agg('max')
# max_time_df = max_time_df.rename(columns={'TIME':'MAX_TIME'})
# new_aki_pt_df = aki_pt_df.merge(aki_df[['ID','AKI_DT']], on=['ID'],how='left').merge(max_time_df[['ID','MAX_TIME']], on=['ID'],how='left')
# new_aki_pt_df[new_aki_pt_df['AKI_DT'] > new_aki_pt_df['MAX_TIME']].drop_duplicates(['ID'])


