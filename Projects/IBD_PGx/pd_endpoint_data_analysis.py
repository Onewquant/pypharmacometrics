## 시간에 따른 모든 사람의 PD 경향성 그려보기

from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## 모델링 데이터셋 로딩

df = pd.read_csv(f'{output_dir}/modeling_df_datacheck(for pda)/infliximab_integrated_datacheck(for pda).csv')
# df['DATETIME'].min()
# df['DATETIME'].max()
id_uid_dict = {row['ID']:row['UID'] for inx, row in df[['ID','UID']].drop_duplicates(['ID']).iterrows()}
df['TIME'] = df['TIME(DAY)'].copy()
uniq_df = df.drop_duplicates(['UID'])[['UID','NAME','START_INDMAINT']].copy()
ind_uids = list(uniq_df[(uniq_df['START_INDMAINT']==0)]['UID'].reset_index(drop=True))
maint_uids = list(uniq_df[(uniq_df['START_INDMAINT']==1)]['UID'].reset_index(drop=True))

## Simulation 데이터셋 로딩

sim_df = pd.read_csv(f"{nonmem_dir}/run/sim86",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
sim_df['ID'] = sim_df['ID'].astype(int)
sim_df['UID'] = sim_df['ID'].map(id_uid_dict)
sim_df['DT_YEAR'] = sim_df['DT_YEAR'].map(lambda x:str(int(x)).zfill(4))
sim_df['DT_MONTH'] = sim_df['DT_MONTH'].map(lambda x:str(int(x)).zfill(2))
sim_df['DT_DAY'] = sim_df['DT_DAY'].map(lambda x:str(int(x)).zfill(2))
sim_df['DATETIME'] = sim_df['DT_YEAR']+'-'+sim_df['DT_MONTH']+'-'+sim_df['DT_DAY']
#
# sim_df['DT_MONTH']
# sim_df['DT_DAY']
# pd_df = pd.read_csv(f'{output_dir}/pd_bsize_df_ori.csv')

## PD 결과 정리한 데이터셋 로딩

pd_res_df = pd.read_csv(f"{output_dir}/PKPD_EDA/PKPD_Corr0/pd_eda_df.csv")
pd_ep_list = list({c[3:] for c in pd_res_df.columns if '_PD_' in c})
for pd_ep in pd_ep_list:
    raise ValueError
    pd_res_df[f'BL_{pd_ep}']
    pd_res_df[f'TG_{pd_ep}']
