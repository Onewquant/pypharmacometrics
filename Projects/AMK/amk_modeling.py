from tools import *
from pynca.tools import *

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
# prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

conc_df = pd.read_csv(f"{output_dir}/final_conc_df(with sampling).csv")
conc_df['DV'] = conc_df['CONC'].copy()
conc_df['MDV'] = 0
conc_df['TIME'] = conc_df['SAMP_DT'].copy()
conc_df['CMT'] = 1
conc_df['AMT'] = '.'
conc_df['RATE'] = '.'
# conc_df = conc_df[conc_df['SAMP_DT'].isna()]
print(f"Deleted / NaN CONC: {len(conc_df[conc_df['SAMP_DT'].isna()])} rows / {len(conc_df[conc_df['SAMP_DT'].isna()]['ID'].unique())} patients")
conc_df = conc_df[~conc_df['SAMP_DT'].isna()].copy()
# conc_df[conc_df['TIME'].map(lambda x:x.split('T')[-1]=='NN:NN')]

dose_df = pd.read_csv(f"{output_dir}/final_dose_df.csv")
dose_df['DV'] = '.'
dose_df['MDV'] = 1
dose_df['TIME'] = dose_df['DATE']+'T'+dose_df['TIME']
dose_df['CMT'] = 1
dose_df['AMT'] = dose_df['DOSE']
dose_df['RATE'] = dose_df['AMT'] / 0.5
print(f"Deleted / Ambiguous dosing time : {len(dose_df[(dose_df['TIME'].map(lambda x:x.split('T')[-1]=='NN:NN'))])} rows / {len(dose_df[(dose_df['TIME'].map(lambda x:x.split('T')[-1]=='NN:NN'))]['ID'].unique())} patients")
dose_df = dose_df[~(dose_df['TIME'].map(lambda x:x.split('T')[-1]=='NN:NN'))]

com_cols = ['ID','NAME','TIME','DV','MDV','CMT','AMT','RATE']

# id_df['TIME'].unique()

modeling_df = pd.concat([conc_df[com_cols], dose_df[com_cols]]).sort_values(['ID','TIME'])
modeling_df = modeling_df[~modeling_df['TIME'].isna()].copy()
# modeling_df = modeling_df[modeling_df['ID']!=10078411].copy()
modeling_df = modeling_df.sort_values(['ID','TIME']).reset_index(drop=True)


val_count_df = modeling_df.groupby('ID')['MDV'].value_counts().unstack(fill_value=0).reset_index(drop=False)
val_count_df.columns.name = None
no_conc_pids = val_count_df[val_count_df[0]==0]['ID']  # CONC Data 가 존재하지 않는 경우
no_dose_pids = val_count_df[val_count_df[1]==0]['ID']  # DOSE Data 가 존재하지 않는 경우
print(f"CONC Data 부재: {len(no_conc_pids)} patients")
print(f"DOSE Data 부재: {len(no_dose_pids)} patients")
no_onside_pids = pd.concat([no_conc_pids, no_dose_pids])
modeling_df = modeling_df[~modeling_df['ID'].isin(no_onside_pids)].copy()

for id, id_df in modeling_df.groupby('ID'): #break
    try:
        first_dose_dt = id_df[id_df['MDV']==1]['TIME'].iloc[0]
    except:
        raise ValueError
    id_df['TIME'] = id_df['TIME'].map(lambda x:(datetime.strptime(x,'%Y-%m-%dT%H:%M')-datetime.strptime(first_dose_dt,'%Y-%m-%dT%H:%M')).total_seconds()/3600)

modeling_df.to_csv(f"{output_dir}/modeling_amk.csv",index=False, encoding='utf-8-sig')
print(f"Total for modeling / {len(modeling_df)} rows / {len(modeling_df['ID'].unique())} patients")

