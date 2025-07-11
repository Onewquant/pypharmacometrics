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


dose_df = pd.read_csv(f"{output_dir}/final_dose_df.csv")
dose_df['DV'] = '.'
dose_df['MDV'] = 1
dose_df['TIME'] = dose_df['DATE']+'T'+dose_df['TIME']
dose_df['CMT'] = 1
dose_df['AMT'] = dose_df['DOSE']
dose_df['RATE'] = dose_df['AMT'] / 0.5

com_cols = ['ID','NAME','TIME','DV','MDV','CMT','AMT','RATE']

# id_df['TIME'].unique()

modeling_df = pd.concat([conc_df[com_cols], dose_df[com_cols]]).sort_values(['ID','TIME'])
modeling_df = modeling_df[~modeling_df['TIME'].isna()].copy()
modeling_df = modeling_df[modeling_df['ID']!=10078411].copy()
modeling_df = modeling_df.sort_values(['ID','TIME']).reset_index(drop=True)
nn_time_df = modeling_df[modeling_df['TIME'].map(lambda x:x.split('T')[-1]=='NN:NN')]

for inx, row in nn_time_df.iterrows():break


for id, id_df in modeling_df.groupby('ID'): #break
    try:
        first_dose_dt = id_df[id_df['MDV']==1]['TIME'].iloc[0]
    except:
        pass
    id_df['TIME'] = id_df['TIME'].map(lambda x:(datetime.strptime(x,'%Y-%m-%dT%H:%M')-datetime.strptime(first_dose_dt,'%Y-%m-%dT%H:%M')).total_seconds()/3600)

modeling_df.to_csv(f"{output_dir}/modeling_amk.csv",index=False, encoding='utf-8-sig')

