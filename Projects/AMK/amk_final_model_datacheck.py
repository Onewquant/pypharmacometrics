from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

# df = pd.read_csv(f"{output_dir}/amk_modeling_covar/AMK_F2.csv",encoding='utf-8-sig')
# df.to_csv(f"{output_dir}/amk_modeling_covar/AMK_F3.csv",encoding='utf-8', index=False)
# set(df['ID'].unique())
# df.columns
# df['TBIL']


cdf = pd.read_csv(f"{output_dir}/amk_modeling_covar/amk_modeling_df_covar.csv",encoding='utf-8-sig')
cdf['DV'] = cdf['DV'].map(lambda x:x if x=='.' else float(x))
cdf['UID'] = cdf['UID'].astype(int)
cdf['TIME'] = cdf['TIME'].map(lambda x:round(x,4))
cdf['AMT'] = cdf['AMT'].map(lambda x:x if x=='.' else float(x))

fdf = pd.read_csv(f"{output_dir}/amk_modeling_covar/AMK_F3.csv",encoding='utf-8-sig')
fdf = fdf.drop(['CRCL'],axis=1)
fdf['DV'] = fdf['DV'].map(lambda x:x if x=='.' else float(x))
fdf['UID'] = fdf['UID'].astype(int)
fdf['TIME'] = fdf['TIME'].map(lambda x:round(x,4))
fdf['AMT'] = fdf['AMT'].map(lambda x:x if x=='.' else float(x))

cdf['CORRECTION'] = 0
fdf['CORRECTION'] = 1
# fdf['DV'].iloc[0]
# fdf['TIME'].iloc[0]

# mdf['ID'].iloc[0]
# mdf['ID'].iloc[1]

exception_uids = cdf[~(cdf['UID']).isin(fdf['UID'].drop_duplicates())].drop_duplicates(['UID'])['UID']
#
# mdf['TIME'].iloc[-3]
# mdf['TIME'].iloc[-4]

# mdf['ID'].iloc[0]
# mdf['ID'].iloc[1]

mdf = pd.concat([cdf, fdf]).sort_values(['UID','TIME','DV'])
# mdf = mdf.drop_duplicates(['ID','TIME','DV','AMT'],keep=False)
mdf = mdf[['ID','TIME','TAD','DV','MDV','CMT','AMT','RATE','UID','CORRECTION']].copy()
mdf = mdf.drop_duplicates(['UID','TIME','DV','AMT'],keep=False)
mdf = mdf.sort_values(['UID','TIME','DV','AMT','CORRECTION'])
# mdf

mdf[mdf['UID'].isin(exception_uids)].to_csv(f"{output_dir}/amk_modeling_covar/deleted_rows.csv",index=False)
mdf[~mdf['UID'].isin(exception_uids)].to_csv(f"{output_dir}/amk_modeling_covar/changed_rows.csv",index=False)