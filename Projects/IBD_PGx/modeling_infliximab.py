from tools import *
from pynca.tools import *

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

induction_df = pd.read_excel(f"{resource_dir}/IBD-PGx_induction_date.xlsx")

totlab_df = pd.read_csv(f"{output_dir}/conc_df(lab).csv")
cumlab_df = pd.read_csv(f"{output_dir}/conc_df(cum_lab).csv")
lab_df = pd.concat([cumlab_df, totlab_df[['오더일','보고일']]], axis=1)
lab_df['DATETIME'] = lab_df['DATETIME'].map(lambda x:x.replace(' ','T'))
lab_df = lab_df.sort_values(['ID','DATETIME'])
lab_df.to_csv(f"{output_dir}/conc_df.csv", encoding='utf-8-sig', index=False)

# lab_df[(lab_df['CONC']!=lab_df['CONC_cumlab'])|(lab_df['DRUG']!=lab_df['DRUG_cumlab'])][['DRUG','CONC','DRUG_cumlab','CONC_cumlab']]

dose_df = pd.read_csv(f"{output_dir}/dose_df.csv")
dose_df.columns
lab_df['DATETIME']
dose_df['DATETIME']

# Induction Phase 만 구분

induction_df
dose_df['DATETIME']