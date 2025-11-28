from tools import *
from pynca.tools import *

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'
# prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

# cm_df = pd.read_csv(f"{output_dir}/final_comorbidity_df.csv")
# comed_df = pd.read_csv(f"{output_dir}/final_comed_df.csv")
# comed_df[comed_df['NEPHTOX_DRUG_YN']==1]['CCM_CATNUM'].unique()
# cm_df
#
# df['']
#
# uids = list(cm_df[(cm_df['CM_CATNUM']==18)]['UID'].unique())
# len(uids)
# df[df['UID'].isin(uids)]
#
# test_pt = df[df['UID']==10182459].iloc[0]
# test_pt = df.iloc[0]
# test_pt_cm.columns
#
# test_pt_cm = cm_df[(cm_df['UID'] == test_pt['UID'])&(cm_df['CM_CATNUM']==18)].copy()
# test_pt_cm[test_pt_cm['DATE']>'2008-05-15']

df = pd.read_csv(f"{output_dir}/final_mlres_data.csv")
set(df.columns)
df['HEPATIC INSUFFICIENCY'].sum()
df['SEPSIS'].sum()
df['SHOCK'].sum()
df['ACIDOSIS'].sum()
df['DIABETES'].sum()
df['ANEMIA'].sum()

df['AT LEAST ONE NEPHROTOXIC AGENT'].sum()
df['LOOP DIURETIC'].sum()