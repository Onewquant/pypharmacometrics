from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

df = pd.read_csv(f"{output_dir}/amk_modeling_covar/AMK_F2.csv",encoding='utf-8-sig')
df.to_csv(f"{output_dir}/amk_modeling_covar/AMK_F3.csv",encoding='utf-8', index=False)
set(df['ID'].unique())
df.columns
df['TBIL']