from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


# pid_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/CCM.csv")
# pid_df = pid_df.rename(columns={'Deidentification_ID':'UID','환자번호':'REQ_ID'}).drop(['Column1', 'Column2'],axis=1)


df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/CM.csv")
cm_index_df = pd.read_excel(f"{resource_dir}/[AMK_AKI_ML_DATA]/CM_index_table.xlsx")
df = df.rename(columns={'ID':'UID'})
df['UID'] = df['UID'].map(lambda x: x.split('-')[0])