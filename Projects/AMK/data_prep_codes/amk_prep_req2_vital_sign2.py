from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr
import msoffcrypto
from io import BytesIO


# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

password = "snubhsnubh"

file_path = f"{resource_dir}/[AMK_AKI_ML_DATA]/VS_NEW.xlsx"

# 암호 해제
decrypted = BytesIO()
with open(file_path, "rb") as f:
    office_file = msoffcrypto.OfficeFile(f)
    office_file.load_key(password=password)
    office_file.decrypt(decrypted)

# pid_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/CCM.csv")
# pid_df = pid_df.rename(columns={'Deidentification_ID':'UID','환자번호':'REQ_ID'}).drop(['Column1', 'Column2'],axis=1)

pid_decode_df1 = pd.read_excel(decrypted, header=1, sheet_name="1-1000")
pid_decode_df2 = pd.read_excel(decrypted, header=1, sheet_name="1001-2000")
pid_decode_df3 = pd.read_excel(decrypted, header=1, sheet_name="2001-3000")
pid_decode_df4 = pd.read_excel(decrypted, header=1, sheet_name="3001-3591")

pid_decode_df = pd.concat([pid_decode_df1,pid_decode_df2,pid_decode_df3,pid_decode_df4], ignore_index=True)


pid_decode_df = pid_decode_df.rename(columns={'환자번호':'UID','Deidentification_ID':'PID'})
pid_decode_df['UID'] = pid_decode_df['UID'].map(lambda x: x.split('-')[0])

df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/VITAL.csv", encoding='euc-kr')
df = df.rename(columns={'ID':'UID'})
df['DATE'] = df['DATETIME'].map(lambda x:x.split(' ')[0])

def x_to_float(x):
    try: return float(x)
    except: return np.nan

vs_type_list = ['SBP','DBP','HR','BT']
for inx, vs_type in enumerate(vs_type_list):
    df[vs_type + '_MAX'] = df[vs_type].map(x_to_float)
    df[vs_type + '_MIN'] = df[vs_type].map(x_to_float)

df = df.drop(vs_type_list+['DATETIME'],axis=1)
min_max_dict = {c+'_MAX':'max' for c in vs_type_list}
min_max_dict.update({c+'_MIN':'min' for c in vs_type_list})
df = df.groupby(['UID','DATE'],as_index=False).agg(min_max_dict)
# df['SBP_MAX'].iloc[0]
df = df.melt(id_vars=['UID','DATE'], value_vars=list(min_max_dict.keys()), var_name='VS_TYPE', value_name="VALUE")
df.to_csv(f"{output_dir}/final_vs_data.csv", index=False, encoding='utf-8-sig')