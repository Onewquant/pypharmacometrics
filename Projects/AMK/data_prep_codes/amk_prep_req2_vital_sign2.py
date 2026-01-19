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
pid_decode_df = pid_decode_df.drop(['순번'],axis=1)
pid_decode_df = pid_decode_df.rename(columns={'환자번호':'UID','기록작성일시':'DATETIME','TEMP':'BT','PR':'HR'})

# pid_decode_df.to_csv(f"{output_dir}/VS_NEW(Parsed).csv", index=False, encoding='utf-8-sig')
pid_decode_df = pd.read_csv(f"{output_dir}/VS_NEW(Parsed).csv")

pid_decode_df['DATE'] = pid_decode_df['DATETIME'].map(lambda x:x.split(' ')[0])

def x_to_float(x):
    try: return float(x)
    except: return np.nan
# pid_decode_df.columns
vs_type_list = ['SBP','DBP','HR','BT','RR']
for inx, vs_type in enumerate(vs_type_list):
    pid_decode_df[vs_type + '_MAX'] = pid_decode_df[vs_type].map(x_to_float)
    pid_decode_df[vs_type + '_MIN'] = pid_decode_df[vs_type].map(x_to_float)

df = pid_decode_df.drop(vs_type_list+['DATETIME'],axis=1)
min_max_dict = {c+'_MAX':'max' for c in vs_type_list}
min_max_dict.update({c+'_MIN':'min' for c in vs_type_list})
df = df.groupby(['UID','DATE'],as_index=False).agg(min_max_dict)
# df['SBP_MAX'].iloc[0]
df = df.melt(id_vars=['UID','DATE'], value_vars=list(min_max_dict.keys()), var_name='VS_TYPE', value_name="VALUE")
# df[(df['VS_TYPE']=='SBP_MAX')&(df['VALUE']<=300)].min()
# df[(df['VS_TYPE']=='SBP_MAX')&(df['VALUE']<=5)].max()
# df[(df['VS_TYPE']=='DBP_MAX')&(df['VALUE']<=300)].min()
# df[(df['VS_TYPE']=='DBP_MAX')&(df['VALUE']<=5)].max()
#
# df[(df['VS_TYPE']=='DBP_MIN')&(df['VALUE']<=40)].max()
# df[(df['VS_TYPE']=='DBP_MIN')&(df['VALUE']<=40)].max()
df['VS_TYPE2'] = df['VS_TYPE'].map(lambda x:x.split('_')[0])
df_list = list()
for vs_type, vs_range in {'SBP':(40,300), 'DBP':(20,180), 'HR':(20,300), 'BT':(28,43)}.items():#break
    vs_type_df = df[(df['VS_TYPE2']==vs_type)&(df['VALUE']>=vs_range[0])&(df['VALUE']<=vs_range[1])].copy()
    df_list.append(vs_type_df.drop(['VS_TYPE2'],axis=1))

df = pd.concat(df_list)

df.to_csv(f"{output_dir}/final_vs_data.csv", index=False, encoding='utf-8-sig')