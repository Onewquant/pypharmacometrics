from tools import *
from pynca.tools import *
import msoffcrypto
import io

# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Dosing 정보 불러오기

pid_df = pd.read_csv(f"{resource_dir}/AMK_REQ_DATA/df_11007_재식별.csv")
pid_df = pid_df.rename(columns={'Deidentification_ID':'UID','환자번호':'REQ_ID'}).drop(['Column1', 'Column2'],axis=1)
# pid_df[pid_df['UID']==11367967]

# vs_df.columns
vs_df = pd.read_csv(f"{resource_dir}/AMK_REQ_DATA/dt_11007_간호기록_321688.csv")
vs_df = vs_df.rename(columns={'환자번호':'REQ_ID','기록종류':'REC_TYPE','간호기록 작성일자':'DATETIME','항목명':'VS','항목값':'VALUE'})

vs_df = vs_df.merge(pid_df, on=['REQ_ID'], how='left').drop(['REQ_ID'],axis=1)
vs_df = vs_df[['UID', 'DATETIME', 'VS', 'VALUE']].copy()
# vs_df['VS'].unique()
vs_df['VS'] = vs_df['VS'].map({'BMI':'BMI','신체계측 > 신장(cm)':'HT','신체계측 > 체중(kg)':'WT','신장(cm)':'HT', '체중(kg)':'WT'})

vs_df.to_csv(f"{output_dir}/final_req_bodysize_data.csv", index=False, encoding='utf-8-sig')