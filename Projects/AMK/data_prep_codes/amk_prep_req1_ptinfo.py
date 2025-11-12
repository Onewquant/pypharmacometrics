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

pid_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/CCM.csv")
pid_df = pid_df.rename(columns={'Deidentification_ID':'UID','환자번호':'REQ_ID'}).drop(['Column1', 'Column2'],axis=1)

ptinfo_result_df = pd.read_csv(f"{resource_dir}/AMK_REQ_DATA/dt_11007_환자정보_3538.csv")
ptinfo_result_df = ptinfo_result_df.rename(columns={'성별':'SEX','환자번호':'REQ_ID','생년월일':'BIRTH_DATE'})

pid_df = pid_df.merge(ptinfo_result_df, on=['REQ_ID'], how='left').drop(['REQ_ID'],axis=1)
pid_df.to_csv(f"{output_dir}/final_req_ptinfo_data.csv", index=False, encoding='utf-8-sig')