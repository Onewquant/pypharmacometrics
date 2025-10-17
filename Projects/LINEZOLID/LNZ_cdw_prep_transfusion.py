from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

prj_name = 'LNZ'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

raw_df = pd.read_csv(f"{resource_dir}/lnz_cdw_transfusion_data.csv")
# raw_df.columns
# raw_df = raw_df.rename(columns={'환자번호':'UID','수행시간':'ACTING','약품 오더일자':'DATE', '[실처방] 용법':'REGIMEN','[실처방] 처방일수':'DAYS', '[함량단위환산] 1일 처방량':'DAILY_DOSE','[함량단위환산] 1회 처방량':'DOSE', '[실처방] 투약위치':'PLACE',"[실처방] 경로":'ROUTE','약품명(성분명)':'DRUG','[실처방] 처방비고':'ETC_INFO'})
raw_df = raw_df.rename(columns={'환자번호':'UID','(처치/수혈) 오더일자':'DATE', '(처치/수혈) 오더명':'TF_TYPE','(처치/수혈) 처방 횟수/수량':'TF_COUNT'})
raw_df['UID'] = raw_df['UID'].map(lambda x:x.split('-')[0])
raw_df['DATE'] = raw_df['DATETIME'].map(lambda x:'T'+x.split(' ')[-1])
raw_df['TF_TYPE'] = raw_df['TF_TYPE'].map(lambda x:x.split('Apheresis ')[-1].split('Washed ')[-1].split(' ')[0]).map({'FFP':'FFP', 'RBC':'RBC', 'Red':'Red', 'PLT':'PLT', 'Cryoprecipitate':'Cryoprecipitate','Whole':'Whole Blood','지정혈액1':'Whole Blood','자가혈액1':'Whole Blood', '자가혈액2':'Whole Blood', 'Plateletpheresis':'PLT'})
raw_df = raw_df.drop(['REC_TYPE','DATETIME'], axis=1)
raw_df = raw_df.dropna(subset=['VALUE'])
raw_df = raw_df[['UID','DATE','TF_TYPE','TF_COUNT']].copy()

if not os.path.exists(f'{output_dir}'):
    os.mkdir(f'{output_dir}')
raw_df.to_csv(f"{output_dir}/lnz_final_bodysize_df.csv", encoding='utf-8-sig', index=False)

# raw_df['기록종류']