from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

prj_name = 'LNZ'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

raw_df = pd.read_csv(f"{resource_dir}/lnz_cdw_bodysize_data.csv")
# raw_df.columns
# raw_df = raw_df.rename(columns={'환자번호':'UID','수행시간':'ACTING','약품 오더일자':'DATE', '[실처방] 용법':'REGIMEN','[실처방] 처방일수':'DAYS', '[함량단위환산] 1일 처방량':'DAILY_DOSE','[함량단위환산] 1회 처방량':'DOSE', '[실처방] 투약위치':'PLACE',"[실처방] 경로":'ROUTE','약품명(성분명)':'DRUG','[실처방] 처방비고':'ETC_INFO'})
raw_df = raw_df.rename(columns={'환자번호':'UID','간호기록 작성일시':'DATETIME','간호기록 작성일자':'DATE', '기록종류':'REC_TYPE', '항목명':'BODYSIZE','항목값':'VALUE'})
raw_df['UID'] = raw_df['UID'].map(lambda x:x.split('-')[0])
raw_df['TIME'] = raw_df['DATETIME'].map(lambda x:'T'+x.split(' ')[-1])
raw_df['BODYSIZE'] = raw_df['BODYSIZE'].map({'몸무게':'WT','키':'HT','BMI':'BMI'})
raw_df = raw_df.drop(['REC_TYPE','DATETIME'], axis=1)
raw_df = raw_df.dropna(subset=['VALUE'])
raw_df = raw_df[['UID','DATE','TIME','BODYSIZE','VALUE']].copy()
raw_df.to_csv(f"{output_dir}/lnz_final_bodysize_df.csv", encoding='utf-8-sig', index=False)

# raw_df['기록종류']