from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

prj_name = 'LNZ'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

raw_df = pd.read_csv(f"{resource_dir}/lnz_cdw_ptinfo_data.csv")
# raw_df.columns
# raw_df = raw_df.rename(columns={'환자번호':'UID','수행시간':'ACTING','약품 오더일자':'DATE', '[실처방] 용법':'REGIMEN','[실처방] 처방일수':'DAYS', '[함량단위환산] 1일 처방량':'DAILY_DOSE','[함량단위환산] 1회 처방량':'DOSE', '[실처방] 투약위치':'PLACE',"[실처방] 경로":'ROUTE','약품명(성분명)':'DRUG','[실처방] 처방비고':'ETC_INFO'})
raw_df = raw_df.rename(columns={'환자번호':'UID','성별':'SEX','생년월일':'BIRTH_DATE', '최종 주소지':'ADDRESS'})
raw_df['UID'] = raw_df['UID'].map(lambda x:x.split('-')[0])

if not os.path.exists(f'{output_dir}'):
    os.mkdir(f'{output_dir}')
raw_df.to_csv(f"{output_dir}/lnz_final_ptinfo_df.csv", encoding='utf-8-sig', index=False)

# raw_df['기록종류']