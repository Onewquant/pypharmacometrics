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

pid_decode_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/재식별 파일.csv")
pid_decode_df = pid_decode_df.rename(columns={'환자번호':'UID','Deidentification_ID':'PID'})
pid_decode_df['UID'] = pid_decode_df['UID'].map(lambda x: x.split('-')[0])

surgery_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/PROCEDURE_AND_OUTCOME/CARDIAC_SURGERY,SURGERY.csv")
surgery_df = surgery_df.rename(columns={'ID':'UID'})
surgery_df['UID'] = surgery_df['UID'].map(lambda x: x.split('-')[0])
surgery_df['PROC_CATNUM'] = surgery_df['CODE'].astype(float).map(lambda x: 1 if x in (35.2, 35.22, 35.24, 35.28, 36.1, 36.11, 36.12, 37.12, 37.31, 37.33, 37.34, 37.5, 37.62, 37.66, 37.8) else 2)
surgery_df['PROC'] = surgery_df['NAME'].copy()
surgery_df = surgery_df[['UID','DATE','PROC_CATNUM','PROC']].copy()

# surgery_df.columns
icu_mortal_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/PROCEDURE_AND_OUTCOME/LOS,MORTALITY,ICU.csv")
icu_mortal_df = icu_mortal_df.rename(columns={'ID':'UID'})
icu_mortal_df['UID'] = icu_mortal_df['UID'].map(lambda x: x.split('-')[0])
# icu_mortal_df.columns
# icu_mortal_df[['UID','LOS','ICU','ADD_DATE','DIS_DATE']]
icu_mortal_df = icu_mortal_df[~icu_mortal_df['ICU'].isna()][['UID','LOS','ICU','ADD_DATE','DIS_DATE']].reset_index(drop=True)
icu_mortal_df['PROC_CATNUM'] = 3
icu_mortal_df['PROC'] = 'ICU'
# DIS_DATE 비어있는 데이터 채우기 (ADD_DATE + ICU 이용 days)
for inx, row in icu_mortal_df.iterrows():
    if type(row['DIS_DATE'])!=str:
        icu_mortal_df.at[inx, 'DIS_DATE'] = (datetime.strptime(row['ADD_DATE'],'%Y-%m-%d')+timedelta(days=int(row['LOS']))).strftime('%Y-%m-%d')

icu_mortal_df['DATE'] = icu_mortal_df.apply(lambda x:pd.date_range(x['ADD_DATE'],x['DIS_DATE']).strftime('%Y-%m-%d').tolist(), axis=1)
icu_mortal_df = icu_mortal_df[['UID','DATE','PROC_CATNUM', 'PROC']].explode('DATE').drop_duplicates(['UID','DATE']).sort_values(['UID','DATE'], ignore_index=True)

# icu_mortal_df[icu_mortal_df['ADD_DATE'].isna()]

# x = {'ADD_DATE':'2008-11-11','DIS_DATE':'2009-05-26'}

mv_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/PROCEDURE_AND_OUTCOME/MV.csv", encoding='euc-kr')
mv_df = mv_df.rename(columns={'ID':'UID'})
mv_df['UID'] = mv_df['UID'].map(lambda x: x.split('-')[0])
mv_df = mv_df[mv_df['MLCODE'] > 0][['UID','DATE']].copy()
mv_df['PROC_CATNUM'] = 4
mv_df['PROC'] = 'Mechanical Ventilation'
mv_df = mv_df[['UID','DATE','PROC_CATNUM', 'PROC']].copy()


proc_df = pd.concat([surgery_df, icu_mortal_df, mv_df])
proc_df = proc_df.drop_duplicates(['UID','DATE','PROC_CATNUM'])

proc_df['UID'] = proc_df['UID'].map(str)
proc_df = proc_df.merge(pid_decode_df, on=['UID'],how='left')
proc_df['UID'] = proc_df['PID'].copy()
proc_df = proc_df.drop(['PID'], axis=1)

proc_index_dict = {1:(365,1),2:(365,1),3:(365,1),4:(14,1)}
proc_df['PROC_PERIOD_FROM'] = proc_df['PROC_CATNUM'].map(lambda x:proc_index_dict[x][0])
proc_df['PROC_PERIOD_UNTIL'] = proc_df['PROC_CATNUM'].map(lambda x:proc_index_dict[x][1])


if not os.path.exists(f'{output_dir}'):
    os.mkdir(f'{output_dir}')
proc_df.to_csv(f"{output_dir}/final_procedure_df.csv", encoding='utf-8-sig', index=False)
