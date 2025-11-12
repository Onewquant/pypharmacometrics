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

# LOS data
los_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/PROCEDURE_AND_OUTCOME/LOS,MORTALITY,ICU.csv")
los_df = los_df.rename(columns={'ID':'UID'})
los_df['UID'] = los_df['UID'].map(lambda x: x.split('-')[0])
# los_df[['LOS','ADD_DATE','DIS_DATE']]
# los_df.columns
# los_df[['UID','LOS','ADD_DATE','DIS_DATE']]
los_df = los_df[~los_df['ADD_DATE'].isna()][['UID','LOS','ADD_DATE','DIS_DATE']].copy()

# DIS_DATE 비어있는 데이터 채우기 (ADD_DATE + ICU 이용 days)
for inx, row in los_df.iterrows():
    if type(row['DIS_DATE'])!=str:
        los_df.at[inx, 'DIS_DATE'] = (datetime.strptime(row['ADD_DATE'],'%Y-%m-%d')+timedelta(days=int(row['LOS']))).strftime('%Y-%m-%d')

los_df['DATE'] = los_df.apply(lambda x:pd.date_range(x['ADD_DATE'],x['DIS_DATE']).strftime('%Y-%m-%d').tolist(), axis=1)
los_df['LOSMOT_CATNUM'] = 1
los_df['LOSMOT'] = 'LOS'
los_df = los_df[['UID','DATE','LOSMOT_CATNUM','LOSMOT']].explode('DATE').drop_duplicates(['UID','DATE'])

# Motality
mot_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/PROCEDURE_AND_OUTCOME/LOS,MORTALITY,ICU.csv")
mot_df = mot_df.rename(columns={'ID':'UID'})
mot_df['UID'] = mot_df['UID'].map(lambda x: x.split('-')[0])
# mot_df[['LOS','ADD_DATE','DIS_DATE']]
# mot_df.columns
# mot_df[['UID','LOS','ADD_DATE','DIS_DATE']]
mot_df = mot_df[~mot_df['DEATH_DATE'].isna()][['UID', 'MORTALITY', 'DEATH_DATE']].copy()

# # DIS_DATE 비어있는 데이터 채우기 (ADD_DATE + ICU 이용 days)
# for inx, row in mot_df.iterrows():
#     if type(row['DIS_DATE'])!=str:
#         mot_df.at[inx, 'DIS_DATE'] = (datetime.strptime(row['ADD_DATE'],'%Y-%m-%d')+timedelta(days=int(row['LOS']))).strftime('%Y-%m-%d')

mot_df['DATE'] = mot_df['DEATH_DATE']
mot_df['LOSMOT_CATNUM'] = 2
mot_df['LOSMOT'] = 'MOTALITY'
mot_df = mot_df[['UID','DATE','LOSMOT_CATNUM','LOSMOT']].drop_duplicates(['UID','DATE'])


losmot_df = pd.concat([los_df, mot_df])
losmot_df = losmot_df.drop_duplicates(['UID','LOSMOT_CATNUM','DATE'])

losmot_df['UID'] = losmot_df['UID'].map(str)
losmot_df = losmot_df.merge(pid_decode_df, on=['UID'],how='left')
losmot_df['UID'] = losmot_df['PID'].copy()
losmot_df = losmot_df.drop(['PID'], axis=1)
losmot_df = losmot_df.sort_values(['LOSMOT_CATNUM','UID','DATE'])


if not os.path.exists(f'{output_dir}'):
    os.mkdir(f'{output_dir}')
losmot_df.to_csv(f"{output_dir}/final_locmotality_df.csv", encoding='utf-8-sig', index=False)
