from tools import *
from pynca.tools import *
from datetime import datetime

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'
# prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

# cm_df = pd.read_csv(f"{output_dir}/final_comorbidity_df.csv")
# comed_df = pd.read_csv(f"{output_dir}/final_comed_df.csv")
# comed_df[comed_df['NEPHTOX_DRUG_YN']==1]['CCM_CATNUM'].unique()
# cm_df
#
# df['']
#
# uids = list(cm_df[(cm_df['CM_CATNUM']==18)]['UID'].unique())
# len(uids)
# df[df['UID'].isin(uids)]
#
# test_pt = df[df['UID']==10182459].iloc[0]
# test_pt = df.iloc[0]
# test_pt_cm.columns
#
# test_pt_cm = cm_df[(cm_df['UID'] == test_pt['UID'])&(cm_df['CM_CATNUM']==18)].copy()
# test_pt_cm[test_pt_cm['DATE']>'2008-05-15']

df = pd.read_csv(f"{output_dir}/final_mlres_data.csv")
# set(df.columns)
# df['HEPATIC INSUFFICIENCY'].sum()
# df['SEPSIS'].sum()
# df['SHOCK'].sum()
# df['ACIDOSIS'].sum()
# df['DIABETES'].sum()
# df['ANEMIA'].sum()
#
# df[df['ACIDOSIS']==1]
#
# df['ACIDOSIS'].unique()
# df['AT LEAST ONE NEPHROTOXIC AGENT'].unique()
# df['AT LEAST ONE NEPHROTOXIC AGENT'].sum()
# df['LOOP DIURETIC'].sum()

rid_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/재식별 파일.csv")
rid_df = rid_df.rename(columns={'Deidentification_ID':'UID','환자번호':'PID'})

mort_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/MORTAL.csv")
mort_df = mort_df.rename(columns={'환자번호':'PID','사망일자':'MOT_DATE'}).merge(rid_df, on=['PID'], how='left').drop(['PID'], axis=1)

los_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/LOS.csv")
los_df = los_df.rename(columns={'환자번호':'PID'}).merge(rid_df, on=['PID'], how='left').drop(['PID'], axis=1)

# los_df[los_df['재원일수']==1]
# los_df[los_df['재원일수'].isna()]

uid_admdate_dict = {uid:set() for uid in los_df[los_df['재원일수']>1]['UID'].drop_duplicates()}
# los_df[los_df['재원일수']>1]['UID']

prev_uid = 0
uid_count = 0
for inx, row in los_df[los_df['재원일수']>1].iterrows():
    if prev_uid!=row['UID']:
        print(f"({uid_count}) {row['UID']}")
        uid_count+=1
    uid = row['UID']
    if type(row['입원일자'])==str:
        if ((row['재원일수']==832) and row['UID']==21151121) or ((row['재원일수']==1889) and row['UID']==33100962) or ((row['재원일수']==143) and row['UID']==35983985):
            row['재원일수'] = 1.0
            row['입원일자'] = row['수진일자']
            row['퇴원일자'] = (datetime.strptime(row['수진일자'],'%Y-%m-%d') + timedelta(days=row['재원일수']-1)).strftime('%Y-%m-%d')

            los_df.at[inx, '입원일자'] = row['입원일자']
            los_df.at[inx, '퇴원일자'] = row['퇴원일자']
            los_df.at[inx, '재원일수'] = row['재원일수']
            continue
        else:
            pass
    else:
        if np.isnan(row['입원일자']):
            row['입원일자'] = row['수진일자']
            row['퇴원일자'] = (datetime.strptime(row['수진일자'],'%Y-%m-%d') + timedelta(days=row['재원일수']-1)).strftime('%Y-%m-%d')

            los_df.at[inx, '입원일자'] = row['입원일자']
            los_df.at[inx, '퇴원일자'] = row['퇴원일자']
        else:
            raise ValueError
    adm_dates = set(pd.date_range(row['입원일자'], row['퇴원일자']).strftime('%Y-%m-%d'))
    uid_admdate_dict[uid] = uid_admdate_dict[uid].union(adm_dates)
    prev_uid = row['UID']

rest_df = pd.DataFrame(columns=['수진일자', '재원일수', '입원일자', '퇴원일자', 'UID'])
count = 0
for uid, adm_dates in uid_admdate_dict.items(): #break
    print(f"({count}) {uid} / 차이 dates 제거 중")
    undefined_uid_df = los_df[(los_df['UID']==uid)&(~(los_df['재원일수']>1))].copy()
    unique_uid_dates = set(undefined_uid_df['수진일자'].drop_duplicates())-(uid_admdate_dict[uid])
    rest_uid_df = undefined_uid_df[undefined_uid_df['수진일자'].isin(unique_uid_dates)].drop_duplicates()
    rest_uid_df['재원일수']=1
    rest_uid_df['입원일자']=rest_uid_df['수진일자']
    rest_uid_df['퇴원일자']=rest_uid_df['수진일자']
    rest_df = rest_df.append(rest_uid_df)
    count+=1

# los_df[(los_df['UID']==21151121)&(los_df['재원일수']==2)]

more_than_one_df = los_df[los_df['재원일수']>1].copy()
gen_los_df=pd.concat([more_than_one_df,rest_df])

# df['LOS'] = np.nan
los_list = list()
mort_list = list()
for inx, row in df.iterrows():
    uid = row['UID']
    aki_date = row['AKI_DATE']
    uid_los_df = gen_los_df[gen_los_df['UID']==uid].copy()
    uid_los_row_df = uid_los_df[(uid_los_df['입원일자'] <= aki_date)&(uid_los_df['퇴원일자'] >= aki_date)].copy()

    uid_mot_df = mort_df[mort_df['UID'] == uid].copy()

    if len(uid_mot_df)==0:
        mort_yn = 0
    else:
        if len(uid_los_row_df)==0:
            mort_yn = 0
        else:
            uid_los_row = uid_los_row_df.iloc[0]
            if (uid_los_row['입원일자'] <= uid_mot_df.iloc[0]['MOT_DATE']) and (uid_los_row['퇴원일자'] >= uid_mot_df.iloc[0]['MOT_DATE']):
                mort_yn = 1
            else:
                mort_yn = 0

    if len(uid_los_row_df)==0:
        print(f'({inx}) {uid} / no row')
        uid_los_days = np.nan
    else:
        print(f'({inx}) {uid} / >= 1 rows')
        uid_los_days = uid_los_row_df.iloc[0]['재원일수']



    mort_list.append(mort_yn)
    los_list.append(uid_los_days)

df['LOS'] = los_list
df['MORTAL'] = mort_list

drop_list = df.loc[:,'BL_LAB_DATE':'LAST_DOSE_DATE'].columns
final_df = df.drop(drop_list, axis=1).rename(columns={'AKI_DATE':'DATE'})
final_df.to_csv(f"{output_dir}/final_mlres_data5.csv", index=False, encoding='utf-8-sig')
# (datetime.strptime('2024-03-07','%Y-%m-%d')-datetime.strptime('2020-08-19','%Y-%m-%d')).days