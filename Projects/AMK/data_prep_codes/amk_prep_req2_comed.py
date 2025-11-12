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

df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/CCM.csv")

# df.columns
df = df.rename(columns={'환자번호': 'UID', '수행시간': 'ACTING', '약품 오더일자': 'DATE', '[실처방] 용법': 'REGIMEN', '[실처방] 처방일수': 'DAYS', '[함량단위환산] 1일 처방량': 'DAILY_DOSE', '[실처방] 1회 처방량': 'DOSE', '[실처방] 투약위치': 'PLACE', "[실처방] 경로": 'ROUTE', '약품명(성분명)': 'DRUG', '[실처방] 처방비고': 'ETC_INFO','ATC코드':'ATC'})
df = df[(~df['약품명(일반명)'].isna())&(~df['DOSE'].isna())].copy()
df['[함량단위환산] 1회 처방량'].unique()
df[~df['[함량단위환산] 1회 처방량'].isna()]['[함량단위환산] 1회 처방량'].map(lambda x:x[-3:]).drop_duplicates()

# try:
#     raw_df['DOSE'] = raw_df.apply(lambda x: float(re.findall(r'\d+[\.]?\d*', x['DOSE'])[0]) if ('mg' in x['DOSE']) else (float(re.findall(r'\d+[\.]?\d*', x['DOSE'])[0]) * float(re.findall(r'\d+[\.]?\d*', x['약품명(일반명)'])[0])), axis=1)
# except:
#     raw_df['DOSE'] = raw_df.apply(lambda x: float(re.findall(r'\d+[\.]?\d*', x['DOSE'])[0]) if ('mg' in x['DOSE']) else (float(re.findall(r'\d+[\.]?\d*', x['[함량단위환산] 1회 처방량'])[0])), axis=1)

df['DOSE'] = df['DOSE'].astype(str)
df['UID'] = df['UID'].map(lambda x: x.split('-')[0])
# df.iloc[0]
# raw_df['UID'] = raw_df['UID'].map(lambda x: x.split('-')[0])
# raw_df['DRUG'] = drugname

    # df = raw_df[~(raw_df['ROUTE'].isin(['Irrigation', 'L-spine']))].copy()
df = df[~(df['DAYS'].isin(['관류하세요', '관류하세요.']))].copy()
df = df[~(df['REGIMEN'].map(lambda x: x.strip()).isin(['관류하세요']))].copy()
df = df[~((df['DRUG'].isna())&(df['주성분명'].isna()))].copy()
    # df = df.rename(columns={'[실처방] 투약위치':'PLACE'})
    # raw_df[raw_df['DRUG']=='lacosamide']
    # df['[함량단위환산] 1회 처방량'].unique()

alter_acting_dict = {
    '1일 1회': ('q24h', '09:00/Y'),
    '1일 2회': ('q12h', '09:00/Y, 21:00/Y'),
    '1일 3회': ('q8h', '08:00/Y, 16:00/Y, 23:59/Y'),
    '1일 4회': ('q6h', '04:00/Y, 10:00/Y, 16:00/Y, 22:00/Y'),
    '1일 5회': ('q4.8h', '04:00/Y, 09:00/Y, 13:00/Y, 18:00/Y, 23:00/Y'),
    '1일 6회': ('q4h', '08:00/Y, 12:00/Y, 16:00/Y, 20:00/Y, 23:59/Y'),
    '1일 7회': ('q3.5h', '02:00/Y, 05:00/Y, 08:00/Y, 11:00/Y, 14:00/Y, 18:00/Y, 22:00/Y'),
    '1일 8회': ('q3h', '02:00/Y, 05:00/Y, 08:00/Y, 11:00/Y, 14:00/Y, 17:00/Y, 20:00/Y, 23:00/Y'),
    '1일 9회': ('q3h', '02:00/Y, 05:00/Y, 08:00/Y, 11:00/Y, 14:00/Y, 17:00/Y, 20:00/Y, 22:00/Y, 23:00/Y'),
    '1일 10회': ('q2.4h', '02:00/Y, 04:00/Y, 05:00/Y, 08:00/Y, 11:00/Y, 14:00/Y, 17:00/Y, 19:00/Y, 21:00/Y, 23:00/Y'),
    '1일 11회': ('q2.2h', '02:00/Y, 04:00/Y, 05:00/Y, 08:00/Y, 10:00/Y, 11:00/Y, 14:00/Y, 17:00/Y, 19:00/Y, 21:00/Y, 23:00/Y'),
    '1일 12회': ('q2h', '02:00/Y, 04:00/Y, 05:00/Y, 08:00/Y, 10:00/Y, 11:00/Y, 14:00/Y, 16:00/Y, 17:00/Y, 19:00/Y, 21:00/Y, 23:00/Y'),
    '필요시 자기전에': ('q[hours]h', '22:00/Y'),
    '필요시 복용하세요': ('q[hours]h', '09:00/Y'),
    '필요시 1시간전에': ('q[hours]h', '09:00/Y'),
    '격일로 자기전': ('q48h', '22:00/Y'),
    '격일로 저녁': ('q48h', '21:00/Y'),
    '격일로 저녁7시에': ('q48h', '19:00/Y'),
    '격일로 아침': ('q48h', '09:00/Y'),
    '격일로 48시간마다': ('q48h', '09:00/Y'),
    '2주 1회': ('q336h', '09:00/Y'),
    '3주 1회': ('q504h', '09:00/Y'),
    '4주 1회': ('q672h', '09:00/Y'),
    '1달 1회': ('q672h', '09:00/Y'),
    '1년 1회': ('q[hours]h', '09:00/Y'),
    '1주 3회': ('q56h', '09:00/Y'),
    '1주 1회': ('q168h', '09:00/Y'),
    '1주 2회': ('q84h', '09:00/Y'),
    '1주 5회': ('q33h', '09:00/Y'),
    '3일 1회': ('q72h', '09:00/Y'),
    '5일 1회': ('q120h', '09:00/Y'),
    '의사지시대로 관류하세요': ('q[hours]h', '09:00/Y'),
    '발열시 복용하세요': ('q24h', '09:00/Y'),
    '필요시 주사하세요': ('q24h', '09:00/Y'),
    '넣으세요': ('q24h', '09:00/Y'),
    '수시로 복용하세요': ('q24h', '09:00/Y'),
    '통증이 있을때': ('q24h', '09:00/Y'),
    '검사전에 지시대로': ('q24h', '09:00/Y'),
}
# df.columns
# df[df['UID']=='352504997098455']
# df['REGIMEN'].unique()
# df['DOSE'].unique()
# df[~df['ATC코드'].isna()]
# df['ATC코드명']
dose_list = list()
alter_acting_list = list()
interval_list = list()
ondemand_list = list()
count = 0
prev_uid = '0000'
# df.iloc[0]
for inx, row in df.iterrows():
    if prev_uid != row['UID']:
        count += 1
        drugname = row['DRUG'] if type(row['DRUG'])==str else row['주성분명'].split(' ')[0]
        drugcode = row['ATC']
        print(f"({count}) {row['DRUG']} / {row['UID']}")
        prev_uid = row['UID']
    ## '250mg-500mg-500mg' 형태로 dose가 쓰여있는 경우
    # if row['']
    seq_dose_patterns = re.findall(r'\d+\.*\d*', row['DOSE'])
    if len(seq_dose_patterns) >= 1:
        dose_val = ('_'.join(seq_dose_patterns))
    else:
        dose_val = seq_dose_patterns[0]
        # dose_val = re.findall(r'\d+',x.split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0]
        # raise ValueError
    dose_list.append(dose_val)

    daily_count_str = ' '.join(row['REGIMEN'].split(' ')[:2])

    if row['REGIMEN'].split(' ')[0] == '필요시':
        # print(f"{row['REGIMEN']} / hour_str: {str(int(24/len(seq_dose_patterns)))} / DOSE: {row['DOSE']} ({row['DAILY_DOSE']})")
        ondemand_list.append(True)
    else:
        ondemand_list.append(False)

    alter_acting_list.append(alter_acting_dict[daily_count_str][-1])
    hour_str = str(int(24 / len(seq_dose_patterns)))
    interval_list.append(alter_acting_dict[daily_count_str][0].replace('[hours]', hour_str))

# df['ROUTE'].unique()
df['DRUG'] = df.apply(lambda row:row['DRUG'] if type(row['DRUG'])==str else row['주성분명'].split(' ')[0], axis=1)
df['DOSE'] = dose_list
df['ALTER_ACTING'] = alter_acting_list
df['INTERVAL'] = interval_list
df['ON_DEMAND'] = ondemand_list
# df['DRUG'] = row['DRUG']

# df['DRUG']
df['ROUTE'] = df['ROUTE'].map(
    {'P.O': 'PO', 'MIV': 'IV', 'PLT': 'PO', 'IVS': 'IV', 'SC': 'SC', 'S.L': "SL", 'IntraPleural': 'IPLEU',
     'IntraPeritoneal': 'IPERI'})  # PLT: ( Par L-tube) / S.L: sublingual
df['ADDL'] = df['DAYS'].copy()
# df = df[df['ROUTE']!='Irrigation'].copy()

df = df.sort_values(['UID', 'DATE']).reset_index(drop=True)
# df[df['UID']==df['UID'].unique()[2]][['DATE','DOSE','ROUTE','INTERVAL','ADDL','ACTING','PLACE']]
# raise ValueError

dose_cond1 = ~(df['ACTING'].isna())  # 본원투약기록 -> 투약 기록대로 추가
dose_cond2 = (df['ADDL'] > 1) & (df['ACTING'].isna())  # 외래자가투약 / 타병원 투약 처방 -> 날짜대로 ALTER_ACTING으로 추가
# dose_cond3 = (df['ACTING'].isna())

df1 = df[dose_cond1].sort_values(['UID', 'DATE']).reset_index(drop=True)
df1['VIRTUALITY'] = False
df2 = df[dose_cond2].sort_values(['UID', 'DATE']).reset_index(drop=True)
df2['ACTING'] = df2['ALTER_ACTING'].copy()
df2['VIRTUALITY'] = True

sdf = pd.concat([df1, df2]).sort_values(['UID', 'DATE', 'ACTING'])[['UID', 'DATE', 'DRUG', 'ATC', 'DOSE', 'ROUTE', 'INTERVAL', 'ADDL', 'ACTING', 'PLACE']].copy()

sdf_row_count = 0
comed_res_df = list()
for uid_drugatc, uid_sdf in sdf.groupby(['UID','ATC']):  #break
    uid = uid_drugatc[0]
    atc_code = uid_drugatc[1]
    drugname = uid_sdf['DRUG'].iloc[0]
    # if len(set(uid_sdf['DRUG'])) > 1:
    #     raise ValueError
    # raise ValueError
    # if uid=='155505674746153':
    #     raise ValueError
    # frag_sdf1 = frag_sdf[frag_sdf['ADDL']==1].copy()
    # frag_sdf1_acting_rows = frag_sdf1[frag_sdf1['ACTING'].map(lambda x: ('/Y' in x) or ('/O' in x))].copy()
    # if len(frag_sdf1_acting_rows)==0:
    #     frag_sdf1_dates = set()
    # else:
    #     frag_sdf1_dates = set(frag_sdf1_acting_rows['DATE'])
    #
    # frag_sdf2 = frag_sdf[frag_sdf['ADDL']>1].copy()
    uid_sdf = uid_sdf[uid_sdf['ACTING'].map(lambda x: False if (('/C' in x) and ('/Y' not in x)) else True)]
    uid_sdf['II_NUM'] = uid_sdf['INTERVAL'].map(lambda x: float(x.replace('q', '').replace('h', '')))
    uid_sdf['II_APPLIED'] = uid_sdf['II_NUM'].map(lambda x: x / 24 if x > 24 else 1.0)

    # frag_sdf1_dates = set(frag_sdf1[frag_sdf1['ACTING'].map(lambda x: ('/Y' in x) or ('/O' in x))]['DATE'])
    # frag_sdf2_dates = set()
    uid_sdf_dates = set()
    for sdf_inx, uid_sdf_row in uid_sdf.iterrows():  # break
        # if uid_sdf_row['ADDL'] == 1:
        #     raise ValueError
        # if uid_sdf_row['II_APPLIED'] > 1:
        #     raise ValueError
        # dates_every_step_days(start_date=uid_sdf_row['DATE'], day_step=uid_sdf_row['II_APPLIED'], addl=uid_sdf_row['ADDL'], fmt="%Y-%m-%d", include_start=True)
        uid_sdf_dates = uid_sdf_dates.union(set(dates_every_step_days(start_date=uid_sdf_row['DATE'], day_step=uid_sdf_row['II_APPLIED'], addl=uid_sdf_row['ADDL'], fmt="%Y-%m-%d", include_start=True)))

    uid_dates_rows = pd.DataFrame(uid_sdf_dates, columns=['DATE']).sort_values(['DATE'], ignore_index=True)
    uid_dates_row = {'UID': uid, 'ATC':atc_code, 'DRUG':drugname, 'DATE_LIST': tuple(uid_dates_rows['DATE'])}
    # drug_res_df.append(uid_dates_row)
    comed_res_df.append(uid_dates_row)

    print(f"({sdf_row_count}) / {drugname} / {uid} / {len(uid_sdf_dates)}")
    sdf_row_count += 1

    # drug_res_df = pd.DataFrame(drug_res_df)
comed_res_df = pd.DataFrame(comed_res_df)

ccm_index_df_ori = pd.read_excel(f"{resource_dir}/[AMK_AKI_ML_DATA]/CCM_index_table.xlsx")
ccm_index_df = ccm_index_df_ori.iloc[:-1,:].copy()
# ccm_index_df.columns
ccm_index_df['CAT_ATC_LIST'] = ccm_index_df['CAT_ATC_LIST'].map(lambda x:x.split(', '))
ccm_index_df = ccm_index_df.explode(column='CAT_ATC_LIST')
ccm_index_df['NEPHTOX_DRUG_YN'] = (ccm_index_df['CAT_NUM'].isin([1, 2, 4, 11, 12, 13, 14, 15]))*1
cat_num_df = ccm_index_df[['CAT_ATC_LIST','CAT_NUM']].copy()
ren_tox_df = ccm_index_df[['CAT_ATC_LIST','NEPHTOX_DRUG_YN']].copy()
cat_num_dict = cat_num_df.set_index(['CAT_ATC_LIST'])['CAT_NUM'].to_dict()
ren_tox_dict = ren_tox_df.set_index(['CAT_ATC_LIST'])['NEPHTOX_DRUG_YN'].to_dict()

comed_res_df['CCM_CATNUM'] = comed_res_df['ATC'].map(cat_num_dict)
comed_res_df['NEPHTOX_DRUG_YN'] = comed_res_df['ATC'].map(ren_tox_dict)
# comed_res_df[~comed_res_df['CCM_CATNUM'].isna()]

# comed_res_df = pd.read_csv(f"{output_dir}/final_comed_df.csv")
comed_res_df['UID'] = comed_res_df['UID'].map(str)
comed_res_df = comed_res_df.merge(pid_decode_df, on=['UID'],how='left')
comed_res_df['UID'] = comed_res_df['PID'].copy()
comed_res_df = comed_res_df.drop(['PID'], axis=1)
comed_res_df = comed_res_df[~((comed_res_df['CCM_CATNUM'].isna())&(comed_res_df['NEPHTOX_DRUG_YN'].isna()))].reset_index(drop=True)
comed_res_df['CCM_CATNUM'] = comed_res_df['CCM_CATNUM'].map(int)
comed_res_df['NEPHTOX_DRUG_YN'] = comed_res_df['NEPHTOX_DRUG_YN'].map(int)
# comed_res_df['UID'].iloc[0]
# pid_decode_df['UID'].iloc[0]

if not os.path.exists(f'{output_dir}'):
    os.mkdir(f'{output_dir}')
comed_res_df.to_csv(f"{output_dir}/final_comed_df.csv", encoding='utf-8-sig', index=False)
# comed_res_df['DRUG'].drop_duplicates()
# raw_df['기록종류']
