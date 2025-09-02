from tools import *
from pynca.tools import *

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'
# prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"



conc_df = pd.read_csv(f"{output_dir}/final_conc_df(with sampling).csv")
# conc_df['ID'].drop_duplicates()
conc_df['DV'] = conc_df['CONC'].copy()
conc_df['MDV'] = 0
conc_df['TIME'] = conc_df['SAMP_DT'].copy()
conc_df['CMT'] = 1
conc_df['AMT'] = '.'
conc_df['RATE'] = '.'
# conc_df = conc_df[conc_df['SAMP_DT'].isna()]
print(f"Deleted / NaN CONC: {len(conc_df[conc_df['SAMP_DT'].isna()])} rows / {len(conc_df[conc_df['SAMP_DT'].isna()]['ID'].unique())} patients")
conc_df = conc_df[~conc_df['SAMP_DT'].isna()].copy()
# conc_df[~conc_df['SAMP_DT'].isna()]['ID'].drop_duplicates()
# conc_df = conc_df[~conc_df['REC_REASON'].isin(['오더비고반영', '결과비고반영(시간AP분)', '결과비고반영(시간_분AP)', '결과비고반영(날짜_시간AP)', '결과비고반영(날짜_시간)'])].copy()
# conc_df[conc_df['TIME'].map(lambda x:x.split('T')[-1]=='NN:NN')]

dose_df = pd.read_csv(f"{output_dir}/final_dose_df.csv")
dose_df['DV'] = '.'
dose_df['MDV'] = 1
dose_df['TIME'] = dose_df['DATE']+'T'+dose_df['TIME']
dose_df['CMT'] = 1
dose_df['AMT'] = dose_df['DOSE']
dose_df['RATE'] = dose_df['AMT'] / 0.5
# dose_df[dose_df['ID']==10036912][['TIME','AMT']]
print(f"Deleted / Ambiguous dosing time : {len(dose_df[(dose_df['TIME'].map(lambda x:x.split('T')[-1]=='NN:NN'))])} rows / {len(dose_df[(dose_df['TIME'].map(lambda x:x.split('T')[-1]=='NN:NN'))]['ID'].unique())} patients")
dose_df = dose_df[~(dose_df['TIME'].map(lambda x:x.split('T')[-1]=='NN:NN'))]
dose_df['REC_REASON'] = ''

com_cols = ['ID','NAME','TIME','DV','MDV','CMT','AMT','RATE']
datacheck_basic_cols = ['ID','NAME','TIME','DV','MDV','CMT','AMT','RATE','REC_REASON']
# id_df['TIME'].unique()

modeling_df = pd.concat([conc_df[datacheck_basic_cols], dose_df[datacheck_basic_cols]]).sort_values(['ID','TIME'])
# modeling_df = pd.concat([conc_df[com_cols], dose_df[com_cols]]).sort_values(['ID','TIME'])
modeling_df = modeling_df[~modeling_df['TIME'].isna()].copy()
modeling_df = modeling_df.sort_values(['ID','TIME']).reset_index(drop=True)

modeling_df['UID'] = modeling_df['ID']
modeling_df['DVfloat'] = modeling_df['DV'].replace('.',np.nan).map(float)
modeling_df['DATETIME'] = modeling_df['TIME']


## 환자 수 필터링 (나이 / HD 필터링)

# pt_info = pd.read_csv(f"{output_dir}/patient_info.csv",encoding='utf-8-sig')
# pt_info['AGE'] = pt_info['AGE'].map(lambda x: float(x.replace('개월',''))/12 if '개월' in x else float(x.replace('세','')))
# adult_pids = pt_info[pt_info['AGE'] >= 19].copy()['ID']
mindt_df = modeling_df.groupby('UID',as_index=False)['DATETIME'].min()
mindt_df['MIN_DATE'] = mindt_df['DATETIME'].map(lambda x:x.split('T')[0])
mindt_df = mindt_df[['UID','MIN_DATE']].copy()

pt_info = pd.read_csv(f"{output_dir}/final_req_ptinfo_data.csv",encoding='utf-8-sig')

agedt_df = mindt_df.merge(pt_info, on=['UID'], how='left')
agedt_df['AGE_AT_1ST_DOSE'] = ((pd.to_datetime(agedt_df['MIN_DATE'])-pd.to_datetime(agedt_df['BIRTH_DATE'])).dt.total_seconds()/(365.25*86400))
adult_pids = agedt_df[agedt_df['AGE_AT_1ST_DOSE'] >= 19]['UID'].drop_duplicates()

hd_df = pd.read_csv(f"{output_dir}/final_hd_df.csv",encoding='utf-8-sig')
hd_df.rename(columns={'ID':'UID'})
hd_pids = hd_df['ID'].drop_duplicates()

modeling_df = modeling_df[(modeling_df['ID'].isin(adult_pids))&(~(modeling_df['ID'].isin(hd_pids)))].copy()
# modeling_df['UID'].iloc[0]
# modeling_df[modeling_df['UID']==24910778]
# modeling_df[modeling_df['AMT']==25]

# modeling_df['UID'] = modeling_df['ID']
# len(modeling_df[modeling_df['TIME'].map(lambda x:x[:4] > '2007')]['ID'].unique())

# val_count_df = modeling_df.groupby('ID')['MDV'].value_counts().unstack(fill_value=0).reset_index(drop=False)
# val_count_df.columns.name = None
# no_conc_pids = val_count_df[val_count_df[0]==0]['ID'].drop_duplicates()  # CONC Data 가 존재하지 않는 경우
# no_dose_pids = val_count_df[val_count_df[1]==0]['ID'].drop_duplicates()  # DOSE Data 가 존재하지 않는 경우
# print(f"CONC Data 부재: {len(no_conc_pids)} patients")
# print(f"DOSE Data 부재: {len(no_dose_pids)} patients")
# no_oneside_pids = pd.concat([no_conc_pids, no_dose_pids])
# print(f"CONC or DOSE Data 부재: {len(no_oneside_pids)} patients")
# modeling_df = modeling_df[~modeling_df['ID'].isin(no_oneside_pids)].copy()


modeling_datacheck_df = list()
no_conc_pids = set()
no_dose_pids = set()
neg_time_pids = set()
short_dosing_interval_pids = set()
mod_zero_non_dv_pids = set()
del_zero_non_dv_pids = set()
abnm_conc_moving_pids = set()
abnm_conc_moving_row_count = 0
for id, id_df in modeling_df.groupby('ID'): #break
    try:
        first_dose_dt = id_df[id_df['MDV']==1]['TIME'].iloc[0]
    except:
        # raise ValueError
        # DOSING DATA가 없는 경우 해당
        continue
    id_df['TIME'] = id_df['TIME'].map(lambda x:(datetime.strptime(x,'%Y-%m-%dT%H:%M')-datetime.strptime(first_dose_dt,'%Y-%m-%dT%H:%M')).total_seconds()/3600)

    ## TIME < 0 필터링
    neg_time_rows = id_df[(id_df['TIME'] < 0)]
    if len(neg_time_rows)>0:
        neg_time_pids.add(id)
        # len(neg_time_pids)
        continue

    ## DV > 100 필터링 (채혈라인/투약라인 일치 or SAMPLING DATA 오류 가능성)
    # id_df = id_df[~(id_df['DVfloat'] > 100)].copy()
    id_df = id_df[~(id_df['DVfloat'] > 80)].copy()

    ## TIME=0 DV=0 추가
    zero_conc_row = id_df.iloc[:1, :].copy()
    zero_conc_row['DV'] = 0.0
    zero_conc_row['DVfloat'] = 0.0
    zero_conc_row['MDV'] = 0
    zero_conc_row['CMT'] = 1
    zero_conc_row['AMT'] = '.'
    zero_conc_row['RATE'] = '.'
    id_df = pd.concat([zero_conc_row, id_df]).sort_values(['ID', 'DATETIME', 'MDV'], ignore_index=True)

    ## 한 ID 내에서 농도 측정 사이에 투약기록 없는데, 농도가 증가하게 측정된 경우
    # if id==10313051:
    #     raise ValueError
    abnm_conc_moving_pids.add(id)
    abnm_conc_moving_row_count += len(id_df[((id_df['DVfloat'] / (id_df['DVfloat'].shift(1)) > 2) & (id_df['MDV'] == 0) & (id_df['MDV'].shift(1) == 0))])

    id_df = id_df[~((id_df['DVfloat']/(id_df['DVfloat'].shift(1)) > 2) & (id_df['MDV']==0) & (id_df['MDV'].shift(1)==0))].copy()

    ## 투약 96시간 초과시 끊는다
    # next_cycle_rows = id_df[id_df['MDV'] == 1].copy()
    next_cycle_rows = id_df[id_df['TIME'].diff() >= 96]
    if len(next_cycle_rows)>0:
        next_cycle_inx = next_cycle_rows.index[0]
        id_df = id_df[id_df.index < next_cycle_inx].reset_index(drop=True)
        # raise ValueError

    ## TIME=0 / DV=0 이외 DV 값 부재시 ID 삭제
    only_zero_dv_rows = id_df[id_df['DVfloat']>0]
    if len(only_zero_dv_rows)==0:
        no_conc_pids.add(id)
        continue
        # raise ValueError

    ## TIME=0 / DV=0 이외 DOSE 기록 부재시 ID 삭제
    only_zero_dosing_rows = id_df[id_df['MDV']==1]
    if len(only_zero_dosing_rows)==0:
        no_dose_pids.add(id)
        continue


    ## TIME=0 일때 DV!=0 기록 존재시 BLQ 이하(0.9)으로 측정되었으면 0으로 추정, BLQ 이상이면 ID 자체 삭제
    zero_non_dv_rows = id_df[(id_df['MDV']==0)&(id_df['TIME']<=2)&(id_df['DVfloat']>0)].copy()
    if len(zero_non_dv_rows)>0:
        if zero_non_dv_rows['DV'].max() >= 1:
            del_zero_non_dv_pids.add(id)
            continue
        else:
            # raise ValueError
            mod_zero_non_dv_pids.add(id)
            id_df = id_df[~(id_df.index).isin(zero_non_dv_rows.index)].reset_index(drop=True)

    # 투약 간격이 너무 좁은 경우 (2h 이하)의 투약 테이터들 중 최신 것들만 남기고 삭제
    # prev_id_df = id_df
    # first_cycle = True
    #
    # while True:
    #     prev_id_df = id_df.copy()
    #
    #     short_dosing_conc_rows = id_df[(id_df['MDV']==0)].copy()
    #     short_dosing_dose_rows = id_df[(id_df['MDV']==1)].reset_index(drop=True)
    #     short_dosing_interval_rows = short_dosing_dose_rows.copy()
    #     short_dosing_interval_rows = short_dosing_interval_rows[(short_dosing_interval_rows['TIME'].diff() <= 2)].copy()
    #
    #     if len(short_dosing_interval_rows) > 0:
    #
    #         short_dosing_interval_pids.add(id)
    #         short_dosing_dose_rows = short_dosing_dose_rows[~(short_dosing_dose_rows.index).isin([sdi_inx - 1 for sdi_inx in short_dosing_interval_rows.index])].copy()
    #         id_df = pd.concat([short_dosing_conc_rows,short_dosing_dose_rows]).sort_values(['ID', 'DATETIME', 'MDV'], ignore_index=True)
    #
    #         # prev_id_df = id_df.copy()
    #         # raise ValueError
    #     if len(prev_id_df)==len(id_df):
    #         # raise ValueError
    #         break
    #
    modeling_datacheck_df.append(id_df)

# print(f"CONC Data 부재: {len(no_conc_pids)} patients")
# print(f"DOSE Data 부재: {len(no_dose_pids)} patients")
no_oneside_pids = no_conc_pids.union(no_dose_pids)


modeling_datacheck_df = pd.concat(modeling_datacheck_df).sort_values(['ID', 'TIME', 'MDV'], ignore_index=True)

id_dict = {uid:id for id, uid in enumerate(modeling_datacheck_df['UID'].drop_duplicates())}
modeling_datacheck_df['ID'] = modeling_datacheck_df['UID'].map(id_dict)
# raise ValueError
# raise ValueError
# pt_info.columns
filt_based_on_date = pt_info.rename(columns={'ID':'UID'})[['UID','TDM_REQ_DATE']].copy()
modeling_datacheck_df = modeling_datacheck_df.merge(filt_based_on_date, on='UID', how='left')
modeling_datacheck_df['TDM_YEAR'] = modeling_datacheck_df['TDM_REQ_DATE'].map(lambda x:x.split('-')[0])

modeling_datacheck_dir = f"{output_dir}/amk_modeling_datacheck"
if not os.path.exists(modeling_datacheck_dir):
    os.mkdir(modeling_datacheck_dir)

modeling_datacheck_df.to_csv(f"{modeling_datacheck_dir}/amk_modeling_datacheck.csv",index=False, encoding='utf-8-sig')

# final_modeling_df.columns

# final_modeling_df = modeling_datacheck_df[com_cols+['UID','TDM_YEAR']].drop(['NAME'],axis=1)
# final_modeling_df = modeling_datacheck_df[datacheck_basic_cols+['UID','DATETIME']].drop(['NAME'],axis=1)
# final_modeling_df.to_csv(f"{modeling_datacheck_dir}/amk_modeling_df_tdmyrs.csv",index=False, encoding='utf-8-sig')
# final_modeling_df = modeling_datacheck_df[datacheck_basic_cols+['UID']].drop(['NAME'],axis=1)
final_modeling_df = modeling_datacheck_df[com_cols+['UID']].drop(['NAME'],axis=1)
final_modeling_df.to_csv(f"{modeling_datacheck_dir}/amk_modeling_df.csv",index=False, encoding='utf-8-sig')

print(f"[Check Point] Modeling Dataset: {len(modeling_df['ID'].unique())} patients")
print(f"CONC or DOSE Data 부재: {len(no_oneside_pids)} patients / NO CONC: {len(no_conc_pids)} + NO DOSE: {len(no_dose_pids)}")
print(f"Short dosing interval: {len(short_dosing_interval_pids)} patients")
print(f"Modified (Non Zero CONC at TIME=0): {len(mod_zero_non_dv_pids)} patients")
print(f"Deleted (Non Zero CONC at TIME=0): {len(del_zero_non_dv_pids)} patients")
print(f"Abnormal conc movement: {len(abnm_conc_moving_pids)} patients / {abnm_conc_moving_row_count} rows")
print(f"[Completed] Final Modeling Dataset: {len(modeling_datacheck_df['ID'].unique())} patients / {len(final_modeling_df['ID'])} rows")


recent_modeling_df = modeling_datacheck_df[modeling_datacheck_df['DATETIME'].map(lambda x:x.split('T')[0][0:4]) > '2020'].copy()
recent_modeling_df.to_csv(f"{modeling_datacheck_dir}/recent_amk_modeling_datacheck.csv",index=False, encoding='utf-8-sig')
recent_modeling_df.drop(['NAME'],axis=1).to_csv(f"{modeling_datacheck_dir}/recent_amk_modeling_df.csv",index=False, encoding='utf-8-sig')

# len(recent_modeling_df.drop(['NAME'],axis=1))

past_modeling_df = modeling_datacheck_df[modeling_datacheck_df['DATETIME'].map(lambda x:x.split('T')[0][0:4]) < '2008'].copy()
past_modeling_df.to_csv(f"{modeling_datacheck_dir}/past_amk_modeling_datacheck.csv",index=False, encoding='utf-8-sig')
past_modeling_df[com_cols+['UID']].drop(['NAME'],axis=1).to_csv(f"{output_dir}/past_amk_modeling_df.csv",index=False, encoding='utf-8-sig')

####### NONMEM SDTAB

# nmsdtab_df = pd.read_csv(f"{nonmem_dir}/run/sdtab008",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
# nmsdtab_df['ID'] = nmsdtab_df['ID'].astype(int)
# nmsdtab_df['TDM_YEAR'] = nmsdtab_df['TDM_YEAR'].astype(int)
# under_pred_df = nmsdtab_df[(nmsdtab_df['DV'] > 10)&(nmsdtab_df['IPRED'] < 7)].copy()
# over_pred_df = nmsdtab_df[(nmsdtab_df['DV'] < 7)&(nmsdtab_df['IPRED'] > 10)].copy()
# mis_pred_df = pd.concat([under_pred_df, over_pred_df])
#
# filt_modeling_df = final_modeling_df[~(final_modeling_df['ID'].isin(mis_pred_df['ID'].drop_duplicates()))]
# filt_modeling_df.to_csv(f"{output_dir}/amk_modeling_df_filt.csv",index=False, encoding='utf-8-sig')






# 연도별

# for year in np.arange(2003,2025):
#     year_str = str(year)
#     recent_modeling_df = modeling_datacheck_df[modeling_datacheck_df['DATETIME'].map(lambda x: x.split('T')[0][0:4]) > year_str].copy()
#     recent_modeling_df.to_csv(f"{output_dir}/recent_amk_modeling_datacheck_{year_str}.csv", index=False, encoding='utf-8-sig')
#     recent_modeling_df[com_cols + ['UID']].drop(['NAME'], axis=1).to_csv(f"{output_dir}/recent_amk_modeling_df_{year_str}.csv",index=False, encoding='utf-8-sig')
# recent_modeling_df = modeling_datacheck_df[modeling_datacheck_df['DATETIME'].map(lambda x:x.split('T')[0][0:4]) > '2020'].copy()


# # final_modeling_df['ID']
# n_split = 30
# # id 기준 unique한 사람 목록 추출 및 셔플
# unique_ids = final_modeling_df['ID'].drop_duplicates().reset_index(drop=True)
#
# # 총 사람 수와 1/5씩 나누기
# n = len(unique_ids)
# split_size = n // n_split
#
# # 5개로 나누어서 저장
# for i in range(n_split):
#     if i < n_split-1:
#         ids_slice = unique_ids[i * split_size: (i + 1) * split_size]
#     else:  # 마지막 조각은 나머지까지 포함
#         ids_slice = unique_ids[i * split_size:]
#
#     df_slice = final_modeling_df[final_modeling_df['ID'].isin(ids_slice)]
#     df_slice.to_csv(f"{output_dir}/amk_modeling_df_{i + 1}.csv", index=False)


# modeling_datacheck_df[(modeling_datacheck_df['DVfloat'] < 0.4)&(modeling_datacheck_df['TIME'] > 900)]

# print(f"Total for modeling / {len(modeling_df)} rows / {len(modeling_df['ID'].unique())} patients")

# raise ValueError

## [해결] DV > 100 필터링(채혈라인/투약라인 일치 or SAMPLING DATA 오류 가능성)
## [해결] TIME=0 DV=0 추가
## 투약 96시간 초과시 끊는다
## TIME=0 / DV=0 이외 DV 값 부재시 ID 삭제
## TIME=0 인데, DV!=0 아닌 경우 존재 (33213604, 25495777)
## 투약 시간 간격이 너무 작은(<2h) 경우 제외

## 도근분(10041604) -> 채혈시간이 기록 안되어 있는데  접수시간으로 들어간 게 DOSING TIME 관련해서 이상하게 들어감
## 남궁진(10011083) -> PT채혈시간이 있는데, dosing time 연관해서 dosing 30분 후로 바로 투약된 것으로 들어감
## 김금덕(10048706) -> 2004-01-28T9P14 / 2004-01-29T8P39


## DOSE 기록 이상: 1005
