from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

ibd_type_df = pd.read_excel(f"{resource_dir}/IBD-PGx_induction_date.xlsx")
ibd_type_df = ibd_type_df.rename(columns={'EMR ID':'ID','name':'NAME','induction_start_date':'IND_START_DATE','IBD type':'IBD_TYPE'})
# ibd_type_df['IND_START_DATE'] = ibd_type_df['IND_START_DATE'].astype(str).replace('NaT','')
#
# induction_df2 = pd.read_csv(f"{resource_dir}/initial_induction_patients.csv")
# induction_df2 = induction_df2.rename(columns={'DATETIME':'START_INDMAINT'})
# induction_df2['START_INDMAINT'] = induction_df2['START_INDMAINT'].astype(str).replace('NaT','')
#
# induction_df = induction_df[['ID','NAME','IND_START_DATE','IBD_TYPE']].merge(induction_df2[['ID','START_INDMAINT']], on=['ID'], how='left')
# induction_df['IND_START_DATE'] = induction_df.apply(lambda x: x['START_INDMAINT'] if type(x['START_INDMAINT'])==str else x['IND_START_DATE'], axis=1)
# induction_df['START_INDMAINT'] = (induction_df['START_INDMAINT'].isna())*1
# inddf = induction_df[~induction_df['START_INDMAINT'].isna()].reset_index(drop=True)
# inddf[inddf['IND_START_DATE']!=inddf['START_INDMAINT']]

totlab_df = pd.read_csv(f"{output_dir}/conc_df(lab).csv")
cumlab_df = pd.read_csv(f"{output_dir}/conc_df(cum_lab).csv")

# totlab_df[totlab_df['ID']==21169146]
# cumlab_df[cumlab_df['ID']==21169146]

# totlab_df[totlab_df['ID']==21911051]
# cumlab_df[cumlab_df['ID']==21911051]
# raise ValueError
lab_df = pd.concat([cumlab_df, totlab_df[['오더일','보고일']]], axis=1)
lab_df['접수일'] = lab_df['DATETIME'].map(lambda x:x.split(' ')[0])
lab_df['접수시간'] = lab_df['DATETIME'].map(lambda x:x.split(' ')[1])
receipt_list = list()
for inx, row in lab_df.iterrows():
# for inx, row in lab_df[lab_df['ID']==21169146].iterrows():

    # row = lab_df[lab_df['ID']==21169146].iloc[1]
    # row = lab_df[lab_df['ID'] == 21169146].iloc[-2]
    if (datetime.strptime(row['보고일'],'%Y-%m-%d')-datetime.strptime(row['오더일'],'%Y-%m-%d')).days > 14:
        date_str = min(row['보고일'],row['접수일'])
        receipt_list.append(date_str+f"T{row['접수시간']}")
    else:
        if np.abs((datetime.strptime(row['오더일'],'%Y-%m-%d')-datetime.strptime(row['접수일'],'%Y-%m-%d')).days) <= 14:
            date_str = min(row['오더일'], min(row['보고일'], row['접수일']))
        else:
            date_str = max(row['보고일'],row['접수일'])
        receipt_list.append(date_str+f"T{row['접수시간']}")
lab_df['DATETIME'] = receipt_list
# lab_df['DATETIME'] = lab_df['DATETIME'].map(lambda x:x.re place(' ','T'))
lab_df = lab_df.sort_values(['ID','DATETIME'])
lab_df.to_csv(f"{output_dir}/conc_df.csv", encoding='utf-8-sig', index=False)

# lab_df[(lab_df['CONC']!=lab_df['CONC_cumlab'])|(lab_df['DRUG']!=lab_df['DRUG_cumlab'])][['DRUG','CONC','DRUG_cumlab','CONC_cumlab']]
lab_df = pd.read_csv(f"{output_dir}/conc_df.csv")
# lab_df[lab_df['ID']==11788526]
# lab_df.columns
# lab_df['DRUG']
lab_df = lab_df.rename(columns={'CONC':'DV'})
lab_df['MDV']=0
lab_df['DUR']='.'
lab_df['AMT']='.'
lab_df['ROUTE']='.'

dose_df = pd.read_csv(f"{output_dir}/dose_df.csv")
# dose_df.columns
# dose_df['DOSE']
dose_df = dose_df.rename(columns={'DOSE':'AMT'})
dose_df['DV'] = '.'
dose_df['MDV'] = 1
dose_df['DUR'] = 1

## Dosing 패턴으로 Induction / Maintenacne 구분
total_indmaint_df = list()
# for drug in ['infliximab', 'adalimumab', 'ustekinumab']: #break
for drug in ['adalimumab']: #break
    # raise ValueError
    # '36959194'
    # indmaint_df[indmaint_df['ID']==36959194]
    indmaint_df = dose_df[dose_df['DRUG']==drug].copy()
    # dose_df[dose_df['DRUG'] == drug].to_csv(f'{output_dir}/modeling_df_datacheck/{drug}_dose_datacheck.csv', index=False, encoding='utf-8-sig')
    init_dosing_dt = indmaint_df.groupby('ID',as_index=False)['DATETIME'].min().rename(columns={"DATETIME":"INIT_DOSE_DT"})
    indmaint_df = indmaint_df.merge(init_dosing_dt, on=['ID'], how='left')
    # indmaint_df.columns
    indmaint_df['DOSING_TIME'] = indmaint_df.apply(lambda x:(datetime.strptime(x['DATETIME'],'%Y-%m-%dT%H:%M') - datetime.strptime(x['INIT_DOSE_DT'],'%Y-%m-%dT%H:%M')), axis=1).dt.total_seconds() / 86400 / 7
    indmaint_df['DOSING_TIME_WKNORM'] = indmaint_df['DOSING_TIME'].map(lambda x:round(x,0))
    # indmaint_df.columns
    indmaint_df[indmaint_df['ID']==(indmaint_df['ID'].drop_duplicates().iloc[50])][['ID','NAME','DATETIME','AMT','DOSING_TIME_WKNORM']]
    indmaint_df[indmaint_df['DOSING_TIME_WKNORM']<3].groupby('ID',as_index=False)['AMT'].sum() ################################
    # raise ValueError
    indmaint_list = list()
    for uid, id_indmaint_df in indmaint_df.groupby('ID'):
        id_indmaint_df['DOSING_INTERVAL'] = id_indmaint_df['DOSING_TIME'].diff()
        id_indmaint_df = id_indmaint_df.dropna(subset=['DOSING_INTERVAL'])
        indmaint_list.append(id_indmaint_df)
    indmaint_df = pd.concat(indmaint_list).sort_values(['ID','DATETIME'], ignore_index=True)
    indmaint_df['DRUG'] = drug

    # Infliximab: 첫 번째 ~ 두 번째 투약 간격이 4 미만 / IV 제제만 Induction으로 인정
    if drug=='infliximab':
        indmaint_df = indmaint_df.groupby('ID', as_index=False).agg({'DRUG': 'first','DOSING_INTERVAL': 'first', 'ROUTE': 'first'})
        indmaint_df['START_INDMAINT'] = (~((indmaint_df['DOSING_INTERVAL'] < 4)&(indmaint_df['ROUTE']=='IV')))*1
    elif drug=='adalimumab':
        indmaint_df = indmaint_df.groupby('ID', as_index=False).agg({'DRUG': 'first','DOSING_INTERVAL': 'first', 'ROUTE': 'first'})
        # raise ValueError
        indmaint_df['START_INDMAINT'] = (~((indmaint_df['DOSING_INTERVAL'] < 4)&(indmaint_df['ROUTE']=='SC')))*1
    elif drug=='ustekinumab':
        indmaint_df = indmaint_df.groupby('ID', as_index=False).agg({'DRUG': 'first','DOSING_INTERVAL': 'first', 'ROUTE': 'first'})
        indmaint_df['START_INDMAINT'] = (~((indmaint_df['DOSING_INTERVAL'] < 4)&(indmaint_df['ROUTE']=='IV')))*1

    total_indmaint_df.append(indmaint_df)
total_indmaint_df = pd.concat(total_indmaint_df).sort_values(['ID','DRUG'],ignore_index=True)
init_dose_dt_df = dose_df.groupby(['ID','DRUG'],as_index=False)['DATETIME'].min().rename(columns={"DATETIME":"INIT_DOSE_DT"})
total_indmaint_df = total_indmaint_df.merge(init_dose_dt_df, on=['ID','DRUG'], how='left')
total_indmaint_df['INIT_DOSE_DATE'] = total_indmaint_df['INIT_DOSE_DT'].map(lambda x:x.split('T')[0])
# total_indmaint_df = total_indmaint_df.rename(columns={'DATETIME':'IND_START_DATE'})
# total_indmaint_df[total_indmaint_df['DRUG']=='adalimumab']
# sns.displot(indmaint_df[(indmaint_df['ROUTE']=='IV')]['DOSING_INTERVAL'])

# indmaint_df[(indmaint_df['DOSING_INTERVAL']<4)&(indmaint_df['DOSING_INTERVAL']=='IV')]
# id_indmaint_df.columns

################################ Clinical Data 받은 것과 Comparison #######################

# inf_initdose_df = total_indmaint_df[total_indmaint_df['DRUG']=='infliximab'].copy()
# inf_induction_df = pd.read_csv(f"{resource_dir}/initial_induction_patients.csv")
# inf_maintenance_df = pd.read_csv(f"{resource_dir}/initial_maintenance_patients.csv")
#
# inf_maint_check_df = inf_initdose_df[inf_initdose_df['ID'].isin(inf_maintenance_df['ID'])].copy()
# inf_maint_check_df = inf_maint_check_df.merge(inf_maintenance_df[['ID','NAME','DRUG_START_DATE']].rename(columns={'DRUG_START_DATE':'CLIN_IND_START_DATE'}))
# filt_inf_maint_check_df = inf_maint_check_df[(inf_maint_check_df['CLIN_IND_START_DATE']==inf_maint_check_df['INIT_DOSE_DATE'])&(inf_maint_check_df['ROUTE']=='IV')].copy()
# # sns.displot(filt_inf_maint_check_df['DOSING_INTERVAL'])

# inf_maint_check_df['FILT'] = ((inf_maint_check_df['CLIN_IND_START_DATE']==inf_maint_check_df['INIT_DOSE_DATE'])&(inf_maint_check_df['ROUTE']=='IV'))
# inf_maint_check_df.to_csv(f"{output_dir}/data comparison/clin_maintenance_patients_comparison.csv", index=False, encoding='utf-8-sig')
# filt_inf_maint_check_df2 = filt_inf_maint_check_df[~(filt_inf_maint_check_df['ID'].isin(inf_induction_df['ID']))]
# filt_inf_maint_check_df2.to_csv(f"{output_dir}/data comparison/clin_maintenance_patients_comparison2.csv", index=False, encoding='utf-8-sig')
# sns.displot(filt_inf_maint_check_df2['DOSING_INTERVAL'])
# # inf_maint_check_df.columns
##########################################################################################

mediator_cols = ['ID','NAME','DRUG','ROUTE','DATETIME','DV','MDV','AMT','DUR']
lab_df = lab_df[mediator_cols].reset_index(drop=True)
dose_df = dose_df[mediator_cols].reset_index(drop=True)
merged_df = pd.concat([lab_df,dose_df]).sort_values(['ID','DATETIME'], ignore_index=True)

merged_df.to_csv(f'{output_dir}/merged_df.csv',index=False, encoding='utf-8-sig')
# merged_df.drop_duplicates(['ID'])

merged_df['DATE'] = merged_df['DATETIME'].map(lambda x:x.split('T')[0])
# merged_df['AZERO'] = 0
merged_df = merged_df.merge(ibd_type_df[['ID', 'IBD_TYPE']], on=['ID'], how='left')

## ULOQ 넘는 농도값은 제거
merged_df = merged_df[~((merged_df['DV']==48)&((merged_df['ID'].isin([21911051,])) & (merged_df['DATE'].isin(['2019-02-22',]))))].copy()
# [21911051, 34019533, 32482522]

# Induction Phase 불일치 환자 구분 (전체 합친 데이터에서)

# min_dose_df = merged_df.groupby(['ID','DRUG']).agg({'NAME':'min','DATETIME':'min'}).reset_index(drop=False)
# min_dose_df['MIN_DOSE_DATE'] = min_dose_df['DATETIME'].map(lambda x:x.split('T')[0])
# comp_df = min_dose_df.merge(total_indmaint_df[['ID','DRUG','START_INDMAINT']], on=['ID','DRUG'], how='left')
# comp_df = min_dose_df.merge(total_indmaint_df[['ID','DRUG','IND_START_DATE']], on=['ID','DRUG'], how='left')

# comp_df


## Adalimumab은 IND/MAINT 구분 또 다르게 해야할 듯!
# ind_pids = induction_df[induction_df['START_INDMAINT']==0]['ID'].reset_index(drop=True)
# maint_pids = induction_df[induction_df['START_INDMAINT']==1]['ID'].reset_index(drop=True)

# maint_pids = induction_df[induction_df['START_INDMAINT']==0]['ID'].reset_index(drop=True)
# ind_pids = induction_df[induction_df['START_INDMAINT']==1]['ID'].reset_index(drop=True)

drug = 'adalimumab'
maint_cons_df = total_indmaint_df[(total_indmaint_df['DRUG']==drug)&(total_indmaint_df['START_INDMAINT']==0)].reset_index(drop=True)
maint_diff_df = total_indmaint_df[(total_indmaint_df['DRUG']==drug)&(total_indmaint_df['START_INDMAINT']==1)].reset_index(drop=True)

# total_indmaint_df



# raise ValueError
# maint_diff_df[maint_diff_df['ID']==36959194]
# maint_cons_df[maint_cons_df['ID']==36959194]
# maint_cons_df = comp_df[comp_df['MIN_DOSE_DATE']==comp_df['IND_START_DATE']].reset_index(drop=True)
# maint_diff_df = comp_df[comp_df['MIN_DOSE_DATE']!=comp_df['IND_START_DATE']].reset_index(drop=True)

# 17439372 -> 23.04.06 infliximab conc: 2.4 이며 타원에서 투약하다 전원됨
# 37625588 -> 22.12.07 부터 infliximab 투약 / 타원에서 투약하다 전원됨
# 35093356 -> 20.11.02 부터 infliximab 투약 / 타원에서 투약하다 전원됨
# 36898756 -> 23.01.03 부터 infliximab 투약 / 타원에서 투약하다 전원됨


# set(uids_from_diff_to_cons)-set(maint_diff_df['ID'])
# set(uids_in_cons)-set(maint_cons_df['ID'])

# uids_in_cons = [25269024, 24590129, 27530683]
# uids_from_diff_to_cons = [32067581, 29336700, 35028484, 36325931]
#
# for moving_uid in uids_from_diff_to_cons:
#     moving_df = maint_diff_df[maint_diff_df['ID']==moving_uid].copy()
#     if len(moving_df)==0:
#         continue
#     maint_diff_df = maint_diff_df[maint_diff_df['ID']!=moving_uid].copy()
#     maint_cons_df = pd.concat([maint_cons_df, moving_df])
#
# maint_cons_df['IND_START_DATE'] = maint_cons_df['MIN_DOSE_DATE']

appended_frag_cols = ['UID', 'NAME', 'DRUG', 'ROUTE', 'TIME', 'DV', 'MDV', 'AMT', 'DUR', 'CMT', 'DATETIME','IBD_TYPE','START_INDMAINT']

ind_df = list()
no_indconc_df = list()
no_indid_df = list()
no_maintid_df = list()
for inx, row in maint_cons_df.iterrows(): #break
    # if row['ID']==36959194:
    #     raise ValueError

    # if inx
    # break
    # start_date = row['IND_START_DATE']
    # end_date = (datetime.strptime(start_date,'%Y-%m-%d') + timedelta(days=127)).strftime('%Y-%m-%d')
    id_df = merged_df[merged_df['ID']==row['ID']].copy()
    if (len(id_df)==0):
        no_indid_df.append(pd.DataFrame([row]))
        continue

    # if row['ID']==10875838:
    #     raise ValueError
    # ind_df_frag = id_df[(id_df['DATE'] >= start_date)].copy()
    ind_df_frag = id_df.copy()

    if (len(ind_df_frag['MDV'].unique())<=1):
        no_indconc_df.append(pd.DataFrame([row]))
        continue

    ind_df_frag = ind_df_frag.sort_values(['ID', 'DATE', 'MDV'], ascending=[True, True, False])
    unique_date_df = ind_df_frag.groupby(['DATE'])['DV'].count().reset_index(drop=False)
    # unique_date_df = unique_date_df[unique_date_df['DV'] >= 2].reset_index(drop=True)
    unique_date_df = unique_date_df[unique_date_df['DV'] >= 2].reset_index(drop=True)
    if len(unique_date_df)>0:
        # if row['ID']==34019533: raise ValueError
        ind_df_frag_list = list()
        for dup_date_inx, dup_date_row in unique_date_df.iterrows():

            dup_date = dup_date_row['DATE']
            ind_date_df_frag = ind_df_frag[ind_df_frag['DATE']==dup_date].copy()
            # raise ValueError
            nodup_df_frag = ind_df_frag[~(ind_df_frag['DATE'].isin(unique_date_df['DATE']))].copy()
            dup_df_frag = ind_date_df_frag[ind_date_df_frag['DATE']==dup_date].copy()

            dupdv_df_frag = dup_df_frag[dup_df_frag['MDV']==0].copy().sort_values(['DV'])
            dupmdv_df_frag = dup_df_frag[dup_df_frag['MDV'] !=0].copy()
            # dupdv_df_frag['DV']
            if len(dupdv_df_frag)>=3:
                print("> 3")
                raise ValueError
            elif (len(dupdv_df_frag)==2) and (len(dupmdv_df_frag)==1):
                # if len(dupdv_df_frag) == 2: raise ValueError
                dupdv_df_frag.at[dupdv_df_frag.index[0],'DATETIME'] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0],'%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
                dupdv_df_frag.at[dupdv_df_frag.index[-1],'DATETIME'] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0], '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
            elif (len(dupdv_df_frag)==1) and (len(dupmdv_df_frag)==2):
                dupdv_df_frag.at[dupdv_df_frag.index[0],'DATETIME'] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0],'%Y-%m-%dT%H:%M')-timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
            elif (len(dupdv_df_frag)==1) and (len(dupmdv_df_frag)==1):
                dupdv_df_frag.at[dupdv_df_frag.index[0],'DATETIME'] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0],'%Y-%m-%dT%H:%M')-timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
            elif (len(dupdv_df_frag)==0) or (len(dupmdv_df_frag)==0):
                # print("No data")
                # raise ValueError
                pass
            else:
                raise ValueError

            ind_df_frag_list.append(pd.concat([nodup_df_frag, dupdv_df_frag, dupmdv_df_frag]))

        ind_df_frag = pd.concat(ind_df_frag_list).drop_duplicates(['ID','DATETIME']).sort_values(['ID','DATETIME'], ignore_index=True)

    for drug_inx, drug in enumerate(ind_df_frag['DRUG'].drop_duplicates()):
        zero_conc_row = ind_df_frag[ind_df_frag['DRUG']==drug].iloc[:1,:].copy()
        zero_conc_row['DV'] = 0.0
        zero_conc_row['MDV'] = 0
        zero_conc_row['AMT'] = '.'
        zero_conc_row['DUR'] = '.'
        zero_conc_row['ROUTE'] = '.'
        ind_df_frag = pd.concat([zero_conc_row, ind_df_frag]).sort_values(['ID','DATETIME'], ignore_index=True)

    # if row['ID']==11788526:
    #     bef_df = ind_df_frag.iloc[:4,:].copy()
    #     aft_df1 = ind_df_frag.iloc[4:5,:].copy()
    #     aft_df2 = ind_df_frag.iloc[5:,:].copy()
    #     aft_df1['DATETIME'] = '2022-05-10T12:18'
    #     aft_df2['DATETIME'] = '2022-05-10T12:16'
    #     ind_df_frag = pd.concat([bef_df, aft_df2, aft_df1]).sort_values(['ID','DATETIME'], ignore_index=True)
    #     # raise ValueError


    ind_date_str = ind_df_frag['DATETIME'].iloc[0]
    ind_df_frag['TIME'] = ind_df_frag['DATETIME'].map(lambda x:(datetime.strptime(x,'%Y-%m-%dT%H:%M') - datetime.strptime(ind_date_str,'%Y-%m-%dT%H:%M'))).dt.total_seconds() * 24 / 86400

    # ind_df_frag['WKTIME'] = ind_df_frag['TIME'] / 7
    # ind_df_frag['DWKTIME'] = ind_df_frag['WKTIME'].diff().fillna(0.0)
    ind_df_frag['CMT'] = ind_df_frag['ROUTE'].map(lambda x: 2 if x!='SC' else 1)
    ind_df_frag['START_INDMAINT'] = row['START_INDMAINT']
    ind_df_frag = ind_df_frag.rename(columns={'ID':'UID'})
    # ind_df_frag['TIME'] = ind_df_frag['TIME'].map(lambda x:x)
    ind_df.append(ind_df_frag[appended_frag_cols])

ind_df = pd.concat(ind_df).reset_index(drop=True)
ind_uniq_df = ind_df.drop_duplicates(['UID'])


maint_df = list()
no_maintconc_df = list()
for inx, row in maint_diff_df.iterrows():

    # if row['ID']==36959194:
    #     raise ValueError

    # if inx
    # break
    # start_date = row['IND_START_DATE']
    # end_date = (datetime.strptime(start_date,'%Y-%m-%d') + timedelta(days=127)).strftime('%Y-%m-%d')
    # id_df['DV']
    id_df = merged_df[merged_df['ID']==row['ID']].copy()
    if (len(id_df)==0):
        no_maintid_df.append(pd.DataFrame([row]))
        continue
    # if (len((id_df['MDV']!=1).sum())<=2):
    #     continue

    # if row['ID']==12681022:
    #     raise ValueError



    # ind_df_frag = id_df[(id_df['DATE'] >= start_date) & (id_df['DATE'] <= end_date)].copy()
    # if len(ind_conc_df_frag)==0: # induction phase에 농도 데이터 부재시 나머지 데이터 중 첫 농도 값을 기준으로 maintenance 시작
    # ind_end_inx = id_df['MDV']
    # maint_df_frag = id_df.loc[ind_end_inx:].copy()
    first_dv_inx = id_df[id_df['MDV']==0].iloc[0].name
    maint_df_frag = id_df.loc[first_dv_inx:].copy()
    # else: # induction phase에 농도 데이터 있으면, induction phase의 마지막 농도 데이터 값을 기준으로 maintenance phase 설정
    #     last_dv_inx = ind_df_frag[ind_df_frag['MDV']=='.'].iloc[-1].name
    #     maint_df_frag = id_df.loc[last_dv_inx:].copy()
    # maint_df_frag.iloc[0]
    # ind_df_frag = id_df[(id_df['DATE'] >= end_date)].sort_values(['DATETIME'], ignore_index=True)
    # first_dv_inx = ind_df_frag[ind_df_frag['MDV']=='.'].iloc[0].name   ###################### 에러남. 확인 요망
    drug_maint_conc_df = maint_df_frag.groupby(['DRUG','MDV'],as_index=False)[['ID']].count().rename(columns={'ID':'COUNT'})
    drug_maint_conc_df = drug_maint_conc_df[(drug_maint_conc_df['MDV']!=1)&(drug_maint_conc_df['COUNT']>1)].copy()

    if len(drug_maint_conc_df)==0:
        no_maintconc_df.append(pd.DataFrame([row]))
        continue
    else:
        maint_df_frag = maint_df_frag[maint_df_frag['DRUG'].isin(list(drug_maint_conc_df['DRUG']))].copy()
    # ind_df_frag = ind_df_frag[ind_df_frag['DATE']=='2023-05-04']['DV']
    # maintenance 기간동안 약을 바꿔서 동시투약 한 사람은 없음
    # if len(list(drug_maint_conc_df['DRUG']))>=2:
    #     raise ValueError
    # else:
    #     continue

    # TROUGH 채혈이 투약보다 먼저 나오도록 설정
    maint_df_frag = maint_df_frag.sort_values(['ID', 'DATE', 'MDV'], ascending=[True, True, False])
    unique_date_df = maint_df_frag.groupby(['DATE'])['DV'].count().reset_index(drop=False)
    unique_date_df = unique_date_df[unique_date_df['DV'] >= 2].reset_index(drop=True)
    if len(unique_date_df)>0:
        # if row['ID']==34019533: raise ValueError
        maint_df_frag_list = list()
        for dup_date_inx, dup_date_row in unique_date_df.iterrows(): # break

            dup_date = dup_date_row['DATE']
            maint_date_df_frag = maint_df_frag[maint_df_frag['DATE']==dup_date].copy()
            # raise ValueError
            nodup_df_frag = maint_df_frag[~(maint_df_frag['DATE'].isin(unique_date_df['DATE']))].copy()
            dup_df_frag = maint_date_df_frag[maint_date_df_frag['DATE']==dup_date].copy()

            dupdv_df_frag = dup_df_frag[dup_df_frag['MDV'] == 0].copy()
            dupmdv_df_frag = dup_df_frag[dup_df_frag['MDV'] != 0].copy()
            if len(dupmdv_df_frag)>0:
                maint_df_frag_list.append(pd.concat([nodup_df_frag, dupdv_df_frag, dupmdv_df_frag]))
                dupdv_df_frag['DATETIME'] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0], '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
            elif len(dupmdv_df_frag)==0:
                maint_df_frag_list.append(pd.concat([nodup_df_frag, dupdv_df_frag]))

        maint_df_frag = pd.concat(maint_df_frag_list).drop_duplicates(['ID','DATETIME']).sort_values(['ID','DATETIME'], ignore_index=True)

    # if row['ID']==29702679:
    #     raise ValueError

    # 약을 여러종류 사용했을경우, 첫 약은 maintenance일 수 있으나, 다음 약부터는 serial 하게 사용하므로 induction부터 시작
    for drug_inx, drug in enumerate(id_df['DRUG'].drop_duplicates()): #break
        if drug_inx==0:
            continue
        else:
            zero_conc_row = maint_df_frag[maint_df_frag['DRUG'] == drug].iloc[:1, :].copy()
            zero_conc_row['DV'] = 0.0
            zero_conc_row['MDV'] = 0
            zero_conc_row['AMT'] = '.'
            zero_conc_row['DUR'] = '.'
            zero_conc_row['ROUTE'] = '.'
            maint_df_frag = pd.concat([zero_conc_row, maint_df_frag]).sort_values(['ID', 'DATETIME'], ignore_index=True)

        # zero_conc_row = ind_df_frag.iloc[:1,:].copy()
        # zero_conc_row['DV'] = 0.0
        # zero_conc_row['MDV'] = '.'
        # zero_conc_row['AMT'] = '.'
        # zero_conc_row['DUR'] = '.'
        # ind_df_frag = pd.concat([zero_conc_row, ind_df_frag]).sort_values(['ID','DATETIME'], ignore_index=True)

    # if row['ID']==11788526:
    #     bef_df = ind_df_frag.iloc[:4,:].copy()
    #     aft_df1 = ind_df_frag.iloc[4:5,:].copy()
    #     aft_df2 = ind_df_frag.iloc[5:,:].copy()
    #     aft_df1['DATETIME'] = '2022-05-10T12:18'
    #     aft_df2['DATETIME'] = '2022-05-10T12:16'
    #     ind_df_frag = pd.concat([bef_df, aft_df2, aft_df1]).sort_values(['ID','DATETIME'], ignore_index=True)
    #     # raise ValueError

    maint_date_str = maint_df_frag['DATETIME'].iloc[0]
    maint_df_frag['TIME'] = maint_df_frag['DATETIME'].map(lambda x:(datetime.strptime(x,'%Y-%m-%dT%H:%M') - datetime.strptime(maint_date_str,'%Y-%m-%dT%H:%M'))).dt.total_seconds() * 24 / 86400

    # maint_df_frag['WKTIME'] = maint_df_frag['TIME'] / 7
    # maint_df_frag['DWKTIME'] = maint_df_frag['WKTIME'].diff().fillna(0.0)
    maint_df_frag['CMT'] = maint_df_frag['ROUTE'].map(lambda x: 2 if x!='SC' else 1)
    maint_df_frag['START_INDMAINT'] = row['START_INDMAINT']
    maint_df_frag = maint_df_frag.rename(columns={'ID':'UID'})
    # ind_df_frag['TIME'] = ind_df_frag['TIME'].map(lambda x:x)
    maint_df.append(maint_df_frag[appended_frag_cols])
    # ind_df_frag.columns
    # raise ValueError

    # IND 중간에 그만 둔 사람?

maint_df = pd.concat(maint_df).reset_index(drop=True)
maint_uniq_df = maint_df.drop_duplicates(['UID'])
if len(no_maintconc_df)==0:
    no_maintconc_df = pd.DataFrame(columns=list(maint_df.columns))
else:
    no_maintconc_df = pd.concat(no_maintconc_df).reset_index(drop=True)


# print(f"# Induction 시작시점 일치: {len(maint_cons_df)} (Infliximab: {len(maint_cons_df[maint_cons_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_cons_df[maint_cons_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_cons_df[maint_cons_df['DRUG']=='ustekinumab'])}) ")
# print(f"# Induction 시작시점 불일치: {len(maint_diff_df)} (Infliximab: {len(maint_diff_df[maint_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_diff_df[maint_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_diff_df[maint_diff_df['DRUG']=='ustekinumab'])}) ")


# print('')
# print(f"# Maintenance 시작시 농도 측정값 존재: {len(maint_uniq_df)} (Infliximab: {len(maint_uniq_df[maint_uniq_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_uniq_df[maint_uniq_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_uniq_df[maint_uniq_df['DRUG']=='ustekinumab'])}) ")
# print(f"# Maintenance 시작시 농도 측정값 부재: {len(no_maintconc_df)} (Infliximab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='infliximab'])} / Adalimumab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='ustekinumab'])}) ")

# raise ValueError
# ind_df[ind_df['UID']==32067581]
# maint_df[maint_df['UID']==32067581]
# 32067581


# inf_ind_df = inf_ind_df
# inf_ind_df


# pd.concat([inf_ind_df, inf_maint_df])['UID'].drop_duplicates().reset_index(drop=True)
# pd.concat([ada_ind_df, ada_maint_df])['UID'].drop_duplicates().reset_index(drop=True)
# maint_df[maint_df['UID']==36959194]['DV']

drug_df_dict = dict()
modeling_cols = ['ID','TIME','DV','MDV','AMT','DUR','CMT','DATETIME','IBD_TYPE','UID','NAME','ROUTE','DRUG','START_INDMAINT']

for drug in set(ind_df['DRUG']).union(set(maint_df['DRUG'])): #break
    # drug_df_dict[drug]
    drug_ind_df = ind_df[ind_df['DRUG']==drug].copy()
    drug_maint_df = maint_df[maint_df['DRUG'] == drug].sort_values(['UID', 'DATETIME'], ignore_index=True)
    drug_df = pd.concat([drug_ind_df, drug_maint_df]).sort_values(['UID','DATETIME'], ignore_index=True)
    drug_df['ID'] = drug_df['UID'].map({uid: uid_inx for uid_inx, uid in enumerate(list(drug_df['UID'].unique()))})
    drug_df = drug_df[modeling_cols].sort_values(['ID','TIME','MDV'])

    ## 약물별로 TIME==0 시작 아닌 경우 0으로 초기화 (e.g. adalimumab -> infliximab -> ustekinumab 으로 사용시 infliximab과 ustekinumab은 TIME!=0 인 상태임)

    drug_nztime_df = drug_df.groupby('UID', as_index=False)['TIME'].min()
    drug_nztime_df = drug_nztime_df[drug_nztime_df['TIME'] != 0].copy()
    for nzt_inx, nzt_row in drug_nztime_df.iterrows():
        uid = nzt_row['UID']
        start_time = nzt_row['TIME']
        drug_rest_of_uids_df = drug_df[~(drug_df['UID'] == uid)].copy()
        drug_uid_df = drug_df[(drug_df['UID'] == uid)].copy()
        drug_uid_df['TIME'] = drug_uid_df['TIME'] - start_time
        drug_df = pd.concat([drug_rest_of_uids_df, drug_uid_df])

    ## SC 자가투약 된 것 투약력 정리

    drug_dose_df = pd.read_csv(f'{output_dir}/dose_df.csv')
    drug_dose_df = drug_dose_df[drug_dose_df['DRUG'] == drug].rename(columns={'ID': 'UID'})
    drug_dose_df['ETC_INFO_TREATED'] = drug_dose_df['ETC_INFO_TREATED'].replace(np.nan, '')
    drug_dose_df['ADDED_ADDL'] = drug_dose_df['ETC_INFO_TREATED'].map(lambda x: 'ADDL반영' in x)
    drug_df = drug_df.merge(drug_dose_df[['UID', 'DATETIME', 'ADDED_ADDL']], on=['UID', 'DATETIME'], how='left')
    drug_df['ADDED_ADDL'] = drug_df['ADDED_ADDL'].replace(np.nan, False)

    drug_del_df = drug_df[(drug_df['MDV'] == 1)].sort_values(['ID', 'TIME'])
    drug_del_df = drug_del_df[(drug_del_df['TIME'].diff() > 0) & (drug_del_df['TIME'].diff() <= 9 * 24) & (drug_del_df['ADDED_ADDL'])].copy()
    drug_df = drug_df[~drug_df.index.isin(drug_del_df.index)].copy()


    # drug_dv_exists_pids = drug_df[(drug_df['TIME']!=0)&(drug_df['DV']!=0)&(drug_df['MDV']!=1)].copy()
    # drug_dv_exists_pids = drug_dv_exists_pids.groupby('UID',as_index=False).agg({'UID':'max','DV':'count'})['UID']
    # drug_no_indconc = set(drug_df[(drug_df['TIME']==0)&(drug_df['DV']==0)&(drug_df['MDV']!=1)].copy()['UID'])-set(drug_dv_exists_pids)

    ## CONC측정이 근처 투여 시간 이후인 경우 투약시간 재조정

    drug_conc_first_rows = drug_df[(drug_df['ID'] == (drug_df['ID'].shift(1))) & (drug_df['TIME'] != 0) & (drug_df['MDV'] == 0) & (drug_df['TIME'].diff(1).map(np.abs) < 120) & (drug_df['TIME'].diff(1).map(np.abs) > 0)].copy()

    for inx, row in drug_conc_first_rows.iterrows():
        near_dinx = inx - 1

        # if (inf_df.at[inx, 'ID']==inf_df.at[near_dinx,'ID']):
        # 근처 Dosing된 포인트가 SC 추가투여로 구성한 데이터일때
        if drug_df.at[near_dinx, 'ADDED_ADDL']:
            drug_df.at[near_dinx, 'TIME'] = row['TIME'] + 0.5
            drug_df.at[near_dinx, 'DATETIME'] = (datetime.strptime(row['DATETIME'], '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
        else:
            # 근처 Dosing 날짜와 같은 날짜이긴 할때
            drug_df.at[inx, 'TIME'] = drug_df.at[near_dinx, 'TIME'] - 0.5
            drug_df.at[inx, 'DATETIME'] = (datetime.strptime(drug_df.at[near_dinx, 'DATETIME'], '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')

    drug_df = drug_df.sort_values(['ID', 'TIME', 'MDV'], ignore_index=True)
    # drug_df[drug_df['TAD']!=0]
    drug_df['TIME(DAY)'] = drug_df['TIME'] / 24
    drug_df['TIME(WEEK)'] = drug_df['TIME'] / (7 * 24)

    ## 저장
    if not os.path.exists(f'{output_dir}/modeling_df_datacheck'):
        os.mkdir(f'{output_dir}/modeling_df_datacheck')
    drug_df.to_csv(f'{output_dir}/modeling_df_datacheck/{drug}_integrated_datacheck.csv', index=False, encoding='utf-8-sig')

    # raise ValueError
    print(f"[{drug}] modeling df / {len(drug_df)} rows")

print(f"\n# 총 환자 수: 140")
print(f"# Induction 시작시점 일치: {len(maint_cons_df)} (Infliximab: {len(maint_cons_df[maint_cons_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_cons_df[maint_cons_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_cons_df[maint_cons_df['DRUG']=='ustekinumab'])}) ")
print(f"# Induction 시작시점 불일치: {len(maint_diff_df)} (Infliximab: {len(maint_diff_df[maint_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_diff_df[maint_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_diff_df[maint_diff_df['DRUG']=='ustekinumab'])}) ")
# print(f"# No conc in Induction 시작시점 일치: {len(no_indconc_df)}")

# drug_df[drug_df['ID']==1]

# print(f"# Induction 시작시점 불일치: {len(maint_diff_df)} (Infliximab: {len(maint_diff_df[maint_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_diff_df[maint_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_diff_df[maint_diff_df['DRUG']=='ustekinumab'])}) ")
# print(f"# No conc in Induction: {len(inf_no_indconc)} (Infliximab: {len(inf_df[(inf_df.isin(inf_no_indconc))&(inf_df['DRUG']=='infliximab')]['UID'].drop_duplicates()} / Adalimumab: {len(inf_no_indconc[inf_no_indconc['DRUG']=='adalimumab'])} / Ustekinumab: {len(inf_no_indconc[inf_no_indconc['DRUG']=='ustekinumab'])}) ")
# print(f"# No conc in Maintenance: {len(no_maintconc_df)} (Infliximab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='infliximab'])} / Adalimumab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='ustekinumab'])}) ")
# raise ValueError
# inf_df = inf_df[inf_df['UID'].isin(inf_dv_exists_pids)].reset_index(drop=True)



# ################ 레거시 코드 ##############
#
# inf_ind_df = ind_df[ind_df['DRUG']=='infliximab'].copy()
# ada_ind_df = ind_df[ind_df['DRUG']=='adalimumab'].copy()
#
# inf_maint_df = maint_df[maint_df['DRUG']=='infliximab'].sort_values(['UID','DATETIME'], ignore_index=True)
# ada_maint_df = maint_df[maint_df['DRUG']=='adalimumab'].sort_values(['UID','DATETIME'], ignore_index=True)
#
# inf_df = pd.concat([inf_ind_df, inf_maint_df]).sort_values(['UID','DATETIME'], ignore_index=True)
# ada_df = pd.concat([ada_ind_df, ada_maint_df]).sort_values(['UID','DATETIME'], ignore_index=True)
#
# inf_df['ID'] = inf_df['UID'].map({uid:uid_inx for uid_inx, uid in enumerate(list(inf_df['UID'].unique()))})
# ada_df['ID'] = ada_df['UID'].map({uid:uid_inx for uid_inx, uid in enumerate(list(ada_df['UID'].unique()))})
#
#
# inf_df = inf_df[modeling_cols].sort_values(['ID','TIME','MDV'])
# ada_df = ada_df[modeling_cols].sort_values(['ID','TIME','MDV'])
#
# ## 약물별로 TIME==0 시작 아닌 경우 0으로 초기화 (e.g. adalimumab -> infliximab -> ustekinumab 으로 사용시 infliximab과 ustekinumab은 TIME!=0 인 상태임)
#
# inf_nztime_df = inf_df.groupby('UID',as_index=False)['TIME'].min()
# inf_nztime_df = inf_nztime_df[inf_nztime_df['TIME']!=0].copy()
# for nzt_inx, nzt_row in inf_nztime_df.iterrows():
#     uid = nzt_row['UID']
#     start_time = nzt_row['TIME']
#     inf_rest_of_uids_df = inf_df[~(inf_df['UID']==uid)].copy()
#     inf_uid_df = inf_df[(inf_df['UID']==uid)].copy()
#     inf_uid_df['TIME']=inf_uid_df['TIME'] - start_time
#     inf_df = pd.concat([inf_rest_of_uids_df , inf_uid_df])
#
# ## SC 자가투약 된 것 투약력 정리
# inf_dose_df = pd.read_csv(f'{output_dir}/dose_df.csv')
# inf_dose_df = inf_dose_df[inf_dose_df['DRUG']=='infliximab'].rename(columns={'ID':'UID'})
# inf_dose_df['ETC_INFO_TREATED'] = inf_dose_df['ETC_INFO_TREATED'].replace(np.nan,'')
# inf_dose_df['ADDED_ADDL'] = inf_dose_df['ETC_INFO_TREATED'].map(lambda x: 'ADDL반영' in x)
# inf_df = inf_df.merge(inf_dose_df[['UID','DATETIME','ADDED_ADDL']], on=['UID','DATETIME'], how='left')
# inf_df['ADDED_ADDL'] = inf_df['ADDED_ADDL'].replace(np.nan, False)
#
#
# inf_del_df = inf_df[(inf_df['MDV']==1)].sort_values(['ID','TIME'])
# inf_del_df = inf_del_df[(inf_del_df['TIME'].diff() > 0)&(inf_del_df['TIME'].diff() <= 9*24) & (inf_del_df['ADDED_ADDL'])].copy()
# inf_df = inf_df[~inf_df.index.isin(inf_del_df.index)].copy()
#
#
#
# ada_dose_df = pd.read_csv(f'{output_dir}/dose_df.csv')
# ada_dose_df = ada_dose_df[ada_dose_df['DRUG']=='adalimumab'].rename(columns={'ID':'UID'})
# ada_dose_df['ETC_INFO_TREATED'] = ada_dose_df['ETC_INFO_TREATED'].replace(np.nan,'')
# ada_dose_df['ADDED_ADDL'] = ada_dose_df['ETC_INFO_TREATED'].map(lambda x: 'ADDL반영' in x)
# ada_df = ada_df.merge(ada_dose_df[['UID','DATETIME','ADDED_ADDL']], on=['UID','DATETIME'], how='left')
# ada_df['ADDED_ADDL'] = ada_df['ADDED_ADDL'].replace(np.nan, False)
#
# ada_del_df = ada_df[(ada_df['MDV']==1)].sort_values(['ID','TIME'])
# ada_del_df = ada_del_df[(ada_del_df['TIME'].diff() > 0)&(ada_del_df['TIME'].diff() <= 9*24) & (ada_del_df['ADDED_ADDL'])].copy()
# ada_df = ada_df[~ada_df.index.isin(ada_del_df.index)].copy()
#
# ## TIME==0 일때 CONC==0 만 DV로 가지고 있는 경우 제외
# # set(inf_df[(inf_df['TIME']==0)&(inf_df['DV']==0)].drop_duplicates('UID')['UID']) - set(inf_df[(inf_df['TIME']==0)&(inf_df['DV']==0)].drop_duplicates('UID',keep=False)['UID'])
#
# # 34665842 -> TIME=0 일때 DV=0 두 번 들어감. 확인필요 !!!!!####################
#
# inf_dv_exists_pids = inf_df[(inf_df['TIME']!=0)&(inf_df['DV']!=0)&(inf_df['MDV']!=1)].copy()
# inf_dv_exists_pids = inf_dv_exists_pids.groupby('UID',as_index=False).agg({'UID':'max','DV':'count'})['UID']
# inf_no_indconc = set(inf_df[(inf_df['TIME']==0)&(inf_df['DV']==0)&(inf_df['MDV']!=1)].copy()['UID'])-set(inf_dv_exists_pids)
#
# print(f"# 총 환자 수: 140")
# print(f"# Induction 시작시점 일치: {len(maint_cons_df)} (Infliximab: {len(maint_cons_df[maint_cons_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_cons_df[maint_cons_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_cons_df[maint_cons_df['DRUG']=='ustekinumab'])}) ")
# print(f"# Induction 시작시점 불일치: {len(maint_diff_df)} (Infliximab: {len(maint_diff_df[maint_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_diff_df[maint_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_diff_df[maint_diff_df['DRUG']=='ustekinumab'])}) ")
# # print(f"# No conc in Induction 시작시점 일치: {len(no_indconc_df)}")
#
# # print(f"# Induction 시작시점 불일치: {len(maint_diff_df)} (Infliximab: {len(maint_diff_df[maint_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_diff_df[maint_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_diff_df[maint_diff_df['DRUG']=='ustekinumab'])}) ")
# # print(f"# No conc in Induction: {len(inf_no_indconc)} (Infliximab: {len(inf_df[(inf_df.isin(inf_no_indconc))&(inf_df['DRUG']=='infliximab')]['UID'].drop_duplicates()} / Adalimumab: {len(inf_no_indconc[inf_no_indconc['DRUG']=='adalimumab'])} / Ustekinumab: {len(inf_no_indconc[inf_no_indconc['DRUG']=='ustekinumab'])}) ")
# print(f"# No conc in Maintenance: {len(no_maintconc_df)} (Infliximab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='infliximab'])} / Adalimumab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='ustekinumab'])}) ")
# # raise ValueError
# # inf_df = inf_df[inf_df['UID'].isin(inf_dv_exists_pids)].reset_index(drop=True)
#
#
#
#
#
# # no_dv.columns
# ## CONC측정이 근처 투여 시간 이후인 경우 투약시간 재조정
#
# inf_conc_first_rows = inf_df[(inf_df['ID']==(inf_df['ID'].shift(1)))&(inf_df['TIME']!=0)&(inf_df['MDV']==0)&(inf_df['TIME'].diff(1).map(np.abs) < 120)&(inf_df['TIME'].diff(1).map(np.abs) > 0)].copy()
# ada_conc_first_rows = ada_df[(ada_df['ID']==(ada_df['ID'].shift(1)))&(ada_df['TIME']!=0)&(ada_df['MDV']==0)&(ada_df['TIME'].diff(1).map(np.abs) < 120)&(ada_df['TIME'].diff(1).map(np.abs) > 0)].copy()
#
# for inx, row in inf_conc_first_rows.iterrows():
#     near_dinx = inx-1
#
#     # if (inf_df.at[inx, 'ID']==inf_df.at[near_dinx,'ID']):
#         # 근처 Dosing된 포인트가 SC 추가투여로 구성한 데이터일때
#     if inf_df.at[near_dinx, 'ADDED_ADDL']:
#         inf_df.at[near_dinx, 'TIME'] = row['TIME'] + 0.5
#         inf_df.at[near_dinx,'DATETIME'] = (datetime.strptime(row['DATETIME'], '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
#     else:
#         # 근처 Dosing 날짜와 같은 날짜이긴 할때
#         inf_df.at[inx, 'TIME'] = inf_df.at[near_dinx, 'TIME'] - 0.5
#         inf_df.at[inx, 'DATETIME'] = (datetime.strptime(inf_df.at[near_dinx, 'DATETIME'], '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
#
# # inf_df['UID'].drop_duplicates()
#
# for inx, row in ada_conc_first_rows.iterrows():
#     near_dinx = inx-1
#
#     # if (inf_df.at[inx, 'ID']==inf_df.at[near_dinx,'ID']):
#         # 근처 Dosing된 포인트가 SC 추가투여로 구성한 데이터일때
#     if ada_df.at[near_dinx, 'ADDED_ADDL']:
#         ada_df.at[near_dinx, 'TIME'] = row['TIME'] + 0.5
#         ada_df.at[near_dinx,'DATETIME'] = (datetime.strptime(row['DATETIME'], '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
#     else:
#         # 근처 Dosing 날짜와 같은 날짜이긴 할때
#         ada_df.at[inx, 'TIME'] = ada_df.at[near_dinx, 'TIME'] - 0.5
#         ada_df.at[inx, 'DATETIME'] = (datetime.strptime(ada_df.at[near_dinx, 'DATETIME'], '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
#
# # inf_df['UID'].drop_duplicates()
#
# inf_df = inf_df.sort_values(['ID','TIME','MDV'],ignore_index=True)
# ada_df = ada_df.sort_values(['ID','TIME','MDV'],ignore_index=True)
#
# inf_df['TIME(DAY)'] = inf_df['TIME']/24
# ada_df['TIME(DAY)'] = ada_df['TIME']/24
#
# inf_df['TIME(WEEK)'] = inf_df['TIME']/(7*24)
# ada_df['TIME(WEEK)'] = ada_df['TIME']/(7*24)
#
# ## 저장
# inf_df.to_csv(f'{output_dir}/infliximab_integrated_datacheck.csv',index=False, encoding='utf-8-sig')
# ada_df.to_csv(f'{output_dir}/adalimumab_integrated_datacheck.csv',index=False, encoding='utf-8-sig')

# ind_modeling_cols = ['ID','TIME','DV','MDV','AMT','DUR','CMT','IBD_TYPE']
# inf_df['IBD_TYPE'] = inf_df['IBD_TYPE'].map({'CD':1,'UC':2})
# ada_df['IBD_TYPE'] = ada_df['IBD_TYPE'].map({'CD':1,'UC':2})
# inf_df = inf_df[ind_modeling_cols].copy()
# ada_df = ada_df[ind_modeling_cols].copy()
#
#
# inf_df.to_csv(f'{output_dir}/infliximab_integrated_df.csv',index=False, encoding='utf-8-sig')
# ada_df.to_csv(f'{output_dir}/adalimumab_integrated_df.csv',index=False, encoding='utf-8-sig')

################################
# len(inf_df['UID'].drop_duplicates())
# ada_maint_df['AMT'] = ada_maint_df['AMT'].map({'1 pen':40,'2 pen':80, '2 pen':160})