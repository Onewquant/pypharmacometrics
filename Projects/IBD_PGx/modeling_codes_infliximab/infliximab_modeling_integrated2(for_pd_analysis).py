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
# lab_df.to_csv(f"{output_dir}/conc_df.csv", encoding='utf-8-sig', index=False)

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
for drug in ['infliximab', 'adalimumab', 'ustekinumab']: #break
    # raise ValueError
    # '36959194'
    # indmaint_df[indmaint_df['ID']==36959194]
    indmaint_df = dose_df[dose_df['DRUG']==drug].copy()
    init_dosing_dt = indmaint_df.groupby('ID',as_index=False)['DATETIME'].min().rename(columns={"DATETIME":"INIT_DOSE_DT"})
    indmaint_df = indmaint_df.merge(init_dosing_dt, on=['ID'], how='left')
    # indmaint_df.columns
    indmaint_df['DOSING_TIME'] = indmaint_df.apply(lambda x:(datetime.strptime(x['DATETIME'],'%Y-%m-%dT%H:%M') - datetime.strptime(x['INIT_DOSE_DT'],'%Y-%m-%dT%H:%M')), axis=1).dt.total_seconds() / 86400 / 7

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
        indmaint_df['START_INDMAINT'] = (~((indmaint_df['DOSING_INTERVAL'] < 4)&(indmaint_df['ROUTE']=='IV')))*1
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
# sns.displot(indmaint_df['DOSING_INTERVAL'])
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

# merged_df.to_csv(f'{output_dir}/merged_df.csv',index=False, encoding='utf-8-sig')
# merged_df.drop_duplicates(['ID'])

merged_df['DATE'] = merged_df['DATETIME'].map(lambda x:x.split('T')[0])
# merged_df['AZERO'] = 0
merged_df = merged_df.merge(ibd_type_df[['ID', 'IBD_TYPE']], on=['ID'], how='left')

## ULOQ 넘는 농도값은 제거
merged_df = merged_df[~((merged_df['DV']==48)&((merged_df['ID'].isin([21911051,])) & (merged_df['DATE'].isin(['2019-02-22',]))))].copy()


drug = 'infliximab'
maint_cons_df = total_indmaint_df[(total_indmaint_df['DRUG']==drug)&(total_indmaint_df['START_INDMAINT']==0)].reset_index(drop=True)
maint_diff_df = total_indmaint_df[(total_indmaint_df['DRUG']==drug)&(total_indmaint_df['START_INDMAINT']==1)].reset_index(drop=True)


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
for inx, row in maint_diff_df.iterrows(): #break


    id_df = merged_df[merged_df['ID']==row['ID']].copy()
    if (len(id_df)==0):
        no_maintid_df.append(pd.DataFrame([row]))
        continue

    # first_dv_inx = id_df[id_df['MDV']==0].iloc[0].name
    # maint_df_frag = id_df.loc[first_dv_inx:].copy()
    maint_df_frag = id_df.copy()

    # drug_maint_conc_df = maint_df_frag.groupby(['DRUG','MDV'],as_index=False)[['ID']].count().rename(columns={'ID':'COUNT'})
    # drug_maint_conc_df = drug_maint_conc_df[(drug_maint_conc_df['MDV']!=1)&(drug_maint_conc_df['COUNT']>1)].copy()
    #
    # if len(drug_maint_conc_df)==0:
    #     no_maintconc_df.append(pd.DataFrame([row]))
    #     continue
    # else:
    #     maint_df_frag = maint_df_frag[maint_df_frag['DRUG'].isin(list(drug_maint_conc_df['DRUG']))].copy()
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


    maint_date_str = maint_df_frag['DATETIME'].iloc[0]
    maint_df_frag['TIME'] = maint_df_frag['DATETIME'].map(lambda x:(datetime.strptime(x,'%Y-%m-%dT%H:%M') - datetime.strptime(maint_date_str,'%Y-%m-%dT%H:%M'))).dt.total_seconds() * 24 / 86400

    # maint_df_frag['WKTIME'] = maint_df_frag['TIME'] / 7
    # maint_df_frag['DWKTIME'] = maint_df_frag['WKTIME'].diff().fillna(0.0)
    maint_df_frag['CMT'] = maint_df_frag['ROUTE'].map(lambda x: 2 if x!='SC' else 1)
    maint_df_frag['START_INDMAINT'] = row['START_INDMAINT']
    maint_df_frag = maint_df_frag.rename(columns={'ID':'UID'})
    # ind_df_frag['TIME'] = ind_df_frag['TIME'].map(lambda x:x)
    maint_df.append(maint_df_frag[appended_frag_cols])


maint_df = pd.concat(maint_df).reset_index(drop=True)
maint_uniq_df = maint_df.drop_duplicates(['UID'])
if len(no_maintconc_df)==0:
    no_maintconc_df = pd.DataFrame(columns=list(maint_df.columns))
else:
    no_maintconc_df = pd.concat(no_maintconc_df).reset_index(drop=True)



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
    drug_df['TIME(DAY)'] = drug_df['TIME'] / 24
    drug_df['TIME(WEEK)'] = drug_df['TIME'] / (7 * 24)

    # drug_df['DT_YEAR'] = drug_df['DATETIME'].map(lambda x: int(x.split('-')[0]))
    # drug_df['DT_MONTH'] = drug_df['DATETIME'].map(lambda x: int(x.split('-')[1]))
    # drug_df['DT_DAY'] = drug_df['DATETIME'].map(lambda x: int(x.split('-')[2].split('T')[0]))
    #
    # drug_df['TIME'] =drug_df['TIME'] / 24
    # drug_df['DUR'] = drug_df['DUR'].map(lambda x: float(x) / 24 if x != '.' else x)
    # # drug_df['']
    # drug_df.drop(['DATETIME','NAME','ROUTE','DRUG','TIME(DAY)','TIME(WEEK)','START_INDMAINT', 'ADDED_ADDL', 'IBD_TYPE','UID'],axis=1).to_csv(f'{output_dir}/modeling_df_covar/{drug}_integrated_presim_df_dayscale(for pda).csv', index=False, encoding='utf-8')
    # drug_df.columns

    ## 저장
    if not os.path.exists(f'{output_dir}/modeling_df_datacheck(for pda)'):
        os.mkdir(f'{output_dir}/modeling_df_datacheck(for pda)')
    drug_df.to_csv(f'{output_dir}/modeling_df_datacheck(for pda)/{drug}_integrated_datacheck(for pda).csv', index=False, encoding='utf-8-sig')

    print(f"[{drug}] modeling df / {len(drug_df)} rows")

print(f"\n# 총 환자 수: 140")
print(f"# 전체 데이터: {len(maint_cons_df)} (Infliximab: {len(maint_cons_df[maint_cons_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_cons_df[maint_cons_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_cons_df[maint_cons_df['DRUG']=='ustekinumab'])}) ")
# print(f"# Induction 시작시점 불일치: {len(maint_diff_df)} (Infliximab: {len(maint_diff_df[maint_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_diff_df[maint_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_diff_df[maint_diff_df['DRUG']=='ustekinumab'])}) ")
