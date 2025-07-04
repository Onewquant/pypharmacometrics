from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

induction_df = pd.read_excel(f"{resource_dir}/IBD-PGx_induction_date.xlsx")
# induction_df.sort_values(['EMR ID']).to_csv(f"{resource_dir}/IBD-PGx_induction_date.csv",encoding='utf-8-sig',index=False)
induction_df = induction_df.rename(columns={'EMR ID':'ID','name':'NAME','induction_start_date':'IND_START_DATE','IBD type':'IBD_TYPE'})
induction_df['IND_START_DATE'] = induction_df['IND_START_DATE'].astype(str).replace('NaT','')

# totlab_df = pd.read_csv(f"{output_dir}/conc_df(lab).csv")
# cumlab_df = pd.read_csv(f"{output_dir}/conc_df(cum_lab).csv")
# lab_df = pd.concat([cumlab_df, totlab_df[['오더일','보고일']]], axis=1)
# lab_df['접수일'] = lab_df['DATETIME'].map(lambda x:x.split(' ')[0])
# lab_df['접수시간'] = lab_df['DATETIME'].map(lambda x:x.split(' ')[1])
# receipt_list = list()
# for inx, row in lab_df.iterrows():
#     if (datetime.strptime(row['보고일'],'%Y-%m-%d')-datetime.strptime(row['오더일'],'%Y-%m-%d')).days > 14:
#         date_str = min(row['보고일'],row['접수일'])
#         receipt_list.append(date_str+f"T{row['접수시간']}")
#     else:
#         if np.abs((datetime.strptime(row['오더일'],'%Y-%m-%d')-datetime.strptime(row['접수일'],'%Y-%m-%d')).days) <= 14:
#             date_str = min(row['오더일'], min(row['보고일'], row['접수일']))
#         else:
#             date_str = max(row['보고일'],row['접수일'])
#         receipt_list.append(date_str+f"T{row['접수시간']}")
# lab_df['DATETIME'] = receipt_list
# # lab_df['DATETIME'] = lab_df['DATETIME'].map(lambda x:x.re place(' ','T'))
# lab_df = lab_df.sort_values(['ID','DATETIME'])
# lab_df.to_csv(f"{output_dir}/conc_df.csv", encoding='utf-8-sig', index=False)

# lab_df[(lab_df['CONC']!=lab_df['CONC_cumlab'])|(lab_df['DRUG']!=lab_df['DRUG_cumlab'])][['DRUG','CONC','DRUG_cumlab','CONC_cumlab']]
lab_df = pd.read_csv(f"{output_dir}/conc_df.csv")
# lab_df[lab_df['ID']==11788526]
# lab_df.columns
# lab_df['DRUG']
lab_df = lab_df.rename(columns={'CONC':'DV'})
lab_df['MDV']='.'
lab_df['RATE']='.'
lab_df['DUR']='.'
lab_df['AMT']='.'
lab_df['ROUTE']='.'

dose_df = pd.read_csv(f"{output_dir}/dose_df.csv")
# dose_df.columns
# dose_df['DOSE']
dose_df = dose_df.rename(columns={'DOSE':'AMT'})
dose_df['DV'] = '.'
dose_df['MDV'] = 1
dose_df['RATE']=dose_df['AMT']
dose_df['DUR'] = '.'

mediator_cols = ['ID','NAME','DRUG','ROUTE','DATETIME','DV','MDV','AMT','RATE','DUR']
lab_df = lab_df[mediator_cols].reset_index(drop=True)
dose_df = dose_df[mediator_cols].reset_index(drop=True)
merged_df = pd.concat([lab_df,dose_df]).sort_values(['ID','DATETIME'], ignore_index=True)

merged_df.to_csv(f'{output_dir}/merged_df.csv',index=False, encoding='utf-8-sig')
# merged_df.drop_duplicates(['ID'])

merged_df['DATE'] = merged_df['DATETIME'].map(lambda x:x.split('T')[0])
merged_df['AZERO'] = 0
merged_df = merged_df.merge(induction_df[['ID','IBD_TYPE']], on=['ID'], how='left')

# Induction Phase 불일치 환자 구분 (전체 합친 데이터에서)

min_dose_df = merged_df.groupby(['ID']).agg({'NAME':'min','DATETIME':'min','DRUG':'first'}).reset_index(drop=False)
min_dose_df['MIN_DOSE_DATE'] = min_dose_df['DATETIME'].map(lambda x:x.split('T')[0])
comp_df = min_dose_df.merge(induction_df[['ID','IND_START_DATE']], on=['ID'], how='left')

maint_cons_df = comp_df.reset_index(drop=True)

maint_cons_df  = comp_df[comp_df['MIN_DOSE_DATE']==comp_df['IND_START_DATE']].reset_index(drop=True)
maint_diff_df  = comp_df[comp_df['MIN_DOSE_DATE']!=comp_df['IND_START_DATE']].reset_index(drop=True)

# print(f"# Induction 시작시점 일치: {len(ind_cons_df)} (Infliximab: {len(ind_cons_df[ind_cons_df['DRUG']=='infliximab'])} / Adalimumab: {len(ind_cons_df[ind_cons_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(ind_cons_df[ind_cons_df['DRUG']=='ustekinumab'])}) ")
# print(f"# Induction 시작시점 불일치: {len(ind_diff_df)} (Infliximab: {len(ind_diff_df[ind_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(ind_diff_df[ind_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(ind_diff_df[ind_diff_df['DRUG']=='ustekinumab'])}) ")


# merged_df[merged_df['DRUG']=='infliximab']


# merged_df.drop_duplicates(['ID'])


# lab_df['DATETIME']
# dose_df['DATETIME']
# dose_df['DRUG']

# # Induction Phase 만 구분
#
# min_dose_df = dose_df.groupby(['ID','NAME']).agg({'DATETIME':'min'}).reset_index(drop=False)
# min_dose_df['MIN_DOSE_DATE'] = min_dose_df['DATETIME'].map(lambda x:x.split('T')[0])
# comp_df = min_dose_df.merge(induction_df[['ID','IND_START_DATE']], on=['ID'], how='left')
#
# ind_cons_df = comp_df[comp_df['MIN_DOSE_DATE']==comp_df['IND_START_DATE']].reset_index(drop=True)
# ind_diff_df = comp_df[comp_df['MIN_DOSE_DATE']!=comp_df['IND_START_DATE']].reset_index(drop=True)

# min_dose_df['ID']
# comp_df['IND_START_DATE'].iloc[0]

appended_frag_cols = ['UID', 'NAME', 'DRUG','ROUTE', 'TIME', 'WKTIME', 'DWKTIME', 'DV', 'MDV', 'AMT','RATE', 'DUR', 'CMT', 'DATETIME','AZERO','IBD_TYPE']

maint_df = list()
no_maintconc_df = list()
for inx, row in maint_cons_df.iterrows():
    # if inx
    # break
    start_date = row['IND_START_DATE']
    end_date = (datetime.strptime(start_date,'%Y-%m-%d') + timedelta(days=127)).strftime('%Y-%m-%d')
    id_df = merged_df[merged_df['ID']==row['ID']].copy()
    if (len(id_df)==0):
        continue
    # if (len((id_df['MDV']!=1).sum())<=2):
    #     continue

    # if row['ID']==12681022:
    #     raise ValueError



    ind_df_frag = id_df[(id_df['DATE'] >= start_date) & (id_df['DATE'] <= end_date)].copy()
    ind_conc_df_frag = ind_df_frag[ind_df_frag['MDV']=='.']
    if len(ind_conc_df_frag)==0: # induction phase에 농도 데이터 부재시 나머지 데이터 중 첫 농도 값을 기준으로 maintenance 시작
        ind_end_inx = ind_df_frag.iloc[-1].name
        maint_df_frag = id_df.loc[ind_end_inx:].copy()
        first_dv_inx = maint_df_frag[maint_df_frag['MDV']=='.'].iloc[0].name
        maint_df_frag = id_df.loc[first_dv_inx:].copy()
    else: # induction phase에 농도 데이터 있으면, induction phase의 마지막 농도 데이터 값을 기준으로 maintenance phase 설정
        last_dv_inx = ind_df_frag[ind_df_frag['MDV']=='.'].iloc[-1].name
        maint_df_frag = id_df.loc[last_dv_inx:].copy()
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

            dupdv_df_frag = dup_df_frag[dup_df_frag['MDV']=='.'].copy()
            dupmdv_df_frag = dup_df_frag[dup_df_frag['MDV'] != '.'].copy()
            if len(dupdv_df_frag)>=3:
                print("> 3")
                raise ValueError
            elif (len(dupdv_df_frag)==2) and (len(dupmdv_df_frag)==1):
                # if len(dupdv_df_frag) == 2: raise ValueError
                dupdv_df_frag.iat[0,3] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0],'%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
                dupdv_df_frag.iat[-1, 3] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0], '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
            elif (len(dupdv_df_frag)==1) and (len(dupmdv_df_frag)==2):
                dupdv_df_frag.iat[0,3] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0],'%Y-%m-%dT%H:%M')-timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
            elif (len(dupdv_df_frag)==1) and (len(dupmdv_df_frag)==1):
                dupdv_df_frag.iat[0,3] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0],'%Y-%m-%dT%H:%M')-timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
            elif (len(dupdv_df_frag)==0) or (len(dupmdv_df_frag)==0):
                # print("No data")
                # raise ValueError
                pass
            else:
                raise ValueError

            maint_df_frag_list.append(pd.concat([nodup_df_frag, dupdv_df_frag, dupmdv_df_frag]))

        maint_df_frag = pd.concat(maint_df_frag_list).drop_duplicates(['ID','DATETIME']).sort_values(['ID','DATETIME'], ignore_index=True)
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

    maint_df_frag['WKTIME'] = maint_df_frag['TIME'] / 7
    maint_df_frag['DWKTIME'] = maint_df_frag['WKTIME'].diff().fillna(0.0)
    maint_df_frag['CMT'] = maint_df_frag['ROUTE'].map(lambda x: 2 if x!='SC' else 1)
    maint_df_frag['RATE'] = maint_df_frag.apply(lambda x: x['RATE'] if (x['ROUTE']!='SC') else -2, axis=1)
    maint_df_frag['RATE'] = maint_df_frag.apply(lambda x: x['RATE'] if (x['ROUTE']!='SC') else -2, axis=1)

    # maint_df_frag['ROUTE']
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

print(f"# 총 환자 수: 140")
print(f"# Induction 시작시점 일치: {len(maint_cons_df)} (Infliximab: {len(maint_cons_df[maint_cons_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_cons_df[maint_cons_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_cons_df[maint_cons_df['DRUG']=='ustekinumab'])}) ")
print(f"# Induction 시작시점 불일치: {len(maint_diff_df)} (Infliximab: {len(maint_diff_df[maint_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_diff_df[maint_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_diff_df[maint_diff_df['DRUG']=='ustekinumab'])}) ")
print('')
print(f"# Maintenance 시작시 농도 측정값 존재: {len(maint_uniq_df)} (Infliximab: {len(maint_uniq_df[maint_uniq_df['DRUG']=='infliximab'])} / Adalimumab: {len(maint_uniq_df[maint_uniq_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(maint_uniq_df[maint_uniq_df['DRUG']=='ustekinumab'])}) ")
print(f"# Maintenance 시작시 농도 측정값 부재: {len(no_maintconc_df)} (Infliximab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='infliximab'])} / Adalimumab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(no_maintconc_df[no_maintconc_df['DRUG']=='ustekinumab'])}) ")



# inf_ind_df = inf_ind_df
# inf_ind_df


inf_maint_df = maint_df[maint_df['DRUG']=='infliximab'].sort_values(['UID','DATETIME'], ignore_index=True)
inf_maint_df['AZERO'] = (inf_maint_df['UID']!=(inf_maint_df['UID'].shift(1).fillna(0)))*1
ada_maint_df = maint_df[maint_df['DRUG']=='adalimumab'].sort_values(['UID','DATETIME'], ignore_index=True)
ada_maint_df['AZERO'] =(ada_maint_df['UID']!=(ada_maint_df['UID'].shift(1).fillna(0)))*1

inf_maint_df['ID'] = inf_maint_df['UID'].map({uid:uid_inx for uid_inx, uid in enumerate(list(inf_maint_df['UID'].unique()))})
ada_maint_df['ID'] = ada_maint_df['UID'].map({uid:uid_inx for uid_inx, uid in enumerate(list(ada_maint_df['UID'].unique()))})

ind_modeling_cols = ['ID','TIME','WKTIME','DWKTIME','DV','MDV','AMT','RATE','DUR','CMT','DATETIME','IBD_TYPE','AZERO','UID','NAME','ROUTE','DRUG']
inf_maint_df = inf_maint_df[ind_modeling_cols].copy()
ada_maint_df = ada_maint_df[ind_modeling_cols].copy()



inf_maint_df.to_csv(f'{output_dir}/infliximab_maintenance_datacheck.csv',index=False, encoding='utf-8-sig')
ada_maint_df.to_csv(f'{output_dir}/adalimumab_maintenance_datacheck.csv',index=False, encoding='utf-8-sig')

ind_modeling_cols = ['ID','TIME','DV','MDV','AMT','RATE','DUR','CMT','IBD_TYPE']
inf_maint_df['IBD_TYPE'] = inf_maint_df['IBD_TYPE'].map({'CD':1,'UC':2})
ada_maint_df['IBD_TYPE'] = ada_maint_df['IBD_TYPE'].map({'CD':1,'UC':2})
inf_maint_df = inf_maint_df[ind_modeling_cols].copy()
ada_maint_df = ada_maint_df[ind_modeling_cols].copy()
# ada_maint_df['AMT'] = ada_maint_df['AMT'].map({'1 pen':40,'2 pen':80, '2 pen':160})

inf_maint_df.to_csv(f'{output_dir}/infliximab_maintenance_df.csv',index=False, encoding='utf-8-sig')
ada_maint_df.to_csv(f'{output_dir}/adalimumab_maintenance_df.csv',index=False, encoding='utf-8-sig')

# len(inf_ind_df['ID'].drop_duplicates())