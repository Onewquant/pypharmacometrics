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

mediator_cols = ['ID','NAME','DRUG','ROUTE','DATETIME','DV','MDV','AMT','DUR']
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

ind_cons_df = comp_df[comp_df['MIN_DOSE_DATE']==comp_df['IND_START_DATE']].reset_index(drop=True)
ind_diff_df = comp_df[comp_df['MIN_DOSE_DATE']!=comp_df['IND_START_DATE']].reset_index(drop=True)

# print(f"# Induction 시작시점 일치: {len(ind_cons_df)} (Infliximab: {len(ind_cons_df[ind_cons_df['DRUG']=='infliximab'])} / Adalimumab: {len(ind_cons_df[ind_cons_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(ind_cons_df[ind_cons_df['DRUG']=='ustekinumab'])}) ")
# print(f"# Induction 시작시점 불일치: {len(ind_diff_df)} (Infliximab: {len(ind_diff_df[ind_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(ind_diff_df[ind_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(ind_diff_df[ind_diff_df['DRUG']=='ustekinumab'])}) ")


# merged_df[merged_df['DRUG']=='infliximab']


merged_df.drop_duplicates(['ID'])


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
appended_frag_cols = ['UID', 'NAME', 'DRUG','ROUTE', 'TIME', 'WKTIME', 'DWKTIME', 'DV', 'MDV', 'AMT', 'DUR', 'CMT', 'DATETIME','AZERO','IBD_TYPE']


ind_df = list()
no_indconc_df = list()
for inx, row in ind_cons_df.iterrows():
    # if inx
    # break
    start_date = row['IND_START_DATE']
    end_date = (datetime.strptime(start_date,'%Y-%m-%d') + timedelta(days=127)).strftime('%Y-%m-%d')
    id_df = merged_df[merged_df['ID']==row['ID']].copy()
    if (len(id_df)==0):
        continue

    # if row['ID']==10875838:
    #     raise ValueError
    ind_df_frag = id_df[(id_df['DATE'] >= start_date)&(id_df['DATE'] <= end_date)].copy()


    if (len(ind_df_frag['MDV'].unique())<=1):
        no_indconc_df.append(pd.DataFrame([row]))
        continue

    ind_df_frag = ind_df_frag.sort_values(['ID', 'DATE', 'MDV'], ascending=[True, True, False])
    unique_date_df = ind_df_frag.groupby(['DATE'])['DV'].count().reset_index(drop=False)
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

            dupdv_df_frag = dup_df_frag[dup_df_frag['MDV']=='.'].copy()
            dupmdv_df_frag = dup_df_frag[dup_df_frag['MDV'] != '.'].copy()
            dupdv_df_frag['DATETIME'] = (datetime.strptime(dupmdv_df_frag['DATETIME'].iloc[0],'%Y-%m-%dT%H:%M')-timedelta(minutes=1)).strftime('%Y-%m-%dT%H:%M')

            ind_df_frag_list.append(pd.concat([nodup_df_frag, dupdv_df_frag, dupmdv_df_frag]))

        ind_df_frag = pd.concat(ind_df_frag_list).drop_duplicates(['ID','DATETIME']).sort_values(['ID','DATETIME'], ignore_index=True)
    zero_conc_row = ind_df_frag.iloc[:1,:].copy()
    zero_conc_row['DV'] = 0.0
    zero_conc_row['MDV'] = '.'
    zero_conc_row['AMT'] = '.'
    zero_conc_row['DUR'] = '.'
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

    ind_df_frag['WKTIME'] = ind_df_frag['TIME'] / 7
    ind_df_frag['DWKTIME'] = ind_df_frag['WKTIME'].diff().fillna(0.0)
    ind_df_frag['CMT'] = 1
    ind_df_frag = ind_df_frag.rename(columns={'ID':'UID'})
    # ind_df_frag['TIME'] = ind_df_frag['TIME'].map(lambda x:x)
    ind_df.append(ind_df_frag[appended_frag_cols])
    # ind_df_frag.columns
    # raise ValueError

    # IND 중간에 그만 둔 사람?

ind_df = pd.concat(ind_df).reset_index(drop=True)
ind_uniq_df = ind_df.drop_duplicates(['UID'])
no_indconc_df = pd.concat(no_indconc_df).reset_index(drop=True)

print(f"# 총 환자 수: 140")
print(f"# Induction 시작시점 일치: {len(ind_cons_df)} (Infliximab: {len(ind_cons_df[ind_cons_df['DRUG']=='infliximab'])} / Adalimumab: {len(ind_cons_df[ind_cons_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(ind_cons_df[ind_cons_df['DRUG']=='ustekinumab'])}) ")
print(f"# Induction 시작시점 불일치: {len(ind_diff_df)} (Infliximab: {len(ind_diff_df[ind_diff_df['DRUG']=='infliximab'])} / Adalimumab: {len(ind_diff_df[ind_diff_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(ind_diff_df[ind_diff_df['DRUG']=='ustekinumab'])}) ")
print('')
print(f"# Induction Trough 농도 측정값 존재(시작시점 일치중): {len(ind_uniq_df)} (Infliximab: {len(ind_uniq_df[ind_uniq_df['DRUG']=='infliximab'])} / Adalimumab: {len(ind_uniq_df[ind_uniq_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(ind_uniq_df[ind_uniq_df['DRUG']=='ustekinumab'])}) ")
print(f"# Induction Trough 농도 측정값 부재(시작시점 일치중): {len(no_indconc_df)} (Infliximab: {len(no_indconc_df[no_indconc_df['DRUG']=='infliximab'])} / Adalimumab: {len(no_indconc_df[no_indconc_df['DRUG']=='adalimumab'])} / Ustekinumab: {len(no_indconc_df[no_indconc_df['DRUG']=='ustekinumab'])}) ")





# inf_ind_df = inf_ind_df
# inf_ind_df


inf_ind_df = ind_df[ind_df['DRUG']=='infliximab'].copy()
ada_ind_df = ind_df[ind_df['DRUG']=='adalimumab'].copy()

inf_ind_df['ID'] = inf_ind_df['UID'].map({uid:uid_inx for uid_inx, uid in enumerate(list(inf_ind_df['UID'].unique()))})
ada_ind_df['ID'] = ada_ind_df['UID'].map({uid:uid_inx for uid_inx, uid in enumerate(list(ada_ind_df['UID'].unique()))})

ind_modeling_cols = ['ID','TIME','WKTIME','DWKTIME','DV','MDV','AMT','DUR','CMT','DATETIME','IBD_TYPE','AZERO','UID','NAME','ROUTE','DRUG']
inf_ind_df = inf_ind_df[ind_modeling_cols].copy()
ada_ind_df = ada_ind_df[ind_modeling_cols].copy()

inf_ind_df.to_csv(f'{output_dir}/infliximab_induction_datacheck.csv',index=False, encoding='utf-8-sig')
ada_ind_df.to_csv(f'{output_dir}/adalimumab_induction_datacheck.csv',index=False, encoding='utf-8-sig')

# Covariates

ind_modeling_cols = ['ID','TIME','DV','MDV','AMT','DUR','CMT','IBD_TYPE']
inf_ind_df['IBD_TYPE'] = inf_ind_df['IBD_TYPE'].map({'CD':1,'UC':2})
ada_ind_df['IBD_TYPE'] = ada_ind_df['IBD_TYPE'].map({'CD':1,'UC':2})





inf_ind_df = inf_ind_df[ind_modeling_cols].copy()
ada_ind_df = ada_ind_df[ind_modeling_cols].copy()
# ada_ind_df['AMT'] = ada_ind_df['AMT'].map({'1 pen':40,'2 pen':80, '2 pen':160})

inf_ind_df.to_csv(f'{output_dir}/infliximab_induction_df.csv',index=False, encoding='utf-8-sig')
ada_ind_df.to_csv(f'{output_dir}/adalimumab_induction_df.csv',index=False, encoding='utf-8-sig')

# len(inf_ind_df['ID'].drop_duplicates())