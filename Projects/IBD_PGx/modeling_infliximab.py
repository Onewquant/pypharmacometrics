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
# lab_df['DATETIME'] = lab_df['DATETIME'].map(lambda x:x.replace(' ','T'))
# lab_df = lab_df.sort_values(['ID','DATETIME'])
# lab_df.to_csv(f"{output_dir}/conc_df.csv", encoding='utf-8-sig', index=False)

# lab_df[(lab_df['CONC']!=lab_df['CONC_cumlab'])|(lab_df['DRUG']!=lab_df['DRUG_cumlab'])][['DRUG','CONC','DRUG_cumlab','CONC_cumlab']]
lab_df = pd.read_csv(f"{output_dir}/conc_df.csv")
# lab_df.columns
# lab_df['DRUG']
lab_df = lab_df.rename(columns={'CONC':'DV'})
lab_df['MDV']='.'
lab_df['DUR']='.'
lab_df['AMT']='.'

dose_df = pd.read_csv(f"{output_dir}/dose_df.csv")
# dose_df.columns
# dose_df['DOSE']
dose_df = dose_df.rename(columns={'DOSE':'AMT'})
dose_df['DV'] = '.'
dose_df['MDV'] = 1
dose_df['DUR'] = 0.5

mediator_cols = ['ID','NAME','DATETIME','DV','MDV','AMT','DUR']
lab_df = lab_df[lab_df['DRUG']=='infliximab'][mediator_cols].reset_index(drop=True)
dose_df = dose_df[dose_df['DRUG']=='infliximab'][mediator_cols].reset_index(drop=True)
merged_df = pd.concat([lab_df,dose_df]).sort_values(['ID','DATETIME'], ignore_index=True)
merged_df.to_csv(f'{output_dir}/merged_df.csv',index=False, encoding='utf-8-sig')

merged_df['DATE'] = merged_df['DATETIME'].map(lambda x:x.split('T')[0])
merged_df = merged_df.merge(induction_df[['ID','IBD_TYPE']], on=['ID'], how='left')
# lab_df['DATETIME']
# dose_df['DATETIME']
# dose_df['DRUG']

# Induction Phase 만 구분

min_dose_df = dose_df.groupby(['ID','NAME']).agg({'DATETIME':'min'}).reset_index(drop=False)
min_dose_df['MIN_DOSE_DATE'] = min_dose_df['DATETIME'].map(lambda x:x.split('T')[0])
comp_df = min_dose_df.merge(induction_df[['ID','IND_START_DATE']], on=['ID'], how='left')

ind_cons_df = comp_df[comp_df['MIN_DOSE_DATE']==comp_df['IND_START_DATE']].reset_index(drop=True)
ind_diff_df = comp_df[comp_df['MIN_DOSE_DATE']!=comp_df['IND_START_DATE']].reset_index(drop=True)

# min_dose_df['ID']
# comp_df['IND_START_DATE'].iloc[0]

inf_ind_df = list()
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

        ind_df_frag = pd.concat(ind_df_frag_list).drop_duplicates(['ID','DATETIME'])

    ind_date_str = ind_df_frag['DATETIME'].iloc[0]
    ind_df_frag['ATIME'] = ind_df_frag['DATETIME'].map(lambda x:(datetime.strptime(x,'%Y-%m-%dT%H:%M') - datetime.strptime(ind_date_str,'%Y-%m-%dT%H:%M'))).dt.total_seconds() / 86400
    ind_df_frag = ind_df_frag.rename(columns={'ID':'UID'})
    # ind_df_frag['ATIME'] = ind_df_frag['ATIME'].map(lambda x:x)
    inf_ind_df.append(ind_df_frag[['UID','NAME','ATIME','DV','MDV','AMT','DUR','DATETIME','IBD_TYPE']])
    # ind_df_frag.columns
    # raise ValueError

    # IND 중간에 그만 둔 사람?

inf_ind_df = pd.concat(inf_ind_df).reset_index(drop=True)



inf_ind_df['ID'] = inf_ind_df['UID'].map({uid:uid_inx for uid_inx, uid in enumerate(list(inf_ind_df['UID'].unique()))})
# inf_ind_df = inf_ind_df
# inf_ind_df


inf_ind_df[['ID','ATIME','DV','MDV','AMT','DUR','DATETIME','IBD_TYPE','UID','NAME']].to_csv(f'{output_dir}/modeling_df_ind_checked.csv',index=False, encoding='utf-8-sig')

# len(inf_ind_df['ID'].drop_duplicates())