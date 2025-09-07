## 시간에 따른 모든 사람의 PD 경향성 그려보기

from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

# drug = 'infliximab'
# mode_str = 'integrated'
# df = pd.read_csv(f'{output_dir}/{drug}_{mode_str}_pdeda_df_dayscale.csv')
df = pd.read_csv(f'{output_dir}/modeling_df_datacheck(for pda)/infliximab_integrated_datacheck(for pda).csv')
df['TIME'] = df['TIME(DAY)'].copy()
# df.columns
# df.drop_duplicates(['UID'])['START_INDMAINT'].sum()
uniq_df = df.drop_duplicates(['UID'])[['UID','START_INDMAINT']].copy()
ind_uids = list(uniq_df[(uniq_df['START_INDMAINT']==0)]['UID'].reset_index(drop=True))
maint_uids = list(uniq_df[(uniq_df['START_INDMAINT']==1)]['UID'].reset_index(drop=True))

pd_df = pd.read_csv(f'{output_dir}/pd_bsize_df_ori.csv')
# pd_df['UID'].drop_duplicates()
df['UID'].drop_duplicates()
pd_cols = [c for c in pd_df.columns if 'PD_' in c]
# pd_df.columns

no_maint_blpd_uids = list()
no_maint_tgpd_uids = list()

res_df = list()
ind_count = 0
maint_count = 0
for uid, uid_df in df.groupby('UID'): #break

    ibd_type = uid_df['IBD_TYPE'].iloc[0]
    uid_pd_df = pd_df[pd_df['UID'] == uid].copy()

    min_date = uid_df['DATETIME'].min().split('T')[0]
    min_date_bf2mo = (datetime.strptime(min_date, '%Y-%m-%d') - timedelta(days=60)).strftime('%Y-%m-%d')

    res_dict = {'UID': uid,
                'IBD_TYPE': ibd_type}

    ## BL 기록하기

    bl_pd_df = uid_pd_df[(uid_pd_df['DATETIME'] <= min_date) & (uid_pd_df['DATETIME'] >= min_date_bf2mo)].copy()
    if len(bl_pd_df) == 0:
        no_maint_blpd_uids.append(uid)
        blpd_date = np.nan
        bl_pd_row = {pdc: np.nan for pdc in pd_cols}

    else:
        # BL_DATE -> PD값들의 NAN이 가장 적은 것 중 가장 나중 날짜
        bl_pd_df['NAN_COUNT'] = (bl_pd_df[pd_cols].isna() * 1).sum(axis=1)
        blpd_date = bl_pd_df[bl_pd_df['NAN_COUNT'] == bl_pd_df['NAN_COUNT'].min()]['DATETIME'].iloc[-1]

        bl_pd_df = bl_pd_df.fillna(method='ffill')
        bl_pd_row = bl_pd_df.iloc[-1]

    ## BL 기록시행

    res_dict['BL_DATE'] = blpd_date
    for pdc in pd_cols:
        res_dict['BL_'+pdc] = bl_pd_row[pdc]

    ## IND / MAINT Phase 구분하여 기록

    if uid in ind_uids:
        ind_count+=1
        # raise ValueError
        print(f'[IND ({ind_count})] / {uid}')
        res_dict['PHASE']='Induction'

        try: dose3rd_date = uid_df[uid_df['MDV']==1].iloc[2]['DATETIME'].split('T')[0]
        except: dose3rd_date = (datetime.strptime(min_date, '%Y-%m-%d') + timedelta(days=90)).strftime('%Y-%m-%d')

        try: dose4th_date = uid_df[uid_df['MDV']==1].iloc[3]['DATETIME'].split('T')[0]
        except: dose4th_date = (datetime.strptime(min_date, '%Y-%m-%d') + timedelta(days=127)).strftime('%Y-%m-%d')

        ## AFTIND 기록하기

        tg_pd_df = uid_pd_df[(uid_pd_df['DATETIME'] <= dose4th_date) & (uid_pd_df['DATETIME'] >= dose3rd_date)].copy()

        if len(tg_pd_df) == 0:
            no_maint_blpd_uids.append(uid)
            tgpd_date = np.nan
            tg_pd_row = {pdc: np.nan for pdc in pd_cols}
        else:
            tg_pd_df['NAN_COUNT'] = (tg_pd_df[pd_cols].isna() * 1).sum(axis=1)
            tgpd_date = tg_pd_df[tg_pd_df['NAN_COUNT'] == tg_pd_df['NAN_COUNT'].min()]['DATETIME'].iloc[0]

            tg_pd_df = tg_pd_df.fillna(method='bfill')
            tg_pd_row = tg_pd_df.iloc[0]

        ## AFTIND 기록시행
        res_dict['TG_DATE'] = tgpd_date
        for pdc in pd_cols:
            res_dict['TG_'+pdc] = tg_pd_row[pdc]

    else:
        maint_count+=1
        print(f'[MAINT ({maint_count})] {uid}')
        res_dict['PHASE']='Maintenance'
        min_date_af11mo = (datetime.strptime(min_date, '%Y-%m-%d') + timedelta(days=335)).strftime('%Y-%m-%d')
        min_date_af13mo = (datetime.strptime(min_date, '%Y-%m-%d') + timedelta(days=395)).strftime('%Y-%m-%d')

        ## 1YR 기록하기

        tg_pd_df = uid_pd_df[(uid_pd_df['DATETIME'] <= min_date_af13mo) & (uid_pd_df['DATETIME'] >= min_date_af11mo)].copy()

        if len(tg_pd_df) == 0:
            no_maint_blpd_uids.append(uid)
            tgpd_date = np.nan
            tg_pd_row = {pdc: np.nan for pdc in pd_cols}
        else:
            tg_pd_df['NAN_COUNT'] = (tg_pd_df[pd_cols].isna() * 1).sum(axis=1)
            tgpd_date = tg_pd_df[tg_pd_df['NAN_COUNT'] == tg_pd_df['NAN_COUNT'].min()]['DATETIME'].iloc[0]

            tg_pd_df = tg_pd_df.fillna(method='bfill')
            tg_pd_row = tg_pd_df.iloc[0]

        ## 1YR 기록시행
        res_dict['TG_DATE'] = tgpd_date
        for pdc in pd_cols:
            res_dict['TG_'+pdc] = tg_pd_row[pdc]

    res_df.append(res_dict)
res_df = pd.DataFrame(res_df)

ind_res_df = res_df[res_df['PHASE']=='Induction'].copy()
ind_res_df[~(ind_res_df['BL_DATE'].isna())]
ind_res_df[~(ind_res_df['TG_DATE'].isna())]
ind_res_df[(~(ind_res_df['BL_DATE'].isna()))&(~(ind_res_df['TG_DATE'].isna()))]