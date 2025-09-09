## 시간에 따른 모든 사람의 PD 경향성 그려보기

from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## 모델링 데이터셋 로딩

df = pd.read_csv(f'{output_dir}/modeling_df_datacheck(for pda)/infliximab_integrated_datacheck(for pda).csv')
id_uid_dict = {row['ID']:row['UID'] for inx, row in df[['ID','UID']].drop_duplicates(['ID']).iterrows()}
df['TIME'] = df['TIME(DAY)'].copy()
uniq_df = df.drop_duplicates(['UID'])[['UID','NAME','START_INDMAINT']].copy()
ind_uids = list(uniq_df[(uniq_df['START_INDMAINT']==0)]['UID'].reset_index(drop=True))
maint_uids = list(uniq_df[(uniq_df['START_INDMAINT']==1)]['UID'].reset_index(drop=True))

## Simulation 데이터셋 로딩

sim_df = pd.read_csv(f"{nonmem_dir}/run/sim86",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
sim_df['ID'] = sim_df['ID'].astype(int)
sim_df['UID'] = sim_df['ID'].map(id_uid_dict)
sim_df['DT_YEAR'] = sim_df['DT_YEAR'].map(lambda x:str(int(x)).zfill(4))
sim_df['DT_MONTH'] = sim_df['DT_MONTH'].map(lambda x:str(int(x)).zfill(2))
sim_df['DT_DAY'] = sim_df['DT_DAY'].map(lambda x:str(int(x)).zfill(2))
sim_df['DATETIME'] = sim_df['DT_YEAR']+'-'+sim_df['DT_MONTH']+'-'+sim_df['DT_DAY']
#
# sim_df['DT_MONTH']
# sim_df['DT_DAY']

## PD 데이터셋 로딩

pd_df = pd.read_csv(f'{output_dir}/pd_bsize_df_ori.csv')

## LAB 중 CRP, fCal 을 PD 데이터셋에 붙이기

totlab_df = pd.read_csv(f"{output_dir}/lab_df.csv")
totlab_df = totlab_df[['UID', 'DATETIME',  'CRP', 'hsCRP',  'Calprotectin (Stool)']].copy()
for c in list(totlab_df.columns)[2:]:  # break
    totlab_df[c] = totlab_df[c].map(lambda x: x if type(x) == float else float(re.findall(r'[\d]+.*[\d]*', str(x))[0]))
totlab_df['PD_FCAL'] = totlab_df[['Calprotectin (Stool)']].max(axis=1)
totlab_df['PD_CRP'] = totlab_df[['CRP', 'hsCRP']].max(axis=1)
totlab_df = totlab_df[['UID','DATETIME','PD_FCAL', 'PD_CRP']].copy()
pd_df = totlab_df.merge(pd_df, on=['UID','DATETIME'], how='outer')
pd_cols = [c for c in pd_df.columns if 'PD_' in c]

no_maint_blpd_uids = list()
no_maint_tgpd_uids = list()

res_df = list()
ind_count = 0
maint_count = 0
for uid, uid_df in df.groupby('UID'): #break
    pt_name = uid_df.iloc[0]['NAME']
    ibd_type = uid_df['IBD_TYPE'].iloc[0]
    uid_pd_df = pd_df[pd_df['UID'] == uid].copy()
    uid_df['DATE'] = uid_df['DATETIME'].map(lambda x: x.split('T')[0])

    min_date = uid_df['DATE'].min()
    min_date_bf2mo = (datetime.strptime(min_date, '%Y-%m-%d') - timedelta(days=60)).strftime('%Y-%m-%d')

    res_dict = {'UID': uid,
                'NAME':pt_name,
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
    maint_ind_type_str = ''
    if uid in ind_uids:
        # raise ValueError

        ## [IND 기간 처리]

        ind_count+=1
        # raise ValueError
        print(f'[IND ({ind_count})] / {uid}')
        res_dict['PHASE']='IND'

        try: dose3rd_date = uid_df[uid_df['MDV']==1].iloc[2]['DATE']
        except: dose3rd_date = (datetime.strptime(min_date, '%Y-%m-%d') + timedelta(days=90)).strftime('%Y-%m-%d')

        try: dose4th_date = uid_df[uid_df['MDV']==1].iloc[3]['DATE']
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
        res_dict['DS_DATE'] = min_date
        res_dict['DE_DATE'] = dose3rd_date
        res_dict['TG_DATE'] = tgpd_date
        for pdc in pd_cols:
            res_dict['TG_'+pdc] = tg_pd_row[pdc]

        res_df.append(res_dict)

        ## [Maint in IND 기간 처리]
        # uid_df['TIME']
        maint_in_ind_uid_df = uid_df[uid_df['DATE'] >= dose4th_date].copy()

        # Induction period만 존재하고 maintenance 기간은 존재하지 않는 경우 -> 비논리적 날짜로 설정하여 안 나오게
        if len(maint_in_ind_uid_df)==0:
            mii_min_date = maint_in_ind_uid_df['DATE'].max()
            mii_max_date = maint_in_ind_uid_df['DATE'].min()
        else:
            mii_min_date = max(maint_in_ind_uid_df['DATE'].min(), uid_pd_df['DATETIME'].min())
            mii_max_date = min(maint_in_ind_uid_df['DATE'].max(), uid_pd_df['DATETIME'].max())
        mii_pd_df = uid_pd_df[(uid_pd_df['DATETIME'] >= mii_min_date)&(uid_pd_df['DATETIME'] <= mii_max_date)].copy()
        mii_df =uid_df[(uid_df['DATE'] >= mii_min_date)&(uid_df['DATE'] <= mii_max_date)].copy()

        res_dict = {'UID': uid,
                    'NAME':pt_name,
                    'IBD_TYPE': ibd_type}

        # Induction period만 존재하고 maintenance 기간은 존재하지 않는 경우 -> 비논리적 날짜로 설정하여 안 나오게

        if len(mii_df)==0:
            min_date = uid_df['DATE'].min()
        else:
            min_date = mii_df['DATE'].min()
        min_date_bf2mo = (datetime.strptime(min_date, '%Y-%m-%d') - timedelta(days=60)).strftime('%Y-%m-%d')

        # BL Maint in IND 기록하기

        bl_pd_df = mii_pd_df[(mii_pd_df['DATETIME'] <= min_date) & (mii_pd_df['DATETIME'] >= min_date_bf2mo)].copy()
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
            res_dict['BL_' + pdc] = bl_pd_row[pdc]

        maint_ind_type_str = '_AFTIND'

    # else:
    maint_count+=1
    print(f'[MAINT ({maint_count})] {uid}')
    res_dict['PHASE']=f'MAINT{maint_ind_type_str}'
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
    res_dict['DS_DATE'] = min_date
    res_dict['DE_DATE'] = uid_df[(uid_df['MDV']==1)&(uid_df['DATE'] <= tgpd_date)]['DATE'].max()
    res_dict['TG_DATE'] = tgpd_date
    for pdc in pd_cols:
        res_dict['TG_'+pdc] = tg_pd_row[pdc]

    res_df.append(res_dict)
res_df = pd.DataFrame(res_df)
res_df['DAY_DELT'] = (pd.to_datetime(res_df['TG_DATE'])-pd.to_datetime(res_df['BL_DATE'])).dt.total_seconds()/86400

if not os.path.exists(f'{output_dir}/PKPD_EDA'):
    os.mkdir(f'{output_dir}/PKPD_EDA')
if not os.path.exists(f'{output_dir}/PKPD_EDA/PKPD_Corr0'):
    os.mkdir(f'{output_dir}/PKPD_EDA/PKPD_Corr0')

index_cols = ['UID','NAME','IBD_TYPE','PHASE','DAY_DELT','DS_DATE','DE_DATE']
res_df = res_df[index_cols + [c for c in res_df if c not in index_cols]]
res_df.to_csv(f"{output_dir}/PKPD_EDA/PKPD_Corr0/pd_eda_df.csv", index=False, encoding='utf-8-sig')


# res_df['PHASE']
# res_df.to_csv(f"{}")
# ind_res_df.columns
# ind_res_df = res_df[res_df['PHASE']=='Induction'].copy()
# ind_res_df[~(ind_res_df['BL_DATE'].isna())]
# ind_res_df[(ind_res_df['BL_DATE'].isna())]
# ind_res_df[~(ind_res_df['TG_DATE'].isna())]

# ind_na_pd_df = ind_res_df[(ind_res_df['TG_DATE'].isna())].copy()
# fdate_df = df.groupby('UID',as_index=False).agg({'NAME':'min','DATETIME':'min'}).rename(columns={'DATETIME':'FIRST_DATA_DATE'})
# ind_na_pd_df = ind_na_pd_df.merge(fdate_df,on=['UID'], how='left')[['UID','NAME','IBD_TYPE','FIRST_DATA_DATE']]
# first_pd_df = pd_df[pd_df['UID'].isin(ind_na_pd_df['UID'])].groupby('UID',as_index=False).agg({'PD_PRO2':'first','DATETIME':'min'}).rename(columns={'DATETIME':'FIRST_PD_DATE'})
# ind_na_pd_df = ind_na_pd_df.merge(first_pd_df,on=['UID'], how='left')
#
# ind_na_pd_df['FIRST_DATA_DATE'] = ind_na_pd_df['FIRST_DATA_DATE'].map(lambda x:x.split('T')[0])
# ind_na_pd_df.to_csv(f"{output_dir}/vacant_indpd_list.csv", index=False, encoding='utf-8-sig')
# ind_res_df[((ind_res_df['BL_DATE'].isna()))&((ind_res_df['TG_DATE'].isna()))]

# ind_res_df[(~(ind_res_df['BL_DATE'].isna()))&(~(ind_res_df['TG_DATE'].isna()))]