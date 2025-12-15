from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

# rdf = pd.read_csv(f'{output_dir}/infliximab_induction_modeling_df.csv')
# rdf['ALB'].median()
# rdf['CRP'].median()
# rdf['CREATININE'].median()

# demo2_df = pd.read_csv(f"{resource_dir}/demo2/IBD_PGx_demo2.csv")
# demo2_df = demo2_df.rename(columns={'EMR ID':'UID','height':'HT','weight':'WT','bmi':'BMI','bsa':'BSA','name':'NAME'})
# demo2_df = demo2_df[['UID','NAME']].copy()
# demo2_df['UID'] = demo2_df['UID'].astype(str)

## WT, HT, ADHERENCE Covariates / PD markers Loading

pd_marker_list = ['PD_PRO2', 'PD_STLFREQ', 'PD_ABDPAIN', 'PD_RECTBLD','PD_TOTALSCORE', 'PD_CR']

cdai_hx_df = pd.read_csv(f"{output_dir}/pdmarker_cdai_df.csv")
cdai_hx_df = cdai_hx_df.rename(columns={'ID':'UID','DATE':'DATETIME','CDAI_TOTALSCORE':'PD_TOTALSCORE','CDAI_ADHERENCEPCT':'ADHERENCE','CDAI_WEIGHT':'WT','CDAI_HEIGHT':'HT', 'CDAI_BMI':'BMI'})
cdai_hx_df['UID'] = cdai_hx_df['UID'].astype(str)
cdai_hx_df['PD_PRO2'] = 2*cdai_hx_df['CDAI_DIARHEACNT'] + 5*cdai_hx_df['CDAI_ABDPAINCNT']
cdai_hx_df['PD_STLFREQ'] = cdai_hx_df['CDAI_DIARHEACNT']
cdai_hx_df['PD_ABDPAIN'] = cdai_hx_df['CDAI_ABDPAINCNT']
cdai_hx_df['PD_RECTBLD'] = np.nan
cdai_hx_df['PD_CR'] = (((cdai_hx_df['PD_STLFREQ'] <= 3)&(cdai_hx_df['PD_ABDPAIN'] <= 1))|(cdai_hx_df['PD_PRO2']<=9))*1
cdai_hx_df['IBD_TYPE'] = 'CD'
# cdai_adh_df.columns

pms_hx_df = pd.read_csv(f"{output_dir}/pdmarker_pms_df.csv")
pms_hx_df = pms_hx_df.rename(columns={'ID':'UID'}).rename(columns={'ID':'UID','DATE':'DATETIME','PMS_TOTALSCORE':'PD_TOTALSCORE','PMS_ADHERENCEPCT':'ADHERENCE','PMS_WEIGHT':'WT','PMS_HEIGHT':'HT', 'PMS_BMI':'BMI'})
pms_hx_df['UID'] = pms_hx_df['UID'].astype(str)
pms_hx_df['PD_PRO2'] = pms_hx_df['PMS_STLCNT'] + pms_hx_df['PMS_HEMATOCHEZIA']
pms_hx_df['PD_STLFREQ'] = pms_hx_df['PMS_STLCNT']
pms_hx_df['PD_ABDPAIN'] = np.nan
pms_hx_df['PD_RECTBLD'] = pms_hx_df['PMS_HEMATOCHEZIA']
pms_hx_df['PD_CR'] = (((pms_hx_df['PD_STLFREQ']<=1)&(pms_hx_df['PD_ABDPAIN']==0))|(pms_hx_df['PD_PRO2']==0))*1
pms_hx_df['IBD_TYPE'] = 'UC'
# pms_adh_df.columns
# pd_df[pd_df['UID']=='29702679'][['DATETIME','PD_PRO2']]
pd_df = pd.concat([cdai_hx_df[['UID','DATETIME','ADHERENCE','IBD_TYPE'] + pd_marker_list], pms_hx_df[['UID','DATETIME','ADHERENCE','IBD_TYPE']+pd_marker_list]]).drop_duplicates(['UID','DATETIME']).sort_values(['UID','DATETIME'])
bsize_df = pd.concat([cdai_hx_df[['UID','DATETIME','WT','HT','BMI']], pms_hx_df[['UID','DATETIME','WT','HT','BMI']]]).drop_duplicates(['UID','DATETIME']).sort_values(['UID','DATETIME'])
# adh_df[adh_df['WT'].isna()]
# adh_df[adh_df['HT'].isna()]

full_result_df = list()
count = 0
for uid, uid_df in pd_df.groupby('UID',as_index=False): #break

    min_lab_date = uid_df['DATETIME'].min()
    max_lab_date = uid_df['DATETIME'].max()

    print(f"({count}) {uid} / Date Range : ({min_lab_date},  {max_lab_date})")

    uid_fulldt_df = pd.DataFrame(columns=['UID','DATETIME'])
    uid_fulldt_df['DATETIME'] = pd.date_range(start=min_lab_date,end=max_lab_date).astype(str)
    uid_fulldt_df['UID'] = uid

    uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')#.fillna(method='bfill')
    # uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')
    # uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')

    if count==0:
        full_result_df = uid_fulldt_df.copy()
    else:
        full_result_df = pd.concat([full_result_df, uid_fulldt_df.copy()])

    count+=1

# raise ValueError
pd_df = full_result_df.copy()
# pdm_med_uc = adh_df[adh_df['IBD_TYPE']=='UC']['PD_MARKER'].median()
# pdm_med_cd = adh_df[adh_df['IBD_TYPE']=='CD']['PD_MARKER'].median()
#
# pdpro2_med_uc = adh_df[adh_df['IBD_TYPE']=='UC']['PD_PRO2'].median()
# pdpro2_med_cd = adh_df[adh_df['IBD_TYPE']=='CD']['PD_PRO2'].median()
#
# med_pd_value = {'UC':pdm_med_uc,'CD':pdm_med_cd}
# full_result_df['PD_MARKER'] = full_result_df.apply(lambda row:row['PD_MARKER'] if not np.isnan(row['PD_MARKER']) else med_pd_value[row['IBD_TYPE']], axis=1)
# full_result_df['PD_PRO2'] = full_result_df.apply(lambda row:row['PD_PRO2'] if not np.isnan(row['PD_PRO2']) else 7777777, axis=1)
# adh_df = full_result_df.fillna(full_result_df.median(numeric_only=True))
# adh_df['PD_PRO2'] = adh_df['PD_PRO2'].replace(7777777,np.nan)

pd_df = pd_df.drop('IBD_TYPE', axis=1)

## bsize

full_result_df = list()
count = 0
for uid, uid_df in bsize_df.groupby('UID',as_index=False): #break

    min_lab_date = uid_df['DATETIME'].min()
    max_lab_date = uid_df['DATETIME'].max()

    print(f"({count}) {uid} / Date Range : ({min_lab_date},  {max_lab_date})")

    uid_fulldt_df = pd.DataFrame(columns=['UID','DATETIME'])
    uid_fulldt_df['DATETIME'] = pd.date_range(start=min_lab_date,end=max_lab_date).astype(str)
    uid_fulldt_df['UID'] = uid

    uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')#.fillna(method='ffill')
    # uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')

    if count==0:
        full_result_df = uid_fulldt_df.copy()
    else:
        full_result_df = pd.concat([full_result_df, uid_fulldt_df.copy()])

    count+=1

# bsize_df = full_result_df.fillna(full_result_df.median(numeric_only=True))


pd_bsize_df = pd_df.merge(bsize_df, on=['UID','DATETIME'], how='left')
pd_bsize_df = pd_bsize_df[(~pd_bsize_df.loc[:,'ADHERENCE':].isna()).sum(axis=1)!=0].copy()
# pd_bsize_df['UID'].drop_duplicates()
# pd_bsize_df = pd_bsize_df.merge(demo2_df[['UID','NAME']].copy(),on=['UID'],how='left')
# pd_bsize_df.to_csv(f"{output_dir}/pd_bsize_df(filt).csv", encoding='utf-8-sig', index=False)
# pd_bsize_df.to_csv(f"{output_dir}/pd_bsize_df_ori.csv", encoding='utf-8-sig', index=False)
