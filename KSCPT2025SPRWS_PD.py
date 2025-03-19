from tools import *
from pynca.tools import *


project_dict = {'KSCPTSPRWS25MULTI':['SGLT2INH'], 'KSCPTSPRWS25SINGLE':['SGLT2INH']}
modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25"
modeling_prepconc_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/modeling_prep_data"
resource_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/resource"

## AUC 구하기

mdprep_conc_files = glob.glob(modeling_prepconc_dir_path+'/MDP_*.csv')
mdprep_conc_df = pd.read_csv(mdprep_conc_files[0]).copy()
mdprep_conc_df = mdprep_conc_df[(mdprep_conc_df['MDV']!=1)&(~(mdprep_conc_df['TAD'].isin([96, 120])))].copy()
mdprep_conc_df = mdprep_conc_df.drop(columns=['MDV','AMT','RATE','DUR','CMT'])
mdprep_conc_df['DAY'] = mdprep_conc_df['TAD'].map(lambda x:1 if x < 144 else 7)

adj_col_df = list()
for inx, df_frag in mdprep_conc_df.groupby(['ID','DAY','UID','GRP']):

    adj_tad_df = df_frag[['ID','UID','GRP', 'TAD', 'TIME']].copy()
    if inx[1]==1:
        adj_tad_df['ADJTAD'] = adj_tad_df['TAD'].copy()
        adj_tad_df['ADJTIME'] = adj_tad_df['TIME'].copy()
        adj_col_df.append(adj_tad_df)
    else:
        min_tad_rows = adj_tad_df[adj_tad_df['TAD']==144]
        if len(min_tad_rows)==0:
            raise ValueError
        else:
            min_tad = min_tad_rows.iloc[0]['TAD']
            min_time = min_tad_rows.iloc[0]['TIME']
        adj_tad_df['ADJTAD'] = adj_tad_df['TAD']-min_tad
        adj_tad_df['ADJTIME'] = adj_tad_df['TIME']-min_time
        adj_col_df.append(adj_tad_df)

adj_col_df = pd.concat(adj_col_df, ignore_index=True)
mdprep_conc_df = mdprep_conc_df.merge(adj_col_df, on=['ID','UID','GRP','TAD','TIME'], how='left')
mdprep_conc_df = mdprep_conc_df[mdprep_conc_df['ADJTAD'] <= 24].copy()
mdprep_conc_df['Short_UID'] = mdprep_conc_df['UID'].map(lambda x:int(x[-2:]))
mdprep_conc_df['UID'] = mdprep_conc_df['UID'].map(lambda x:x[1:])

# mdprep_conc_df.to_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/modeling_prep_data/calculate_iAUC.csv", index=False, encoding='utf-8-sig')

# mdprep_conc_df['TADADJ'] =

# iAUC_df=pd.DataFrame({'Name':['T0-6','T6-12','T12-24'], 'Start':[0,6,12], 'End':[6,12,24]})
# result = tblNCA(concData=mdprep_conc_df,key=['ID','DAY'],colTime="A_TIME",colConc='CONC',down = "Linear",iAUC=iAUC_df,dose=0.5,slopeMode='BEST',colStyle='pw')
# iauc_result = tblNCA(concData=mdprep_conc_df,key=['ID','DAY'],colTime="ADJTIME",colConc='DV',down = "Log",iAUC=iAUC_df,dose=0.5,slopeMode='BEST',colStyle='pw')


## PD와 연결

iauc_result = pd.read_csv(modeling_prepconc_dir_path+'/sglt2i_iauc.csv')
iauc_result = iauc_result.rename(columns={'Short_UID':'UID'})
iauc_result['UID'] = iauc_result['UID'].map(lambda x:'S'+str(x).zfill(3))
iauc_res_df = iauc_result.melt(id_vars=list(iauc_result.columns[:-4]), var_name='AUCTIMEINT', value_name='IAUC')
iauc_res_df['AUCTIMEINT'] = iauc_res_df['AUCTIMEINT'].map(lambda x:x.split('[')[-1].split('h]')[0])

## Blood PD dataset 정리

bpd_df = pd.read_excel(f"{resource_dir_path}/PD prep_Final_0707_BASE group.xlsx")
bpd_df = bpd_df.iloc[1:]
bpd_df = bpd_df.replace('.',np.nan)
bpd_df['DAY'] = bpd_df['N_DAY'].astype(float)
bpd_df['N_TIME'] = bpd_df['N_TIME'].astype(float)
bpd_df['C_Serum GLU'] = bpd_df['C_Serum GLU'].astype(float)
bpd_df = bpd_df[(bpd_df['DAY'].isin([1,7]))&(bpd_df['N_TIME'] <= 24)].copy()

pg_zero_df = bpd_df[bpd_df['N_TIME']==0][['ID','DAY','C_Serum GLU']].rename(columns={'C_Serum GLU':'PG_ZERO'})

glu_nca_result = tblNCA(concData=bpd_df,key=['ID','DAY'],colTime="N_TIME",colConc='C_Serum GLU',down = "Log",dose=0.1,slopeMode='BEST',colStyle='pw')
# glu_nca_result.columns
s_glu_df = glu_nca_result[['ID','DAY','AUClast']].iloc[1:]
s_glu_df['PG_AVG'] = s_glu_df['AUClast']/24
s_glu_df = s_glu_df.merge(pg_zero_df, on=['ID','DAY'], how='left')
s_glu_df = s_glu_df.rename(columns={'AUClast':'AUC_GLU','DAY':'N_Day'})

# bpd_df[bpd_df['ID']=='S031']

## Urine PD dataset 정리

upd_df = pd.read_excel(f"{resource_dir_path}/PD urine_prep_Final_0707_prep_BASE group.xlsx")
# upd_df = upd_df.iloc[1:].copy()
# upd_df.columns
upd_df = upd_df[upd_df['N_Time'].isin(['0-24','0-6','6-12','12-24'])]
upd_df['UID'] = upd_df['ID']
upd_df['DAY'] = upd_df['N_Day']
upd_df['GRP'] = upd_df['GRP'].map(int)
upd_df['AUCTIMEINT'] = upd_df['N_Time']
upd_df['Amount_GLU'] = upd_df['Amount_GLU'].astype(float)
upd_df['Volume'] = upd_df['Volume'].astype(float)
upd_addcol_df = upd_df.groupby(['ID','DAY']).agg({'Amount_GLU':'sum','Volume':'sum'}).reset_index(drop=False).rename(columns={'Amount_GLU':'UGE24','Volume':'VOLUME24'})
upd_baseline_df = upd_addcol_df[upd_addcol_df['DAY']==-1].copy().rename(columns={'UGE24':'BASELINE'}).drop(['DAY','VOLUME24'],axis=1)
upd_addcol_df = upd_addcol_df.merge(upd_baseline_df, on=['ID'],how='left')
upd_df = upd_df.merge(upd_addcol_df, on=['ID','DAY'],how='left')

# daily auc만으로 정리
upd_daily_df = upd_df.copy()
upd_daily_df = upd_daily_df[upd_daily_df['Accumulated Amount_GLU']==upd_daily_df['UGE24']].copy()
upd_daily_df['N_Time'] = '0-24'
upd_daily_df['Volume'] = upd_daily_df['VOLUME24']
upd_daily_df['Amount_GLU'] = upd_daily_df['UGE24']
upd_daily_df['AUCTIMEINT'] = upd_daily_df['N_Time']
upd_daily_df = upd_daily_df.drop_duplicates(['UID','N_Day','N_Time'])

# 공변량 붙이기
covar_df = mdprep_conc_df.rename(columns={'ID':'NID','UID':'ID'})[['NID','ID']+list(mdprep_conc_df.columns)[7:15]].drop_duplicates(['ID'])
upd_daily_df = upd_daily_df.merge(covar_df, on=['ID'],how='left')
upd_daily_df = upd_daily_df.merge(s_glu_df, on=['ID','N_Day'],how='left')

# EDA
pre_df = upd_daily_df[upd_daily_df['DAY']==-1].copy()
nss_df = upd_daily_df[upd_daily_df['DAY']==1].copy()
ss_df = upd_daily_df[upd_daily_df['DAY']==7].copy()

nss_df['UGE24']
nss_df['UGE24']