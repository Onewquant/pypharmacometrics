from tools import *
from pynca.tools import *
from statsmodels.stats.outliers_influence import variance_inflation_factor


# project_dict = {'KSCPTSPRWS25MULTI':['SGLT2INH'], 'KSCPTSPRWS25SINGLE':['SGLT2INH']}
# modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25"
# modeling_prepconc_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/modeling_prep_data"
resource_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/resource"
r_dataset_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/sglt2i_dataset"
results_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/results_r"


pk_res_df = pd.read_csv(results_dir_path + "/[WSCT] NCARes_PK.csv")

## Blood PD dataset 정리

# bpd_df = pd.read_excel(f"{resource_dir_path}/PD prep_Final_0707_BASE group.xlsx")
# bpd_df = bpd_df.iloc[1:]
# bpd_df = bpd_df.replace('.',np.nan)
# bpd_df['N_Day'] = bpd_df['N_DAY'].astype(float)
# bpd_df['N_TIME'] = bpd_df['N_TIME'].astype(float)
# bpd_df['C_Serum GLU'] = bpd_df['C_Serum GLU'].astype(float)
# bpd_df = bpd_df[(bpd_df['N_Day'].isin([1]))&(bpd_df['N_TIME'] <= 24)].copy()
# bpd_df = bpd_df.drop(columns=['HOMA-IR','C_Insulin','N_DAY'])
#
# bpd_df.to_csv(f"{r_dataset_dir_path}/KSCPTSPRWS25_SGLT2i_Blood_PD.csv", index=False)

bpd_df = pd.read_csv(f"{r_dataset_dir_path}/KSCPTSPRWS25_SGLT2i_Blood_PD.csv")
bpd_baseline_df = bpd_df[bpd_df['N_TIME']==0][['ID','N_Day','C_Serum GLU','HbA1c']].rename(columns={'C_Serum GLU':'PG_base', 'HbA1c':'HbA1c_base'})

glu_nca_result = tblNCA(concData=bpd_df,key=['ID','N_Day'],colTime="N_TIME",colConc='C_Serum GLU',down = "Log",dose=0.1,slopeMode='BEST',colStyle='pw')
s_glu_df = glu_nca_result[['ID','N_Day','AUClast']].iloc[1:]
s_glu_df['PG_avg'] = s_glu_df['AUClast'].astype(float)/24
s_glu_df = s_glu_df.merge(bpd_baseline_df, on=['ID','N_Day'], how='left')
s_glu_df = s_glu_df.rename(columns={'AUClast':'AUC_glu'})

## Urine PD dataset 정리


upd_df = pd.read_excel(f"{resource_dir_path}/PD urine_prep_Final_0707_prep_BASE group.xlsx")
# upd_df = upd_df.iloc[1:].copy()
# upd_df.columns
upd_df = upd_df[upd_df['N_Time'].isin(['0-24','0-6','6-12','12-24'])]
upd_df = upd_df.rename(columns={'ID':'UID'})
upd_df['DAY'] = upd_df['N_Day']
upd_df['GRP'] = upd_df['GRP'].map(int)
upd_df['AUCTIMEINT'] = upd_df['N_Time']
upd_df['Amount_GLU'] = upd_df['Amount_GLU'].astype(float)
upd_df['Volume'] = upd_df['Volume'].astype(float)
upd_addcol_df = upd_df.groupby(['UID','DAY']).agg({'Amount_GLU':'sum','Volume':'sum'}).reset_index(drop=False).rename(columns={'Amount_GLU':'UGE24','Volume':'VOLUME24'})
upd_baseline_df = upd_addcol_df[upd_addcol_df['DAY']==-1].copy().rename(columns={'UGE24':'UGEbase'}).drop(['DAY','VOLUME24'],axis=1)
upd_addcol_df = upd_addcol_df.merge(upd_baseline_df, on=['UID'],how='left')
upd_df = upd_df.merge(upd_addcol_df, on=['UID','DAY'],how='left')
# upd_df[['UID','UGE24']]
# upd_df[(upd_df['DAY']==1)&(upd_df['UID']=='S003')]

# daily auc만으로 정리
upd_daily_df = upd_df.copy()
upd_daily_df = upd_daily_df[upd_daily_df['Accumulated Amount_GLU']==upd_daily_df['UGE24']].copy()
upd_daily_df['N_Time'] = '0-24'
upd_daily_df['Volume'] = upd_daily_df['VOLUME24']
upd_daily_df['Amount_GLU'] = upd_daily_df['UGE24']
upd_daily_df['AUCTIMEINT'] = upd_daily_df['N_Time']
upd_daily_df = upd_daily_df.drop_duplicates(['UID','N_Day','N_Time'])
upd_daily_df = upd_daily_df[upd_daily_df['N_Day']==1]
len(upd_daily_df['UID'].unique())