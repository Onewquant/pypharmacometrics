from tools import *
from pynca.tools import *
from statsmodels.stats.outliers_influence import variance_inflation_factor


# project_dict = {'KSCPTSPRWS25MULTI':['SGLT2INH'], 'KSCPTSPRWS25SINGLE':['SGLT2INH']}
# modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25"
# modeling_prepconc_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/modeling_prep_data"
resource_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/resource"
r_dataset_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/sglt2i_dataset"
results_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/results_r"

pk_res_df = pd.read_csv(r_dataset_dir_path + "/KSCPTSPRWS25_SGLT2i_PK.csv")

## Blood PD dataset 정리

# bpd_df = pd.read_excel(f"{resource_dir_path}/PD prep_Final_0707_BASE group.xlsx")
# bpd_df = bpd_df.iloc[1:]
# bpd_df = bpd_df.replace('.',np.nan)
# bpd_df['N_Day'] = bpd_df['N_DAY'].astype(float)
# bpd_df['N_TIME'] = bpd_df['N_TIME'].astype(float)
# bpd_df['C_Serum GLU'] = bpd_df['C_Serum GLU'].astype(float)
# bpd_df = bpd_df[(bpd_df['N_Day'].isin([1]))&(bpd_df['N_TIME'] <= 24)].copy()
# bpd_df = bpd_df.drop(columns=['HOMA-IR','C_Insulin','N_DAY','N_Day'])
# # bpd_df = bpd_df.rename(columns={'N_DAY':'N_DAY'})
#
# bpd_df.to_csv(f"{r_dataset_dir_path}/KSCPTSPRWS25_SGLT2i_Blood_PD.csv", index=False)
# bpd_df['ID'].unique()
bpd_df = pd.read_csv(f"{r_dataset_dir_path}/KSCPTSPRWS25_SGLT2i_Blood_PD.csv")
bpd_baseline_df = bpd_df[bpd_df['N_TIME']==0][['ID','N_TIME','C_Serum GLU','HbA1c']].rename(columns={'C_Serum GLU':'PG_base', 'HbA1c':'HbA1c_base'})

glu_nca_result = tblNCA(concData=bpd_df,key=['ID'],colTime="N_TIME",colConc='C_Serum GLU',down = "Log",dose=0.1,slopeMode='BEST',colStyle='pw')
glu_auc_df = glu_nca_result[['ID', 'AUClast']].rename(columns={'AUClast':'AUC_glu'}).iloc[1:]
glu_auc_df['PG_avg'] = glu_auc_df['AUC_glu'].astype(float)/24
glu_auc_df = glu_auc_df.merge(bpd_baseline_df[['ID','PG_base','HbA1c_base']], on=['ID'], how='left')

## Urine PD dataset 정리

# upd_df = pd.read_excel(f"{resource_dir_path}/PD urine_prep_Final_0707_prep_BASE group.xlsx")
# upd_df = upd_df[upd_df['N_Day'].isin([-1,1])].rename(columns={'N_Day':'N_DAY','N_Time':'N_TIME'}).drop(columns=['Accumulated Amount_GLU'])
# upd_df['Amount_GLU'] = upd_df['Amount_GLU'].astype(float)
# upd_df['Volume'] = upd_df['Volume'].astype(float)
# upd_df.to_csv(f'{r_dataset_dir_path}/KSCPTSPRWS25_SGLT2i_Urine_PD.csv',index=False)


upd_df = pd.read_csv(f'{r_dataset_dir_path}/KSCPTSPRWS25_SGLT2i_Urine_PD.csv')
uge24_df = upd_df.groupby(['ID','N_DAY']).agg({'Amount_GLU':'sum','Volume':'sum'}).reset_index(drop=False).rename(columns={'Amount_GLU':'UGE24','Volume':'VOLUME24'})
ugebase_df = uge24_df[uge24_df['N_DAY']==-1].copy().rename(columns={'UGE24':'UGEbase','VOLUME24':'VOLUMEbase'}).drop(['N_DAY'],axis=1)
uge24_df = uge24_df[(uge24_df['N_DAY']==1)].merge(ugebase_df, on=['ID'],how='left')

## COVAR dataset 정리

# id_grp_dict = {row['ID']:row['GRP'] for inx, row in pk_res_df[['ID','GRP']].drop_duplicates().iterrows()}
# covar_df['GRP'] = covar_df['ID'].map(id_grp_dict)
# covar_df = covar_df[['ID','GRP','AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'GFR', 'TBIL', 'SODIUM', 'AST', 'ALT']].rename(columns={'GFR':'eGFR'})
# covar_df.to_csv(f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/sglt2i_dataset/KSCPTSPRWS25_SGLT2i_COVAR.csv',index=False)

covar_df = pd.read_csv(f'{r_dataset_dir_path}/KSCPTSPRWS25_SGLT2i_COVAR.csv')

## Merge & Generate final prep-PD df

pd_res_df = glu_auc_df.merge(uge24_df, on=['ID'], how='left').merge(covar_df, on=['ID'], how='left')
pd_res_df['pd_res_df']