from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr



prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## DEMO Covariates Loading

simulation_df = pd.read_csv(f"{output_dir}/infliximab_integrated_pdeda_df_dayscale.csv")
final_sim_df = pd.read_csv(f"{nonmem_dir}/run/sim57",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')

# sim_conc_df = final_sim_df[final_sim_df['MDV']==1].copy()
# realworld_df = simulation_df[simulation_df['MDV']==0].copy()

ibd_type_dict = {1:'CD',2:'UC'}

# final_sim_df.columns
# final_sim_df
sim_conc_df = final_sim_df[(final_sim_df['MDV']==0)&(final_sim_df['REALDATA']==1)].copy()
# final_sim_df[['ID','TIME','DV','MDV','PD_PRO2','REALDATA']]
# sim_conc_df = sim_conc_df[['ID','IBD_TYPE','TIME','DV','PD_PRO2']].copy()



# 변수 간 spearman 상관계수 및 p-value 계산
corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['PD_PRO2'])
print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")


corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['CRP'])
print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")


corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['CALPRTSTL'])
print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")


corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['SEX'])
print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")


corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['AGE'])
print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")


corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['HT'])
print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")