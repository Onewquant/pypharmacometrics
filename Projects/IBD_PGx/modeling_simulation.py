from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## DEMO Covariates Loading

simulation_df = pd.read_csv(f"{output_dir}/infliximab_integrated_modeling_df_dayscale.csv")
# simulation_df['TIME_DIFF'] = simulation_df['TIME'].diff()
# simulation_df[simulation_df['TIME_DIFF']>0]['TIME_DIFF'].sort_values()
# interval = 3/24
interval = 1
final_sim_df = get_model_population_sim_df(df=simulation_df, interval=interval, add_on_period=0)
final_sim_df.to_csv(f"{output_dir}/infliximab_integrated_simulation_df.csv",index=False, encoding='utf-8-sig')





# uniq_ids = list(final_sim_df['ID'].unique())[:30]
# final_sim_df[final_sim_df['ID'].isin(uniq_ids)].to_csv(f"{output_dir}/amk_simulation_partial_df.csv",index=False, encoding='utf-8-sig')

# final_sim_df['ID'].drop_duplicates()
# final_sim_df['TAD'] = add_time_after_dosing_column(df=final_sim_df)



