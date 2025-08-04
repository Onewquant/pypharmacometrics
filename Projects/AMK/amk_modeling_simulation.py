from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## DEMO Covariates Loading

simulation_df = pd.read_csv(f"{output_dir}/amk_modeling_df_covar_filt.csv")
# simulation_df['TIME_DIFF'] = simulation_df['TIME'].diff()
# simulation_df[simulation_df['TIME_DIFF']>0]['TIME_DIFF'].sort_values()
min_interval = 10/60

final_sim_df = list()
for id, id_sim_df in simulation_df.groupby('ID'): #break
    # raise ValueError
    print(f"({id}) {id_sim_df['UID'].iloc[0]}")
    id_conc_rows = id_sim_df[id_sim_df['MDV']==0].copy()
    id_conc_rows['DV'] = np.nan
    id_dose_rows = id_sim_df[id_sim_df['MDV']==1].copy()

    internal_start_time = id_sim_df['TIME'].min()
    internal_end_time = id_sim_df['TIME'].max()
    external_end_time = internal_end_time + 4*24

    time_arr = np.arange(internal_start_time, external_end_time + min_interval, min_interval)  # stop보다 큰 수까지 포함하려면 stop + step
    time_arr = time_arr[time_arr <= external_end_time]  # 정확히 stop보다 작거나 같은 값만 남기기
    uniq_time_arr = list(set(time_arr)-set(id_conc_rows['TIME']))
    uniq_time_arr.sort()

    added_sim_df = pd.DataFrame(columns=id_conc_rows.columns)
    added_sim_df['TIME'] = uniq_time_arr
    added_sim_df['ID'] = id
    # added_sim_df['DV'] = np.nan
    # added_sim_df['DV'] = np.nan
    added_sim_df['MDV'] = 0
    added_sim_df['CMT'] = np.int64(id_conc_rows['CMT'].median())
    added_sim_df['AMT'] = '.'
    added_sim_df['RATE'] = '.'

    final_id_sim_df = pd.concat([id_conc_rows, id_dose_rows, added_sim_df]).sort_values(['TIME','MDV'])
    final_id_sim_df = final_id_sim_df.fillna(method='ffill').fillna(method='bfill')
    # id_conc_rows['RATE'].iloc[0]
    # id_conc_rows['CMT']

    # added_sim_df.columns
    final_sim_df.append(final_id_sim_df)

final_sim_df = pd.concat(final_sim_df, ignore_index=True)
final_sim_df['TAD'] = final_sim_df['TIME']
final_sim_df.to_csv(f"{output_dir}/amk_simulation_df.csv",index=False, encoding='utf-8-sig')

# final_sim_df[final_sim_df['ID'].isin([1,2,4,5])].to_csv(f"{output_dir}/amk_simulation_partial_df.csv",index=False, encoding='utf-8-sig')

# final_sim_df['ID'].drop_duplicates()
# final_sim_df['TAD'] = add_time_after_dosing_column(df=final_sim_df)



