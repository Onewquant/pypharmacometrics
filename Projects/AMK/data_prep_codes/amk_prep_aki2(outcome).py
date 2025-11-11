import pandas as pd
import numpy as np
from datetime import datetime, timedelta

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

# sim_df = pd.read_csv(f'{output_dir}/sim114.csv')
sim_df = pd.read_csv(f'{nonmem_dir}/amk_modeling_df_covar.csv')
sim_df['ID'] = sim_df['ID'].astype(int)
dose_df = sim_df[sim_df['MDV']==1].copy()
conc_df = sim_df[sim_df['MDV']==0].copy()
# dose_df = pd.read_csv('./prepdata/amikacin_dose.csv')
# lab_df = pd.read_csv('./prepdata/amikacin_lab_Cr.csv')

"""
# (1) 48h 이내에 Cr 0.3 증가
# (2) 7일 이내에 기저치 * 1.5 이상으로 상승 
"""

cycle_df = list()
for inx, frag_df in dose_df.groupby('ID'):# break
    frag_df = frag_df.sort_values(['TIME'])
    cycle_df.append({'ID':inx, 'SDOSE_TIME':frag_df['TIME'].iloc[0], 'LDOSE_TIME':frag_df['TIME'].iloc[-1],'DOSING_COUNT':len(frag_df), 'TOTAL_DAYS':(frag_df['TIME'].iloc[-1]-frag_df['TIME'].iloc[0])/24, 'TOTAL_DOSE':frag_df['AMT'].sum()})

cycle_df = pd.DataFrame(cycle_df)
# dose_df[dose_df['ID']==148484881559882]

# def dtstr_calculation(dtstr, days, format='%Y-%m-%dT%H:%M'):
#     return (datetime.strptime(dtstr, format) + timedelta(days=days)).strftime(format)

aki_df = list()
result_cols = ['ID','COND_TYPE','BLCr_DT','BLCr_RSLT','AKI_DT','AKI_RSLT']
result_df = list()
for inx, row in cycle_df.iterrows(): #break
    # if inx[0]==148484876560382: raise ValueError
    # if row['ID']==9:
    #     raise ValueError
    id_cr_df = sim_df[sim_df['ID']==row['ID']].sort_values(['TIME'])

    print(f"({inx}) {row['ID']}")

    # before_month_dt = dtstr_calculation(dtstr=row['SDOSE_TIME'], days=-30)
    # after_2days_dt = dtstr_calculation(dtstr=row['LDOSE_TIME'], days=2)

    # before_month_dt = row['SDOSE_TIME'] -
    # after_2days_dt = (datetime.strptime(row['LDOSE_TIME'], '%Y-%m-%dT%H:%M') + timedelta(days=2)).strftime('%Y-%m-%dT%H:%M')

    # afterdose_df = id_lab_df[(id_lab_df['TIME'] > row['SDOSE_TIME'])].copy()
    # recentbl_df = id_lab_df.iloc[:1,:].copy()

    # id_cr_df.columns
    # id_cr_df = pd.concat([recentbl_df, afterdose_df])
    for cr_inx, cr_row in id_cr_df.iterrows():#break

        bl_dt = cr_row['TIME']
        bl_cr = cr_row['CREATININE']

        aki_48h_dt = cr_row['TIME'] + 24*2
        aki_7days_dt = cr_row['TIME'] + 24*7

        # aki_48h_dt = dtstr_calculation(dtstr=cr_row['DT'], days=2)
        # aki_7days_dt = dtstr_calculation(dtstr=cr_row['DT'], days=7)

        # (1) 48h 이내에 Cr 0.3 증가
        # (2) 7일 이내에 기저치 * 1.5 이상으로 상승

        if (((id_cr_df['TIME'] < aki_48h_dt) & (id_cr_df['TIME'] > bl_dt))*1).sum() > 0:
            cond_exist_df = id_cr_df[(id_cr_df['TIME'] < aki_48h_dt) & (id_cr_df['TIME'] > bl_dt)].copy()
            aki_exist_df = cond_exist_df[cond_exist_df['CREATININE'] >= bl_cr+0.3].copy()
            if len(aki_exist_df)>0:
                aki_exist_df = aki_exist_df.rename(columns={'TIME':'AKI_DT','CREATININE':'AKI_RSLT'})
                aki_exist_df['COND_TYPE'] = '48h'
                aki_exist_df['BLCr_DT'] = bl_dt
                aki_exist_df['BLCr_RSLT'] = bl_cr
                result_df.append(aki_exist_df[result_cols].copy())


        if (((id_cr_df['TIME'] < aki_7days_dt) & (id_cr_df['TIME'] > bl_dt)) * 1).sum() > 0:
            cond_exist_df = id_cr_df[(id_cr_df['TIME'] < aki_7days_dt) & (id_cr_df['TIME'] > bl_dt)].copy()
            aki_exist_df = cond_exist_df[cond_exist_df['CREATININE'] >= bl_cr*1.5].copy()
            if len(aki_exist_df)>0:
                aki_exist_df = aki_exist_df.rename(columns={'TIME':'AKI_DT','CREATININE':'AKI_RSLT'})
                aki_exist_df['COND_TYPE'] = '7days'
                aki_exist_df['BLCr_DT'] = bl_dt
                aki_exist_df['BLCr_RSLT'] = bl_cr
                result_df.append(aki_exist_df[result_cols].copy())

result_df = pd.concat(result_df)

filt_res_df = result_df[result_df['BLCr_RSLT']<=1.2].copy()
filt_res_df_new = filt_res_df[filt_res_df['AKI_RSLT']>1.2].copy()

result_df.to_csv(f'{output_dir}/amk_aki.csv', index=False, encoding='utf-8')
print(f"AKI cases: {len(result_df['ID'].drop_duplicates())}")
# filt_res_df.to_csv(f'{output_dir}/amk_aki(filtered).csv', index=False, encoding='utf-8')
# print(f"AKI cases: {len(filt_res_df['ID'].drop_duplicates())}")
# filt_res_df_new.to_csv(f'{output_dir}/amk_aki(filtered).csv', index=False, encoding='utf-8')
# print(f"AKI cases: {len(filt_res_df_new['ID'].drop_duplicates())}")


# result_df['ID'].drop_duplicates()
# filt_res_df_new['ID'].drop_duplicates()
# filt_res_df['ID'].drop_duplicates()
# sim_df['ID'].drop_duplicates()

#(datetime.strptime(row['SDOSE_DT'], '%Y-%m-%dT%H:%M') - timedelta(days=30)).strftime('%Y-%m-%dT%H:%M')


# cycle_df = pd.DataFrame(cycle_df)

# cycle_df.sort_values(['ID'])
# lab_df.drop_duplicates(['ID']).sort_values('ID')
# 148484876560382
# 148484930717343