import pandas as pd
import numpy as np
from typing import Dict, Any
import glob

output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'
# comed_df = pd.read_csv(f"{output_dir}/lnz_final_comed_df.csv")
# transfusion_df = pd.read_csv(f"{output_dir}/lnz_final_transfusion_df.csv")
# surv_res_df = pd.read_csv(f"{output_dir}/b1da/b1da_lnz_surv_res_df(365).csv")
#
# multivar_res_df = list()
#
# demo_cols = ['SUBTOTAL_N','EV = 1','SEX = 1','ELD = 1','AGE','HT','WT','BMI']
# total_dose_cols = ['CUM_DOSE(TOTAL)','DOSE_PERIOD(TOTAL)','DOSE24(TOTAL_ACTIVE)','DOSE24(TOTAL_PERIOD)','DOSE24PERWT(TOTAL_ACTIVE)','DOSE24PERWT(TOTAL_PERIOD)','DOSE_INTERVAL(TOTAL)']
# ev_dose_cols = ['CUM_DOSE', 'DOSE_PERIOD','DOSE_PERIOD(ACTIVE)', 'DOSE24','DOSE24(ACTIVE)', 'DOSE24PERWT','DOSE24PERWT(ACTIVE)', 'DOSE_INTERVAL']
#
# comed_cols = list(comed_df['DRUG'].drop_duplicates())
# transfusion_cols = list(transfusion_df['TF_TYPE'].drop_duplicates())

incidence_files = glob.glob(f"{output_dir}/b1da/incidence_output/b1da_lnz_incidence_table(*).csv")
incidence_files.sort(reverse=True)
inc_period_list = list()
crude_df = pd.DataFrame()
km_df = pd.DataFrame()
for inc_inx,inc_file in enumerate(incidence_files): #break
    inc_period = int(inc_file.split('(')[-1].split(')')[0])
    inc_period_list.append(inc_period)

    fdf = pd.read_csv(inc_file)
    # fdf.columns

    crude_fdf = fdf[['ENDPOINT','incidence (crude analysis)']].rename(columns={'incidence (crude analysis)':inc_period})
    km_fdf = fdf[['ENDPOINT','cumulative incidence (KM analysis)']].rename(columns={'cumulative incidence (KM analysis)':inc_period})

    if inc_inx==0:
        crude_df = crude_fdf.copy()
        km_df = km_fdf.copy()
    else:
        crude_df = crude_df.merge(crude_fdf, on=['ENDPOINT'], how='left')
        km_df = km_df.merge(km_fdf, on=['ENDPOINT'], how='left')

sort_day_cols = list(crude_df.columns)[1:]
sort_day_cols.sort()

crude_df = crude_df[['ENDPOINT']+sort_day_cols].copy()
km_df = km_df[['ENDPOINT']+sort_day_cols].copy()


vacant_df = pd.DataFrame([{c:'' for c in crude_df.columns}])
crude_header_df = vacant_df.copy()
crude_header_df.at[0,'ENDPOINT'] = 'Cumulative incidence with Crude analysis - Percent value (event/total)'

km_header_df = vacant_df.copy()
km_header_df.at[0,'ENDPOINT'] = 'Cumulative incidence with K-M survival analysis - Percent estimate (95%CI)'

incidence_table = pd.concat([crude_header_df, crude_df, km_header_df, km_df])
incidence_table.to_csv(f"{output_dir}/B1DA/b1da_lnz_cummulative_incidence_table.csv", index=False)

