import pandas as pd
import numpy as np
from datetime import datetime, timedelta

conc_df = pd.read_csv('./prepdata/amikacin_conc.csv')
dose_df = pd.read_csv('./prepdata/amikacin_dose.csv')
lab_df = pd.read_csv('./prepdata/amikacin_lab_Cr.csv')

"""
# (1) 48h 이내에 Cr 0.3 증가
# (2) 7일 이내에 기저치 * 1.5 이상으로 상승 
"""

cycle_df = list()
for inx, frag_df in dose_df.groupby(['ID']):# break
    # if inx[0]==148484876560382: raise ValueError
    cycle_period = (pd.to_datetime(frag_df['DT'], format='%Y-%m-%dT%H:%M').diff().dt.total_seconds()/86400) >= 4
    if len(cycle_period[cycle_period])>0:
        cycle1_end_inx = cycle_period[cycle_period].index[0]
        rem_dose_df = frag_df.loc[:cycle1_end_inx-1].sort_values(['DT'])
    else:
        rem_dose_df = frag_df.sort_values(['DT'])

    cycle_df.append({'ID':inx[0], 'SDOSE_DT':rem_dose_df['DT'].iloc[0], 'LDOSE_DT':rem_dose_df['DT'].iloc[-1],'DOSING_COUNT':len(rem_dose_df), 'TOTAL_DAYS':(pd.to_datetime(rem_dose_df['DT'], format='%Y-%m-%dT%H:%M').diff().dt.total_seconds()/86400).sum(), 'TOTAL_DOSE':rem_dose_df['DOSE'].sum()})

cycle_df = pd.DataFrame(cycle_df)
# dose_df[dose_df['ID']==148484881559882]

def dtstr_calculation(dtstr, days, format='%Y-%m-%dT%H:%M'):
    return (datetime.strptime(dtstr, format) + timedelta(days=days)).strftime(format)

aki_df = list()
result_cols = ['ID','COND_TYPE','BLCr_DT','BLCr_RSLT','AKI_DT','AKI_RSLT']
result_df = list()
for inx, row in cycle_df.iterrows(): #break
    # if inx[0]==148484876560382: raise ValueError
    id_lab_df = lab_df[lab_df['ID']==row['ID']].sort_values(['DT'])

    print(f"({inx}) {row['ID']}")

    before_month_dt = dtstr_calculation(dtstr=row['SDOSE_DT'], days=-30)
    after_2days_dt = dtstr_calculation(dtstr=row['LDOSE_DT'], days=2)

    before_month_dt = (datetime.strptime(row['SDOSE_DT'], '%Y-%m-%dT%H:%M') - timedelta(days=30)).strftime('%Y-%m-%dT%H:%M')
    after_2days_dt = (datetime.strptime(row['LDOSE_DT'], '%Y-%m-%dT%H:%M') + timedelta(days=2)).strftime('%Y-%m-%dT%H:%M')

    baseline_df = id_lab_df[(id_lab_df['DT'] >= before_month_dt)&(id_lab_df['DT'] <= row['SDOSE_DT'])]
    afterdose_df = id_lab_df[(id_lab_df['DT'] > row['SDOSE_DT'])&(id_lab_df['DT'] <= after_2days_dt)]
    recentbl_df = baseline_df.iloc[-1:,:].copy()

    id_cr_df = pd.concat([recentbl_df, afterdose_df])
    for cr_inx, cr_row in id_cr_df.iterrows():

        bl_dt = cr_row['DT']
        bl_cr = cr_row['RSLT']
        aki_48h_dt = dtstr_calculation(dtstr=cr_row['DT'], days=2)
        aki_7days_dt = dtstr_calculation(dtstr=cr_row['DT'], days=7)

        # (1) 48h 이내에 Cr 0.3 증가
        # (2) 7일 이내에 기저치 * 1.5 이상으로 상승

        if (((id_cr_df['DT'] < aki_48h_dt) & (id_cr_df['DT'] != cr_row['DT']))*1).sum() > 0:
            cond_exist_df = id_cr_df[(id_cr_df['DT'] < aki_48h_dt) & (id_cr_df['DT'] != cr_row['DT'])]
            aki_exist_df = cond_exist_df[cond_exist_df['RSLT'] >= bl_cr+0.3]
            if len(aki_exist_df)>0:
                aki_exist_df = aki_exist_df.rename(columns={'DT':'AKI_DT','RSLT':'AKI_RSLT'})
                aki_exist_df['COND_TYPE'] = '48h'
                aki_exist_df['BLCr_DT'] = bl_dt
                aki_exist_df['BLCr_RSLT'] = bl_cr
                result_df.append(aki_exist_df[result_cols].copy())


        if (((id_cr_df['DT'] < aki_7days_dt) & (id_cr_df['DT'] != cr_row['DT'])) * 1).sum() > 0:
            cond_exist_df = id_cr_df[(id_cr_df['DT'] < aki_7days_dt) & (id_cr_df['DT'] != cr_row['DT'])]
            aki_exist_df = cond_exist_df[cond_exist_df['RSLT'] >= bl_cr*1.5]
            if len(aki_exist_df)>0:
                aki_exist_df = aki_exist_df.rename(columns={'DT':'AKI_DT','RSLT':'AKI_RSLT'})
                aki_exist_df['COND_TYPE'] = '7days'
                aki_exist_df['BLCr_DT'] = bl_dt
                aki_exist_df['BLCr_RSLT'] = bl_cr
                result_df.append(aki_exist_df[result_cols].copy())

result_df = pd.concat(result_df)
result_df['ID'].drop_duplicates(ignore_index=True)

result_df.to_csv('./prepdata/amikacin_aki.csv', index=False, encoding='utf-8')


filt_res_df = result_df[result_df['BLCr_RSLT']<=1.2]
filt_res_df.to_csv('./prepdata/amikacin_aki(filtered).csv', index=False, encoding='utf-8')

    #(datetime.strptime(row['SDOSE_DT'], '%Y-%m-%dT%H:%M') - timedelta(days=30)).strftime('%Y-%m-%dT%H:%M')
filt_res_df['ID'].drop_duplicates(ignore_index=True)

cycle_df = pd.DataFrame(cycle_df)

# cycle_df.sort_values(['ID'])
# lab_df.drop_duplicates(['ID']).sort_values('ID')
# 148484876560382
# 148484930717343