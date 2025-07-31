# from emr_scraper.tools import *
import os
import pandas as pd
from datetime import datetime, timedelta

drug_name = 'AMK'
prj_dir = f"C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{drug_name}"
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

df = pd.read_csv(f"{resource_dir}/AMK_tdm.csv",encoding='euc-kr')
df = df.rename(columns={'등록번호':'EMR ID', '검사일':'TDM_DATE','환자명':'name', '담당부서':'DEP'})
df = df.drop_duplicates(['EMR ID'], ignore_index=True)


req_df = list()
for inx, row in df.iterrows():

    row['name'] = row['name'].split('(')[0]
    pid = row['EMR ID']

    tdm_dt = datetime.strptime(row['TDM_DATE'], '%Y-%m-%d %H:%M')
    before_date = (tdm_dt - timedelta(days=60)).strftime('%Y-%m-%d')
    after_date = (tdm_dt + timedelta(days=100)).strftime('%Y-%m-%d')

    req_df.append({'EMR ID': pid, 'NAME': row['name'], 'START_DATE':before_date, 'END_DATE':after_date})

req_df = pd.DataFrame(req_df)
req_df.to_csv(f"{output_dir}/req_df_for_conc.csv", encoding='utf-8-sig', index=False)