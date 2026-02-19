import os
import glob
import numpy as np
import pandas as pd
# from exam_VS import count

prj_name = 'VPA'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'

resource_dir = f'{prj_dir}/resource/ysy_req'
output_dir = f"{prj_dir}/results"


uid_df_dict = {}

for i in range(1,5):
    res_filename = f'251202_재식별 파일_DILI_{i}.csv'
    uid_df_dict[i] = pd.read_csv(f'{resource_dir}/re_identification/{res_filename}')
    uid_df_dict[i] = uid_df_dict[i].rename(columns={'Deidentification_ID':'UID','환자번호':'ID'})
    uid_df_dict[i]['ID'] = uid_df_dict[i]['ID'].map(lambda x:x.split('-')[0])
    uid_df_dict[i] = uid_df_dict[i][['ID','UID']].copy()

lab_files = glob.glob(f"{resource_dir}/LAB_exam/VPA/VPA_LAB_*.csv")
for finx, lab_file in enumerate(lab_files):#break
    # lab_name = 'Albumin'
    lab_df = pd.read_csv(lab_file)
    lab_df['ID'] = lab_df['ID'].astype(str)

    for i in range(1,5):
        # i = 2
        lab_df_test = lab_df.merge(uid_df_dict[i], on=['ID'], how='left')
        uid_pct = ((~(lab_df_test['UID'].isna())).sum()) / len(lab_df_test)
        print(uid_pct)
        if uid_pct > 0.7:
            # raise ValueError
            # lab_df_test
            lab_df_test['ID'] = lab_df_test['UID']
            lab_df_test = lab_df_test.drop(['UID'], axis=1)
            lab_df_test.to_csv(lab_file, index=False, encoding='utf-8-sig')
            print(f"({finx}) 파일 치환 완료 | {lab_file}")
            break
        else:
            continue