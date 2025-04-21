import numpy as np
import pandas as pd

med_df = pd.read_csv('Projects/amikacin_tdm/resource/dt_9083_약품_50207.csv')
conc_df = pd.read_csv('Projects/amikacin_tdm/resource/dt_9084_검사_1090473.csv')

## 약품 오더데이터 편집
result_cols = ['ID', 'DATE', 'DOSE','ACTION']
med_df = med_df[~med_df['간호사 수행 여부'].isna()].copy()
med_df['ID'] = med_df['환자번호'].map(lambda x:x.split('-')[0])
med_df['DATE'] = med_df['약품 오더일자']
med_df['DOSE'] = med_df['[함량단위환산] 1회 처방량'].map(lambda x:float(x.replace('mg','').split('-')[0]))
med_df['ACTION'] = med_df['수행시간'].map(lambda x:x.split(', '))

med_df = med_df[result_cols].explode(column='ACTION').copy()
med_df['ACTION'] = med_df.apply(lambda x:x['ACTION'] if len(x['ACTION'].split(' '))>1 else f"{x['DATE']} {x['ACTION']}", axis=1)
med_df['ACTION'] = med_df['ACTION'].map(lambda x:x.split('/')[0].replace(' ','T'))
med_df['DT'] = med_df['ACTION'].copy()

med_df = med_df[['ID','DT','DOSE']].sort_values(['ID','DT'])
med_df = med_df.groupby(by=['ID','DT']).agg({'DOSE':'sum'})
med_df = med_df.reset_index(drop=False)
med_df.to_csv('Projects/amikacin_tdm/resource/amikacin_dose.csv', index=False)

##



# for inx, df_frag in med_df.groupby(by=['ID','DT']):
#     if len(df_frag)>1:
#         print('stop')
#         break
# med_df



# med_df[med_df['ID']=='148484876560382']
# med_df[med_df['[함량단위환산] 1회 처방량'].map(lambda x:x.replace('mg',''))=='430-430']['[함량단위환산] 1회 처방량']


