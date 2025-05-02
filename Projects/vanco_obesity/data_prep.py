from tools import *
from pynca.tools import *

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/vanco_obesity"

res_cols = ['ID','DATETIME','DV','MDV','CMT','AMT','RATE']

dose_df = pd.read_excel(f"{prj_dir}/resource/Dose(2020-2024).xlsx")
conc_df = pd.read_excel(f"{prj_dir}/resource/VCMconc(2020-2024).xlsx")
# lab_df1 = pd.read_excel(f"{prj_dir}/resource/Lab(2020-2024).xlsx", sheet_name="약제")
# lab_df2 = pd.read_excel(f"{prj_dir}/resource/Lab(2020-2024).xlsx", sheet_name="진검약리")
# lab_df = pd.concat([lab_df1, lab_df2], ignore_index=True)


# dose_df.columns
dose_df = dose_df[['등록번호','약품 오더일자','[함량단위환산] 1회 처방량','[실처방] 경로','수행시간']].copy()
dose_df = dose_df[dose_df['[실처방] 경로'].isin(['IV','MIV','IVS'])].copy()
dose_df = dose_df[~dose_df['수행시간'].isna()].copy()
dose_df['수행시간'] = dose_df['수행시간'].map(lambda x:x.split(', '))
dose_df = dose_df.explode('수행시간')
dose_df['약품 오더일자'] = dose_df['약품 오더일자'].astype(str)
dose_df['수행시간'] = dose_df.apply(lambda x:x['약품 오더일자']+' '+x['수행시간'] if len(x['수행시간'].split(' '))<2 else x['수행시간'], axis=1)

dose_df = dose_df.rename(columns={'등록번호':'ID','수행시간':'DATETIME','[함량단위환산] 1회 처방량':'AMT'})
dose_df = dose_df[['ID','DATETIME','AMT']].copy()
# dose_df.to_csv(f"{prj_dir}/resource/dose_prep.csv", index=False)
dose_df = dose_df[dose_df['DATETIME'].map(lambda x:x.split('/')[-1]).isin(['Y','DR','Z','O'])].copy()
dose_df['DATETIME'] = dose_df['DATETIME'].map(lambda x:x.split('/')[0])
dose_df = dose_df.sort_values(['ID','DATETIME'])
dose_df['DV'] = '.'
dose_df['MDV'] = 1
dose_df['AMT'] = dose_df['AMT'].map(lambda x: x.replace('mg',''))
dose_df['RATE'] = dose_df['AMT']
dose_df['CMT'] = 1
dose_df = dose_df[res_cols].copy()

conc_df = conc_df[['등록번호','채혈일시','검사결과']].copy()
conc_df = conc_df.rename(columns={'등록번호':'ID', '채혈일시':'DATETIME','검사결과':'DV'})
conc_df = conc_df.sort_values(['ID','DATETIME'])
conc_df['DATETIME'] = conc_df['DATETIME'].astype(str).map(lambda x:x[:-3])
conc_df['MDV'] = '.'
conc_df['AMT'] = '.'
conc_df['RATE'] = '.'
conc_df['CMT'] = 1
conc_df = conc_df[res_cols].copy()

md_df = pd.concat([dose_df, conc_df], ignore_index=True)
md_df = md_df.sort_values(['ID','DATETIME'], ignore_index=True)

