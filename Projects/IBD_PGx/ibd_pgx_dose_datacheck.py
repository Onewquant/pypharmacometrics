from tools import *
from pynca.tools import *

prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx'
outcome_dir = f"C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx/results"

info_df = pd.read_csv(f"{prj_dir}/resource/anti_TNFa_patients.csv").sort_values(['EMR ID'])
info_df = info_df.rename(columns={'EMR ID':'ID','name':'NAME'})[['ID','NAME']].copy()

dose_df = pd.read_csv(f"{outcome_dir}/dose_df.csv")
u_dose_df = (dose_df.groupby(['ID'], as_index=False).agg({'DRUG': set})).rename(columns={'DRUG':'U_DRUG'})
# u_dose_df['U_DRUG'] = u_dose_df['U_DRUG'].replace(np.nan, '')
u_dose_df[['ID','U_DRUG']].drop_duplicates(['ID'], ignore_index=True)
# u_dose_df.to_csv(f"{outcome_dir}/dose_check_df.csv", encoding='utf-8-sig', index=False)


merge_df = info_df.merge(u_dose_df, on=['ID'],how='left')
# merge_df.to_csv(f"{outcome_dir}/dose_check_df.csv", encoding='utf-8-sig', index=False)
# merge_df
# merge_df = merge_df.rename(columns={'DRUG':'DOSE_EXISTANCE'})
merge_df['U_DRUG'] = merge_df['U_DRUG'].replace(np.nan,'')
merge_df['EXISTANCE'] = merge_df['U_DRUG'].map(lambda x:'O' if len(x)>0 else 'X')
merge_df[['ID','NAME','EXISTANCE','U_DRUG']].to_csv(f"{outcome_dir}/dose_check_df.csv", index=False, encoding='utf-8-sig')