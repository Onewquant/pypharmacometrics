from tools import *
from pynca.tools import *

prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx'

info_df = pd.read_csv(f"{prj_dir}/resource/anti_TNFa_patients.csv").sort_values(['EMR ID'])
info_df = info_df.rename(columns={'EMR ID':'ID','name':'NAME'})[['ID','NAME']].copy()


dose_df = pd.read_csv(f"{prj_dir}/resource/results/dose_df.csv")

merge_df = info_df.merge(dose_df.drop_duplicates(['ID'])[['ID','DRUG']], on=['ID'],how='left')
merge_df = merge_df.rename(columns={'DRUG':'DOSE_EXISTANCE'})
merge_df['DOSE_EXISTANCE'] = merge_df['DOSE_EXISTANCE'].replace(np.nan,'X').map({'adalimumab':'O','infliximab':'O'})
merge_df.to_csv(f"{prj_dir}/resource/dose_check_df.csv", index=False, encoding='utf-8-sig')