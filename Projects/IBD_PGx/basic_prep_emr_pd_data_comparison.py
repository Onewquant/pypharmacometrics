from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

modeling_df = pd.read_csv(f"{output_dir}/modeling_df_covar/infliximab_integrated_datacheck(covar).csv")
dosing_df = modeling_df[modeling_df['MDV']==1].copy()
start_dosing_df = dosing_df.groupby('UID', as_index=False)['DATETIME'].min().rename(columns={'UID':'ID'})
cdai_df = pd.read_csv(f"{output_dir}/pdmarker_cdai_df.csv")
pms_df = pd.read_csv(f"{output_dir}/pdmarker_pms_df.csv")
# cdai_df.columns

clin_ind_df = pd.read_excel(f"{output_dir}/PKPD_EDA/GI_IM_Clinical Data/IBD WGS TDM_IFX induction TL patient list 250802.xlsx")
clin_maint_df = pd.read_excel(f"{output_dir}/PKPD_EDA/GI_IM_Clinical Data/IBD WGS TDM_IFX 1y TL_patient list 250802.xlsx")

# clin_ind_df.columns
# clin_maint_df.columns

rename_indcols_dict = {'EMR ID':'ID','name':'NAME','sex':'SEX','IBD type':'IBD_TYPE','age_at_dx':'AGE_AT_DX','anti_tnfa_type':'DRUG','anti_tnfa_date':'DRUG_START_DATE', '1_IFX TL_date':'CONC_DATE', '1_IFX TL_result':'CONC', 'IFX_IND_CD_diarrhea':'CDAI_DIARHEACNT', 'IFX_IND_CD_abd':'CDAI_ABDPAINCNT', 'IFX_IND_CD_PRO2':'CD_PRO2', 'IFX_IND_UC_SF':'PMS_STLCNT', 'IFX_IND_UC_RB':'PMS_HEMATOCHEZIA', 'IFX_IND_UC_PRO2':'UC_PRO2', 'IFX_IND_UC_MES':'PMS_TOTALSCORE', 'IFX_IND_UC_modifiedMayo':'MMS_TOTALSCORE', 'IFX_IND_CRP':'PD_CRP', 'IFX_IND_FCP':'PD_FCAL', '1ADA_date':'ADA_DATE', '1ADA_result':'ADA'}
clin_ind_df = clin_ind_df.rename(columns=rename_indcols_dict)
clin_ind_df = clin_ind_df[list(rename_indcols_dict.values())].copy()
clin_ind_df['NAME'] = clin_ind_df['NAME'].map(lambda x:x.split(':')[-1].replace(')',""))
drug_sdate_ind_df = clin_ind_df[['ID','NAME','DRUG_START_DATE']].copy()

ind_comp_df = drug_sdate_ind_df.merge(start_dosing_df, on=['ID'], how='left')
ind_comp_df[ind_comp_df['ID']!=36898756].to_csv(f'{resource_dir}/initial_induction_patients.csv', index=False, encoding='utf-8-sig')
# ind_diff_df = ind_comp_df[ind_comp_df['DRUG_START_DATE'] < ind_comp_df['DATETIME']].copy()
ind_diff_df = ind_comp_df[ind_comp_df['DRUG_START_DATE'] != ind_comp_df['DATETIME']].copy()

rename_maintcols_dict = {'EMR ID':'ID','name':'NAME','sex':'SEX','IBD type':'IBD_TYPE','age_at_dx':'AGE_AT_DX','anti_tnfa_type':'DRUG','anti_tnfa_date':'DRUG_START_DATE', '1Y_IFX TL_date':'CONC_DATE', '1Y IFX TL':'CONC', 'IFX_1y_CD_diarrhea':'CDAI_DIARHEACNT', 'IFX_1y_CD_abd':'CDAI_ABDPAINCNT', 'IFX_1y_CD_PRO2':'CD_PRO2', 'IFX_1y_SF':'PMS_STLCNT', 'IFX_1y_RB':'PMS_HEMATOCHEZIA', 'IFX_1y_PRO2':'UC_PRO2', 'IFX_1y_MES':'PMS_TOTALSCORE', 'IFX_1y_modifiedMayo':'MMS_TOTALSCORE', 'IFX_1y_CRP':'PD_CRP', 'IFX_1y_FCP':'PD_FCAL', '1Y ADA_date':'ADA_DATE', '1Y ADA':'ADA'}
clin_maint_df = clin_maint_df.rename(columns=rename_maintcols_dict)
clin_maint_df = clin_maint_df[list(rename_maintcols_dict.values())].copy()
clin_maint_df['NAME'] = clin_maint_df['NAME'].map(lambda x:x.split(':')[-1].replace(')',""))
drug_sdate_maint_df = clin_maint_df[['ID','NAME','DRUG_START_DATE']].copy()

maint_comp_df = drug_sdate_maint_df.merge(start_dosing_df, on=['ID'], how='left')
# maint_diff_df = maint_comp_df[maint_comp_df['DRUG_START_DATE'] < maint_comp_df['DATETIME']].copy()
maint_diff_df = maint_comp_df[maint_comp_df['DRUG_START_DATE'] != maint_comp_df['DATETIME']].copy()

integrated_diff_df = pd.concat([ind_diff_df, maint_diff_df]).drop_duplicates(['ID'], ignore_index=True)

# drug_sdate_df = pd.concat([drug_sdate_ind_df, drug_sdate_maint_df])

# clin_maint_df.columns

# len(clin_ind_df)
# len(clin_maint_df)

# Induction에 해당하는 사람 수 세기
ind_list = list()
first_dv_df = modeling_df[(modeling_df['MDV']==0)].copy()
first_dv_df['DV'] = first_dv_df['DV'].map(float)
first_dv_df = first_dv_df[(first_dv_df['DV']!=0)].drop_duplicates(['UID'], ignore_index=True)
for inx, row in first_dv_df.iterrows(): #break
    ind_disc_df = modeling_df[(modeling_df['ID']==row['ID'])&(modeling_df['TIME'] < row['TIME'])&(modeling_df['MDV']!=0)].copy()
    if len(ind_disc_df) == 3:
        ind_list.append(row['ID'])
len(ind_list)

## Modeling Data
# INDUCTION:

## Clinical Data
# INDUCTION: 3번째 투약 (atfter 14 weeks) 후 TL 측정된 경우 -> Induction Phase 속함 -> 35명
# Maintenance: 첫 투약 이후 1년 기준 ± 6개월에 TL 측정된 경우 -> Maintenance Phase 속함 -> 53명
# Total: 57명

# drug_sdate_df.drop_duplicates(['ID'], ignore_index=True)






