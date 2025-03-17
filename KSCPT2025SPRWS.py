from tools import *

# df = pd.read_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/prep_data/KSCPT2025SPRWS_ConcPrep_SGLT2INH(R).csv")
# df[df['GRP'].isin([1,3])].reset_index(drop=True).to_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/prep_data/KSCPTSPRWS25MULTI_ConcPrep_SGLT2INH(R).csv", index=False, encoding='utf-8-sig')
# df[df['GRP'].isin([2,4])].reset_index(drop=True).to_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/prep_data/KSCPTSPRWS25SINGLE_ConcPrep_SGLT2INH(R).csv", index=False, encoding='utf-8-sig')

## 기본정보입력

project_dict = {'KSCPTSPRWS25MULTI':['SGLT2INH'], 'KSCPTSPRWS25SINGLE':['SGLT2INH']}
modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25"
prepconc_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/prep_data"

## Modeling and Dosing Policy 파일 불러오기

mdpolicy_file_path = f'{modeling_dir_path}/Modeling_and_Dosing_Policy - MDP.csv'
dspol_df = modeling_dosing_policy(mdpolicy_file_path, selected_models=[], model_colname='MODEL')

## Project, Drug 별로 Conc data 모으기

drugconc_dict=get_drug_conc_data_dict_of_multiple_projects(project_dict, prepconc_dir_path=prepconc_dir_path, conc_filename_format="[project]_ConcPrep_[drug](R).csv")
# drugconc_dict['SGLT2INH'].columns

## NONMEM 형식으로 Dataprep

uid_cols = []
covar_cols=['GRP','SITENM', 'AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'GFR', 'TBIL']
add_covar_df=pd.DataFrame(columns=['UID'])
uid_on=True
term_dict={'TIME': 'ATIME', 'TAD': 'NTIME', 'DV': 'CONC', 'ID': 'ID'}
formatting_data_nca_to_nonmem(drugconc_dict=drugconc_dict, dspol_df=dspol_df, uid_cols=uid_cols, modeling_dir_path=modeling_dir_path, covar_cols=covar_cols, add_covar_df=pd.DataFrame(columns=['UID']), uid_on=True, term_dict=term_dict)
# formatting_data_nca_to_nonmem(drugconc_dict=drugconc_dict, dspol_df=dspol_df[dspol_df['ABS_POLICY'].map(lambda x:x[0])=='E'], uid_cols=uid_cols, modeling_dir_path=modeling_dir_path, covar_cols=covar_cols, add_covar_df=pd.DataFrame(columns=['UID']), uid_on=True, term_dict=term_dict)


