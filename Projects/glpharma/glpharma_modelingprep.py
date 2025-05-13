from tools import *

## 기본정보입력

project_dict = {'GLPHARMA':['W2406']}
modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/glpharma/resource"
prepconc_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/glpharma/resource/prep_data"
if not os.path.exists(prepconc_dir_path):
    os.mkdir(prepconc_dir_path)

## Modeling and Dosing Policy 파일 불러오기

mdpolicy_file_path = f'{modeling_dir_path}/Modeling_and_Dosing_Policy - MDP_GLPHARMA_W2406.csv'
dspol_df = modeling_dosing_policy(mdpolicy_file_path, selected_models=[], model_colname='MODEL')

## Project, Drug 별로 Conc data 모으기

drugconc_dict=get_drug_conc_data_dict_of_multiple_projects(project_dict, prepconc_dir_path=prepconc_dir_path, conc_filename_format="[project]_ConcPrep_[drug]_R.csv")

## NONMEM 형식으로 Dataprep

uid_cols = ['DRUG']
covar_cols=['DRUG']
add_covar_df=pd.DataFrame(columns=['UID'])
uid_on=True
term_dict={'TIME': 'ATIME', 'TAD': 'NTIME', 'DV': 'CONC', 'ID': 'ID'}
formatting_data_nca_to_nonmem(drugconc_dict=drugconc_dict, dspol_df=dspol_df, uid_cols=uid_cols, modeling_dir_path=modeling_dir_path, covar_cols=covar_cols, add_covar_df=add_covar_df, uid_on=True)


modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/glpharma/resource"
df = pd.read_csv(f"{modeling_dir_path}/modeling_prep_data/MDP_C2A1E1_W2406.csv")
df['']

