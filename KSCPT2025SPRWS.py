from tools import *

## 기본정보입력

project_dict = {'CKD379':['Empagliflozin','Sitagliptin','Metformin'], 'CKD383':['Empagliflozin','Lobeglitazone','Metformin']}
modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/metformin_xr"

## Modeling and Dosing Policy 파일 불러오기

mdpolicy_file_path = f'{modeling_dir_path}/Modeling_and_Dosing_Policy - MDP.csv'
dspol_df = modeling_dosing_policy(mdpolicy_file_path, selected_models=[], model_colname='MODEL')

## Project, Drug 별로 Conc data 모으기

drugconc_dict=get_drug_conc_data_dict_of_multiple_projects(project_dict, modeling_dir_path=modeling_dir_path, conc_filename_format="[project]_ConcPrep_[drug](R).csv")

## NONMEM 형식으로 Dataprep

uid_inf_cols = ['PROJECT','PERIOD','SEQUENCE']
formatting_data_nca_to_nonmem(drugconc_dict=drugconc_dict, dspol_df=dspol_df, uid_cols=uid_inf_cols, modeling_dir_path=modeling_dir_path, covar_cols=['FEEDING'], add_covar_df=pd.DataFrame(columns=['UID']), uid_on=True)


