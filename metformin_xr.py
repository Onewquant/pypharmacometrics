from tools import *

## 기본정보입력

project_dict = {'CKD379':['Empagliflozin','Sitagliptin','Metformin'], 'CKD383':['Empagliflozin','Lobeglitazone','Metformin']}
modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/metformin_xr"

## Modeling and Dosing Policy 파일 불러오기

selected_models = []
# selected_models = ['C1A0E1']
dspol_df = pd.read_csv(f'{modeling_dir_path}/Modeling_and_Dosing_Policy - MDP.csv')
if len(selected_models)!=0:
    dspol_df = dspol_df[dspol_df['MODEL'].isin(selected_models)].reset_index(drop=True)
dspol_df = dspol_df.applymap(lambda x:convert_to_numeric_value(x))


## Project, Drug 별로 Conc data 모으기

project_drugconc_dict = dict()
drugconc_dict = dict()
for projectname, drug_list in project_dict.items():
    filename_format = f"{projectname}_ConcPrep_[drug](R).csv"

    ct_project_dir_path = modeling_dir_path + '/' + projectname

    project_drugconc_dict[projectname] = load_data_dict(drug_list, filename_format=filename_format, input_file_dir_path=ct_project_dir_path)

    ## Drug별 데이터 따로 모으기

    for drug, drugconc_df in project_drugconc_dict[projectname].items():
        drugconc_df['PROJECT'] = projectname
        if drug not in drugconc_dict:
            drugconc_dict[drug] = drugconc_df.copy()
        else:
            drugconc_dict[drug] = pd.concat([drugconc_dict[drug], drugconc_df], ignore_index=True)

        # drugconc_df = drugconc_df.drop(['PROJECT'],axis=1)
        # drugconc_df['SEQUENCE'] = drugconc_df['ID'].map(lambda x:x[0]).map({'A':1,'B':2})
        # drugconc_df.to_csv(f'{modeling_dir_path}/{projectname}/{projectname}_ConcPrep_{drug}(R).csv', index=False, encoding='utf-8-sig')

## NONMEM 형식으로 Dataprep

uid_inf_cols = ['PROJECT','PERIOD','SEQUENCE']


def generate_nonmem_subject_id(df, sid_col, uid_cols):
    ## ID Column

    dcdf = df[[sid_col] + uid_cols].copy()
    dcdf['NONMEM_ID'] = '1'  # 처음에 0이 안나오도록 설정: 사람이면 1, 동물이면 2로 하면 될듯

    # uid (Project 내의 고유 대상자 번호를 포함한 ID) 생성
    # uid에 구분할 수 있도록 포함시키려는 Column의 Value 값들에 숫자로 구분하여 배치표(dictionary)를 생성해두고, Nonmem에서 사용할 ID 생성
    uid_inf_dict = dict()
    for uicol in uid_cols:
        if len(dcdf[uicol].unique()) > 10:
            raise ValueError(f"{uicol} 컬럼의 Unique Value가 10개 이상입니다. 새로운 구분법이 필요합니다.")
        else:
            if str in dcdf[uicol].map(lambda x: type(x)).unique():
                uid_inf_dict[uicol] = {c: i for i, c in enumerate(dcdf[uicol].unique())}
            else:
                uid_inf_dict[uicol] = {c: c for c in dcdf[uicol].unique()}

        dcdf['NONMEM_ID'] += dcdf[uicol].map(lambda x: str(uid_inf_dict[uicol][x]))

    dcdf['NONMEM_ID'] += dcdf[sid_col]
    # dcdf['NONMEM_ID'] = dcdf['NONMEM_ID'].map(int)
    # dcdf['NONMEM_UID'] = dcdf['NONMEM_ID'].copy()
    return dcdf['NONMEM_ID'].copy()

def formatting_data_nca_to_nonmem(drugconc_dict, dspol_df, uid_cols, modeling_dir_path, covar_cols=[], add_covar_df=pd.DataFrame(columns=['UID']), uid_on=True, term_dict={'TIME':'ATIME','TAD':'NTIME','DV':'CONC','ID':'ID'}):

    """
    # add_covar_df 를 추가할때는 uid_cols를 반드시 포함하는 primary key로서 사용
    """

    unique_models = dspol_df['MODEL'].unique()
    for model_name in unique_models: #print(model_name)
        for drug in drugconc_dict.keys():

            model_drug_dspol = dspol_df[(dspol_df['MODEL'] == model_name) & (dspol_df['DRUG'] == drug)].copy()
            # if drug=='Metformin': break
            dcdf = drugconc_dict[drug].copy()

            ## conc data에 이미 covar이 존재하는 경우 covar col 남김

            drug_covar_cols = list()
            for covar_c in covar_cols:
                if covar_c in dcdf.columns:
                    dcdf[f'NONMEM_{covar_c}'] = dcdf[covar_c].copy()
                    drug_covar_cols.append(covar_c)
                else:
                    print(f"{model_name} / {drug} / Raw data에 {covar_c}라는 컬럼이 존재하지 않습니다.")
                    continue

            ## Result columns 설정

            result_cols = ['ID'] + ['UID'] if uid_on else [] + ['TIME','TAD','DV','MDV','AMT','RATE','DUR','CMT'] + drug_covar_cols

            ## add_covar_df 추가하는 경우 (uid_cols 제외하고)

            add_covar_df_copy = add_covar_df.copy()
            new_add_covar_df_cols = list()
            for add_var_df_col in add_covar_df_copy.columns:
                if add_var_df_col in uid_cols:
                    new_add_covar_df_cols.append(add_var_df_col)
                else:
                    new_add_covar_df_cols.append(f"NONMEM_{add_var_df_col}")
            add_covar_df_copy.columns = new_add_covar_df_cols

            dcdf = dcdf.merge(add_covar_df_copy, on=uid_cols, how='left')

            ## nonmem_subject_id 생성

            dcdf['NONMEM_ID'] = generate_nonmem_subject_id(df=dcdf, sid_col=term_dict['ID'], uid_cols=uid_cols)
            if uid_on:
                dcdf['NONMEM_UID'] = dcdf['NONMEM_ID'].copy()

            # ## ID Column
            #
            # dcdf['NONMEM_ID']='1'  # 처음에 0이 안나오도록 설정: 사람이면 1, 동물이면 2로 하면 될듯
            #
            # # uid (Project 내의 고유 대상자 번호를 포함한 ID) 생성
            # # uid에 구분할 수 있도록 포함시키려는 Column의 Value 값들에 숫자로 구분하여 배치표(dictionary)를 생성해두고, Nonmem에서 사용할 ID 생서
            # for uicol in uid_inf_cols:
            #     if len(dcdf[uicol].unique())>10:
            #         raise ValueError(f"{uicol} 컬럼의 Unique Value가 10개 이상입니다. 새로운 구분법이 필요합니다.")
            #     else:
            #         if str in dcdf[uicol].map(lambda x:type(x)).unique():
            #             uid_inf_dict[uicol] = {c:i for i, c in enumerate(dcdf[uicol].unique())}
            #         else:
            #             uid_inf_dict[uicol] = {c:c for c in dcdf[uicol].unique()}
            #
            #     dcdf['NONMEM_ID']+=dcdf[uicol].map(lambda x:str(uid_inf_dict[uicol][x]))
            # dcdf['NONMEM_ID']+=dcdf['ID'].map(lambda x:re.findall(r'[\d]+', x)[0])
            # dcdf['NONMEM_ID']=dcdf['NONMEM_ID'].map(int)
            # dcdf['NONMEM_UID']=dcdf['NONMEM_ID'].copy()

            ## TIME Columns

            dcdf['NONMEM_TIME'] = dcdf[term_dict['TIME']].copy()
            dcdf['NONMEM_TAD'] = dcdf[term_dict['TAD']].copy()

            ## DV, MDV Column

            dcdf['NONMEM_DV'] = dcdf[term_dict['DV']].copy()
            dcdf['NONMEM_MDV'] = 0

            ## AMT, RATE, DUR

            dcdf['NONMEM_AMT'] = '.'

            ## CMT / DUR / RATE

            dcdf['NONMEM_RATE'] = '.'
            dcdf['NONMEM_DUR'] = '.'
            dcdf['NONMEM_CMT'] = sampling_cmt_for_specific_advan_type(advan=model_drug_dspol['ADVAN'].iloc[0], forced_sampling_cmt=np.nan)

            """
            # 주로 농도측정이 Systemic compartment에서 진행되므로 2로 표시하는게 나을듯. 나중에 수정필요할수도
            # dosing policy 파일에 이 정보가 반영되도록 만들어야. 
            # (dosing policy / modeling policy 가 모두 영향을 주는 듯)
            """

            # Dosing row 추가 : Project / Drug 종류에 따라 Dosing Policy 다를 수 있다고 가정. (Drug은 위에서 dcdf로 이미 구분되어있음)

            drug_nmprep_df = list()
            for projectname, prj_dcdf in dcdf.groupby(['PROJECT']):

                # 해당 dosing policy
                prj_dspol = model_drug_dspol[(dspol_df['PROJECT'] == projectname)].reset_index(drop=True)

                for dspol_inx, dspol_row in prj_dspol.iterrows(): #break

                    for nmid, nmid_df in prj_dcdf.groupby(['NONMEM_ID']): #break

                        nmid_df = nmid_df.reset_index(drop=True)

                        # 추가할 Dosing row

                        add_row = nmid_df.iloc[0:1, :].copy()
                        add_row['NONMEM_TIME'] = dspol_row['RELTIME']
                        add_row['NONMEM_TAD'] = dspol_row['RELTIME']
                        add_row['NONMEM_DV'] = '.'
                        add_row['NONMEM_MDV'] = 1
                        add_row['NONMEM_CMT'] = dosing_cmt_for_advan_type(advan=dspol_row['ADVAN'], route=dspol_row['ROUTE'], forced_dosing_cmt=np.nan)
                        add_row['NONMEM_AMT'] = dspol_row['DOSE']

                        # add_row.iloc[0]
                        """
                        # RATE=-1 : Rate 추정 코드 추가해야
                        # RATE=-2 : Dur 추정 코드 추가해야
                        
                        # [ADVAN 1 - 1구획 iv] : 1CMT (TRANS1 - K / TRANS2 - CL, V)
                        # [ADVAN 2 - 1구획 PO] : 1CMT-depot, 2CMT-central (TRANS1 - K, KA / TRANS2 - CL, V, KA)
                        # [ADVAN 3 - 2구획 iv] : 1CMT-central, 2CMT-peripheral (TRANS1 - K, K12, K21 / TRANS3 - CL, V, Q, VSS / TRANS4 - CL, V1, Q, V2 / TRANS5 - AOB, ALPHA, BETA / TRANS6 - ALPHA, BETA, K21)
                        # [ADVAN 4 - 2구획 PO] : 1CMT-depot, 2CMT-central, 3CMT-peripheral (TRANS1 - K, K23, K32, KA / TRANS3 - CL, V, Q, VSS, KA / TRANS4 - CL, V2, Q, V3, KA / TRANS5 - AOB, ALPHA, BETA, KA / TRANS6 - ALPHA, BETA, K32, KA)
                        # [ADVAN 11 - 3구획 iv] : 1CMT-central, 2CMT-peripheral1, 3CMT-peripheral2
                        # [ADVAN 12 - 3구획 PO] : 1CMT-depot, 2CMT-central, 3CMT-peripheral1, 4CMT-peripheral2
                        """

                        """
                        [ADVAN1 TRANS2 (1-comp IV세팅, dosing : CMT 1)] 
                        - AMT(O), RATE( _ or point ), DUR( _ or point ) => 안 돌아감
                        - AMT(O), RATE(-2), DUR(point) => 돌아감
                        
                        [ADVAN2 TRANS2 (1-comp PO세팅, dosing : CMT 2)] 
                        - AMT(O), RATE( _ or point ), DUR( _ or point ) => 돌아감
                        """

                        add_row['NONMEM_RATE'] = dspol_row['RATE']
                        add_row['NONMEM_DUR'] = dspol_row['DUR']
                        add_row['NONMEM_CMT'] = dosing_cmt_for_advan_type(advan=dspol_row['ADVAN'],route=dspol_row['ROUTE'])

                        # DUR값을 음수인 경우 -> 에러처리
                        if (dspol_row['DUR'] not in ('.', 0)):
                            if dspol_row['DUR'] < 0:
                                raise ValueError(f"Dosing Policy에서 {dspol_row['ROUTE'].upper()} 투여인데 DUR 값이 음수로 기록됨.")
                            else:
                                pass

                        if (dspol_row['DOSE'] < 0): ### 여기
                            raise ValueError(f"DOSE가 음수로 기록되었습니다.")
                        else:
                            if (dspol_row['DOSE'] > 0) and (dspol_row['RATE'] in (0,'.')):
                                print('Bolus와 같이 추정합니다. / DUR이 필요한 경우일지도. / 아직 어떻게 진행될지 모름 -> 해보자')
                                add_row['NONMEM_RATE'] = dspol_row['RATE']
                                add_row['NONMEM_DUR'] = dspol_row['DUR']
                            else:
                                if (dspol_row['RATE']==-1):
                                    # print('Absorption의 Rate을 추정합니다. / NONMEM 코드에 Rn = theta(m) 구문을 추가하세요')
                                    add_row['NONMEM_RATE'] = -1
                                    add_row['NONMEM_DUR'] = '.'
                                elif (dspol_row['RATE']==-2):
                                    # print('Absorption의 Duration을 추정합니다. / NONMEM 코드에 Dn = theta(m) 구문을 추가하세요')
                                    add_row['NONMEM_DUR'] = '.'
                                elif (dspol_row['RATE']>0):
                                    # print('양수인 RATE을 가지고 있습니다. / 코드에 어떻게 추가해야하는지 아직 모름 -> 해보자')
                                    add_row['NONMEM_DUR'] = '.'
                                else:
                                    raise ValueError(f"RATE, DURATION 값을 확인하세요.")

                        # 추가할 Dosing row를 위치할 곳 설정

                        if dspol_row['RELPOSITION'] not in (0, 1):
                            raise ValueError(f"Dosing Policy에서 RELPOSITION 은 0 또는 1이어야함. 해당 row 기준으로 (0: 직전 / 1: 직후)에 추가 row 삽입")

                        rt_inx = nmid_df[nmid_df[f"NONMEM_{dspol_row['RELTIMECOL']}"] == dspol_row['RELTIME']].iloc[0].name
                        add_pos_inx = rt_inx + dspol_row['RELPOSITION']

                        if add_pos_inx < 0:
                            raise ValueError(f"Dosing row를 추가하려는 위치가 0번째 row보다 작은값입니다. / {projectname}, {drug}, {nmid}")
                        elif add_pos_inx == 0:
                            nmid_df = pd.concat([add_row, nmid_df], ignore_index=True)
                        elif (add_pos_inx > 0) and (add_pos_inx < len(nmid_df)):
                            nmid_df = pd.concat(
                                [nmid_df.loc[:add_pos_inx - 1, :].copy(), add_row, nmid_df.loc[add_pos_inx:, :].copy()],
                                ignore_index=True)
                        elif (add_pos_inx == len(nmid_df)):
                            nmid_df = pd.concat([nmid_df, add_row], ignore_index=True)
                        else:
                            raise ValueError(f"Dosing row를 추가하려는 위치가 마지막 row보다 큰 값입니다. / {projectname}, {drug}, {nmid}")

                        drug_nmprep_df.append(nmid_df)


            # Prep Data 생성

            drug_nmprep_df = pd.concat(drug_nmprep_df, ignore_index=True)

            # Prep Data 컬럼편집

            drug_nonmem_cols = [c for c in drug_nmprep_df.columns if c.split('_')[0]=='NONMEM']
            drug_nmprep_df=drug_nmprep_df[drug_nonmem_cols].copy()
            drug_nmprep_df.columns=[c[7:] for c in drug_nonmem_cols]

            # ID 값을 작은 값으로 재조정 및 dictionary 생성

            short_id_dict = {id:short_id+1 for short_id, id in enumerate(drug_nmprep_df['ID'].unique())}
            drug_nmprep_df['ID'] = drug_nmprep_df['ID'].map(short_id_dict)

            # Prep Data 에 Covariate정보 추가 (추후 구현)

            # Prep Data 저장

            if not os.path.exists(f"{modeling_dir_path}/prep_data"):
                os.mkdir(f"{modeling_dir_path}/prep_data")
            drug_nmprep_df[result_cols].to_csv(f"{modeling_dir_path}/prep_data/MDP_{model_name}_{drug}.csv", index=False, encoding='utf-8-sig')

            print(f"{model_name} / {drug} / data-formatting completed")

            # CMT, RATE, DUR, EVID, SS, ADDL

            ## Covariates



formatting_data_nca_to_nonmem(drugconc_dict, dspol_df, uid_cols, modeling_dir_path=, covar_cols=[], add_covar_df=pd.DataFrame(columns=['UID']), uid_on=True)


