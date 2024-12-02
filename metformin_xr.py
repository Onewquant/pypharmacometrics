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
unique_models = dspol_df['MODEL'].unique()

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
uid_inf_dict = dict()
for model_name in unique_models:
    for drug in drugconc_dict.keys():
        dcdf = drugconc_dict[drug].copy()

        ## ID Column

        dcdf['NONMEM_ID']='1'  # 처음에 0이 안나오도록 설정: 사람이면 1, 동물이면 2로 하면 될듯

        # uid (Project 내의 고유 대상자 번호를 포함한 ID) 생성
        # uid에 구분할 수 있도록 포함시키려는 Column의 Value 값들에 숫자로 구분하여 배치표(dictionary)를 생성해두고, Nonmem에서 사용할 ID 생서
        for uicol in uid_inf_cols:
            if len(dcdf[uicol].unique())>10:
                raise ValueError(f"{uicol} 컬럼의 Unique Value가 10개 이상입니다. 새로운 구분법이 필요합니다.")
            else:
                if str in dcdf[uicol].map(lambda x:type(x)).unique():
                    uid_inf_dict[uicol] = {c:i for i, c in enumerate(dcdf[uicol].unique())}
                else:
                    uid_inf_dict[uicol] = {c:c for c in dcdf[uicol].unique()}

            dcdf['NONMEM_ID']+=dcdf[uicol].map(lambda x:str(uid_inf_dict[uicol][x]))
        dcdf['NONMEM_ID']+=dcdf['ID'].map(lambda x:re.findall(r'[\d]+', x)[0])
        dcdf['NONMEM_ID']=dcdf['NONMEM_ID'].map(int)

        ## TIME Columns

        dcdf['NONMEM_TIME']=dcdf['ATIME'].copy()
        dcdf['NONMEM_TAD'] = dcdf['NTIME'].copy()


        ## DV, MDV Column

        dcdf['NONMEM_DV'] = dcdf['CONC'].copy()
        dcdf['NONMEM_MDV'] = 0

        ## AMT, RATE, DUR

        dcdf['NONMEM_AMT'] = '.'

        ## CMT / DUR / RATE

        dcdf['NONMEM_RATE'] = '.'
        dcdf['NONMEM_DUR'] = '.'
        dcdf['NONMEM_CMT'] = 0

        """
        # 주로 농도측정이 Systemic compartment에서 진행되므로 2로 표시하는게 나을듯. 나중에 수정필요할수도
        # dosing policy 파일에 이 정보가 반영되도록 만들어야. 
        # (dosing policy / modeling policy 가 모두 영향을 주는 듯)
        """

        # Dosing row 추가 : Project / Drug 종류에 따라 Dosing Policy 다를 수 있다고 가정. (Drug은 위에서 dcdf로 이미 구분되어있음)

        drug_nmprep_df = list()
        for projectlist, prj_dcdf in dcdf.groupby(['PROJECT']):
            projectname = projectlist[0]
            prj_dspol = dspol_df[(dspol_df['MODEL'] == model_name) & (dspol_df['DRUG'] == drug) & (dspol_df['PROJECT'] == projectname)].reset_index(drop=True)

            for dspol_inx, dspol_row in prj_dspol.iterrows():

                for nmid, nmid_df in dcdf.groupby(['NONMEM_ID']):

                    nmid_df = nmid_df.reset_index(drop=True)

                    # 추가할 Dosing row

                    add_row = nmid_df.iloc[0:1, :].copy()
                    add_row['NONMEM_TIME'] = dspol_row['RELTIME']
                    add_row['NONMEM_TAD'] = dspol_row['RELTIME']
                    add_row['NONMEM_DV'] = '.'
                    add_row['NONMEM_MDV'] = 1
                    add_row['NONMEM_CMT'] = dspol_row['DOSING_CMT']  # Compartment 는 PO인 경우 CMT=1 / IV인 경우 CMT=2 로 설정 / 기타 IM 등은 모델에 따라 CMT 결정해서 Dosing_Policy에 정확히 기입요망
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

                    add_row['NONMEM_RATE'] = dspol_row['RATE']
                    add_row['NONMEM_DUR'] = dspol_row['DUR']
                    add_row['NONMEM_CMT'] = dosing_cmt_for_advan_type(advan=dspol_row['ADVAN'],route=dspol_row['ROUTE'])

                    # DUR값을 음수인 경우 -> 에러처리
                    if (dspol_row['DUR'] not in ('.', 0)):
                        if dspol_row['DUR'] < 0:
                            raise ValueError(f"Dosing Policy에서 {dspol_row['ROUTE'].upper()} 투여인데 DUR 값이 음수로 기록됨.")
                        else:
                            pass

                    if (['AMT'] < 0): ### 여기
                        raise ValueError(f"AMT가 음수로 기록되었습니다.")
                    else:
                        if (dspol_row['AMT'] > 0) and (dspol_row['RATE'] in (0,'.')):
                            print('Bolus와 같이 추정합니다. / DUR이 필요한 경우일지도. / 아직 어떻게 진행될지 모름 -> 해보자')
                            add_row['NONMEM_RATE'] = dspol_row['RATE']
                            add_row['NONMEM_DUR'] = dspol_row['DUR']
                        else:
                            if (dspol_row['RATE']==-1):
                                print('Absorption의 Rate을 추정합니다. / NONMEM 코드에 Rn = theta(m) 구문을 추가하세요')
                                add_row['NONMEM_RATE'] = -1
                                add_row['NONMEM_DUR'] = '.'
                            elif (dspol_row['RATE']==-2):
                                print('Absorption의 Duration을 추정합니다. / NONMEM 코드에 Dn = theta(m) 구문을 추가하세요')
                                add_row['NONMEM_DUR'] = '.'
                            elif (dspol_row['RATE']>0):
                                print('양수인 RATE을 가지고 있습니다. / 코드에 어떻게 추가해야하는지 아직 모름 -> 해보자')
                                add_row['NONMEM_DUR'] = '.'
                            else:
                                raise ValueError(f"RATE, DURATION 값을 확인하세요.")

                    # if (dspol_row['ROUTE'] != 'IV') and (dspol_row['ABS_POLICY'] == 1):
                    #     add_row['NONMEM_RATE'] = '.'
                    #     add_row['NONMEM_DUR'] = '.'
                    #     add_row['NONMEM_CMT'] = dosing_cmt_for_advan_type(advan=dspol_row['ADVAN'],route=dspol_row['ROUTE'])
                    # else:
                    #
                    #
                    #     # PO with zero-order absorption
                    #     if dspol_row['ROUTE'].upper() == 'PO':
                    #         if (dspol_row['ABS_ORD'] == 0):
                    #             if (dspol_row['RATE'] == -2):
                    #                 if dspol_row['DUR'] in (0, '.'):
                    #                     add_row['NONMEM_RATE'] = dspol_row['RATE']
                    #                     add_row['NONMEM_DUR'] = dspol_row['DUR']
                    #                     print(
                    #                         f'PO로 쓰여있으며 Constant Absorption 모델을 의도하였는데, 입력한 DUR은 bolus와 같습니다. / {projectname}, {drug}, {nmid}')
                    #                 else:
                    #                     add_row['NONMEM_RATE'] = dspol_row['RATE']
                    #                     add_row['NONMEM_DUR'] = dspol_row['DUR']
                    #             else:
                    #                 raise ValueError(
                    #                     f"Dosing Policy에서 {dspol_row['ROUTE'].upper()} 투여 및 Constant Absorption인데 RATE가 -2가 아닙니다.")
                    #         else:
                    #             print(
                    #                 f'PO로 쓰여있으며 0, 1차 absorption 이외의 다른 모델을 의도한 것 같습니다. / {projectname}, {drug}, {nmid}')
                    #
                    #     # iV 투여시
                    #     elif dspol_row['ROUTE'].upper() == 'IV_BOLUS':
                    #         if (dspol_row['RATE'] in (-1, -2)) and (dspol_row['DUR'] in ('.', 0)):
                    #             add_row['NONMEM_RATE'] = dspol_row['RATE']
                    #             add_row['NONMEM_DUR'] = dspol_row['DUR']
                    #         else:
                    #             raise ValueError(
                    #                 f"Dosing Policy에서 {dspol_row['ROUTE'].upper()} 투여인데 RATE와 DUR 값이 바르지 않게 기록됩.")
                    #     elif dspol_row['ROUTE'].upper() == 'IV_INFUSION':
                    #         # 잘못 써 놓은 경우
                    #         if (dspol_row['RATE'] in ('.', 0)):
                    #             raise ValueError(f"Dosing Policy에서 {dspol_row['ROUTE'].upper()} 투여인데 RATE 값이 0으로 기록됨.")
                    #
                    #         # Infusion이라고 서놨지만 iv bolus 투여와 같이 써놓은 경우
                    #         elif (dspol_row['RATE'] == -1) and (type(dspol_row['DUR']) in (int, float)) or (
                    #                 dspol_row['DUR'] == '.'):
                    #             add_row['NONMEM_RATE'] = dspol_row['RATE']
                    #             add_row['NONMEM_DUR'] = dspol_row['DUR']
                    #             print(
                    #                 f'iv infusion으로 쓰여있으나, 입력한 RATE와 DUR은 iv bolus와 같습니다. / {projectname}, {drug}, {nmid}')
                    #         elif (dspol_row['RATE'] == -2) and (dspol_row['DUR'] in ('.', 0)):
                    #             add_row['NONMEM_RATE'] = dspol_row['RATE']
                    #             add_row['NONMEM_DUR'] = dspol_row['DUR']
                    #             print(
                    #                 f'iv infusion으로 쓰여있으나, 입력한 RATE와 DUR은 iv bolus와 같습니다. / {projectname}, {drug}, {nmid}')
                    #
                    #         # Infusion에서 RATE를 모델에서 추정하는 경우
                    #         elif (dspol_row['RATE'] == -2) and (type(dspol_row['DUR']) in (int, float)):
                    #             add_row['NONMEM_RATE'] = dspol_row['RATE']
                    #             add_row['NONMEM_DUR'] = dspol_row['DUR']
                    #
                    #         # Infusion에서 RATE와 DUR를 모두 숫자로 잘 기입해 둔 경우
                    #         elif (type(dspol_row['RATE']) in (int, float)) and (type(dspol_row['DUR']) in (int, float)):
                    #             if (dspol_row['RATE'] > 0):
                    #                 add_row['NONMEM_RATE'] = dspol_row['RATE']
                    #                 add_row['NONMEM_DUR'] = dspol_row['DUR']
                    #             elif (dspol_row['RATE'] == -2) and (dspol_row['DUR'] not in ('.', 0)):
                    #                 add_row['NONMEM_RATE'] = dspol_row['RATE']
                    #                 add_row['NONMEM_DUR'] = dspol_row['DUR']
                    #             elif (dspol_row['RATE'] == -1) and (dspol_row['DUR'] not in ('.', 0)):
                    #                 raise ValueError(
                    #                     f"Dosing Policy에서 {dspol_row['ROUTE'].upper()} 투여인데 RATE=-1 / DUR=양수로 논리적으로 맞지 않습니다 (Bolus 인데, Duration 존재)")
                    #             else:
                    #                 raise ValueError(
                    #                     f"Dosing Policy에서 {dspol_row['ROUTE'].upper()} 투여인데 기타 다른 케이스임. 확인필요")

                    # 추가할 Dosing row를 위치할 곳 설정

                    if dspol_row['RELPOSITION'] not in (0, 1):
                        raise ValueError(
                            f"Dosing Policy에서 RELPOSITION 은 0 또는 1이어야함. 해당 row 기준으로 (0: 직전 / 1: 직후)에 추가 row 삽입")

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

        # Prep Data 에 Covariate정보 추가 (추후 구현)



        # Prep Data 저장

        if not os.path.exists(f"{modeling_dir_path}/prep_data"):
            os.mkdir(f"{modeling_dir_path}/prep_data")
        drug_nmprep_df.to_csv(f"{modeling_dir_path}/prep_data/MDP_{model_name}_({drug}).csv", index=False, encoding='utf-8-sig')

        # CMT, RATE, DUR, EVID, SS, ADDL


        ## Covariates





