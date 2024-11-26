from tools import *

## 기본정보입력

project_dict = {'CKD379':['Empagliflozin','Sitagliptin','Metformin'], 'CKD383':['Empagliflozin','Lobeglitazone','Metformin']}

modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/metformin_xr"

## Dosing_Policy 파일 불러오기

dspol_df = pd.read_csv(f'{modeling_dir_path}/Dosing_Policy.csv')

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
for drug in drugconc_dict.keys(): break
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
    dcdf['NONMEM_RATE'] = '.'
    dcdf['NONMEM_DUR'] = '.'

    ## CMT

    """
    # 주로 농도측정이 Systemic compartment에서 진행되므로 2로 표시하는게 나을듯. 
    # dosing policy 파일에 이 정보가 반영되도록 만들어야. 
    # (dosing policy / modeling policy 가 모두 영향을 주는 듯)
    """

    dcdf['NONMEM_CMT'] = 2

    # Project/Drug별로 Dosing Policy 다르다고 가정, Dosing Policy에 따라 Dosing row 추가 (drug은 위에서 이미 구분되어있음 : dcdf)

    for projectname, prj_dcdf in dcdf.groupby(['PROJECT']):

        prj_dspol = dspol_df[(dspol_df['DRUG']==drug)&(dspol_df['PROJECT']==projectname)]

        for dspol_inx, dspol_row in prj_dspol.iterrows():

            for nmid, nmid_df in dcdf.groupby(['NONMEM_ID']):
                
                # 추가할 Dosing row

                add_row = nmid_df.iloc[0:1, :]
                # add_row['NONMEM_ID']
                add_row['NONMEM_TIME']=dspol_row['RELTIME']
                add_row['NONMEM_TAD']=dspol_row['RELTIME']
                add_row['NONMEM_DV'] = '.'
                add_row['NONMEM_MDV'] = 1

                add_row['NONMEM_AMT'] = dspol_row['DOSE']

                if dspol_row['ROUTE']=='IV':
                    add_row['NONMEM_RATE'] = dspol_row['RATE']
                    add_row['NONMEM_DUR'] = dspol_row['DUR']
                elif (dspol_row['ROUTE']=='PO') and (dspol_row['ABS_ORD']==0):
                    add_row['NONMEM_RATE'] = -2
                    add_row['NONMEM_DUR'] = np.nan
                else:
                    add_row['NONMEM_RATE'] = np.nan
                    add_row['NONMEM_DUR'] = np.nan

                # 추가할 Dosing row를 위치할 곳 설정

                rt_inx = nmid_df[nmid_df[f"NONMEM_{dspol_row['RELTIMECOL']}"]==dspol_row['RELTIME']].iloc[0].name
                add_pos_inx = rt_inx + dspol_row['RELPOSITION']
                if add_pos_inx < 0:
                    raise ValueError(f"Dosing row를 추가하려는 위치가 0번째 row보다 작은값입니다. / {projectname}, {drug}, {nmid}")
                elif add_pos_inx==0:
                    pass
                else:
                    
        break

    # CMT, RATE, DUR, EVID, SS, ADDL


    ## Covariates





