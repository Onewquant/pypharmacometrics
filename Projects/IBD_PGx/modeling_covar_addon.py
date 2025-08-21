from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

## DEMO Covariates Loading

demo_df = pd.read_csv(f"{resource_dir}/demo/IBD_PGx_demo.csv")
demo_df = demo_df.rename(columns={'EMR ID':'UID','birthdate':'AGE','sex':'SEX','name':'NAME'})
demo_df = demo_df[['UID','AGE', 'SEX']].copy()
demo_df['UID'] = demo_df['UID'].astype(str)

## LAB Covariates Loading

totlab_df = pd.read_csv(f"{output_dir}/lab_df.csv")
# [x for x in totlab_df.columns if 'infliximab' in x.lower()]
totlab_df = totlab_df[['UID', 'DATETIME', 'Albumin', 'AST', 'AST(GOT)', 'ALT', 'ALT(GPT)', 'CRP', 'hsCRP', 'Calprotectin (Serum)', 'Calprotectin (Stool)', 'eGFR-CKD-EPI', 'Cr (S)', 'Creatinine','Anti-Infliximab Ab [정밀면역검사] (정량)']].copy()
# totlab_df = totlab_df[['UID', 'DATETIME', 'Albumin', 'AST', 'AST(GOT)', 'ALT', 'ALT(GPT)', 'CRP', 'hsCRP', 'Calprotectin (Serum)', 'Calprotectin (Stool)', 'Fecal Calprotectin', 'eGFR-CKD-EPI', 'Cr (S)', 'Creatinine','Anti-Infliximab Ab [정밀면역검사] (정량)']].copy()
totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'] = totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'].map(lambda x: float(re.findall(r'\d*\.\d+|\d+',str(x))[0]) if str(x)!='nan' else np.nan)
for c in list(totlab_df.columns)[2:]:  # break
    totlab_df[c] = totlab_df[c].map(lambda x: x if type(x) == float else float(re.findall(r'[\d]+.*[\d]*', str(x))[0]))
totlab_df['UID'] = totlab_df['UID'].astype(str)
totlab_df['AST'] = totlab_df[['AST', 'AST(GOT)']].max(axis=1)
totlab_df['ALT'] = totlab_df[['ALT', 'ALT(GPT)']].max(axis=1)
totlab_df['CREATININE'] = totlab_df[['Cr (S)', 'Creatinine']].max(axis=1)
# totlab_df['FCAL'] = totlab_df[['Calprotectin (Stool)', 'Fecal Calprotectin']].max(axis=1)
totlab_df['FCAL'] = totlab_df[['Calprotectin (Stool)']].max(axis=1)
totlab_df['CRP'] = totlab_df[['CRP', 'hsCRP']].max(axis=1)
# set(totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'].unique())
totlab_df = totlab_df.rename(columns={'Albumin': 'ALB', 'Calprotectin (Serum)': 'CALPRTSER', 'Anti-Infliximab Ab [정밀면역검사] (정량)':'ADA'})
totlab_df = totlab_df[['UID', 'DATETIME', 'ALB', 'AST', 'ALT', 'CRP', 'FCAL', 'CREATININE','ADA']].copy()
totlab_df['ADA'] = (totlab_df['ADA'] > 9)*1
# totlab_df = totlab_df[['UID', 'DATETIME', 'ALB', 'AST', 'ALT', 'CRP', 'FCAL', 'CREATININE']].copy()
for c in list(totlab_df.columns)[2:]:  # break
    totlab_df[c] = totlab_df[c].map(lambda x: x if type(x) == float else float(re.findall(r'[\d]+.*[\d]*', str(x))[0]))

totlab_df = totlab_df.drop_duplicates(['UID','DATETIME'])

## PD_BODY SIZE Lading

pd_bsize_df = pd.read_csv(f"{output_dir}/pd_bsize_df.csv")
pd_bsize_df['UID'] = pd_bsize_df['UID'].astype(str)

## Lab 중 CRP, FCAL 도 PD marker로 사용

pd_addon_df = totlab_df[['UID','DATETIME','CRP','FCAL']].rename(columns={'CRP':'PD_CRP','FCAL':'PD_FCAL'})
pd_bsize_df = pd_bsize_df.merge(pd_addon_df, on=['UID','DATETIME'], how='left')

## Modeling Data Loading
for drug in ['infliximab','adalimumab']:
    for mode_str in ['integrated',]:
        raw_modeling_df = pd.read_csv(f'{output_dir}/modeling_df_datacheck/{drug}_{mode_str}_datacheck.csv')
        raw_modeling_df['UID']= raw_modeling_df['UID'].astype(str)
        raw_modeling_df['DATETIME'] = raw_modeling_df['DATETIME'].map(lambda x:x.split('T')[0])

        # PD Baseline 설정
        uid_start_date_df = raw_modeling_df.groupby('UID',as_index=False)['DATETIME'].min()
        # raise ValueError
        uid_indphase_date_df = raw_modeling_df[(raw_modeling_df['TIME']==0)&(raw_modeling_df['MDV']==0)][['UID','DV']].copy()
        uid_indphase_date_df['PD_INDEXISTS'] = (uid_indphase_date_df['DV'].astype(float)==0)*1
        pd_basicinfo_df = uid_start_date_df.merge(uid_indphase_date_df[['UID','PD_INDEXISTS']], on=['UID'], how='left')

        new_pd_df = list()
        no_basepd_uids = list()
        no_1yrpd_uids = list()
        no_4mopd_uids = list()
        # pd_bsize df에서 'PD_' prefix가 들어가 있으면 PD marker로 인식함
        nobase_pd_marker_list = pd.DataFrame(pd_bsize_df.columns)[0]
        nobase_pd_marker_list = list(nobase_pd_marker_list[nobase_pd_marker_list.map(lambda x: 'PD_' in x)])
        for uid_inx, uid_pdbinfo_row in pd_basicinfo_df.iterrows(): #break
            uid = uid_pdbinfo_row['UID']
            ind_phase_existance = uid_pdbinfo_row['PD_INDEXISTS']
            uid_pd_df = pd_bsize_df[pd_bsize_df['UID']==uid].copy()

            # Baseline 기간에 해당하는 Data 보기
            uid_pdbase_df = uid_pd_df[uid_pd_df['DATETIME'] < uid_pdbinfo_row['DATETIME']].copy()
            if len(uid_pdbase_df)==0:
                print(f'({uid_inx}) / {uid} / Baseline 기간에 해당하는 Data 없음')
                uid_pdbase_row = dict()
                no_basepd_uids.append(uid)
            else:
                uid_pdbase_row = uid_pdbase_df.iloc[-1]

            # 1년 후 기간에 해당하는 Data 보기
            yr1frombaseline_datetime = (datetime.strptime(uid_pdbinfo_row['DATETIME'],'%Y-%m-%d')+timedelta(days=365)).strftime('%Y-%m-%d')
            uid_pd1yr_df = uid_pd_df[uid_pd_df['DATETIME'] >= yr1frombaseline_datetime].copy()
            if len(uid_pd1yr_df)==0:
                print(f'({uid_inx}) / {uid} / 1년 후 기간에 해당하는 Data 없음')
                uid_pd1yr_row = dict()
                no_1yrpd_uids.append(uid)
            else:
                uid_pd1yr_row = uid_pd1yr_df.iloc[0]

            # 4개월 후(127일) 기간에 해당하는 Data 보기
            mo4frombaseline_datetime = (datetime.strptime(uid_pdbinfo_row['DATETIME'],'%Y-%m-%d')+timedelta(days=127)).strftime('%Y-%m-%d')
            uid_pd4mo_df = uid_pd_df[uid_pd_df['DATETIME'] >= mo4frombaseline_datetime].copy()
            if len(uid_pd4mo_df)==0:
                print(f'({uid_inx}) / {uid} / 1년 후 기간에 해당하는 Data 없음')
                uid_pd4mo_row = dict()
                no_4mopd_uids.append(uid)
            else:
                uid_pd4mo_row = uid_pd4mo_df.iloc[0]

            uid_pd_df[f'PD_INDEXISTS'] = ind_phase_existance

            for pd_marker in nobase_pd_marker_list:
                uid_pd_df[f'{pd_marker}_BL'] = uid_pdbase_row[pd_marker] if len(uid_pdbase_row)!=0 else np.nan
                uid_pd_df[f'{pd_marker}_DELT'] = uid_pd_df[f'{pd_marker}']-uid_pd_df[f'{pd_marker}_BL']
                uid_pd_df[f'{pd_marker}_A4MO'] = uid_pd4mo_row[pd_marker] if len(uid_pd4mo_row)!=0 else np.nan
                uid_pd_df[f'{pd_marker}_A1YR'] = uid_pd1yr_row[pd_marker] if len(uid_pd1yr_row)!=0 else np.nan
                uid_pd_df[f'{pd_marker}_DELT1YR'] = uid_pd_df[f'{pd_marker}_A1YR']-uid_pd_df[f'{pd_marker}_BL']

            new_pd_df.append(uid_pd_df)

        new_pd_df = pd.concat(new_pd_df)

        print(f"No PD Baseline: {len(no_basepd_uids)}")
        print(f"No PD 1yr-after: {len(no_1yrpd_uids)}")

        raw_modeling_df = raw_modeling_df.merge(totlab_df, on=['UID','DATETIME'], how='left')
        raw_modeling_df = raw_modeling_df.merge(demo_df, on=['UID'], how='left')
        raw_modeling_df = raw_modeling_df.merge(new_pd_df, on=['UID','DATETIME'], how='left')
        # new_pd_df.columns
        # raw_modeling_df['UID'].iloc[0]

        # raise ValueError

        pd_marker_list = pd.DataFrame(new_pd_df.columns)[0]
        pd_marker_list = list(pd_marker_list[pd_marker_list.map(lambda x: 'PD_' in x)])

        ## Covariates의 NA value 처리 (ffill 먼저 시도, 없으면 bfill, 그것도 없으면 전체의 median 값)

        md_df_list = list()
        pd_marker_df = raw_modeling_df[pd_marker_list].copy()
        for md_inx, md_df in raw_modeling_df.drop(pd_marker_list, axis=1).groupby(['UID']):
            md_df = md_df.sort_values(['DATETIME']).fillna(method='ffill').fillna(method='bfill')
            md_df_list.append(md_df)
        modeling_df = pd.concat(md_df_list)
        # modeling_df.columns
        # PD마커들은 NAN 유지 (PD_PRO2, PD_MARKER)
        modeling_df.fillna(modeling_df.median(numeric_only=True), inplace=True)

        ## PD Marker들을 NAN 유지하면서 나머지들을 fillna 하는 작업중 -> 모델링, 시뮬레이션때 PD 마커들 여러 개는 원본 그대로 유지하도록 만들어보려함

        # raise ValueError
        # modeling_df.columns
        modeling_df = pd.concat([modeling_df,pd_marker_df], axis=1).reset_index(drop=True)

        ## Covariates의 NA value 처리 (ffill 먼저 시도, 없으면 bfill, 그것도 없으면 전체의 median 값)

        # md_df_list = list()
        # for md_inx, md_df in modeling_df.groupby(['UID']):
        #     # md_df = md_df.sort_values(['DATETIME']).fillna(method='ffill')
        #     md_df = md_df.sort_values(['DATETIME'])
        #     md_df_list.append(md_df)
        # modeling_df = pd.concat(md_df_list).reset_index(drop=True)

        # modeling_df.fillna('.', inplace=True)
        # modeling_df['UID'].drop_duplicates()
        # raise ValueError

        ## Modeling Data 생성
        right_covar_col = 'TIME(WEEK)'
        datacheck_cols = ['ID',	'UID', 'NAME', 'DATETIME','TIME(WEEK)','TIME(DAY)','TIME', 'DV', 'MDV', 'AMT', 'DUR', 'CMT', 'IBD_TYPE', 'ROUTE','DRUG', 'ADDED_ADDL'] + list(modeling_df.loc[:,right_covar_col:].iloc[:,1:].columns)

        if not os.path.exists(f'{output_dir}/modeling_df_covar'):
            os.mkdir(f'{output_dir}/modeling_df_covar')

        ## Modeling Data Check용 파일 저장
        modeling_df[datacheck_cols].to_csv(f'{output_dir}/modeling_df_covar/{drug}_{mode_str}_datacheck(covar).csv', index=False, encoding='utf-8-sig')

        modeling_cols = ['ID','TIME','DV','MDV','AMT','DUR','CMT','ROUTE','IBD_TYPE'] + list(modeling_df.loc[:,right_covar_col:].iloc[:,1:].columns)

        modeling_df['IBD_TYPE'] = modeling_df['IBD_TYPE'].map({'CD':1,'UC':2})
        modeling_df['AGE'] = modeling_df.apply(lambda x: int((datetime.strptime(x['DATETIME'],'%Y-%m-%d') - datetime.strptime(x['AGE'],'%Y-%m-%d')).days/365.25), axis=1)
        modeling_df['SEX'] = modeling_df['SEX'].map({'남':1,'여':2})
        modeling_df['ROUTE'] = modeling_df['ROUTE'].map({'IV':1,'SC':2,'.':'.'})

        # raise ValueError
        # modeling_df[modeling_df['AGE'] <= 18][['UID','NAME','AGE','IBD_TYPE']].drop_duplicates(['UID'], ignore_index=True)


        modeling_df = modeling_df[modeling_cols].sort_values(['ID','TIME'], ignore_index=True)
        modeling_input_line = str(list(modeling_df.columns)).replace("', '"," ")

        print(f"Mode: {mode_str} / {modeling_input_line}")

        ## 실제 NONMEM Modeling용 파일 저장
        modeling_df.drop(columns=pd_marker_list, axis=1).to_csv(f'{output_dir}/modeling_df_covar/{drug}_{mode_str}_modeling_df.csv',index=False, encoding='utf-8-sig')
        modeling_df['TIME']= modeling_df['TIME']/24
        modeling_df['DUR'] = modeling_df['DUR'].map(lambda x: float(x)/24 if x!='.' else x)

        ## PD 분석용 파일 저장
        modeling_df.drop(columns=pd_marker_list, axis=1).to_csv(f'{output_dir}/modeling_df_covar/{drug}_{mode_str}_modeling_df_dayscale.csv',index=False, encoding='utf-8')
        modeling_df.to_csv(f'{output_dir}/modeling_df_covar/{drug}_{mode_str}_pdeda_df_dayscale.csv',index=False, encoding='utf-8')



""" 
# 결과에서 Covariate 일부분이 비어있는 이유
# - basic_prep_lab_covar에서 각 사람마다의 lab 수치가 존재하는 날짜 기준으로 date range 를 생성했는데, 
# - 이게 order data (dosing date range)와 안 맞을 수 있다 
"""

# modeling_df[modeling_df['AST'].isna()]['UID'].unique()

####### NONMEM SDTAB
# nonmem_dir = 'C:/Users/ilma0/NONMEMProjects/IBDPGx/run'
#
# nmsdtab_df = pd.read_csv(f"{nonmem_dir}/sdtab39",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
# nmsdtab_df['ID'] = nmsdtab_df['ID'].astype(int)
# # nmsdtab_df['TDM_YEAR'] = nmsdtab_df['TDM_YEAR'].astype(int)
# under_pred_df = nmsdtab_df[(nmsdtab_df['MDV'] != 1)&(nmsdtab_df['DV'] > 5)&(nmsdtab_df['IPRED'] < 3)].copy()
# over_pred_df = nmsdtab_df[(nmsdtab_df['MDV'] != 1)&(nmsdtab_df['DV'] < 3)&(nmsdtab_df['IPRED'] > 5)].copy()
# mis_pred_df = pd.concat([under_pred_df, over_pred_df])
#
# modeling_df