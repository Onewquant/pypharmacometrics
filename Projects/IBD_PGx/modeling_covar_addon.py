from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

# rdf = pd.read_csv(f'{output_dir}/infliximab_induction_modeling_df.csv')
# rdf['ALB'].median()
# rdf['CRP'].median()
# rdf['CREATININE'].median()

## DEMO Covariates Loading

demo_df = pd.read_csv(f"{resource_dir}/demo/IBD_PGx_demo.csv")
demo_df = demo_df.rename(columns={'EMR ID':'UID','birthdate':'AGE','sex':'SEX','name':'NAME'})
demo_df = demo_df[['UID','AGE', 'SEX']].copy()
demo_df['UID'] = demo_df['UID'].astype(str)

# demo2_df = pd.read_csv(f"{resource_dir}/demo2/IBD_PGx_demo2.csv")
# demo2_df = demo2_df.rename(columns={'EMR ID':'UID','height':'HT','weight':'WT','bmi':'BMI','bsa':'BSA'})
# demo2_df = demo2_df[['UID','HT', 'WT', 'BMI','BSA']].copy()
# demo2_df['UID'] = demo2_df['UID'].astype(str)

# ## WT, HT, ADHERENCE Covariates / PD markers Loading
#
# pd_marker_list = ['PD_PRO2', 'PD_STLFREQ', 'PD_ABDPAIN', 'PD_RECTBLD','PD_TOTALSCORE']
#
# cdai_hx_df = pd.read_csv(f"{output_dir}/pdmarker_cdai_df.csv")
# cdai_hx_df = cdai_hx_df.rename(columns={'ID':'UID','DATE':'DATETIME','CDAI_TOTALSCORE':'PD_TOTALSCORE','CDAI_ADHERENCEPCT':'ADHERENCE','CDAI_WEIGHT':'WT','CDAI_HEIGHT':'HT', 'CDAI_BMI':'BMI'})
# cdai_hx_df['UID'] = cdai_hx_df['UID'].astype(str)
# cdai_hx_df['PD_PRO2'] = 2*cdai_hx_df['CDAI_DIARHEACNT'] + 5*cdai_hx_df['CDAI_ABDPAINCNT']
# cdai_hx_df['PD_STLFREQ'] = cdai_hx_df['CDAI_DIARHEACNT']
# cdai_hx_df['PD_ABDPAIN'] = cdai_hx_df['CDAI_ABDPAINCNT']
# cdai_hx_df['PD_RECTBLD'] = np.nan
# cdai_hx_df['IBD_TYPE'] = 'CD'
# # cdai_adh_df.columns
#
# pms_hx_df = pd.read_csv(f"{output_dir}/pdmarker_pms_df.csv")
# pms_hx_df = pms_hx_df.rename(columns={'ID':'UID'}).rename(columns={'ID':'UID','DATE':'DATETIME','PMS_TOTALSCORE':'PD_TOTALSCORE','PMS_ADHERENCEPCT':'ADHERENCE','PMS_WEIGHT':'WT','PMS_HEIGHT':'HT', 'PMS_BMI':'BMI'})
# pms_hx_df['UID'] = pms_hx_df['UID'].astype(str)
# pms_hx_df['PD_PRO2'] = pms_hx_df['PMS_STLCNT'] + pms_hx_df['PMS_HEMATOCHEZIA']
# pms_hx_df['PD_STLFREQ'] = pms_hx_df['PMS_STLCNT']
# pms_hx_df['PD_ABDPAIN'] = np.nan
# pms_hx_df['PD_RECTBLD'] = pms_hx_df['PMS_HEMATOCHEZIA']
# pms_hx_df['IBD_TYPE'] = 'UC'
# # pms_adh_df.columns
#
# adh_df = pd.concat([cdai_hx_df[['UID','DATETIME','ADHERENCE','IBD_TYPE'] + pd_marker_list], pms_hx_df[['UID','DATETIME','ADHERENCE','IBD_TYPE']+pd_marker_list]]).drop_duplicates(['UID','DATETIME']).sort_values(['UID','DATETIME'])
# bsize_df = pd.concat([cdai_hx_df[['UID','DATETIME','WT','HT','BMI']], pms_hx_df[['UID','DATETIME','WT','HT','BMI']]]).drop_duplicates(['UID','DATETIME']).sort_values(['UID','DATETIME'])
# # adh_df[adh_df['WT'].isna()]
# # adh_df[adh_df['HT'].isna()]
#
# full_result_df = list()
# count = 0
# for uid, uid_df in adh_df.groupby('UID',as_index=False): #break
#
#     min_lab_date = uid_df['DATETIME'].min()
#     max_lab_date = uid_df['DATETIME'].max()
#
#     print(f"({count}) {uid} / Date Range : ({min_lab_date},  {max_lab_date})")
#
#     uid_fulldt_df = pd.DataFrame(columns=['UID','DATETIME'])
#     uid_fulldt_df['DATETIME'] = pd.date_range(start=min_lab_date,end=max_lab_date).astype(str)
#     uid_fulldt_df['UID'] = uid
#
#     # uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left').fillna(method='bfill')
#     uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')
#     # uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')
#
#     if count==0:
#         full_result_df = uid_fulldt_df.copy()
#     else:
#         full_result_df = pd.concat([full_result_df, uid_fulldt_df.copy()])
#
#     count+=1
#
# # raise ValueError
# adh_df = full_result_df.copy()
# # pdm_med_uc = adh_df[adh_df['IBD_TYPE']=='UC']['PD_MARKER'].median()
# # pdm_med_cd = adh_df[adh_df['IBD_TYPE']=='CD']['PD_MARKER'].median()
# #
# # pdpro2_med_uc = adh_df[adh_df['IBD_TYPE']=='UC']['PD_PRO2'].median()
# # pdpro2_med_cd = adh_df[adh_df['IBD_TYPE']=='CD']['PD_PRO2'].median()
# #
# # med_pd_value = {'UC':pdm_med_uc,'CD':pdm_med_cd}
# # full_result_df['PD_MARKER'] = full_result_df.apply(lambda row:row['PD_MARKER'] if not np.isnan(row['PD_MARKER']) else med_pd_value[row['IBD_TYPE']], axis=1)
# # full_result_df['PD_PRO2'] = full_result_df.apply(lambda row:row['PD_PRO2'] if not np.isnan(row['PD_PRO2']) else 7777777, axis=1)
# # adh_df = full_result_df.fillna(full_result_df.median(numeric_only=True))
# # adh_df['PD_PRO2'] = adh_df['PD_PRO2'].replace(7777777,np.nan)
#
# adh_df = adh_df.drop('IBD_TYPE', axis=1)
#
#
# full_result_df = list()
# count = 0
# for uid, uid_df in bsize_df.groupby('UID',as_index=False): #break
#
#     min_lab_date = uid_df['DATETIME'].min()
#     max_lab_date = uid_df['DATETIME'].max()
#
#     print(f"({count}) {uid} / Date Range : ({min_lab_date},  {max_lab_date})")
#
#     uid_fulldt_df = pd.DataFrame(columns=['UID','DATETIME'])
#     uid_fulldt_df['DATETIME'] = pd.date_range(start=min_lab_date,end=max_lab_date).astype(str)
#     uid_fulldt_df['UID'] = uid
#
#     uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left').fillna(method='ffill')
#     # uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')
#
#     if count==0:
#         full_result_df = uid_fulldt_df.copy()
#     else:
#         full_result_df = pd.concat([full_result_df, uid_fulldt_df.copy()])
#
#     count+=1
# bsize_df = full_result_df.fillna(full_result_df.median(numeric_only=True))

# adh_df['ADHERENCE'].drop_duplicates()
# adh_df['UID'].drop_duplicates()
# adh_df.to_csv(f"{output_dir}/adherence_df.csv", encoding='utf-8-sig', index=False)

## LAB Covariates Loading

totlab_df = pd.read_csv(f"{output_dir}/lab_df.csv")
totlab_df = totlab_df[['UID', 'DATETIME', 'Albumin', 'AST', 'AST(GOT)', 'ALT', 'ALT(GPT)', 'CRP', 'Calprotectin (Serum)', 'Calprotectin (Stool)', 'eGFR-CKD-EPI', 'Cr (S)', 'Creatinine','Anti-Infliximab Ab [정밀면역검사] (정량)']].copy()
totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'] = totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'].map(lambda x: float(re.findall(r'\d*\.\d+|\d+',str(x))[0]) if str(x)!='nan' else np.nan)
for c in list(totlab_df.columns)[2:]:  # break
    totlab_df[c] = totlab_df[c].map(lambda x: x if type(x) == float else float(re.findall(r'[\d]+.*[\d]*', str(x))[0]))
totlab_df['UID'] = totlab_df['UID'].astype(str)
totlab_df['AST'] = totlab_df[['AST', 'AST(GOT)']].max(axis=1)
totlab_df['ALT'] = totlab_df[['ALT', 'ALT(GPT)']].max(axis=1)
totlab_df['CREATININE'] = totlab_df[['Cr (S)', 'Creatinine']].max(axis=1)
# set(totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'].unique())
totlab_df = totlab_df.rename(columns={'Albumin': 'ALB', 'Calprotectin (Stool)': 'CALPRTSTL', 'Calprotectin (Serum)': 'CALPRTSER', 'Anti-Infliximab Ab [정밀면역검사] (정량)':'ADA'})
totlab_df = totlab_df[['UID', 'DATETIME', 'ALB', 'AST', 'ALT', 'CRP', 'CALPRTSTL', 'CREATININE','ADA']].copy()
totlab_df['ADA'] = (totlab_df['ADA'] > 9)*1
# totlab_df = totlab_df[['UID', 'DATETIME', 'ALB', 'AST', 'ALT', 'CRP', 'CALPRTSTL', 'CREATININE']].copy()
for c in list(totlab_df.columns)[2:]:  # break
    totlab_df[c] = totlab_df[c].map(lambda x: x if type(x) == float else float(re.findall(r'[\d]+.*[\d]*', str(x))[0]))

totlab_df = totlab_df.drop_duplicates(['UID','DATETIME'])

## PD_BODY SIZE Lading

pd_bsize_df = pd.read_csv(f"{output_dir}/pd_bsize_df.csv")
pd_bsize_df['UID'] = pd_bsize_df['UID'].astype(str)
pd_bsize_df = pd_bsize_df.drop(['IBD_TYPE'],axis=1)

## Lab 중 CRP, CALPRTSTL 도 PD marker로 사용

pd_addon_df = totlab_df[['UID','DATETIME','CRP','CALPRTSTL']].rename(columns={'CRP':'PD_CRP','CALPRTSTL':'PD_CALPRTSTL'})
pd_bsize_df = pd_bsize_df.merge(pd_addon_df, on=['UID','DATETIME'], how='left')

# [c for c in totlab_df.columns.unique() if 'ada' in c.lower()]

# totlab_df[~totlab_df['ALT'].isna()]['ALT']

## Modeling Data Loading
# for drug in ['infliximab',]:
for drug in ['infliximab','adalimumab']:
    for mode_str in ['integrated',]:
        raw_modeling_df = pd.read_csv(f'{output_dir}/modeling_df_datacheck/{drug}_{mode_str}_datacheck.csv')
        raw_modeling_df['UID']= raw_modeling_df['UID'].astype(str)
        raw_modeling_df['DATETIME'] = raw_modeling_df['DATETIME'].map(lambda x:x.split('T')[0])

        # PD Baseline 설정
        uid_start_date_df = raw_modeling_df.groupby('UID',as_index=False)['DATETIME'].min()
        uid_indphase_date_df = raw_modeling_df[(raw_modeling_df['TIME']==0)&(raw_modeling_df['MDV']==0)][['UID','DV']].copy()
        uid_indphase_date_df['PD_INDEXISTS'] = (uid_indphase_date_df['DV'].astype(float)==0)*1
        pd_basicinfo_df = uid_start_date_df.merge(uid_indphase_date_df[['UID','PD_INDEXISTS']], on=['UID'], how='left')

        new_pd_df = list()
        no_basepd_uids = list()
        no_1yrpd_uids = list()

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

            uid_pd_df[f'PD_INDEXISTS'] = ind_phase_existance

            for pd_marker in nobase_pd_marker_list:
                uid_pd_df[f'{pd_marker}_BL'] = uid_pdbase_row[pd_marker] if len(uid_pdbase_row)!=0 else np.nan
                uid_pd_df[f'{pd_marker}_DELT'] = uid_pd_df[f'{pd_marker}']-uid_pd_df[f'{pd_marker}_BL']
                uid_pd_df[f'{pd_marker}_A1YR'] = uid_pd1yr_row[pd_marker] if len(uid_pd1yr_row)!=0 else np.nan
                uid_pd_df[f'{pd_marker}_DELT1YR'] = uid_pd_df[f'{pd_marker}_A1YR']-uid_pd_df[f'{pd_marker}_BL']

            new_pd_df.append(uid_pd_df)

        new_pd_df = pd.concat(new_pd_df)
        # new_pd_df['UID'] = new_pd_df['UID'].astype(str)

        print(f"No PD Baseline: {len(no_basepd_uids)}")
        print(f"No PD 1yr-after: {len(no_1yrpd_uids)}")

        # new_pd_df

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

        ## Modeling Data Saving
        # data_check_cols = ['ID','UID','NAME','DATETIME','TIME','DV','MDV','AMT','DUR','CMT','IBD_TYPE','ADDED_ADDL'] + list(modeling_df.loc[:,'DRUG':].iloc[:,1:].columns)
        # modeling_df[data_check_cols].to_csv(f'{output_dir}/{drug}_{mode_str}_datacheck_covar.csv', index=False, encoding='utf-8-sig')
        right_covar_col = 'TIME(WEEK)'
        datacheck_cols = ['ID',	'UID', 'NAME', 'DATETIME','TIME(WEEK)','TIME(DAY)','TIME', 'DV', 'MDV', 'AMT', 'DUR', 'CMT', 'IBD_TYPE', 'ROUTE','DRUG', 'ADDED_ADDL'] + list(modeling_df.loc[:,right_covar_col:].iloc[:,1:].columns)

        if not os.path.exists(f'{output_dir}/modeling_df_covar'):
            os.mkdir(f'{output_dir}/modeling_df_covar')

        modeling_df[datacheck_cols].to_csv(f'{output_dir}/modeling_df_covar/{drug}_{mode_str}_datacheck(covar).csv', index=False, encoding='utf-8-sig')
        # dcheck_df = modeling_df[~(modeling_df['DV'].isin(['.','0.0']))][['ID',	'UID', 'NAME', 'DATETIME', 'DV']].copy()
        # dcheck_df[dcheck_df['DV'].map(float) > 40]

        modeling_cols = ['ID','TIME','DV','MDV','AMT','DUR','CMT','ROUTE','IBD_TYPE'] + list(modeling_df.loc[:,right_covar_col:].iloc[:,1:].columns)

        modeling_df['IBD_TYPE'] = modeling_df['IBD_TYPE'].map({'CD':1,'UC':2})
        modeling_df['AGE'] = modeling_df.apply(lambda x: int((datetime.strptime(x['DATETIME'],'%Y-%m-%d') - datetime.strptime(x['AGE'],'%Y-%m-%d')).days/365.25), axis=1)
        modeling_df['SEX'] = modeling_df['SEX'].map({'남':1,'여':2})
        modeling_df['ROUTE'] = modeling_df['ROUTE'].map({'IV':1,'SC':2,'.':'.'})


        modeling_df = modeling_df[modeling_cols].sort_values(['ID','TIME'], ignore_index=True)

        modeling_input_line = str(list(modeling_df.columns)).replace("', '"," ")

        print(f"Mode: {mode_str} / {modeling_input_line}")

        # if mode_str=='maintenance':  # Time decrease가 생김... 확인 !
        #     modeling_df["prev_time"] = modeling_df.groupby("ID")["TIME"].shift(1)
        #     modeling_df["time_decrease"] = (modeling_df["TIME"] < modeling_df["prev_time"])
        #     problem_rows = modeling_df[modeling_df["time_decrease"]]
        #     print(problem_rows)

        # modeling_df['CALPRTSTL'].median()
        # modeling_df['CREATININE'].median()
        # modeling_df['ALB'].median()
        # modeling_df['CRP'].median()
        # raise ValueError
        # modeling_df['A_0FLG'] = (modeling_df['ID'].shift(1)!=modeling_df['ID'])*1

        # raise ValueError
        modeling_df.drop(columns=pd_marker_list, axis=1).to_csv(f'{output_dir}/modeling_df_covar/{drug}_{mode_str}_modeling_df.csv',index=False, encoding='utf-8-sig')
        # modeling_df[~((modeling_df['TIME'] == 0) & (modeling_df['DV'] == '0.0'))].copy().to_csv(f'{output_dir}/{drug}_{mode_str}_modeling_df_without_zero_dv.csv',index=False, encoding='utf-8-sig')
        modeling_df['TIME']= modeling_df['TIME']/24
        modeling_df['DUR'] = modeling_df['DUR'].map(lambda x: float(x)/24 if x!='.' else x)

        modeling_df.drop(columns=pd_marker_list, axis=1).to_csv(f'{output_dir}/modeling_df_covar/{drug}_{mode_str}_modeling_df_dayscale.csv',index=False, encoding='utf-8')
        modeling_df.to_csv(f'{output_dir}/modeling_df_covar/{drug}_{mode_str}_pdeda_df_dayscale.csv',index=False, encoding='utf-8')


        # modeling_df['AMT'] = (modeling_df['AMT'].replace('.',0).map(float)) * (np.abs((modeling_df['ADHERENCE'].replace(0,0.0001)) / 100)).replace(0,'.')
        # modeling_df.to_csv(f'{output_dir}/{drug}_{mode_str}_modeling_df_dayscale_adherence.csv', index=False, encoding='utf-8')

        # raise ValueError
        # modeling_df['AMT'] = ((modeling_df['AMT'].replace('.', 0).map(float)) * ((((modeling_df['CMT']==1)&(modeling_df['MDV']==1))*1.5) + ((modeling_df['CMT']==2)&(modeling_df['MDV']==1))*1)).replace(0, '.')
        # modeling_df.to_csv(f'{output_dir}/{drug}_{mode_str}_modeling_df_dayscale_scdouble.csv', index=False, encoding='utf-8')

        # modeling_df['CALPRTSTL'].median()
        # modeling_df[modeling_df['ROUTE'] == 2]['AMT'].unique()
        # len(modeling_df['ID'].unique())
        # len(modeling_df[(modeling_df['TIME']==0)&(modeling_df['MDV']==0)]['ID'].unique())
        # len(modeling_df[(modeling_df['TIME'] == 0)&(modeling_df['MDV'] == 1)]['ID'].unique())
        # raise ValueError
        # len(modeling_df[(modeling_df['TIME']==0)&(modeling_df['DV']=='0.0')]['ID'].unique())
        # len(modeling_df[(modeling_df['TIME']==0)&(modeling_df['DV']!='0.0')&(modeling_df['MDV']!=1)]['ID'].unique())
        # raise ValueError

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