from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

prj_name = 'VPA'
prj_dir = './Projects/VPA'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

raw_df = pd.read_csv(f"{resource_dir}/vpa_cdw_order_data.csv")
# df.columns
raw_df = raw_df.rename(columns={'환자번호':'UID','수행시간':'ACTING','약품 오더일자':'DATE', '[실처방] 용법':'REGIMEN','[실처방] 처방일수':'DAYS', '[함량단위환산] 1일 처방량':'DAILY_DOSE','[함량단위환산] 1회 처방량':'DOSE', '[실처방] 투약위치':'PLACE',"[실처방] 경로":'ROUTE','약품명(성분명)':'DRUG','[실처방] 처방비고':'ETC_INFO'})
raw_df['UID'] = raw_df['UID'].map(lambda x:x.split('-')[0])
# df.columns
# df['DRUG'].unique()
drugname_dict = {'valproic acid':'valproate','lacosamide':'lacosamide'}


drug_adm_res = dict()
for raw_drugname, drugname in drugname_dict.items():

    drug_res_df = list()

    df = raw_df[raw_df['DRUG']==raw_drugname].copy()
    # df = df.rename(columns={'[실처방] 투약위치':'PLACE'})
# raw_df[raw_df['DRUG']=='lacosamide']
# df['[함량단위환산] 1회 처방량'].unique()

    alter_acting_dict = {
        '1일 1회':('q24h','09:00/Y'),
        '1일 2회':('q12h','09:00/Y, 21:00/Y'),
        '1일 3회':('q8h','08:00/Y, 16:00/Y, 23:59/Y'),
        '1일 4회':('q6h','04:00/Y, 10:00/Y, 16:00/Y, 22:00/Y'),
        '1일 5회':('q4.8h','04:00/Y, 09:00/Y, 13:00/Y, 18:00/Y, 23:00/Y'),
        '1일 6회':('q4h','08:00/Y, 12:00/Y, 16:00/Y, 20:00/Y, 23:59/Y'),
        '1일 8회':('q3h','02:00/Y, 05:00/Y, 08:00/Y, 11:00/Y, 14:00/Y, 17:00/Y, 20:00/Y, 23:00/Y'),
        '필요시 자기전에':('q[hours]h','22:00/Y'),
        '필요시 복용하세요':('q[hours]h','09:00/Y'),
        '필요시 1시간전에': ('q[hours]h', '09:00/Y'),
        '격일로 자기전':('q48h','22:00/Y'),
        '격일로 저녁':('q48h','21:00/Y'),
        '격일로 아침':('q48h','09:00/Y'),
        '1년 1회':('q[hours]h','09:00/Y'),
        '1주 3회':('q56h','09:00/Y'),
        '의사지시대로 관류하세요':('q[hours]h','09:00/Y'),
     }

    # df[df['UID']=='352504997098455']
    # df['REGIMEN'].unique()
    # df['DOSE'].unique()
    dose_list = list()
    alter_acting_list = list()
    interval_list = list()
    ondemand_list = list()
    count = 0
    prev_uid = '0000'
    for inx, row in df.iterrows():
        if prev_uid!=row['UID']:
            count+=1
            print(f"({count}) {drugname} / {row['UID']}")
            prev_uid = row['UID']
        ## '250mg-500mg-500mg' 형태로 dose가 쓰여있는 경우
        # if row['']
        seq_dose_patterns = re.findall(r'\d+\.*\d*', row['DOSE'])
        if len(seq_dose_patterns) >= 1:
            dose_val = ('_'.join(seq_dose_patterns))
        else:
            dose_val = seq_dose_patterns[0]
            # dose_val = re.findall(r'\d+',x.split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0]
            # raise ValueError
        dose_list.append(dose_val)

        daily_count_str = ' '.join(row['REGIMEN'].split(' ')[:2])


        if row['REGIMEN'].split(' ')[0]=='필요시':
            # print(f"{row['REGIMEN']} / hour_str: {str(int(24/len(seq_dose_patterns)))} / DOSE: {row['DOSE']} ({row['DAILY_DOSE']})")
            ondemand_list.append(True)
        else:
            ondemand_list.append(False)

        alter_acting_list.append(alter_acting_dict[daily_count_str][-1])
        hour_str = str(int(24/len(seq_dose_patterns)))
        interval_list.append(alter_acting_dict[daily_count_str][0].replace('[hours]',hour_str))

    # df['ROUTE'].unique()
    df['DOSE'] = dose_list
    df['ALTER_ACTING'] = alter_acting_list
    df['INTERVAL'] = interval_list
    df['ON_DEMAND'] = ondemand_list
    df['DRUG'] = drugname

    # df['DRUG']
    df['ROUTE'] = df['ROUTE'].map({'P.O':'PO','MIV':'IV','PLT':'PO','IVS':'IV','SC':'SC','S.L':"SL", 'IntraPleural':'IP'}) # PLT: ( Par L-tube) / S.L: sublingual
    df['ADDL'] = df['DAYS'].copy()
    # df = df[df['ROUTE']!='Irrigation'].copy()

    df = df.sort_values(['UID','DATE']).reset_index(drop=True)
    # df[df['UID']==df['UID'].unique()[2]][['DATE','DOSE','ROUTE','INTERVAL','ADDL','ACTING','PLACE']]
    # raise ValueError


    dose_cond1 = ~(df['ACTING'].isna())                         # 본원투약기록 -> 투약 기록대로 추가
    dose_cond2 = (df['ADDL']>1)&(df['ACTING'].isna())           # 외래자가투약 / 타병원 투약 처방 -> 날짜대로 ALTER_ACTING으로 추가
    # dose_cond3 = (df['ACTING'].isna())

    df1 = df[dose_cond1].sort_values(['UID','DATE']).reset_index(drop=True)
    df1['VIRTUALITY'] = False
    df2 = df[dose_cond2].sort_values(['UID','DATE']).reset_index(drop=True)
    df2['ACTING'] = df2['ALTER_ACTING'].copy()
    df2['VIRTUALITY'] = True

    sdf = pd.concat([df1, df2]).sort_values(['UID', 'DATE', 'ACTING'])[['UID','DATE','DOSE','ROUTE','INTERVAL','ADDL','ACTING','PLACE']].copy()

    sdf_row_count = 0
    for uid, uid_sdf in sdf.groupby('UID'): #break
        # frag_sdf1 = frag_sdf[frag_sdf['ADDL']==1].copy()
        # frag_sdf1_acting_rows = frag_sdf1[frag_sdf1['ACTING'].map(lambda x: ('/Y' in x) or ('/O' in x))].copy()
        # if len(frag_sdf1_acting_rows)==0:
        #     frag_sdf1_dates = set()
        # else:
        #     frag_sdf1_dates = set(frag_sdf1_acting_rows['DATE'])
        #
        # frag_sdf2 = frag_sdf[frag_sdf['ADDL']>1].copy()
        uid_sdf['II_NUM'] = uid_sdf['INTERVAL'].map(lambda x: float(x.replace('q','').replace('h','')))
        uid_sdf['II_APPLIED'] = uid_sdf['II_NUM'].map(lambda x: x/24 if x>24 else 1.0)

        # frag_sdf1_dates = set(frag_sdf1[frag_sdf1['ACTING'].map(lambda x: ('/Y' in x) or ('/O' in x))]['DATE'])
        # frag_sdf2_dates = set()
        uid_sdf_dates = set()
        for sdf_inx, uid_sdf_row in uid_sdf.iterrows(): #break
            # if uid_sdf_row['ADDL'] == 1:
            #     raise ValueError
            # if uid_sdf_row['II_APPLIED'] > 1:
            #     raise ValueError
            # dates_every_step_days(start_date=uid_sdf_row['DATE'], day_step=uid_sdf_row['II_APPLIED'], addl=uid_sdf_row['ADDL'], fmt="%Y-%m-%d", include_start=True)
            uid_sdf_dates = uid_sdf_dates.union(set(dates_every_step_days(start_date=uid_sdf_row['DATE'], day_step=uid_sdf_row['II_APPLIED'], addl=uid_sdf_row['ADDL'], fmt="%Y-%m-%d", include_start=True)))

        uid_dates_rows = pd.DataFrame(uid_sdf_dates, columns=['DATE']).sort_values(['DATE'], ignore_index=True)
        uid_dates_row = {'UID':uid,'DRUG':drugname,'DATE_LIST':tuple(uid_dates_rows['DATE'])}
        drug_res_df.append(uid_dates_row)

        print(f"({sdf_row_count}) / {drugname} / {uid} / {len(uid_sdf_dates)}")
        sdf_row_count+=1

    drug_res_df = pd.DataFrame(drug_res_df)
    # drug_res_df.to_csv(f"{output_dir}/{drugname}_adm_dates.csv", encoding='utf-8-sig',index=False)
    drug_res_df.to_excel(f"{output_dir}/{drugname}_adm_dates.xlsx", encoding='utf-8-sig', index=False)
    # drug_adm_res['valproate'].to_excel(f"{output_dir}/valproate_adm_dates.xlsx", encoding='utf-8-sig', index=False)
    #
    drug_adm_res[drugname] = drug_res_df.copy()

        # uid_dates_rows['DATE_DIFF'] = pd.to_datetime(uid_dates_rows['DATE']).diff().dt.days.fillna(1)
        # uid_dates_rows.loc[uid_dates_rows['DATE_DIFF'] >= 5].index.min()

        # raise ValueError


    # df2['ADDL'].mean() * len(df2)

    # df[dose_cond3]

    # df = pd.concat([df1, df2]).sort_values(['UID','DATE','ACTING'])
    # df1[df1['DAYS']>1][['UID','DATE','INTERVAL','ADDL','ACTING']]

    # ADDL 먼저 추가 작업해야함

    # df2["DATE"] = pd.to_datetime(df2["DATE"], errors="coerce")
    #
    # raise ValueError
    #
    # addl_applied_df2 = list()
    # for inx, row in df2.iterrows():#break
    #     interval_num = float(row['INTERVAL'].replace('q','').replace('h',''))
    #     print(f"({inx}) {row['UID']} / {row['ADDL']} ADDL")
    #     if interval_num > 24:
    #         day_step = interval_num/24
    #         # raise ValueError
    #         # pass
    #     else:
    #         day_step = 1
    #
    #     # end_date = day_step * row["ADDL"]
    #
    #     for i in range(row["ADDL"]):  # ADDL 포함해서 +1
    #         new_row = row.copy()
    #         new_row["DATE"] = (row["DATE"] + pd.Timedelta(days=day_step*i)).strftime("%Y-%m-%d")
    #         new_row["ADDL"] = 1
    #         addl_applied_df2.append(new_row)
    #
    # df2 = pd.DataFrame(addl_applied_df2)
    #
    # #######################################
    #
    # # df = pd.concat([df1,df2]).sort_values(['UID','DATE']).reset_index(drop=True)
    #
    # # DATETIME/DOSE 연결 작업
    #
    # dt_dose_series = list()
    # vacant_data = '0000-00-00TNN:NN'
    # for inx, row in df.iterrows():
    #     # raise ValueError
    #
    #     # if (row['ID']=='10023985')&(row['DATETIME']=='2004-04-19T06:51'):
    #     #     raise ValueError
    #     # else:
    #     #     continue
    #
    #     # if (row['DATETIME']=='2005-06-27T20:00'):
    #     #     raise ValueError
    #     # else:
    #     #     continue
    #
    #     rep_acting_str = row['ACTING'].replace('O.C,','').strip()
    #
    #     new_actval_str=''
    #     for actval in rep_acting_str.split(','): #break
    #
    #         ## 빈칸일때
    #         if actval=='':
    #             # new_actval_str += f'_{vacant_data}'
    #             continue
    #
    #
    #         ## 날짜 시간/Y일때
    #         y_val = re.findall(r"\d\d\d\d-\d\d-\d\d \d\d:\d\d/Y",actval)
    #         if len(y_val) > 0:
    #             new_actval_str+='_' + y_val[0].replace("/Y","").replace(" ","T")
    #             continue
    #         else:
    #             pass
    #
    #         ## 시간/Y 일떄
    #         y_val = re.findall(r"\d\d:\d\d/Y",actval)
    #         if len(y_val) > 0:
    #             new_actval_str+=f'_{row["DATE"]}T{y_val[0].replace("/Y","")}'
    #             continue
    #         else:
    #             rest_y_val = actval.split('/')[-1]
    #
    #             """
    #             # 확인 필요한 부분: D나 N은 투약이 된것인지?
    #             """
    #             if (('Y ' in rest_y_val) or ('O ' in rest_y_val)):
    #                 raise ValueError
    #             else:
    #                 new_actval_str += f'_{vacant_data}'
    #                 continue
    #
    #         raise ValueError
    #
    #     new_actval_str = new_actval_str[1:]
    #     # print(22)
    #     # raise ValueError
    #     # DOSE 붙이기 작업
    #     new_actval_split = new_actval_str.split('_')
    #
    #
    #     new_actval_str = '_'.join([f"{nav}DOSE{row['DOSE']}" for nav in new_actval_split])
    #     dt_dose_series.append(new_actval_str)
    # df['DT_DOSE'] = dt_dose_series
    #
    # df.to_csv(f"{output_dir}/{raw_drugname}_dt_dose_df.csv", encoding='utf-8-sig', index=False)
    # # df[['UID', 'DATE', 'DRUG', 'DOSE', 'ROUTE','REGIMEN','DAYS', 'ETC_INFO', 'ACTING', 'ALTER_ACTING', 'INTERVAL', 'VIRTUALITY', 'DT_DOSE']]
    # # df['INTERVAL']
    # # dose_result_df[['DATE','DOSE','ACTING','DT_DOSE']]
    #
    # #
    # # dose_result_df = pd.read_csv(f"{output_dir}/dt_dose_df.csv")
    # # vacant_data = '0000-00-00TNN:NN'
    # ## ACTING 기록 개별 분리작업
    #
    # final_dose_df = list()
    # cur_id = ''
    # # df.columns
    # for inx, row in df.iterrows(): #break
    #     if cur_id!=row['UID']:
    #         print(f"({inx} / {len(df)}) {row['UID']} / ACTING 기록 개별 분리작업")
    #         cur_id=row['UID']
    #     row_df = pd.DataFrame(columns=['UID','DRUG','INTERVAL','DT_DOSE','ETC_INFO','VIRTUALITY'])
    #     row_df['DT_DOSE'] = row['DT_DOSE'].split('_')
    #     for c in ['UID','DRUG','INTERVAL','ETC_INFO','VIRTUALITY']:
    #         row_df[c] = row[c]
    #     final_dose_df.append(row_df)
    # final_dose_df = pd.concat(final_dose_df, ignore_index=True)
    # final_dose_df['DATE'] = final_dose_df['DT_DOSE'].map(lambda x:x.split('T')[0])
    # final_dose_df['TIME'] = final_dose_df['DT_DOSE'].map(lambda x:'T'+x.split('T')[-1].split('DOSE')[0])
    # final_dose_df['DOSE'] = final_dose_df['DT_DOSE'].map(lambda x:float(x.split('DOSE')[-1]))
    # final_dose_df = final_dose_df[(final_dose_df['DATE']!=vacant_data.split('T')[0])]
    #
    # # 비품용 제외
    # final_dose_df = final_dose_df[~(final_dose_df['ETC_INFO'].map(lambda x:'비품' in x if type(x)!=float else False))].copy()
    #
    #
    # final_dose_df.to_csv(f"{output_dir}/{raw_drugname}_final_dose_df.csv", encoding='utf-8-sig', index=False)

## READ drug adm dates

drug_adm_res = dict()
for drugname in drugname_dict.values(): #break
    df = pd.read_excel(f"{output_dir}/{drugname}_adm_dates.xlsx")
    df['DATE_LIST'] = df['DATE_LIST'].map(lambda x:x.replace("('",'').replace("',)",'').replace("')",'').split("', '"))
    drug_adm_res[drugname] = df.copy()

## Intersection
vpa_df = drug_adm_res['valproate'].copy()
lcm_df = drug_adm_res['lacosamide'].copy()

# vpa_df[vpa_df['UID'].isin(lcm_df['UID'])]
vpa_df = vpa_df[['UID','DATE_LIST']].rename(columns={'DATE_LIST':'VPA_DATES'})
lcm_df = lcm_df[['UID','DATE_LIST']].rename(columns={'DATE_LIST':'LCM_DATES'})

mdf = vpa_df.merge(lcm_df, on=['UID'], how='left')
mdf = mdf[~mdf['LCM_DATES'].isna()].reset_index(drop=True)
# mdf['UID'].drop_duplicates()
mres_df = list()
np_intsec_adm = list()
no_intsec_adm_count = 0
no_single_adm = list()
no_single_adm_count = 0
for inx, row in mdf.iterrows(): #break

    intsec_dates = pd.DataFrame((set(row['VPA_DATES']).intersection(set(row['LCM_DATES']))), columns=['INTSEC_DATE']).sort_values(['INTSEC_DATE'],ignore_index=True)

    print(f"({inx}) / {row['UID']} / {len(intsec_dates)}")

    intsec_dates['DATE_DIFF'] = pd.to_datetime(intsec_dates['INTSEC_DATE']).diff().dt.days.fillna(1)
    wo_inx = intsec_dates.loc[intsec_dates['DATE_DIFF'] >= 5].index.min()
    
    # 5일 이상 간격 벌어진 index point가 없는 경우 -> 전체 데이터 사용
    if np.isnan(wo_inx):
        pass
    else:
        # raise ValueError
        intsec_dates= intsec_dates[intsec_dates.index < wo_inx].copy()

    if len(intsec_dates)==0:
        no_intsec_adm_count += 1
        np_intsec_adm.append(row['UID'])
        print(f"({no_intsec_adm_count}) / {row['UID']} / No single administration")
        continue

    min_vpa_date = min(row['VPA_DATES'])
    min_lcm_date = min(row['LCM_DATES'])
    min_intsec_date = intsec_dates['INTSEC_DATE'].min()
    max_intsec_date = intsec_dates['INTSEC_DATE'].max()

    if min_vpa_date==min_intsec_date:
        no_single_adm_count+=1
        single_drug = 'none'
        # ds_dates = pd.Series([None])
        single_drug_dates = list()
        no_single_adm.append(row['UID'])
        print(f"({no_single_adm_count}) / {row['UID']} / No single administration")
        continue
    elif min_vpa_date < min_intsec_date:
        single_drug = 'valproate'
        ds_dates = pd.Series(row['VPA_DATES'])
        single_drug_dates = ds_dates[ds_dates < min_intsec_date]
    elif min_lcm_date < min_intsec_date:
        single_drug = 'lacosamide'
        ds_dates = pd.Series(row['LCM_DATES'])
        single_drug_dates = ds_dates[ds_dates < min_intsec_date]

    single_drug_dates_diff = pd.to_datetime(single_drug_dates).diff().dt.days.fillna(1)
    # single_drug_dates.index <
    swo_inx = single_drug_dates_diff[single_drug_dates_diff >= 5].index.min()
    if np.isnan(swo_inx):
        pass
    else:
        # raise ValueError
        single_drug_dates = single_drug_dates[single_drug_dates.index < swo_inx].copy()

    # wo_inx = single_drug_dates.loc[single_drug_dates['DATE_DIFF'] >= 5].index.min()
    mres_dict = {'UID':row['UID'],'SINGLE_DRUG':single_drug,'SINGLE_DAYS':len(single_drug_dates),'INTSEC_DAYS':len(intsec_dates['INTSEC_DATE']),"MIN_SINGLE_DATE":min(single_drug_dates),"MAX_SINGLE_DATE":max(single_drug_dates),'SINGLE_DATES':list(single_drug_dates),'MIN_INTSEC_DATE':min_intsec_date,'MAX_INTSEC_DATE':max_intsec_date,'INTSEC_DATES':list(intsec_dates['INTSEC_DATE'])}
    mres_df.append(mres_dict)

mres_df = pd.DataFrame(mres_df)

# mres_df['SINGLE_DAYS'].mean()
# mres_df['INTSEC_DAYS'].mean()

# mres_df['SINGLE_DAYS'].hist()
# mres_df['INTSEC_DAYS'].hist()

print(f"Mean Days (Single: {round(mres_df['SINGLE_DAYS'].mean(),1)},Intsec: {round(mres_df['INTSEC_DAYS'].mean(),1)}) / No Intsec: {no_intsec_adm_count} / No Single: {no_single_adm_count}")
mres_df.to_excel(f"{output_dir}/merged_adm_dates.xlsx", encoding='utf-8-sig', index=False)
