from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

prj_name = 'CFZ'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/CFZ'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

raw_df = pd.read_csv(f"{resource_dir}/cfz_cdw_order_data.csv")
# df.columns
raw_df = raw_df.rename(columns={'환자번호':'UID','수행시간':'ACTING','약품 오더일자':'DATE', '[실처방] 용법':'REGIMEN','[실처방] 처방일수':'DAYS', '[함량단위환산] 1일 처방량':'DAILY_DOSE','[함량단위환산] 1회 처방량':'DOSE', '[실처방] 투약위치':'PLACE',"[실처방] 경로":'ROUTE','약품명(성분명)':'DRUG','[실처방] 처방비고':'ETC_INFO'})
raw_df['UID'] = raw_df['UID'].map(lambda x:x.split('-')[0])
# df.columns
# df['DRUG'].unique()
drugname_dict = {'cefazolin':'cefazolin','furosemide':'furosemide'}
# drugname_dict = {'furosemide':'furosemide'}


drug_adm_res = dict()
for raw_drugname, drugname in drugname_dict.items():

    drug_res_df = list()

    df = raw_df[raw_df['DRUG']==raw_drugname].copy()
    df = df[~(df['ROUTE'].isin(['Irrigation','L-spine']))].copy()
    # df = df[~(df['DAYS'].isin(['관류하세요', '관류하세요.']))].copy()
    df = df[~(df['REGIMEN'].map(lambda x:x.strip()).isin(['관류하세요']))].copy()
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
        '1일 7회':('q3.5h','02:00/Y, 05:00/Y, 08:00/Y, 11:00/Y, 14:00/Y, 18:00/Y, 22:00/Y'),
        '1일 8회':('q3h','02:00/Y, 05:00/Y, 08:00/Y, 11:00/Y, 14:00/Y, 17:00/Y, 20:00/Y, 23:00/Y'),
        '1일 9회':('q3h','02:00/Y, 05:00/Y, 08:00/Y, 11:00/Y, 14:00/Y, 17:00/Y, 20:00/Y, 22:00/Y, 23:00/Y'),
        '1일 10회':('q2.4h','02:00/Y, 04:00/Y, 05:00/Y, 08:00/Y, 11:00/Y, 14:00/Y, 17:00/Y, 19:00/Y, 21:00/Y, 23:00/Y'),
        '1일 12회':('q2.4h','02:00/Y, 04:00/Y, 05:00/Y, 08:00/Y, 10:00/Y, 11:00/Y, 14:00/Y, 16:00/Y, 17:00/Y, 19:00/Y, 21:00/Y, 23:00/Y'),
        '필요시 자기전에':('q[hours]h','22:00/Y'),
        '필요시 복용하세요':('q[hours]h','09:00/Y'),
        '필요시 1시간전에': ('q[hours]h', '09:00/Y'),
        '격일로 자기전':('q48h','22:00/Y'),
        '격일로 저녁':('q48h','21:00/Y'),
        '격일로 아침':('q48h','09:00/Y'),
        '격일로 48시간마다':('q48h','09:00/Y'),
        '1년 1회':('q[hours]h','09:00/Y'),
        '1주 3회':('q56h','09:00/Y'),
        '1주 1회':('q168h','09:00/Y'),
        '1주 2회':('q84h','09:00/Y'),
        '3일 1회':('q72h','09:00/Y'),
        '의사지시대로 관류하세요':('q[hours]h','09:00/Y'),
        '필요시 주사하세요':('q24h','09:00/Y'),
        '넣으세요':('q24h','09:00/Y'),
        '수시로 복용하세요':('q24h','09:00/Y'),
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
    df['ROUTE'] = df['ROUTE'].map({'P.O':'PO','MIV':'IV','PLT':'PO','IVS':'IV','SC':'SC','S.L':"SL", 'IntraPleural':'IPLEU', 'IntraPeritoneal':'IPERI'}) # PLT: ( Par L-tube) / S.L: sublingual
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




## READ drug adm dates

drug_adm_res = dict()
for drugname in drugname_dict.values(): #break
    df = pd.read_excel(f"{output_dir}/{drugname}_adm_dates.xlsx")
    df['DATE_LIST'] = df['DATE_LIST'].map(lambda x:x.replace("('",'').replace("',)",'').replace("')",'').split("', '"))
    drug_adm_res[drugname] = df.copy()

## Intersection
cfz_df = drug_adm_res['cefazolin'].copy()
fur_df = drug_adm_res['furosemide'].copy()

# vpa_df[vpa_df['UID'].isin(lcm_df['UID'])]
cfz_df = cfz_df[['UID','DATE_LIST']].rename(columns={'DATE_LIST':'CFZ_DATES'})
fur_df = fur_df[['UID','DATE_LIST']].rename(columns={'DATE_LIST':'FUR_DATES'})

mdf = cfz_df.merge(fur_df, on=['UID'], how='left')
mdf = mdf[~mdf['FUR_DATES'].isna()].reset_index(drop=True)
# mdf['UID'].drop_duplicates()
mres_df = list()
np_intsec_adm = list()
no_intsec_adm_count = 0
no_single_adm = list()
no_single_adm_count = 0
for inx, row in mdf.iterrows(): #break

    intsec_dates = pd.DataFrame((set(row['CFZ_DATES']).intersection(set(row['FUR_DATES']))), columns=['INTSEC_DATE']).sort_values(['INTSEC_DATE'],ignore_index=True)

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

    min_vpa_date = min(row['CFZ_DATES'])
    min_lcm_date = min(row['FUR_DATES'])
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
        single_drug = 'cefazolin'
        ds_dates = pd.Series(row['CFZ_DATES'])
        single_drug_dates = ds_dates[ds_dates < min_intsec_date]
    elif min_lcm_date < min_intsec_date:
        single_drug = 'furosemide'
        ds_dates = pd.Series(row['FUR_DATES'])
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
