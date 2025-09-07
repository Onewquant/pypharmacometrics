from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr



prj_name = 'LNZ'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

raw_df = pd.read_csv(f"{resource_dir}/lnz_cdw_dose_data.csv")
# df.columns
raw_df = raw_df.rename(columns={'환자번호':'UID','수행시간':'ACTING','약품 오더일자':'DATE', '[실처방] 용법':'REGIMEN','[실처방] 처방일수':'DAYS', '[함량단위환산] 1일 처방량':'DAILY_DOSE','[실처방] 1회 처방량':'DOSE', '[실처방] 투약위치':'PLACE',"[실처방] 경로":'ROUTE','약품명(성분명)':'DRUG','[실처방] 처방비고':'ETC_INFO'})
raw_df['DOSE'] = raw_df.apply(lambda x:float(re.findall(r'\d+[\.]?\d*',x['DOSE'])[0]) if ('mg' in x['DOSE']) else (float(re.findall(r'\d+[\.]?\d*',x['DOSE'])[0]) * float(re.findall(r'\d+[\.]?\d*',x['약품명(일반명)'])[0])), axis=1)
raw_df['DOSE'] = raw_df['DOSE'].astype(str)

raw_df['UID'] = raw_df['UID'].map(lambda x:x.split('-')[0])
# df.columns
# df['DRUG'].unique()
drugname_dict = {'linezolid':'linezolid'}
for raw_drugname, drugname in drugname_dict.items():

    df = raw_df[raw_df['DRUG']==raw_drugname].copy()
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
            print(f"({count}) {row['UID']}")
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


    dose_cond1 = ~(df['ACTING'].isna())                         # 본원투약기록 -> 투약 기록대로 추가
    dose_cond2 = (df['ADDL']>1)&(df['ACTING'].isna())           # 외래자가투약 / 타병원 투약 처방 -> 날짜대로 ALTER_ACTING으로 추가
    # dose_cond3 = (df['ACTING'].isna())

    df1 = df[dose_cond1].sort_values(['UID','DATE']).reset_index(drop=True)
    df1['VIRTUALITY'] = False
    df2 = df[dose_cond2].sort_values(['UID','DATE']).reset_index(drop=True)
    df2['ACTING'] = df2['ALTER_ACTING'].copy()
    df2['VIRTUALITY'] = True

    # df[dose_cond3]

    # df = pd.concat([df1, df2]).sort_values(['UID','DATE','ACTING'])
    # df1[df1['DAYS']>1][['UID','DATE','INTERVAL','ADDL','ACTING']]

    # ADDL 먼저 추가 작업해야함

    df2["DATE"] = pd.to_datetime(df2["DATE"], errors="coerce")

    # raise ValueError

    addl_applied_df2 = list()
    for inx, row in df2.iterrows():#break
        interval_num = float(row['INTERVAL'].replace('q','').replace('h',''))
        print(f"({inx}) {row['UID']} / {row['ADDL']} ADDL")
        if interval_num > 24:
            day_step = interval_num/24
            # raise ValueError
            # pass
        else:
            day_step = 1

        # end_date = day_step * row["ADDL"]

        for i in range(row["ADDL"]):  # ADDL 포함해서 +1
            new_row = row.copy()
            new_row["DATE"] = (row["DATE"] + pd.Timedelta(days=day_step*i)).strftime("%Y-%m-%d")
            new_row["ADDL"] = 1
            addl_applied_df2.append(new_row)

    df2 = pd.DataFrame(addl_applied_df2)

    #######################################

    df = pd.concat([df1,df2]).sort_values(['UID','DATE']).reset_index(drop=True)

    # DATETIME/DOSE 연결 작업

    dt_dose_series = list()
    vacant_data = '0000-00-00TNN:NN'
    for inx, row in df.iterrows():
        # raise ValueError

        # if (row['ID']=='10023985')&(row['DATETIME']=='2004-04-19T06:51'):
        #     raise ValueError
        # else:
        #     continue

        # if (row['DATETIME']=='2005-06-27T20:00'):
        #     raise ValueError
        # else:
        #     continue

        rep_acting_str = row['ACTING'].replace('O.C,','').strip()

        new_actval_str=''
        for actval in rep_acting_str.split(','): #break

            ## 빈칸일때
            if actval=='':
                # new_actval_str += f'_{vacant_data}'
                continue


            ## 날짜 시간/Y일때
            y_val = re.findall(r"\d\d\d\d-\d\d-\d\d \d\d:\d\d/Y",actval)
            if len(y_val) > 0:
                new_actval_str+='_' + y_val[0].replace("/Y","").replace(" ","T")
                continue
            else:
                pass

            ## 시간/Y 일떄
            y_val = re.findall(r"\d\d:\d\d/Y",actval)
            if len(y_val) > 0:
                new_actval_str+=f'_{row["DATE"]}T{y_val[0].replace("/Y","")}'
                continue
            else:
                rest_y_val = actval.split('/')[-1]

                """
                # 확인 필요한 부분: D나 N은 투약이 된것인지?
                """
                if (('Y ' in rest_y_val) or ('O ' in rest_y_val)):
                    raise ValueError
                else:
                    new_actval_str += f'_{vacant_data}'
                    continue

            raise ValueError

        new_actval_str = new_actval_str[1:]
        # print(22)
        # raise ValueError
        # DOSE 붙이기 작업
        new_actval_split = new_actval_str.split('_')


        new_actval_str = '_'.join([f"{nav}DOSE{row['DOSE']}" for nav in new_actval_split])
        dt_dose_series.append(new_actval_str)
    df['DT_DOSE'] = dt_dose_series

    df.to_csv(f"{output_dir}/{raw_drugname}_dt_dose_df.csv", encoding='utf-8-sig', index=False)
    # df[['UID', 'DATE', 'DRUG', 'DOSE', 'ROUTE','REGIMEN','DAYS', 'ETC_INFO', 'ACTING', 'ALTER_ACTING', 'INTERVAL', 'VIRTUALITY', 'DT_DOSE']]
    # df['INTERVAL']
    # dose_result_df[['DATE','DOSE','ACTING','DT_DOSE']]

    #
    # dose_result_df = pd.read_csv(f"{output_dir}/dt_dose_df.csv")
    # vacant_data = '0000-00-00TNN:NN'
    ## ACTING 기록 개별 분리작업

    final_dose_df = list()
    cur_id = ''
    # df.columns
    for inx, row in df.iterrows(): #break
        if cur_id!=row['UID']:
            print(f"({inx} / {len(df)}) {row['UID']} / ACTING 기록 개별 분리작업")
            cur_id=row['UID']
        row_df = pd.DataFrame(columns=['UID','DRUG','INTERVAL','DT_DOSE','ETC_INFO','VIRTUALITY'])
        row_df['DT_DOSE'] = row['DT_DOSE'].split('_')
        for c in ['UID','DRUG','INTERVAL','ETC_INFO','VIRTUALITY']:
            row_df[c] = row[c]
        final_dose_df.append(row_df)
    final_dose_df = pd.concat(final_dose_df, ignore_index=True)
    final_dose_df['DATE'] = final_dose_df['DT_DOSE'].map(lambda x:x.split('T')[0])
    final_dose_df['TIME'] = final_dose_df['DT_DOSE'].map(lambda x:'T'+x.split('T')[-1].split('DOSE')[0])
    final_dose_df['DOSE'] = final_dose_df['DT_DOSE'].map(lambda x:float(x.split('DOSE')[-1]))
    final_dose_df = final_dose_df[(final_dose_df['DATE']!=vacant_data.split('T')[0])]

    # 비품용 제외
    final_dose_df = final_dose_df[~(final_dose_df['ETC_INFO'].map(lambda x:'비품' in x if type(x)!=float else False))].copy()


    final_dose_df.to_csv(f"{output_dir}/lnz_final_dose_df.csv", encoding='utf-8-sig', index=False)

#

# df2
# df['REGIMEN'].unique()
# df['ALTER_ACTING'] = ''




# odf['ROUTE'].unique()
# odf[odf['ROUTE']=='IV'][['UID','오더일','DAYS','ROUTE','[실처방] 투약위치','ACTING']]
# df['ROUTE'].unique()
# df['ETC_INFO'].unique()






