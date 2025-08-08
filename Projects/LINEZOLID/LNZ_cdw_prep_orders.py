from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr



prj_name = 'LINEZOLID'
prj_dir = './Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

df = pd.read_csv(f"{resource_dir}/lnz_cdw_order_data.csv")
# df.columns
df = df.rename(columns={'환자번호':'UID','수행시간':'ACTING','약품 오더일자':'DATE', '[실처방] 용법':'REGIMEN','[실처방] 처방일수':'DAYS', '[실처방] 1회 처방량':'DOSE', "[실처방] 경로":'ROUTE',"약품명(일반명)":'DRUG','[실처방] 처방비고':'ETC_INFO'})
df['UID'] = df['UID'].map(lambda x:x.split('-')[0])


alter_acting_dict = {
    '1일 2회 12시간마다 주사하세요':('q12h','09:00/Y, 21:00/Y'),
 '1일 1회 24시간마다 주사하세요':('q24h','09:00/Y'),
 '1일 2회 아침,저녁 식후에 복용하세요':('q12h','09:00/Y, 21:00/Y'),
 '1일 2회 12시간마다 복용하세요':('q12h','09:00/Y, 21:00/Y'),
 '1일 1회 아침 식후에 복용하세요':('q24h','09:00/Y'),
 '1일 2회 의사지시대로 복용하세요':('q12h','09:00/Y, 21:00/Y'),
 '1일 1회 24시간마다 복용하세요':('q24h','09:00/Y'),
 '1일 1회 저녁 식후에 복용하세요':('q24h','21:00/Y'),
 '1일 2회 아침,저녁 식사직후에 복용하세요':('q12h','09:00/Y, 21:00/Y'),
 '1일 6회 의사지시대로 주사하세요':('q4h','08:00/Y, 12:00/Y, 16:00/Y, 20:00/Y, 23:59/Y'),
 '1일 1회 아침 식전30분에 복용하세요':('q24h','08:00/Y'),
 '1일 1회 의사지시대로 주사하세요':('q24h','09:00/Y'),
 '1일 1회 자기전에 복용하세요':('q24h','08:00/Y'),
 '1일 4회 6시간마다 주사하세요':('q6h','04:00/Y, 10:00/Y, 16:00/Y, 22:00/Y'),
 '1일 3회 8시간마다 주사하세요':('q8h','08:00/Y, 16:00/Y, 23:59/Y'),
 '1일 1회 아침 식사직후에 복용하세요':('q24h','09:00/Y'),
 '1일 1회 저녁 식사직후에 복용하세요':('q24h','21:00/Y'),
 '1일 3회 8시간마다 복용하세요':('q8h','08:00/Y, 16:00/Y, 23:59/Y')
 }


dose_list = list()
alter_acting_list = list()
interval_list = list()
for inx, row in df.iterrows():

    # if row['DOSE']=='0.5tab':
        # raise ValueError

    if 'mg' in row['DOSE']:
        adm_dose = float(row['DOSE'].split('mg')[0])
    else:
        unit_dose = float(row['DRUG'].split('mg')[0].replace('Linezolid ',''))
        unit_cnt = float(re.findall(r'\d+\.?\d*',row['DOSE'])[0])
        adm_dose = unit_dose * unit_cnt
    dose_list.append(adm_dose)

    alter_acting_list.append(alter_acting_dict[row['REGIMEN']][-1])
    interval_list.append(alter_acting_dict[row['REGIMEN']][0])


df['DOSE'] = dose_list
df['ALTER_ACTING'] = alter_acting_list
# df['REGIMEN_ORI'] = df['REGIMEN'].copy()
df['INTERVAL']= interval_list
# df['DATETIME']

df['DRUG'] = 'Linezolid'
df['ROUTE'] = df['ROUTE'].map({'P.O':'PO','MIV':'IV','PLT':'PO','IVS':'IV'}) # PLT: ( Par L-tube)
df = df[df['ROUTE']!='Irrigation'].copy()
dose_cond1 = ~(df['ACTING'].isna())                         # 본원투약기록 -> 투약 기록대로 추가
dose_cond2 = (df['DAYS']>1)&(df['ACTING'].isna())           # 외래자가투약 / 타병원 투약 처방 -> 날짜대로 ALTER_ACTING으로 추가


df1 = df[dose_cond1].copy()
df1['VIRTUALITY'] = False
df2 = df[dose_cond2].copy()
df2['ACTING'] = df2['ALTER_ACTING'].copy()
df2['VIRTUALITY'] = True
df = pd.concat([df1, df2]).sort_values(['UID','DATE','ACTING'])
# df1[df1['DAYS']>1]['ACTING']

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

df.to_csv(f"{output_dir}/lnz_dt_dose_df.csv", encoding='utf-8-sig', index=False)
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






