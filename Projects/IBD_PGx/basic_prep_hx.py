from tools import *
from pynca.tools import *

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'IBD_PGX'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Orders
drug_order_set = set()

hx_files = glob.glob(f'{resource_dir}/hx/IBD_PGx_hx(*).txt')
dose_result_df = list()
wierd_result_df = list()
# result_cols = ['ID','NAME','DATETIME','DRUG','DOSE','ROUTE','ACTING','PERIOD','ETC_INFO','PLACE']
# no_dup_cols = [c for c in result_cols if c!='NAME']
for finx, fpath in enumerate(hx_files): break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    print(f"({finx}) {pname} / {pid}")
    # if pname=='이학준':
    #     raise ValueError

    # if pid in ("15322168", "19739357", "34835292", "37366865", "21618097", "36898756", "36975211", "37858047"):       # lab, order 파일 다시 수집 필요
    #     continue

    with open(fpath, "r", encoding="utf-8") as f:
        f_content = f.read()
        # print(content)
    for fct_inx, fct_str in enumerate([fct.split('\n작성자') for fct in f_content.split('작성과: ')[1:]]):
        'Partial Mayo Score for Ulcerative Colitis'
        'Mayo Score'
        break



    fdf = pd.read_excel(fpath)

    # fdf.columns
    # fdf.to_csv(f"{outcome_dir}/error_dose_df.csv", encoding='utf-8-sig', index=False)

    dose_df = fdf[fdf['처방지시'].map(lambda x: (('adalimumab' in x.lower()) or ('infliximab' in x.lower()) or ('ustekinumab' in x.lower())) and ('quantification' not in x.lower()))].copy()
    # dose_df.iloc[0]['처방지시']
    # dose_df.iloc[1]['처방지시']
    # dose_df.iloc[2]

    dose_df = dose_df[(~dose_df['Acting'].isna())].copy()

    dose_df['처방지시비고'] = dose_df['처방지시'].map(lambda x:x.split(' : ')[-1] if len(x.split(' : '))>1 else '')

    dose_df['ID'] = pid
    dose_df['NAME'] = pname

    # if pname=='김옥순': raise ValueError

    #### 약국_검사가 NA인 것은 제외하고 진행함 (추후 필요시 추가)

    # dose_df[dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test.xlsx",index=False)
    # dose_df[~dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test2.xlsx",index=False)
    dose_df = dose_df[(~dose_df['약국_검사'].isna())].copy()

    #### 성분명이 IBD Biologics 인 경우만으로 필터링 (Infliximab, Adalimumab)

    regex_pattern = r'\(infliximab|\(adalimumab|\(ustekinumab'
    dose_df = dose_df[dose_df['처방지시'].map(lambda x: bool(re.search(regex_pattern, x, flags=re.IGNORECASE)))].copy()

    dose_df['DT1'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[0].split('[')[-1].replace(' ','T'))
    dose_df['DT2'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[-1].split('[')[-1].replace(' ','T')[:-3])
    dose_df['DATETIME'] = dose_df[['DT1','DT2']].min(axis=1)
    dose_df['ETC_INFO'] = dose_df['처방지시비고'].copy()
    dose_df['DRUG'] = dose_df['처방지시'].map(lambda x: re.search(regex_pattern, x, flags=re.IGNORECASE).group().lower().replace('(','').replace(')',''))
    dose_df['ROUTE'] = dose_df['처방지시'].map(lambda x: 'IV' if " [SC] " not in x else 'SC')
    dose_df['PLACE'] = dose_df['주사시행처']
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+',x.split('▣')[-1].split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split(' ')[-1].strip())[0] if " [SC] " not in x else x.split('(Infliximab)')[-1].split('(Adalimumab)')[-1].split('▣')[-1].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+', x.split('mg')[1].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0] if (" [SC] " not in x) else x.split('(Infliximab')[-1].split('(Adalimumab')[-1].split('(Ustekinumab')[-1].split('▣')[-1].split('srg')[0].split('via')[0].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # mg_inx_dict = {'ustekinumab':0,'infliximab':1,'adalimumab':1}
    dose_series = list()
    for dose_inx, dose_row in dose_df.iterrows():
        x = dose_row['처방지시']
        if (" [SC] " not in x):
            # raise ValueError
            try: dose_val = int(re.findall(r'\d+mg',x)[1].replace('mg','').strip())
            except:
                dose_val = int(re.findall(r'\d+mg', x)[0].replace('mg', '').strip()) * int(re.findall(r'\d+ via', x)[0].replace('via', '').strip())
            # dose_val = re.findall(r'\d+',x.split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0]
            # raise ValueError
            dose_series.append(dose_val)
        else:
        #     if ('infliximab' in x.lower()) and (not ('Remsima 120mg pen (Infliximab) ▣ 120mg [SC]' in x)) and (not ('Remsima 120mg pen (Infliximab) 120mg [SC]' in x)):
        #         raise ValueError
        #     else:continue
        # continue
        # pid
            try: dose_val = int(re.findall(r'\d+mg',x)[1].replace('mg','').strip())
            except:
                if '(Ustekinumab' in x:
                    # raise ValueError
                    dose_val = int(re.findall(r'\d+mg',x)[0].replace('mg','').strip()) * int(re.findall(r'\d+ srg',x)[0].replace('srg','').strip())
                elif '(Adalimumab' in x:
                    # raise ValueError
                    # if 'ml pen' in x:
                    #     x = x.replace('ml pen','ml 1 pen')
                    try: pen_count = int(re.findall(r'\d+ pen',x)[0].replace('pen','').strip())
                    except: pen_count = int(re.findall(r'\d+ srg',x)[0].replace('srg','').strip())
                    dose_val = int(re.findall(r'\d+mg',x)[0].replace('mg','').strip()) * pen_count
                    # raise ValueError
                else:
                    raise ValueError
                    dose_val = x.split('(Infliximab')[-1].split('(Adalimumab')[-1].split('(Ustekinumab')[-1].split('▣')[-1].split('srg')[0].split('via')[0].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip()
            dose_series.append(dose_val)
    # raise ValueError
    dose_df['DOSE'] = dose_series
    # dose_df['DOSE'] = dose_df['DOSE'].map({'1 pen': 40, '2 pen': 80, '2 pen': 160})
    # (1) [원내] Remsima 100mg inj (Infliximab Korea) ...
    dose_df['ACTING'] = dose_df['Acting']

    dose_df['PERIOD'] = dose_df['처방지시'].map(lambda x: 'x1' if " [SC] " not in x else x.split(' [SC] ')[-1].split(':')[0].strip())
    # dose_df.to_csv(f"{output_dir}/dose_df_lhj.csv", encoding='utf-8-sig', index=False)
    # dose_df.loc[3203,'처방지시']

    dose_result_df.append(dose_df[result_cols].drop_duplicates(no_dup_cols))

    # drug_order_set = drug_order_set.union(set(dose_df['처방지시'].map(lambda x:''.join(x.split(':')[0].replace('  ',' ').split(') ')[1:]).replace('[원내]','').replace('[D/C]','').replace('[보류]','').replace('[반납]','').replace('[Em] ','').strip()).drop_duplicates()))


dose_result_df = pd.concat(dose_result_df, ignore_index=True).sort_values(['ID','DATETIME'], ascending=[True,False])


## DATETIME 추가 조정
dt_series = list()
for inx, row in dose_result_df.iterrows():
    if 'Y' in row['ACTING']:
        ymdt_patterns = re.findall(r"[\d][\d][\d][\d]-[\d][\d]-[\d][\d] [\d][\d]:[\d][\d]",row['ACTING'])
        if len(ymdt_patterns)>0:
            new_dt = ymdt_patterns[0].replace(' ','T')
        else:
            time_patterns = re.findall(r"[\d][\d]:[\d][\d]",row['ACTING'])
            new_dt = row['DATETIME'].split('T')[0] + 'T' + time_patterns[0]
        dt_series.append(new_dt)
    else:
        dt_series.append(row['DATETIME'])

dose_result_df['DATETIME'] = dt_series

## PERIOD 추가 조정

dose_result_df['PERIOD'] = dose_result_df['PERIOD'].map(lambda x:x.replace('2주 간격 유지용량 ','').replace('#1 ','').replace('3/20 ','').replace('매주 ','').replace(' X1 Day','').replace('1회 ','').replace('/wk ','/1wks').replace('x1 ','1'))


## H, C, Z는 제거

dose_result_df = dose_result_df[dose_result_df['ACTING'].map(lambda x: False if (('H' in x) or ('Z' in x) or ('C' in x) or ('N' in x)) else True)].reset_index(drop=True)
dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)
dose_result_df['ETC_INFO_TREATED'] = ''

## 처방비고 반영 및 기타 정리

# 2주 뒤에 맞기 반영

ETC_INFO_cond = (dose_result_df['ETC_INFO']=='2주 뒤에 자가로 2개 맞으세요.')
PERIOD_cond = list()
PERIOD_added_list1 = ['2주 뒤에 맞으세요 ','2주 뒤에 자가 접종 ','2주 뒤 ','2주 뒤에 자가로 2개 맞으세요. '] + ["오늘 2amp, 내일 2mp 맞으세요 "]
for inx, period_str in enumerate(dose_result_df['PERIOD']):
    period_tf = False
    for comment in PERIOD_added_list1:
        if comment in period_str:
            period_tf=True
    PERIOD_cond.append(period_tf)
PERIOD_cond = pd.Series(PERIOD_cond)

change_df = dose_result_df[ETC_INFO_cond|PERIOD_cond].copy()
for inx, row in change_df.iterrows():
    # dose_result_df.iloc[inx]
    dose_result_df.at[inx,'DATETIME'] = (datetime.strptime(dose_result_df.at[inx,'DATETIME'],'%Y-%m-%dT%H:%M') + timedelta(days=14)).strftime('%Y-%m-%dT%H:%M')
    # dose_result_df.at[inx,'ETC_INFO'] = '처방비고 반영완료'
    dose_result_df.at[inx,'ETC_INFO_TREATED'] = '처방비고 반영완료'
dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)

# Dose 나눠서 이틀에 맞기

ETC_INFO_cond = dose_result_df['ETC_INFO'].isin(['금일 주사실에서 2개 맞고 가고 그 다음날 자가로 2개 더 맞으세요.','금일 두 개, 내일 두 개 맞으세요.','하루에 2개씩 이틀에 나눠 맞으세요','금일 2개, 내일 2개 맞으세요.','금일 2개, 내일 2개','금일 주사실에서 2개 맞고 가고 그 다음날 자가로 2개 더 맞으세요. '])
PERIOD_cond = list()
PERIOD_added_list2 = ['오늘 2개, 내일 2개 맞으세요 ', '금일 2개, 내일 2개 맞으세요 ', '금일 2개, 내일 2개 맞으세요. ','주사실에서 2개 맞고 가고 그 다음날 자가로 2개 더 맞으세요. ']
for inx, period_str in enumerate(dose_result_df['PERIOD']):
    period_tf = False
    for comment in PERIOD_added_list2:
        if comment in period_str:
            period_tf=True
    PERIOD_cond.append(period_tf)
PERIOD_cond = pd.Series(PERIOD_cond)

change_df = dose_result_df[ETC_INFO_cond|PERIOD_cond].copy()
change_df['DATETIME'] = change_df['DATETIME'].map(lambda x:(datetime.strptime(x, '%Y-%m-%dT%H:%M') + timedelta(days=1)).strftime('%Y-%m-%dT%H:%M'))
change_df['DOSE'] = (change_df['DOSE']/2).map(int)
change_df['ETC_INFO_TREATED'] = '처방비고 반영완료'
for inx, row in change_df.iterrows():
    # dose_result_df.iloc[inx]
    dose_result_df.at[inx, 'DOSE'] = int(dose_result_df.at[inx, 'DOSE']/2)
    # dose_result_df.at[inx,'ETC_INFO'] = '처방비고 반영완료'
    dose_result_df.at[inx,'ETC_INFO_TREATED'] = '처방비고 반영완료'
dose_result_df = pd.concat([dose_result_df,change_df]).reset_index(drop=True)

# Period에 들어있는 추가 글자 제거

for c in PERIOD_added_list1 + PERIOD_added_list2 + [' 컨펌 후 투여',' 맞던 날짜대로',' [D]','금일 ', '매 주','3주 간격 ']:
    dose_result_df['PERIOD'] = dose_result_df['PERIOD'].map(lambda x:x.replace(c,''))
for c in ['1주 1회','ut dict ','ut dict','prn ']:
    dose_result_df['PERIOD'] = dose_result_df['PERIOD'].replace(c,'x1')

# ROUTE가 IV인데 ACTING에 Y가 안 들어 있는 것 제거

dose_result_df = dose_result_df[~((dose_result_df['ROUTE']=='IV')&(dose_result_df['ACTING'].map(lambda x:'Y' not in x)))].copy()

# ADDL, II 추가 및 다음 DOSE 보다 넘어가는 ADDL 확인

dose_result_df['PERIOD'] = dose_result_df['PERIOD'].map(lambda x:x.replace('  ',' ').replace('X',' X').replace('  ',' '))
dose_result_df['ADDL'] = dose_result_df['PERIOD'].map(lambda x: int(x.split('/')[0].strip()) * int(x.split('X')[-1].split('Week')[0].strip())/ int(x.split('wks')[0].split('/')[-1].strip())-1 if len(re.findall(r'\d+wks',x)) > 0 else 0).map(int)
dose_result_df['II'] = dose_result_df['PERIOD'].map(lambda x: int(x.split('wks')[0].split('/')[-1].strip()) / int(x.split('/')[0].strip())*7*24 if len(re.findall(r'\d+wks',x)) > 0 else 2*7*24).map(int)
dose_result_df['NXTLST_DOSE_DT'] = dose_result_df.apply(lambda x: (datetime.strptime(x['DATETIME'],'%Y-%m-%dT%H:%M')+timedelta(x['ADDL']*x['II']/24)).strftime('%Y-%m-%dT%H:%M'), axis=1)

# ADDL을 하나씩 다 넣어보자
addl_dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True).copy()
addl_added_df = list()
for inx, row in addl_dose_result_df[addl_dose_result_df['ADDL']>0].iterrows():
    # dose_result_df[dose_result_df.index==inx]
    # if row['DATETIME']=='2021-07-19T12:05':
        # raise ValueError
    max_addl = row['ADDL']
    init_datetime = row['DATETIME']
    addl_frag_df = list()
    for addl_num in range(1,row['ADDL']+1):
        row['ETC_INFO_TREATED'] = f'ADDL반영_{addl_num}'
        row['PERIOD'] = 'x1'
        row['ADDL']=0
        row['DATETIME'] = (datetime.strptime(init_datetime,'%Y-%m-%dT%H:%M')+timedelta(addl_num*row['II']/24)).strftime('%Y-%m-%dT%H:%M')
        addl_frag_df.append(pd.DataFrame([row]))

    addl_added_df.append(pd.concat(addl_frag_df))
    addl_dose_result_df.at[inx,'ADDL']=0
    addl_dose_result_df.at[inx,'ETC_INFO_TREATED'] = f"ADDL존재_{max_addl}"
    # dose_result_df.at[inx,'PERIOD'] = 'x1'

addl_added_df = pd.concat(addl_added_df)


# 요일 정리

dose_result_df['DAY_OF_WK'] = dose_result_df['DATETIME'].map(lambda x: (datetime.strptime(x,'%Y-%m-%dT%H:%M').strftime('%w'))).map({'0':'일','1':'월','2':'화','3':'수','4':'목','5':'금','6':'토'})
addl_dose_result_df['DAY_OF_WK'] = addl_dose_result_df['DATETIME'].map(lambda x: (datetime.strptime(x,'%Y-%m-%dT%H:%M').strftime('%w'))).map({'0':'일','1':'월','2':'화','3':'수','4':'목','5':'금','6':'토'})


# 정리 및 저장

dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True)
dose_result_df.to_csv(f"{output_dir}/dose_df(addl_col).csv", encoding='utf-8-sig', index=False)

addl_dose_result_df = pd.concat([addl_dose_result_df, addl_added_df]).sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True)
addl_dose_result_df.to_csv(f"{output_dir}/dose_df.csv", encoding='utf-8-sig', index=False)

# ADDL시 다음 투약 날짜보다 큰 사람 확인

addl_overflow_df = list()
for inx, frag_df in dose_result_df[['ID','NAME','DATETIME', 'NXTLST_DOSE_DT', 'DRUG','PLACE']].groupby('ID'): #break
    frag_df['SHIFT_NXTLST_DT'] = frag_df['NXTLST_DOSE_DT'].shift(1).fillna('0001-01-01T00:00')
    # frag_df['SHIFT_NXTLST_DT']
    frag_df['DELT_DT'] = frag_df.apply(lambda x: (datetime.strptime(x['DATETIME'],'%Y-%m-%dT%H:%M')-datetime.strptime(x['SHIFT_NXTLST_DT'],'%Y-%m-%dT%H:%M')).days, axis=1)
    addl_overflow_frag = frag_df[frag_df['DELT_DT'] <= 0]
    if len(addl_overflow_frag)>0:
        # raise ValueError
        # addl_overflow_df.append(frag_df[frag_df.index==(addl_overflow_frag.iloc[0].name-1)])
        addl_overflow_df.append(addl_overflow_frag)
addl_overflow_df = pd.concat(addl_overflow_df, ignore_index=True)
addl_overflow_df.to_csv(f"{output_dir}/addl_overflow_dose_df.csv", encoding='utf-8-sig', index=False)

"""
# (1) 보통 아래 예시처럼 외래에서 처방을 받게 되면, 당일 주사는 병원에서 맞고 가는지? 아니면 약만 처방받아 가는 것인지? 
처방 받는 당일 자가로 투약했다고 볼 수 있을지?

18836363	장--	2023-10-02T10:07	infliximab	120	SC	09:00/, 	1/2wks X14 Weeks		자가

# (2) 아래 예시에서, 첫 4회째 투약 시점(2021-05-05) 부근에 외래에서 추가 4회 처방을 받게 되는데, 각 오더 낸 시점에 투여 받고, 이후 
각 투여 간격별로 투약했다고(즉, 아래 예시에서 첫 번째 row의 오더는 03/24, 04/07, 04/21, 05/05 에 투약했다고) 볼 수 있을지? 그런데,
첫 번째 오더의 마지막 투약(05/05) 시점 부근에 다시 외래를 방문하시는데, 이때 또 투약이 된 것인지? 

17677819	김--	2021-03-24T15:18	adalimumab	40	SC	09:00/, 	1/2wks X8 Weeks	2021-05-05T15:18	자가
17677819	김--	2021-05-06T11:04	adalimumab	40	SC	09:00/, 	1/1wks X4 Weeks	2021-05-27T11:04	자가

# (3) 아래 예시에서 22/12/19 이후 대체로 월요일날 처방 오더를 받았는데, 24/05/02 부터는 목요일에 처방 오더를 받습니다. 이때는 
새로 받은 약이라도 투약은 계속 월요일마다 진행되었다고 볼 수 있는 것일지?

10933347	장--	2022-10-29T19:34	토	infliximab	300	IV	19:34/Y, 	x1	2022-10-29T19:34
10933347	장--	2022-12-19T11:25	월	infliximab	120	SC	11:25/Y, 	1/2wks X12 Weeks	2023-02-27T11:25
10933347	장--	2023-03-13T09:35	월	infliximab	120	SC	09:00/, 	1/2wks X12 Weeks	2023-05-22T09:35
10933347	장--	2023-06-05T09:14	월	infliximab	120	SC	09:00/, 	1/2wks X12 Weeks	2023-08-14T09:14
10933347	장--	2023-08-28T09:17	월	infliximab	120	SC	09:00/, 	1/2wks 월요일 X12 Weeks	2023-11-06T09:17
10933347	장--	2023-11-20T09:29	월	infliximab	120	SC	09:00/, 	1/2wks 월요일 X12 Weeks	2024-01-29T09:29
10933347	장--	2024-02-05T10:13	월	infliximab	120	SC	09:00/, 	1/2wks 월요일 X12 Weeks	2024-04-15T10:13
10933347	장--	2024-05-02T09:25	목	infliximab	120	SC	09:00/, 	1/2wks 월요일 X12 Weeks	2024-07-11T09:25
10933347	장--	2024-07-25T09:19	목	infliximab	120	SC	09:00/, 	1/2wks 월요일 X12 Weeks	2024-10-03T09:19
10933347	장--	2024-10-17T09:37	목	infliximab	120	SC	09:00/, 	1/2wks 월요일 X12 Weeks	2024-12-26T09:37
10933347	장--	2025-01-11T09:29	토	infliximab	120	SC	09:00/, 	1/2wks 월요일 X12 Weeks	2025-03-22T09:29
10933347	장--	2025-03-11T14:03	화	infliximab	120	SC	09:00/, 	1/2wks 월요일 X12 Weeks	2025-05-20T14:03


"""


    # break
# len(addl_overflow_df)
# dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)



# 남은 것: 월요일 투약 처리하기 / #2 ~ / 각 사람들 투약 간격 잘 맞는지 눈으로 확인(당일 맞는것 포함해서 약을 주는 건지)

# dose_result_df.drop_duplicates(['ID'], ignore_index=True)
# ot_list = list()
# for inx_ot, order_text in enumerate(drug_order_set):
#     if '[SC]' not in order_text:
#         ot_list.append(order_text)
#         continue
#
# for inx_ot, order_text in enumerate(ot_list):
#     print(f"({inx_ot}) {order_text}")
# len(order_text)

# ['비고','처방지시', '발행처', '발행의', '수납', '약국/검사', '주사시행처', 'Acting', '변경의']
#
# len(order_files)
#
#
# df = pd.read_csv(f'{resource_dir}/glpharma_CONC.csv')
# seq_df = pd.read_csv(f'{prj_dir}/glpharma_SEQUENCE.csv')

# {'adlimumab SC 1 pen': 40, }