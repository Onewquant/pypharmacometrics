from tools import *
from pynca.tools import *

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Orders
drug_order_set = set()

order_files = glob.glob(f'{resource_dir}/order/IBD_PGx_order(*).xlsx')
dose_result_df = list()
wierd_result_df = list()
result_cols = ['ID','NAME','DATETIME','DRUG','DOSE','ROUTE','ACTING','PERIOD','ETC_INFO']
no_dup_cols = [c for c in result_cols if c!='NAME']
for finx, fpath in enumerate(order_files): #break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    print(f"({finx}) {pname} / {pid}")
    # if pname=='이학준':
    #     raise ValueError

    # if pid in ("15322168", "19739357", "34835292", "37366865", "21618097", "36898756", "36975211", "37858047"):       # lab, order 파일 다시 수집 필요
    #     continue

    # if pid in ("34835292","37366865"):       # lab, order 파일 다시 수집 필요 (EMR에 infliximab 혹은 adalimumab 자체가 안 보임)
    #     continue

    # if pid in ("34835292", "37366865", "26309684","26675590", "26875965","27141223","30013487", "35322270", "37590846", "37858047"):       # lab, order 파일 다시 수집 필요
    #     if pid=="37366865":
    #         raise ValueError

    if pid in ('37366865',):       # lab, order 파일 다시 수집 필요
        continue

    # if pid not in ('17638960',):
    #     continue
    #     if pid in ('17638960',):
    #         raise ValueError

    # ['37366865']  # infliximab quantification은 있는데, dose는 없음. 확인 요망
    # ['12541876']  # 2020년도 쪽 누락 (EMR 기록과 다름)  - 재수집 필요
    # ['17457541']  # 2019년도 쪽 누락 (EMR 기록과 다름) - 재수집 필요
    # ['21640049']  # 2020년도 쪽 누락 (EMR 기록과 다름)  - 재수집 필요
    # ['23807949']  # 2016~2023년도 누락 (EMR 기록과 다름)  - 재수집 필요
    # ['30665275']  # 2020년도 쪽 누락 (EMR 기록과 다름)  - 재수집 필요
    # ['34734236']  # 2022.03.29~2022.08.01 누락 (EMR 기록과 다름)  - 재수집 필요

    # ['18037407']  # 5월 가장 최신 dose 누락
    # ['18839115']  # 5월 가장 최신 dose 누락
    # ['22967071']  # 5월 가장 최신 dose 누락
    # ['27141223']  # 5월 가장 최신 dose 누락
    # ['32482522']  # 5월 가장 최신 dose 누락
    # ['33352150']  # 5월 가장 최신 dose 누락
    # ['36898756']  # 5월 가장 최신 dose 누락
    # ['37002659']  # 5월 가장 최신 dose 누락

    # ['14642249']  # 2019.11.26 으로 기록되었는데, EMR 기록에는 2019.11.27 오더로 나와있음 (참고)
    # ['17619635']  # [self] remsima 로 되어있고, Acting 시간은 존재. 투약날짜는 없는데, 2025.10.10로 되어 있음. 따로 추가코드 삽입.
    # ['17638960']  # [self] humira 로 되어있고, Acting 시간은 존재. 투약날짜는 없는데, 2019.08.26 및 2024.01.31 비어있음. 따로 추가코드 삽입.
    # ['19599395']  # [self] humira 는 Acting 시간은 존재 / 거기 비고에, '가지고 오신걸로 맞으시도록' 이라고 써있고, 그 밑에 Acting이 따로 없는 1 pen 오더도 하나 있음.
    # ['17677819']  # Acting은 비어있어서, Filtering 되었는데, 처방지시에 [SC] x1 3/6 X1 Day 이렇게 날짜 써있는 경우 있음(2014.02.22 adalimumab 투약).... 이건 어떻게 처리되는건지 확인
    # ['19605074']  # 2016.03.17~2016.03.21 에 H나 N(NPO) 로 되어 있고, 비고에는 '맞던 날짜대로'라고 적혀있음
    # ['21267158']  # 2028.07.18 에 self도 아니고 Acting에도 Y로는 안 되어 있다. 이건 투약이 된 것인지?
    # ['23298660']  # 2019.04.27 오더에 Humira 40mg - 수납상태: 미수납, 약국_검사: 비어있음(날짜 안찍혀 있음), 주사시행처: 자가
    # ['25464470']  # 수납상태: 미수납, 주사시행처: 주사처치실로 되어 있고, Acting은 비어있음 그리고 처방지시 비고에 '4/4에 자가로 2개 맞으세요.'
    # ['25498341']  # 2021.08.24T15:37 으로 기록되었는데, EMR 기록에는 2021.08.25 로 나와있음 (참고)  // 2021.08.19, 2021.08.18 에는 [Self]가 있긴하나, ER에서 투여되었다고 써있고 Acting에 H 및 Z로 되어있음
    # ['29277034']  # 임상 SEAVUE_Ustekinumab/Adlimumab/위약(SC) - 2019.09.20 ~ 2020.10.21 동안 임상시험 참여한듯. 이 임상약은 현재 dose로 기록 안 되어 있음. 반영해야할지 결정해야.
    #

    # (1) [self] / 자가 인것 남겨두기? (어떻게 처리된 오더인지?)
    # (2) Acting에 H 나 C, Z 인것은 지우기. N 인것은 물어보기
    # (3) 당일 SC의 경우 pen이 4pen 처방하면, 당일 맞는것도 포함되어 있는것?
    # (4) 목록에 있는 환자 모두가 induction을 이 병원에서 한 것? (37554763 이런분은 SC부터 기록이 있는데...)
    # (5) ID (NAME 제외), 나머지 값들로 drop_duplicates (중간에 개명한 사람 중복됨)

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

dose_result_df = dose_result_df[dose_result_df['ACTING'].map(lambda x: False if (('H' in x) or ('Z' in x) or ('C' in x)) else True)].reset_index(drop=True)
dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)

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
    dose_result_df.at[inx,'ETC_INFO'] = '처방비고 반영완료'
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
change_df['ETC_INFO'] = '처방비고 반영완료'
for inx, row in change_df.iterrows():
    # dose_result_df.iloc[inx]
    dose_result_df.at[inx, 'DOSE'] = int(dose_result_df.at[inx, 'DOSE']/2)
    dose_result_df.at[inx, 'ETC_INFO'] = '처방비고 반영완료'
dose_result_df = pd.concat([dose_result_df,change_df]).reset_index(drop=True)

# Period에 들어있는 추가 글자 제거

for c in PERIOD_added_list1 + PERIOD_added_list2 + [' 컨펌 후 투여',' 맞던 날짜대로',' [D]','금일 ', '매 주','3주 간격 ']:
    dose_result_df['PERIOD'] = dose_result_df['PERIOD'].map(lambda x:x.replace(c,''))
for c in ['1주 1회','ut dict ','ut dict','prn ']:
    dose_result_df['PERIOD'] = dose_result_df['PERIOD'].replace(c,'x1')

# ADDL, II 추가 및 다음 DOSE 보다 넘어가는 ADDL 확인

dose_result_df['PERIOD'] = dose_result_df['PERIOD'].map(lambda x:x.replace('  ',' ').replace('X',' X').replace('  ',' '))
dose_result_df['ADDL'] = dose_result_df['PERIOD'].map(lambda x: int(x.split('/')[0].strip()) * int(x.split('X')[-1].split('Week')[0].strip())/ int(x.split('wks')[0].split('/')[-1].strip())-1 if len(re.findall(r'\d+wks',x)) > 0 else 0).map(int)
dose_result_df['II'] = dose_result_df['PERIOD'].map(lambda x: int(x.split('wks')[0].split('/')[-1].strip()) / int(x.split('/')[0].strip())*7*24 if len(re.findall(r'\d+wks',x)) > 0 else 2*7*24).map(int)
dose_result_df['NXTLST_DOSE_DT'] = dose_result_df.apply(lambda x: (datetime.strptime(x['DATETIME'],'%Y-%m-%dT%H:%M')+timedelta(x['ADDL']*x['II']/24)).strftime('%Y-%m-%dT%H:%M'), axis=1)

# ROUTE가 IV인데 ACTING에 Y가 안 들어 있는 것 제거

dose_result_df = dose_result_df[~((dose_result_df['ROUTE']=='IV')&(dose_result_df['ACTING'].map(lambda x:'Y' not in x)))].copy()

dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True)
dose_result_df.to_csv(f"{output_dir}/dose_df.csv", encoding='utf-8-sig', index=False)

# ADDL시 다음 투약 날짜보다 큰 사람 확인

addl_overflow_df = list()
for inx, frag_df in dose_result_df[['ID','NAME','DATETIME', 'NXTLST_DOSE_DT']].groupby('ID'): #break
    frag_df['SHIFT_NXTLST_DT'] = frag_df['NXTLST_DOSE_DT'].shift(1).fillna('0001-01-01T00:00')
    # frag_df['SHIFT_NXTLST_DT']
    frag_df['DELT_DT'] = frag_df.apply(lambda x: (datetime.strptime(x['DATETIME'],'%Y-%m-%dT%H:%M')-datetime.strptime(x['SHIFT_NXTLST_DT'],'%Y-%m-%dT%H:%M')).days, axis=1)
    addl_overflow_frag = frag_df[frag_df['DELT_DT'] <= 0]
    if len(addl_overflow_frag)>0:
        addl_overflow_df.append(addl_overflow_frag)
addl_overflow_df = pd.concat(addl_overflow_df, ignore_index=True)
addl_overflow_df.to_csv(f"{output_dir}/addl_overflow_dose_df.csv", encoding='utf-8-sig', index=False)


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