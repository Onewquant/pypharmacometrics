from tools import *
from pynca.tools import *

# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Orders
drug_order_set = set()

order_files = glob.glob(f'{resource_dir}/order/{prj_name}_order(*).xlsx')
dose_result_df = list()
no_data_pid_list = list()
result_cols = ['ID','NAME','DT1','DT2','DATETIME','DRUG','DOSE','ACTING','PERIOD','PLACE']
no_dup_cols = [c for c in result_cols if c!='NAME']
for finx, fpath in enumerate(order_files): #break



    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    # if finx <= 3137:
    #     print(f"({finx}) {pname} / {pid} / 확인 완료")
    #     continue



    # if pid in ['10078411', '10125999', '10218916', '10350663', '11216216', '11381622', '11573801', '11719638', '12038983', '12178975', '12273678', '12654044', '13210151', '13599924', '13669366', '13707655', '13999322', '14771293', '14874899', '14935073', '15792433', '19844189', '20141664', '20229672', '20403449', '21010679', '21232138', '21439650', '22110109', '22633695', '23082625', '23844892', '24166533', '24902412', '25122721', '25474011', '25566415', '26129958', '26381198', '26432351', '28026675', '28062981', '28181060', '28258092', '28323927', '29329889', '29428762', '30280977', '30520101', '30825589', '31263386', '31744575', '32364181', '32450020', '32647406', '32668852', '33015987', '33199885', '33991012', '34090905', '34200788', '35207054', '36009912', '36131701', '36623534', '37411570', '37769077', '38780318', '38946884']:       # lab, order 파일 다시 수집 필요
    #     continue
    # if pid in ('10078411','10218916'):       # lab, order 파일 다시 확인 필요 (amikacin 7/10-7/12 이런 식으로 쓰여있음)
    #     continue

    fdf = pd.read_excel(fpath)
    # raise ValueError

    # fdf.columns
    # fdf.to_csv(f"{outcome_dir}/error_dose_df.csv", encoding='utf-8-sig', index=False)

    # fdf[fdf['처방지시'].map(lambda x: ('amikacin' in x.lower()))]
    # fdf.to_csv(f"{output_dir}/error_dose_df.csv", encoding='utf-8-sig', index=False)
    dose_df = fdf[fdf['처방지시'].map(lambda x: ('amikacin' in str(x).lower()) and ('농도' not in str(x).lower()) and ('TDM (자문)' not in str(x).lower()))].copy()
    # dose_df.iloc[0]['처방지시']
    # dose_df.iloc[1]['처방지시']
    # dose_df.iloc[2]

    dose_df = dose_df[(~dose_df['Acting'].isna())].copy()
    if len(dose_df)==0:
        no_data_pid_list.append(pid)
        print(f"({finx}) {pname} / {pid} / No dosing data")
        continue

    print(f"({finx}) {pname} / {pid}")

    # dose_df['처방지시비고'] = dose_df['처방지시'].map(lambda x:x.split(' : ')[-1] if len(x.split(' : '))>1 else '')

    dose_df['ID'] = pid
    dose_df['NAME'] = pname

    # if pname=='김옥순': raise ValueError

    #### 약국_검사가 NA인 것은 제외하고 진행함 (추후 필요시 추가)

    # dose_df[dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test.xlsx",index=False)
    # dose_df[~dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test2.xlsx",index=False)
    dose_df = dose_df[(~dose_df['약국_검사'].isna())].copy()

    #### 성분명이 IBD Biologics 인 경우만으로 필터링 (Infliximab, Adalimumab)

    regex_pattern = r'\(amikacin'
    dose_df = dose_df[dose_df['처방지시'].map(lambda x: bool(re.search(regex_pattern, x, flags=re.IGNORECASE)))].copy()
    # dose_df.to_csv(f"{output_dir}/test_dose_df.csv", encoding='utf-8-sig', index=False)

    dose_df['DT1'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[0].split('[')[-1].replace(' ','T'))
    dose_df['DT2'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[-1].split('[')[-1].replace(' ','T')[:-3])
    dose_df['DATETIME'] = dose_df[['DT1','DT2']].min(axis=1)
    # dose_df['ETC_INFO'] = dose_df['처방지시비고'].copy()
    dose_df['DRUG'] = dose_df['처방지시'].map(lambda x: re.search(regex_pattern, x, flags=re.IGNORECASE).group().lower().replace('(','').replace(')',''))
    # dose_df['ROUTE'] = dose_df['처방지시'].map(lambda x: 'IV' if " [SC] " not in x else 'SC')
    dose_df['PLACE'] = dose_df['주사시행처']
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+',x.split('▣')[-1].split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split(' ')[-1].strip())[0] if " [SC] " not in x else x.split('(Infliximab)')[-1].split('(Adalimumab)')[-1].split('▣')[-1].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+', x.split('mg')[1].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0] if (" [SC] " not in x) else x.split('(Infliximab')[-1].split('(Adalimumab')[-1].split('(Ustekinumab')[-1].split('▣')[-1].split('srg')[0].split('via')[0].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # mg_inx_dict = {'ustekinumab':0,'infliximab':1,'adalimumab':1}
    dose_series = list()
    for dose_inx, dose_row in dose_df.iterrows():
        x = dose_row['처방지시']
        # raise ValueError

        ## '(430-430)mg' 형태로 dose가 쓰여있는 경우
        seq_dose_patterns = re.findall(r'\(\d+[^\)]*\)mg',x)
        if len(seq_dose_patterns)>=1:
            dose_val = '_'.join(re.findall(r'\d+',seq_dose_patterns[0]))
        else:
            try: dose_val = int(re.findall(r'\d+mg',x)[1].replace('mg','').strip())
            except:
                """
                (12) [D/C] Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] x1 
                (6)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (6)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (6)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (6)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (5)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (11)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                이런 모양 처리 고민
                """
                if 'via' in x: via_count = int(re.findall(r'\d+ via', x)[0].replace('via', '').strip())
                if 'amp' in x: via_count = int(re.findall(r'\d+ amp', x)[0].replace('amp', '').strip())

                dose_val = int(re.findall(r'\d+mg', x)[0].replace('mg', '').strip()) * via_count
            # dose_val = re.findall(r'\d+',x.split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0]
            # raise ValueError
        dose_series.append(dose_val)

    # raise ValueError
    dose_df['DOSE'] = dose_series
    # dose_df['DOSE'] = dose_df['DOSE'].map({'1 pen': 40, '2 pen': 80, '2 pen': 160})
    # (1) [원내] Remsima 100mg inj (Infliximab Korea) ...
    dose_df['ACTING'] = dose_df['Acting']

    dose_df['PERIOD'] = dose_df['처방지시'].map(lambda x: re.findall(r'q[\d]+h',x)[0] if len(re.findall(r'q[\d]+h',x))>=1 else x.strip())
    # dose_df.to_csv(f"{output_dir}/dose_df_lhj.csv", encoding='utf-8-sig', index=False)
    # dose_df.loc[3203,'처방지시']

    dose_result_df.append(dose_df[result_cols].drop_duplicates(no_dup_cols))

    # drug_order_set = drug_order_set.union(set(dose_df['처방지시'].map(lambda x:''.join(x.split(':')[0].replace('  ',' ').split(') ')[1:]).replace('[원내]','').replace('[D/C]','').replace('[보류]','').replace('[반납]','').replace('[Em] ','').strip()).drop_duplicates()))


dose_result_df = pd.concat(dose_result_df, ignore_index=True).sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)
# dose_result_df.to_csv(f"{output_dir}/part_dose_df.csv", encoding='utf-8-sig', index=False)

dose_result_df = pd.read_csv(f"{output_dir}/part_dose_df.csv")


## 날짜T시간D도스 형태로 정리

dose_result_df['DATE'] = dose_result_df['DATETIME'].map(lambda x:x.split('T')[0])
dose_result_df['ACTING'] = dose_result_df['ACTING'].map(lambda x:x.replace('MG','mg').replace('G','mg').replace('MD','12A').replace('MN',' 12P')) # # MG 이나 MICU 등 보정 (M 들어가 있어 아래 코드에서 잘못인식)

dt_dose_series = list()
vacant_data = '0000-00-00TNN:NN'
for inx, row in dose_result_df.iterrows():
    # raise ValueError

    new_actval_str=''
    for actval in row['ACTING'].strip().split(','): #break


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

        ## 나머지 처리
        else:
            vacant_y_val = re.findall(r"\d\d:\d\d/", actval)
            if len(vacant_y_val)>0:
                new_actval_str += f'_{vacant_data}'
                continue
            elif ('X' in actval) or ('M' in actval) or ('H' in actval) or ('C' in actval):
                new_actval_str+= f'_{vacant_data}'
                continue
            else:
                rest_y_val = actval.split('/')[-1]
                ## Y가 있긴 한데 다음과 같은 format들
                """
                (344) /Y 2005-05-17(2P ER) / ['2005-05-17'] / ['2P']
                (922) /Y 2005-02-05(ER 4P30) / ['2005-02-05'] / ['4P30']
                (1284)  /Y 2005-04-16(외출9A) / ['2005-04-16'] / ['9A']
                (1932) /Y 2004-01-06(2P30 ER) / ['2004-01-06'] / ['2P30']
                (2214) /Y 청구용 2003-12-11(10A-청구용) / ['2003-12-11'] / ['10A']
                (2506) /Y 2005-06-27(반납) / ['2005-06-27'] / []
                (2671)  /Y 2004-05-06(OA) / ['2004-05-06'] / []
                (3125) /Y 2005-01-28(PRE 적용) / ['2005-01-28'] / []
                (3813) /Y 2004-06-18(10P-청구용) / ['2004-06-18'] / ['10P']
                (4002) /Y 2004-07-21(청구용-6A) / ['2004-07-21'] / ['6A']
                (4012) /Y 2004-06-12(7A-청구용) / ['2004-06-12'] / ['7A']
                (4394) /Y 2004-07-06(비품) / ['2004-07-06'] / []
                (4733) /Y 2004-07-21(청구용) / ['2004-07-21'] / []
                (4850) /Y 2004-08-02(비품) / ['2004-08-02'] / []
                (4944) /Y 2005-05-03(7P-ER) / ['2005-05-03'] / ['7P']
                (5085) /Y 2004-08-05(4P-ER) / ['2004-08-05'] / ['4P']
                (5794) /Y 2005-04-04(청구) / ['2005-04-04'] / []
                (5805) /Y 2004-10-27(청구용) / ['2004-10-27'] / []
                (5887)  /Y 2004-11-09(OA) / ['2004-11-09'] / []
                (6358) /Y 2005-01-02(청구용) / ['2005-01-02'] / []
                (6456)  /Y 2005-01-22(ER5P) / ['2005-01-22'] / ['5P']
                (6664) /Y 2005-02-25(청구용) / ['2005-02-25'] / []
                (6974) /Y 2005-04-22(청) / ['2005-04-22'] / []
                (7213) /Y 2005-09-12(청구용) / ['2005-09-12'] / []
                (7461) /Y 2005-06-01(비품용) / ['2005-06-01'] / []
                (7504) /Y 2005-06-14(청구) / ['2005-06-14'] / []
                (7763)  /Y 2005-08-06(ER 1P) / ['2005-08-06'] / ['1P']
                (7972) /Y 2005-09-04(N) / ['2005-09-04'] / []
                (8074) /Y 2005-09-13(600 9P) / ['2005-09-13'] / ['9P']
                
                # 확인 필요한 부분: D나 N은 투약이 된것인지?
                """
                if ('Y ' in rest_y_val) or ('O ' in rest_y_val) or ('N ' in rest_y_val):
                    date_pattern = re.findall(r"\d\d\d\d-\d\d-\d\d",actval)
                    time_pattern = re.findall(r"\d+[P|A]\d*",actval.upper())
                    if (len(date_pattern) > 0) and (len(time_pattern) > 0):
                        # 시간 처리
                        if 'P' in time_pattern[0]:
                            time_pattern_split = time_pattern[0].split('P')
                            hour_str = str(int(time_pattern_split[0]) + 12)
                        else:
                            time_pattern_split = time_pattern[0].split('A')
                            hour_str = str(int(time_pattern_split[0]))
                        # 분 처리
                        if time_pattern_split[-1]=='':
                            minute_str = '00'
                        else:
                            minute_str = time_pattern_split[-1].zfill(2)
                        num_time_str = f"{hour_str}:{minute_str}"
                        new_actval_str += f'_{date_pattern[0]}T{num_time_str}'
                        continue
                        
                    elif (len(date_pattern) > 0) and (len(time_pattern) == 0):
                        new_actval_str += f'_{date_pattern[0]}TNN:NN'
                        continue
                    elif (len(date_pattern) == 0) and (len(time_pattern) > 0):
                        new_actval_str += f'_{row["DATE"]}T{time_pattern[0]}'
                        continue
                    # print(f"({inx}) {actval} / {date_pattern} / {time_pattern}")
                    else:
                        new_actval_str+=f'_{vacant_data}'
                        continue
                #         
                # elif ('O ' in rest_y_val) or ('N ' in rest_y_val):
                #     date_pattern = re.findall(r"\d\d\d\d-\d\d-\d\d",actval)
                #     # print(f"({inx}) {actval}")
                #     new_actval_str+='_'
                #     continue
                else:
                    new_actval_str+=f'_{vacant_data}'
                    continue
    new_actval_str = new_actval_str[1:]


    # DOSE 붙이기 작업
    dose_split = row['DOSE'].strip().replace(')','').replace('(','').split('_')
    dose_split_set = set(dose_split)
    new_actval_split = new_actval_str.split('_')
    if (len(dose_split_set)>1):
        if len(dose_split)==len(new_actval_split):
            print(f'({inx}) C / {new_actval_str}')
            new_actval_str = '_'.join([f"{new_actval_split[nav_dose_inx]}DOSE{dose_split[nav_dose_inx]}" for nav_dose_inx in range(len(dose_split))])
            # raise ValueError
        else:
            if len(new_actval_split)==1:
                print(f'({inx}) N / {new_actval_str}')
                new_actval_str += f"DOSE{dose_split[-1]}"
            else:
                if len(dose_split) > len(new_actval_split):
                    # 500, 340 이면 뒤의 340만으로 Acting에 적용
                    print(f'({inx}) N / {new_actval_str}')
                    new_actval_str = '_'.join([f"{nav}DOSE{dose_split[-len(new_actval_split):][nav_inx]}" for nav_inx,nav in enumerate(new_actval_split)])
                #elif len(dose_split) < len(new_actval_split):
                else:
                    # 500, 340, 340 이면 뒤의 340을 연장해서 Acting에 적용
                    print(f'({inx}) N / {new_actval_str}')
                    new_actval_str = '_'.join([f"{nav}DOSE{(dose_split+[dose_split[-1]]*(len(new_actval_split)-len(dose_split)))[nav_inx]}" for nav_inx,nav in enumerate(new_actval_split)])
                # else:
                #     pass
                # '500_340_340' 이렇게 되어 있는 분 (11836850) 중 ACTING에 2개만 있는 row는 340_340으로 일단 처리
                # '500_340_340' 인데, ['2007-11-08T02:18', '2007-11-08T11:18', '2007-11-08T16:55', '2007-11-08T23:26']
    else:
        print(f'({inx}) O / {new_actval_str}')
        # print('one dose')
        # Acting에 dose가 들어있는 경우도 생각해서 넣어줘야.
        new_actval_str = '_'.join([f"{nav}DOSE{list(dose_split_set)[0]}" for nav in new_actval_split])

        # raise ValueError
        
        # Dose가 한 자리수인 경우도 존재 / dose_split 3개, acting_split 2개 인 경우 있음

    dt_dose_series.append(new_actval_str)
dose_result_df['DT_DOSE'] = dt_dose_series
dose_result_df.to_csv(f"{output_dir}/dt_dose_df.csv", encoding='utf-8-sig', index=False)




dose_result_df = pd.read_csv(f"{output_dir}/dt_dose_df.csv")
"""
# 각 actval에 Dose 추가하는 작업
# 각 actval을 dt_dose_series 에 append하여 한 컬럼 구성
"""
                # print(f"({inx}) {actval}")
                # continue
                # raise ValueError

    # if 'Y' in row['ACTING']:
    #     ymdt_patterns = re.findall(r"[\d][\d][\d][\d]-[\d][\d]-[\d][\d] [\d][\d]:[\d][\d]",row['ACTING'])
    #     if len(ymdt_patterns)>0:
    #         new_dt = ymdt_patterns[0].replace(' ','T')
    #     else:
    #         time_patterns = re.findall(r"[\d][\d]:[\d][\d]",row['ACTING'])
    #         new_dt = row['DATETIME'].split('T')[0] + 'T' + time_patterns[0]
    #     dt_series.append(new_dt)
    # else:
    #     dt_series.append(row['DATETIME'])

dose_result_df['DATETIME'] = dt_series

## PERIOD 추가 조정

dose_result_df['PERIOD'] = dose_result_df['PERIOD'].map(lambda x:x.replace('2주 간격 유지용량 ','').replace('#1 ','').replace('3/20 ','').replace('매주 ','').replace(' X1 Day','').replace('1회 ','').replace('/wk ','/1wks').replace('x1 ','1'))


## H, C, Z, M는 제거

dose_result_df = dose_result_df[dose_result_df['ACTING'].map(lambda x: False if (('H' in x) or ('Z' in x) or ('C' in x)) else True)].reset_index(drop=True)
dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)
dose_result_df['ETC_INFO_TREATED'] = ''

## 처방비고 반영 및 기타 정리