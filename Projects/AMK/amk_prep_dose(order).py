from tools import *
from pynca.tools import *

# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Orders
drug_order_set = set()

order_files = glob.glob(f'{resource_dir}/order/{prj_name}_order(*).xlsx')
dose_result_df = list()
no_data_pid_list = list()
no_dosing_dt_data_list = list()
no_drug_keyword_pid_list = list()
result_cols = ['ID','NAME','DT1','DT2','DATETIME','DRUG','DOSE','ACTING','PERIOD','PLACE','ETC_INFO']
no_dup_cols = [c for c in result_cols if c not in ['NAME','ETC_INFO']]
# len(order_files)
# len(no_data_pid_list)
for finx, fpath in enumerate(order_files): #break


    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    # if finx <= 3137:
    #     print(f"({finx}) {pname} / {pid} / 확인 완료")
    #     continue

    # if pid=='11114219':
    #     raise ValueError
    # else:
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

    # dose_df['처방지시비고'] = dose_df['처방지시'].map(lambda x:x.split(' : ')[-1] if len(x.split(' : '))>1 else '')

    dose_df['ID'] = pid
    dose_df['NAME'] = pname

    # if pname=='김옥순': raise ValueError

    #### 약국_검사가 NA인 것은 제외하고 진행함 (추후 필요시 추가)

    # dose_df[dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test.xlsx",index=False)
    # dose_df[~dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test2.xlsx",index=False)
    dose_df = dose_df[(~dose_df['약국_검사'].isna())].copy()
    if len(dose_df)==0:
        no_dosing_dt_data_list.append(pid)
        print(f"({finx}) {pname} / {pid} / No dosing datetime data")
        continue

    #### 성분명이 IBD Biologics 인 경우만으로 필터링 (Infliximab, Adalimumab)

    regex_pattern = r'amikacin'
    # regex_pattern = r'\(amikacin'
    dose_df = dose_df[dose_df['처방지시'].map(lambda x: bool(re.search(regex_pattern, x, flags=re.IGNORECASE)) and not re.search(r'normal\s+saline', x, re.IGNORECASE) and not re.search(r'dextrose', x, re.IGNORECASE))].copy()

    if len(dose_df)==0:
        no_drug_keyword_pid_list.append(pid)
        print(f"({finx}) {pname} / {pid} / Amikacin Dosing")
        continue

    print(f"({finx}) {pname} / {pid}")

    if pid == '28340528':  # 이 사람 날짜가 제대로 입력되어 있지 않은 row 존재
        dose_df.at[dose_df.index[-1], '약국_검사'] = dose_df.at[dose_df.index[-2], '약국_검사']

    dose_df['DT1'] = dose_df['약국_검사'].map(lambda x:re.findall(r'\d\d\d\d-\d\d-\d\d \d\d:\d\d',x)[0].replace(' ','T'))
    dose_df['DT2'] = dose_df['약국_검사'].map(lambda x:re.findall(r'\d\d\d\d-\d\d-\d\d \d\d:\d\d',x)[-1].replace(' ','T'))
    dose_df['DATETIME'] = dose_df[['DT1','DT2']].min(axis=1)
    # dose_df[['DT1', 'DT2','Acting']]
    # dose_df[['약국_검사', 'Acting']]
    # dose_df[['DATETIME', 'Acting']]
    # dose_df[['DATETIME','처방지시']]
    # dose_df[['처방지시']]
    # dose_df.to_csv(f"{output_dir}/test_dose_df.csv", encoding='utf-8-sig', index=False)

    if pid == '11114219':
        dose_df.at[855,'DATETIME'] = dose_df.at[855,'DATETIME'].replace('2013-05-19','2013-05-20')
        dose_df.at[856,'DATETIME'] = dose_df.at[856,'DATETIME'].replace('2013-05-19','2013-05-20')
        dose_df.at[965,'DATETIME'] = dose_df.at[965,'DATETIME'].replace('2013-05-17','2013-05-18')
        dose_df.at[1279,'처방지시'] = dose_df.at[1279,'처방지시'].replace('500mg [IVS]','550mg [IVS]')
        dose_df.at[1453,'처방지시'] = dose_df.at[1453,'처방지시'].replace('500mg [IVS]','550mg [IVS]')
    elif pid == '17963774':  # 2010-01-10것이 2010-01-09로 잘못 들어감
        # raise ValueError
        change_inx = 939
        dose_df.at[change_inx ,'DATETIME'] = dose_df.at[change_inx ,'DATETIME'].replace('2010-01-09','2010-01-10')
    elif pid=='24184973':
        change_inx = 284
        dose_df.at[change_inx,'DATETIME'] = dose_df.at[change_inx,'DATETIME'].replace('2023-01-29','2023-01-30')
    elif pid=='17169200':
        dose_df.at[2372,'DATETIME'] = dose_df.at[2372,'DATETIME'].replace('2013-12-13','2013-12-14')
        dose_df.at[2732,'DATETIME'] = dose_df.at[2732,'DATETIME'].replace('2013-12-10','2013-12-11')
        dose_df.at[2828,'DATETIME'] = dose_df.at[2828,'DATETIME'].replace('2013-12-09','2013-12-10')
        dose_df.at[2934,'DATETIME'] = dose_df.at[2934,'DATETIME'].replace('2013-12-08','2013-12-09')
        dose_df.at[3043,'DATETIME'] = dose_df.at[3043,'DATETIME'].replace('2013-12-07','2013-12-08')
        dose_df.at[3151,'DATETIME'] = dose_df.at[3151,'DATETIME'].replace('2013-12-06','2013-12-07')
        dose_df.at[3252,'DATETIME'] = dose_df.at[3252,'DATETIME'].replace('2013-12-05','2013-12-06')
        dose_df.at[3362,'DATETIME'] = dose_df.at[3362,'DATETIME'].replace('2013-12-04','2013-12-05')
        dose_df.at[3476,'DATETIME'] = dose_df.at[3476,'DATETIME'].replace('2013-12-03','2013-12-04')
    elif pid=='14188310':
        dose_df.at[533,'DATETIME'] = dose_df.at[533,'DATETIME'].replace('2015-04-25','2015-04-26')
        dose_df.at[534,'DATETIME'] = dose_df.at[534,'DATETIME'].replace('2015-04-25','2015-04-26')

        # raise ValueError
    # if pid == '28340528':  # 2010-01-10것이 2010-01-09로 잘못 들어감
    #     dose_df.at[792,'DATETIME'] = dose_df.at[792,'DATETIME'].replace('2010-01-09','2010-01-10')
    #     # raise ValueError

    # for
    # dose_df[dose_df['DATETIME']=='']
    # dose_df['ETC_INFO'] = dose_df['처방지시비고'].copy()
    dose_df['DRUG'] = dose_df['처방지시'].map(lambda x: re.search(regex_pattern, x, flags=re.IGNORECASE).group().lower().replace('(','').replace(')',''))
    # dose_df['ROUTE'] = dose_df['처방지시'].map(lambda x: 'IV' if " [SC] " not in x else 'SC')
    dose_df['PLACE'] = dose_df['주사시행처']
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+',x.split('▣')[-1].split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split(' ')[-1].strip())[0] if " [SC] " not in x else x.split('(Infliximab)')[-1].split('(Adalimumab)')[-1].split('▣')[-1].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+', x.split('mg')[1].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0] if (" [SC] " not in x) else x.split('(Infliximab')[-1].split('(Adalimumab')[-1].split('(Ustekinumab')[-1].split('▣')[-1].split('srg')[0].split('via')[0].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # mg_inx_dict = {'ustekinumab':0,'infliximab':1,'adalimumab':1}
    dose_series = list()
    for dose_inx, dose_row in dose_df.iterrows(): #break
        # if dose_inx==1413:
        #     raise ValueError
        x = dose_row['처방지시'].replace('inj신풍 2ml','inj신풍 500mg').split('h ')[0].split('ut dict')[0]  # 이걸해야 dose 2인 것 안 나오긴 함. 근데 그냥 아래서 수정
        # x = dose_row['처방지시']
        # raise ValueError

        ## '(430-430)mg' 형태로 dose가 쓰여있는 경우
        seq_dose_patterns = re.findall(r'\(\d+[^\)]*\)mg',x)
        if len(seq_dose_patterns)>=1:
            dose_val = '_'.join(re.findall(r'\d+',seq_dose_patterns[0]))
        else:
            try: dose_val = float(re.findall(r'\d+\.*\d*m',x)[1].replace('m','').strip())
            except:
                """
                (12) [D/C] Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] x1 
                (6)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (6)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (6)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (6)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (5)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (11)  Amikacin soln 500mg 신풍(Amikacin) (430-430)mg [IVS] q12h 
                (26)  Amikacin 500mg inj신풍 2ml [MIV] q12h 

                이런 모양 처리 고민
                """
                if 'via' in x:
                    via_count = float(re.findall(r'\d*\.?\d+ via', x)[0].replace('via', '').strip())
                elif 'amp' in x:
                    via_count = float(re.findall(r'\d*\.?\d+ amp', x)[0].replace('amp', '').strip())
                else:
                    via_count = 1
                    print(f"({finx}) {pname} / {pid} / (No vial count) / {x}")
                    # continue
                # re.findall(r'\d+\.?\d*m', '500mg')
                mg_num_list = re.findall(r'\d+\.?\d*m', x)
                if len(mg_num_list)==0:
                    mg_num = np.nan
                    print(f"({finx}) {pname} / {pid} / (No mg number) / {x}")
                    # continue
                else:
                    mg_num = int(mg_num_list[0].replace('m', '').strip())
                dose_val = mg_num * via_count
            # dose_val = re.findall(r'\d+',x.split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0]
            # raise ValueError

        # 잘못 기록된 것 있는 듯 (250 -> 25로 기록 해둠 / 10 배 해줘야)
        # if type(dose_val)!=str:
        #     if dose_val <= 140:
        #         raise ValueError
        #         dose_val = dose_val*10
        dose_series.append(dose_val)

    # raise ValueError
    dose_df['DOSE'] = dose_series
    # dose_df[['DATETIME','DOSE','Acting']]
    # dose_df['DOSE'] = dose_df['DOSE'].map({'1 pen': 40, '2 pen': 80, '2 pen': 160})
    # (1) [원내] Remsima 100mg inj (Infliximab Korea) ...
    dose_df['ACTING'] = dose_df['Acting']

    dose_df['PERIOD'] = dose_df['처방지시'].map(lambda x: re.findall(r'q[\d]+h',x)[0] if len(re.findall(r'q[\d]+h',x))>=1 else x.split(':')[0].split(']')[-1].strip())
    dose_df['ETC_INFO'] = dose_df['처방지시'].map(lambda x:x.split(':')[-1].strip() if ':' in x else '')
    dose_df = dose_df[~dose_df['DOSE'].isna()].copy()

    # dose_df[['DATETIME','DOSE','ACTING']]

    # dose_df.to_csv(f"{output_dir}/dose_df_lhj.csv", encoding='utf-8-sig', index=False)
    # dose_df.loc[3203,'처방지시']

    dose_result_df.append(dose_df[result_cols])

    # drug_order_set = drug_order_set.union(set(dose_df['처방지시'].map(lambda x:''.join(x.split(':')[0].replace('  ',' ').split(') ')[1:]).replace('[원내]','').replace('[D/C]','').replace('[보류]','').replace('[반납]','').replace('[Em] ','').strip()).drop_duplicates()))

# dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True)
# dose_result_df = pd.concat(dose_result_df, ignore_index=True).sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)
dose_result_df = pd.concat(dose_result_df, ignore_index=True).sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True)
dose_result_df.to_csv(f"{output_dir}/part_dose_df.csv", encoding='utf-8-sig', index=False)
# dose_result_df.columns

print(f"TOTAL 오더 파일 수 : {len(order_files)} 명")
print(f"ACTING 정보 존재 : {len(order_files)-len(no_data_pid_list)} 명 (-{len(no_data_pid_list)} 명)")
print(f"약국_검사에 Datetime 정보 존재 : {len(order_files)-len(no_data_pid_list)-len(no_dosing_dt_data_list)} 명 (-{len(no_dosing_dt_data_list)} 명)")
print(f"처방지시에 Amikacin 키워드 존재 : {len(order_files)-len(no_data_pid_list)-len(no_dosing_dt_data_list)-len(no_drug_keyword_pid_list)} 명 (-{len(no_drug_keyword_pid_list)} 명)")
print(f"모든 DOSING 정보 존재 : {len(dose_result_df.drop_duplicates(['ID']))} 명")
# dose_result_df.drop_duplicates(['ID'])

# dose_result_df = pd.read_csv(f"{output_dir}/part_dose_df.csv")

# dose_result_df[dose_result_df['NAME']=='박원필'][['DATETIME','ACTING']]

## 날짜T시간DOSE도스 형태로 정리

dose_result_df['ID'] = dose_result_df['ID'].astype(str)
dose_result_df['DATE'] = dose_result_df['DATETIME'].map(lambda x:x.split('T')[0])
dose_result_df['ACTING'] = dose_result_df['ACTING'].map(lambda x:x.replace('MG','mg').replace('G','mg').replace('MD','12A').replace('MN',' 12P')) # # MG 이나 MICU 등 보정 (M 들어가 있어 아래 코드에서 잘못인식)
dose_result_df['DOSE'] = dose_result_df['DOSE'].astype(str)

dt_dose_series = list()
vacant_data = '0000-00-00TNN:NN'
for inx, row in dose_result_df.iterrows():
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

        ## 날짜 시간/Y일때


        ## 나머지 처리
        else:
            # actval = actval
            vacant_y_val = re.findall(r"\d\d:\d\d/", actval)
            if len(vacant_y_val)>0:
                new_actval_str += f'_{vacant_data.replace("NN:NN",vacant_y_val[0][:-1])}'
                continue
            elif ('X' in actval) or ('M' in actval) or ('H' in actval) or ('C' in actval) or ('N' in actval):
                new_actval_str+= f'_{vacant_data}'
                continue
            else:
                rest_y_val = actval.split('/')[-1]

                """
                # 확인 필요한 부분: D나 N은 투약이 된것인지?
                """
                if (('Y ' in rest_y_val) or ('O ' in rest_y_val)):
                # if ('Y ' in rest_y_val) or ('O ' in rest_y_val) or ('N ' in rest_y_val):
                    date_pattern = re.findall(r"\d\d\d\d-\d\d-\d\d",actval)
                    time_pattern = re.findall(r"\d+[P|A]\d*",actval.upper())
                    if (len(date_pattern) > 0) and (len(time_pattern) > 0):
                        # 시간 처리
                        if 'P' in time_pattern[0]:
                            time_pattern_split = time_pattern[0].split('P')
                            hour_str = str(int(time_pattern_split[0]) + 12).zfill(2)
                        else:
                            time_pattern_split = time_pattern[0].split('A')
                            hour_str = str(int(time_pattern_split[0])).zfill(2)
                        # 분 처리
                        if time_pattern_split[-1]=='':
                            minute_str = '00'
                        else:
                            minute_str = time_pattern_split[-1].zfill(2)
                        num_time_str = f"{hour_str}:{minute_str}"
                        new_actval_str += f'_{date_pattern[0]}T{num_time_str}'
                        continue

                    elif (len(date_pattern) > 0) and (len(time_pattern) == 0):
                        # ACTING에 날짜만 기록되어 있고 시간 정보는 없을떄 -> 같은 날짜에 시간 기록 있는지 확인해서 그 기록으로 중복시킨 후 아래서 중복처리에 포함시킴
                        # rep_acting_str.split(',')
                        # dose_result_df
                        same_date_other_rows = dose_result_df[(dose_result_df['ID']==row['ID'])&(dose_result_df.index!=inx)&(dose_result_df['DATE']==date_pattern[0])].copy()
                        if len(same_date_other_rows)==1:

                            same_date_rec_times = [x.replace('/Y','') for x in re.findall(r"\d\d:\d\d/Y", same_date_other_rows['ACTING'].iloc[0])]
                            if same_date_other_rows.index[0] < inx:
                                reversed_sorting = True
                            else:
                                reversed_sorting = False
                            same_date_rec_times.sort(reverse=reversed_sorting)
                            # new_actval_str = ''
                            add_actval_completion = False
                            for sdr_time in same_date_rec_times:
                                add_actval_str = f'{date_pattern[0]}T{sdr_time}'
                                if add_actval_str not in new_actval_str:
                                    new_actval_str += f'_{add_actval_str}'
                                    add_actval_completion = True
                                    break
                                else:
                                    continue

                            if add_actval_completion==False:
                                print('날짜만 존재하는 Acting 정보에 시간 정보 추가 어려움')
                                raise ValueError

                            """
                            ID         /  약국_검사                                           Acting
                            10908086   / 접수 [2005-06-27 20:00]   2005-06-27 20:00:53		/Y 2005-06-27(반납), 
                            10948899  / 접수 [2004-05-06 06:00]   2004-05-06 06:00:34		12:00/ 2004-05-06(), /Y 2004-05-06(OA), 
                            """
                        elif len(same_date_other_rows)==0:
                            print('같은 날짜 기록 0 개')

                        else:
                            print('같은 날짜 기록 여러개')
                            if ('(비품용)' in rest_y_val):
                                new_actval_str += f'_{vacant_data}'
                            else:
                                raise ValueError
                            # dose_result_df[(dose_result_df['ID'] == row['ID'])][['DATETIME','ACTING','PERIOD','DOSE']]
                            # raise ValueError
                        # pid


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

# dose_result_df[['DATE','DOSE','ACTING','DT_DOSE']]

#
# dose_result_df = pd.read_csv(f"{output_dir}/dt_dose_df.csv")
# vacant_data = '0000-00-00TNN:NN'
## ACTING 기록 개별 분리작업

final_dose_df = list()
cur_id = ''
for inx, row in dose_result_df.iterrows(): #break
    if cur_id!=row['ID']:
        print(f"({inx} / {len(dose_result_df)}) {row['ID']} / ACTING 기록 개별 분리작업")
        cur_id=row['ID']
    row_df = pd.DataFrame(columns=['ID','NAME','DRUG','PERIOD','DT_DOSE','ETC_INFO'])
    row_df['DT_DOSE'] = row['DT_DOSE'].split('_')
    for c in ['ID','NAME','DRUG','PERIOD','ETC_INFO']:
        row_df[c] = row[c]
    final_dose_df.append(row_df)
final_dose_df = pd.concat(final_dose_df, ignore_index=True)
final_dose_df['DATE'] = final_dose_df['DT_DOSE'].map(lambda x:x.split('T')[0])
final_dose_df['TIME'] = final_dose_df['DT_DOSE'].map(lambda x:x.split('T')[-1].split('DOSE')[0])
final_dose_df['DOSE'] = final_dose_df['DT_DOSE'].map(lambda x:float(x.split('DOSE')[-1]))
final_dose_df = final_dose_df[(final_dose_df['DATE']!=vacant_data.split('T')[0])]

# final_dose_df['DOSE'] = final_dose_df['DOSE'].replace(2.0,500.0)
# sns.displot(final_dose_df['DOSE'].min())
dose_unique_list = list(final_dose_df['DOSE'].unique())
dose_unique_list.sort()
# final_dose_df[final_dose_df['DOSE']==6]

## 같은 DT에 다른 DOSE 값 가지고 있으면, SUM 하여 사용
# final_dose_df = pd.read_csv(f"{output_dir}/final_dose_df.csv")
final_dose_df['ETC_INFO'] = (final_dose_df['ETC_INFO'].replace(np.nan, '')+'||').replace('||', '')
final_dose_df['DATETIME'] = final_dose_df['DATE'] + 'T' + final_dose_df['TIME']
final_dose_df = final_dose_df.sort_values(['ID','DATETIME'], ignore_index=True)
# final_dose_df[final_dose_df['ID']=='10023985']

## 중복 오더 처리
print('중복 오더 처리 작업 시작')
final_result_df = list()
cur_id==''
cur_count = 1
for pid, frag_df in final_dose_df.groupby('ID'):
    # if pid=='13415653':
    if cur_id!=pid:
        print(f"({cur_count}) {pid} / {len(frag_df)} rows / 중복 오더 처리 작업 시작")
        cur_count+=1
        cur_id=pid


    ## 2시간 이내로 중복된 같은 값의 중복된 오더 있고, 그 앞뒤로 비어있는 날 있으면, 비어있는 날의 order로 설정
    dose_pct_diff = np.abs(frag_df['DOSE'].diff() / frag_df['DOSE'])
    dt = pd.to_datetime(frag_df['DATETIME'])
    dt_diff = dt.diff().dt.total_seconds() / 3600

    short_dosing_interval_rows = frag_df[(dt_diff < 2)&(dose_pct_diff==0)].copy()
    if len(short_dosing_interval_rows) > 0:
        min_date = frag_df['DATE'].min()
        max_date = frag_df['DATE'].max()
        dates_in_range = pd.date_range(start=min_date, end=max_date).astype(str)
        dose_vacant_dates = set(dates_in_range) - set(frag_df['DATE'])


        # if str(pid) not in ['13415653','10015489','10020669','10022083','10020669','10023985','10024128']:
        #     raise ValueError

        if len(dose_vacant_dates)==0:
            pass
        else:
            for sdi_inx, sdi_row in short_dosing_interval_rows.iterrows():
                prev_date = (datetime.strptime(sdi_row['DATE'],'%Y-%m-%d') - timedelta(days=1)).strftime('%Y-%m-%d')
                next_date = (datetime.strptime(sdi_row['DATE'],'%Y-%m-%d') + timedelta(days=1)).strftime('%Y-%m-%d')
                if prev_date in dose_vacant_dates:
                    frag_df.at[sdi_inx, 'DATE'] = prev_date
                    frag_df.at[sdi_inx, 'DATETIME'] = prev_date + 'T' + frag_df.at[sdi_inx, 'DATETIME'].split('T')[-1]
                elif next_date in dose_vacant_dates:
                    frag_df.at[sdi_inx, 'DATE'] = next_date
                    frag_df.at[sdi_inx, 'DATETIME'] = next_date + 'T' + frag_df.at[sdi_inx, 'DATETIME'].split('T')[-1]
                else:
                    pass
                    # raise ValueError
        frag_df = frag_df.sort_values(['ID','DATE','TIME'], ignore_index=True)
    else:
        pass

    ## DOSE 투여 기록들의 DATETIME이 2시간 이내이고 DOSE 크기 차이가 40% 이상 나면 DOSE합치고, 크기 차이가 30% 미만이면 최근 것 하나만 선택
    dose_pct_diff = np.abs(frag_df['DOSE'].diff() / frag_df['DOSE'])
    dt = pd.to_datetime(frag_df['DATETIME'])
    dt_diff = dt.diff().dt.total_seconds() / 3600

    short_dosing_interval_rows = frag_df[(dt_diff < 2)].copy()
    if len(short_dosing_interval_rows) > 0:
        small_dose_diff_rows = frag_df[(dt_diff < 2)&(dose_pct_diff <= 0.3)].copy()
        large_dose_diff_rows = frag_df[(dt_diff < 2)&(dose_pct_diff > 0.3)].copy()
        for sdd_inx,sdd_row in small_dose_diff_rows.iterrows():
            frag_df = frag_df[frag_df.index != (sdd_inx-1)].copy()

        for sdd_inx,sdd_row in large_dose_diff_rows.iterrows():
            modi_frag_df = frag_df[frag_df.index.isin([sdd_inx-1,sdd_inx])].copy()
            modi_frag_df.at[modi_frag_df.index[-1],'DOSE'] = modi_frag_df['DOSE'].sum()
            modi_frag_df = modi_frag_df.iloc[-1:,:]
            keep_frag_df = frag_df[~frag_df.index.isin([sdd_inx-1,sdd_inx])].copy()
            frag_df = pd.concat([modi_frag_df, keep_frag_df]).sort_values(['ID','DATE','TIME'])
        frag_df = frag_df.sort_values(['ID','DATE','TIME']).reset_index(drop=True)

        # raise ValueError

    final_result_df.append(frag_df)
    # break
# print(final_dose_df[final_dose_df['ID']==13415653])
# final_dose_df= final_dose_df.groupby(['ID','NAME','DRUG','PERIOD','DATE','TIME'],as_index=False).agg({'DOSE':'sum','ETC_INFO':'sum'})

# final_dose_df[final_dose_df['ID']=='13415653']
# final_dose_df.to_csv(f"{output_dir}/final_dose_datacheck.csv", encoding='utf-8-sig', index=False)

# final_dose_df.to_csv(f"{output_dir}/final_dose_datacheck.csv", encoding='utf-8-sig', index=False)
# raise ValueError
# final_dose_df = final_dose_df.groupby(['ID','NAME','DRUG','DATE','TIME'],as_index=False).agg({'DOSE':'sum','ETC_INFO':'sum','PERIOD':'min'})
# print(final_dose_df[final_dose_df['ID']==18115888])
# print(len(final_dose_df['ETC_INFO'].unique()))
# final_dose_df = final_dose_df[(final_dose_df['TIME']!=vacant_data.split('T')[-1])]
# raise ValueError

# final_dose_df = final_dose_df[['ID','NAME','DRUG','PERIOD','DATE','TIME','DOSE','ETC_INFO']].sort_values(['ID','DATE','TIME'], ignore_index=True)
# final_dose_df.to_csv(f"{output_dir}/final_dose_df.csv", encoding='utf-8-sig', index=False)

final_dose_df = pd.concat(final_result_df).sort_values(['ID','DATE','TIME'], ignore_index=True)
# final_dose_df.to_csv(f"{output_dir}/final_dose_datacheck.csv", encoding='utf-8-sig', index=False)
# final_dose_df = final_dose_df[(final_dose_df['TIME']!=vacant_data.split('T')[-1])]
final_dose_df.to_csv(f"{output_dir}/final_dose_df.csv", encoding='utf-8-sig', index=False)

# final_dose_df = pd.read_csv(f"{output_dir}/final_dose_df.csv")
# final_dose_df[final_dose_df['DOSE']==8]
# sns.displot(final_dose_df['DOSE'])


no_dosing_data_df = pd.DataFrame(columns=['ID'])
no_dosing_data_df['ID'] = no_data_pid_list
no_dosing_data_df.to_csv(f"{output_dir}/err_no_dose_df.csv", encoding='utf-8-sig', index=False)

no_dosing_dt_df = pd.DataFrame(columns=['ID'])
no_dosing_dt_df['ID'] = no_dosing_dt_data_list
no_dosing_dt_df.to_csv(f"{output_dir}/err_no_dosing_dt_df.csv", encoding='utf-8-sig', index=False)


no_no_drug_keyword_df = pd.DataFrame(columns=['ID'])
no_no_drug_keyword_df['ID'] = no_drug_keyword_pid_list
no_no_drug_keyword_df.to_csv(f"{output_dir}/err_no_drug_keyword_df.csv", encoding='utf-8-sig', index=False)


# no_data_pid_list = list()
# no_dosing_dt_data_list = list()
# no_drug_keyword_pid_list

# final_dose_df[final_dose_df['PERIOD'].map(lambda x:len(x)>20)]

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

# dose_result_df['DATETIME'] = dt_series

## PERIOD 추가 조정

# dose_result_df['PERIOD'] = dose_result_df['PERIOD'].map(lambda x:x.replace('2주 간격 유지용량 ','').replace('#1 ','').replace('3/20 ','').replace('매주 ','').replace(' X1 Day','').replace('1회 ','').replace('/wk ','/1wks').replace('x1 ','1'))


## H, C, Z, M는 제거

# dose_result_df = dose_result_df[dose_result_df['ACTING'].map(lambda x: False if (('H' in x) or ('Z' in x) or ('C' in x)) else True)].reset_index(drop=True)
# dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)
# dose_result_df['ETC_INFO_TREATED'] = ''

## 처방비고 반영 및 기타 정리