from tools import *
from pynca.tools import *

# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


## 첫 TDM 요청 이전에 HD 하고 있는 환자 목록 수집
pt_info = pd.read_csv(f"{output_dir}/patient_info.csv", encoding='utf-8-sig')

## Orders
drug_order_set = set()

order_files = glob.glob(f'{resource_dir}/order/{prj_name}_order(*).xlsx')
hd_result_df = list()
no_hd_data_pid_list = list()
no_hd_dt_data_list = list()
no_drug_keyword_pid_list = list()
result_cols = ['ID','NAME','DATETIME','HD_TYPE','ACTING','PLACE']
no_dup_cols = [c for c in result_cols if c not in ['NAME']]
# len(order_files)
# len(no_data_pid_list)
for finx, fpath in enumerate(order_files): #break


    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    # if finx <= 3137:
    #     print(f"({finx}) {pname} / {pid} / 확인 완료")
    #     continue

    # if pid=='10908086':
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
    hd_df = fdf[fdf['처방지시'].map(lambda x: ('renal replacement' in str(x).lower()) or ('신대체요법' in str(x).lower()) or ('crrt' in str(x).lower()) or ('투석' in str(x).lower()) or ('hemodialysis' in str(x).lower()))].copy()
    # dose_df.iloc[0]['처방지시']
    # dose_df.iloc[1]['처방지시']
    # dose_df.iloc[2]


    # if len(hd_df)>0:
    #     raise ValueError
    # else:
    #     continue

    hd_df = hd_df[(~hd_df['Acting'].isna())].copy()


    if len(hd_df)==0:
        no_hd_data_pid_list.append(pid)
        print(f"({finx}) {pname} / {pid} / No hemodialysis data")
        continue
    else:
        print(f"({finx}) {pname} / {pid} / Yes hemodialysis data")

    raise ValueError

    hd_df['ID'] = pid
    hd_df['NAME'] = pname

    # dose_df['HD_TYPE'] = dose_df['처방지시'].map(lambda x: re.search(regex_pattern, x, flags=re.IGNORECASE).group().lower().replace('(','').replace(')',''))
    hd_df['PLACE'] = hd_df['주사시행처']
    hd_df['ACTING'] = hd_df['Acting']
    # hd_df.columns

    hd_df['DATETIME'] = '0000-00-00TNN:NN'
    hd_df['HD_TYPE'] = hd_df['처방지시']
    # hd_df['변경의']

    # dose_df[['DATETIME','DOSE','ACTING']]

    # dose_df.to_csv(f"{output_dir}/dose_df_lhj.csv", encoding='utf-8-sig', index=False)
    # dose_df.loc[3203,'처방지시']

    hd_result_df.append(hd_df[result_cols])

    # drug_order_set = drug_order_set.union(set(dose_df['처방지시'].map(lambda x:''.join(x.split(':')[0].replace('  ',' ').split(') ')[1:]).replace('[원내]','').replace('[D/C]','').replace('[보류]','').replace('[반납]','').replace('[Em] ','').strip()).drop_duplicates()))

# dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True)
# dose_result_df = pd.concat(dose_result_df, ignore_index=True).sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)
hd_result_df = pd.concat(hd_result_df, ignore_index=True).sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True)


hd_result_df.to_csv(f"{output_dir}/hd_df_datacheck.csv", encoding='utf-8-sig', index=False)
# len(no_hd_data_pid_list)
# hd_result_df['']

# dose_result_df.columns

# print(f"TOTAL 오더 파일 수 : {len(order_files)} 명")
# print(f"ACTING 정보 존재 : {len(order_files)-len(no_data_pid_list)} 명 (-{len(no_data_pid_list)} 명)")
# print(f"약국_검사에 Datetime 정보 존재 : {len(order_files)-len(no_data_pid_list)-len(no_dosing_dt_data_list)} 명 (-{len(no_dosing_dt_data_list)} 명)")
# print(f"처방지시에 Amikacin 키워드 존재 : {len(order_files)-len(no_data_pid_list)-len(no_dosing_dt_data_list)-len(no_drug_keyword_pid_list)} 명 (-{len(no_drug_keyword_pid_list)} 명)")
# print(f"모든 DOSING 정보 존재 : {len(dose_result_df.drop_duplicates(['ID']))} 명")
# dose_result_df.drop_duplicates(['ID'])

# dose_result_df = pd.read_csv(f"{output_dir}/part_dose_df.csv")

# dose_result_df[dose_result_df['NAME']=='박원필'][['DATETIME','ACTING']]

## 날짜T시간DOSE도스 형태로 정리

hd_result_df['ID'] = hd_result_df['ID'].astype(str)
hd_result_df['DATE'] = hd_result_df['DATETIME'].map(lambda x:x.split('T')[0])
# hd_result_df['ACTING'] = hd_result_df['ACTING'].map(lambda x:x.replace('MG','mg').replace('G','mg').replace('MD','12A').replace('MN',' 12P')) # # MG 이나 MICU 등 보정 (M 들어가 있어 아래 코드에서 잘못인식)


dt_dose_series = list()
vacant_data = '0000-00-00TNN:NN'
for inx, row in hd_result_df.iterrows():
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
                        same_date_other_rows = hd_result_df[(hd_result_df['ID']==row['ID'])&(hd_result_df.index!=inx)&(hd_result_df['DATE']==date_pattern[0])].copy()
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


    dt_dose_series.append(new_actval_str)
hd_result_df['DT'] = dt_dose_series
# hd_result_df.to_csv(f"{output_dir}/final_hd_df.csv", encoding='utf-8-sig', index=False)

# dose_result_df[['DATE','DOSE','ACTING','DT_DOSE']]

#
# dose_result_df = pd.read_csv(f"{output_dir}/dt_dose_df.csv")
# vacant_data = '0000-00-00TNN:NN'
## ACTING 기록 개별 분리작업

final_hd_df = list()
cur_id = ''
for inx, row in hd_result_df.iterrows(): #break
    if cur_id!=row['ID']:
        print(f"({inx} / {len(hd_result_df)}) {row['ID']} / ACTING 기록 개별 분리작업")
        cur_id=row['ID']
    row_df = pd.DataFrame(columns=['ID','NAME','DRUG','PERIOD','DT_DOSE','ETC_INFO'])
    row_df['DT_DOSE'] = row['DT_DOSE'].split('_')
    for c in ['ID','NAME','DRUG','PERIOD','ETC_INFO']:
        row_df[c] = row[c]
    final_hd_df.append(row_df)
final_hd_df = pd.concat(final_hd_df, ignore_index=True)
final_hd_df['DATE'] = final_hd_df['DT'].map(lambda x:x.split('T')[0])
final_hd_df['TIME'] = final_hd_df['DT'].map(lambda x:x.split('T')[-1])
final_hd_df = final_hd_df[(final_hd_df['DATE']!=vacant_data.split('T')[0])]
final_hd_df['DATETIME'] = final_hd_df['DATE'] + 'T' + final_hd_df['TIME']
final_hd_df = final_hd_df.sort_values(['ID','DATETIME'], ignore_index=True)
final_hd_df.to_csv(f"{output_dir}/final_dose_df.csv", encoding='utf-8-sig', index=False)


