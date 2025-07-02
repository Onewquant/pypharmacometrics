from tools import *
from pynca.tools import *

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
# prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Orders
lab_set = set()

conc_files = glob.glob(f'{resource_dir}/lab/{prj_name}_lab(*).xlsx')
conc_result_df = list()
for finx, fpath in enumerate(conc_files): #break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    print(f"({finx}) {pname} / {pid}")

    # if pid in ("15322168", "19739357", "34835292", "37366865", "21618097", "36898756", "36975211", "37858047"):       # lab, order 파일 다시 수집 필요
    #     continue

    fdf = pd.read_excel(fpath)
    # lab_set = lab_set.union(set(fdf[fdf['검사명'].map(lambda x:'amikacin' in x.lower())]['검사명'].unique()))

    # fdf['검사명'].unique()
    # break

    # fdf.columns

    # cumlab_set = cumlab_set.union(set(fdf['Lab'].unique()))
    # [c for c in cumlab_set if ('infliximab' in c.lower()) or ('adalimumab' in c.lower())]
    # lab_set = lab_set.union(set(fdf['검사명'].unique()))
    # [c for c in lab_set if ('infliximab' in c.lower()) or ('adalimumab' in c.lower())]

    # 'Infliximab 정량: 재검한 결과입니다.'
    # 'Infliximab Quantification'
    # 'Adalimumab Quantification'

    # conc_df.columns
    conc_df = fdf[fdf['검사명'].isin(['Amikacin'])].copy()     # Amikacin 농도 측정 항목만 추출
    """
    try:
        conc_df['검사결과'].astype(float)
    except:
        print(f"({finx}) {pname} / {pid} / float 아님 추가")
        conc_result_df.append(conc_df[['검사결과']])
        continue        
    print(f"({finx}) {pname} / {pid}")
    continue
    conc_result_df = pd.concat(conc_result_df, ignore_index=True)
    
    strconc_inx_list = list()
    for inx, res in enumerate(conc_result_df['검사결과']): #break
        try: float(res.replace('<','').replace('>',''))
        except:strconc_inx_list.append(inx)

    strconc_df = conc_result_df.iloc[strconc_inx_list].copy()
    for inx, row in strconc_df.reset_index(drop=True).iterrows():
        print(f"({inx}) {row['검사결과']}")

    strconc_df.to_csv(f"{output_dir}/str_conc_df(lab).csv", encoding='utf-8-sig', index=False)
    """


    # fdf = fdf[~fdf['Acting'].isna()].copy()
    # dose_df = fdf[fdf['처방지시'].map(lambda x: (('adalimumab' in x.lower()) or ('infliximab' in x.lower())) and ('quantification' not in x.lower()))].copy()
    # dose_df['처방지시비고'] = dose_df['처방지시'].map(lambda x:x.split(' : ')[-1] if len(x.split(' : '))>1 else '')
    #
    #
    conc_df['ID'] = pid
    conc_df['NAME'] = pname
    conc_df['DRUG'] = conc_df['검사명'].map(lambda x:x.split(' ')[0].lower())
    conc_df['검사결과'] = conc_df['검사결과'].replace('중복 오더임',np.nan).astype(str)
    conc_df['CONC'] = conc_df['검사결과'].map(lambda x:float(x.split('(')[0].replace('<','').replace('>','').strip()) if x.lower()!='nan' else np.nan)
    conc_df['POTENTIAL_SAMPLING_INFO'] = conc_df['검사결과'].map(lambda x:x.split('(')[-1].strip() if x.lower() != 'nan' else np.nan)


    conc_df = conc_df[~conc_df['CONC'].isna()].copy()
    # conc_df[] = conc_df[['보고일','오더일']].max(axis=1)

    conc_result_df.append(conc_df[['ID','NAME','보고일','오더일','DRUG','CONC', 'POTENTIAL_SAMPLING_INFO']])
    # conc_result_df['ID'].drop_duplicates()
    # # if pname=='김옥순': raise ValueError
    #
    # #### 약국_검사가 NA인 것은 제외하고 진행함 (추후 필요시 추가)
    #
    # # dose_df[dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test.xlsx",index=False)
    # # dose_df[~dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test2.xlsx",index=False)
    # dose_df = dose_df[~dose_df['약국_검사'].isna()].copy()
    #
    # #### 성분명이 IBD Biologics 인 경우만으로 필터링 (Infliximab, Adalimumab)
    #
    # regex_pattern = r'\(infliximab|adalimumab\)'
    # dose_df = dose_df[dose_df['처방지시'].map(lambda x: bool(re.search(regex_pattern, x, flags=re.IGNORECASE)))].copy()
    #
    # dose_df['DT1'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[0].split('[')[-1].replace(' ','T'))
    # dose_df['DT2'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[-1].split('[')[-1].replace(' ','T')[:-3])
    # dose_df['DATETIME'] = dose_df[['DT1','DT2']].min(axis=1)
    # dose_df['ETC_INFO'] = dose_df['처방지시비고'].copy()
    # dose_df['DRUG'] = dose_df['처방지시'].map(lambda x: re.search(regex_pattern, x, flags=re.IGNORECASE).group().lower().replace('(','').replace(')',''))
    # dose_df['ROUTE'] = dose_df['처방지시'].map(lambda x: 'IV' if " [SC] " not in x else 'SC')
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+',x.split('▣')[-1].split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split(' ')[-1].strip())[0] if " [SC] " not in x else x.split('(Infliximab)')[-1].split('(Adalimumab)')[-1].split('▣')[-1].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # dose_df['PERIOD'] = dose_df['처방지시'].map(lambda x: 'x1' if " [SC] " not in x else x.split(' [SC] ')[-1].split(':')[0].strip())
    #
    # # dose_df.loc[3203,'처방지시']
    # dose_result_df.append(dose_df[['ID','NAME','DATETIME','DRUG','DOSE','ROUTE','PERIOD','ETC_INFO']].copy())

    # drug_order_set = drug_order_set.union(set(dose_df['처방지시'].map(lambda x:''.join(x.split(':')[0].replace('  ',' ').split(') ')[1:]).replace('[원내]','').replace('[D/C]','').replace('[보류]','').replace('[반납]','').replace('[Em] ','').strip()).drop_duplicates()))

conc_result_df = pd.concat(conc_result_df, ignore_index=True).drop_duplicates(['ID','보고일','오더일','DRUG','CONC'])

etc_info_list = list()
for x in conc_result_df['POTENTIAL_SAMPLING_INFO']:
    try:
        float(x)
        etc_info_list.append(np.nan)
    except:
        etc_info_list.append(x)
conc_result_df['POTENTIAL_SAMPLING_INFO'] = etc_info_list

conc_result_df.to_csv(f"{output_dir}/conc_df(lab).csv", encoding='utf-8-sig', index=False)




conc_result_df = pd.read_csv(f'{output_dir}/conc_df(lab).csv')

## 랩 자체에 sampling 시간 같이 쓰여 있는 경우 parsing

pot_time_list = list()
pot_date_list = list()
for inx, row in conc_result_df.iterrows():

    pot_date_str = ''
    year = max(row['보고일'],row['오더일'])[:4]

    ## 날짜 Parsing
    samp_md_patterns = re.findall(r'\d{1,2}\/\d+', str(row['POTENTIAL_SAMPLING_INFO']))
    if len(samp_md_patterns) > 0:
        sp_md = samp_md_patterns[0]
        sp_md_split = sp_md.split('/')
        month = sp_md_split[0]
        day = sp_md_split[-1]
        if (int(month) > 12) or (int(day) > 31):
            print('(month > 12) or (day > 31)')

            str_day = day
            if len(str_day) == 4:
                day = str_day[:2]
            elif len(str_day) == 3:
                sp_first = sp_md_split[-1]
                if ((int((sp_first[:2])) <= 31) and (int((sp_first[2:])) <= 12)) and ((int((sp_first[:1])) <= 31) and (int((sp_first[-2:])) <= 12)):
                    raise ValueError
                else:
                    if (int(sp_first[:2]) <= 31) and (int(sp_first[2:]) <= 12):
                        day = sp_first[:2]
                    elif (int(sp_first[:1]) <= 31) and (int(sp_first[-2:]) <= 12):
                        day = int(sp_first[:1])
            # raise ValueError

        sp_md = f"{year}-{month.zfill(2)}-{day.zfill(2)}"
        pot_date_str = sp_md
    pot_date_list.append(pot_date_str)


    ## 시간 Parsing

    search_samp_str = str(row['POTENTIAL_SAMPLING_INFO']).upper().replace('MD', '12:').replace('::', ':')

    # NN:NN 타입
    pot_samp_str = ''

    samp_patterns = re.findall(r'\d+:\d+',search_samp_str)
    if len(samp_patterns)>0:
        sp = samp_patterns[0]
        hour = int(sp.split(':')[0])
        minute = int(sp.split(':')[-1])
        if (hour >= 24) or (minute >= 60):
            print('(hour >= 24) or (minute >= 60)')
            raise ValueError
        else:
            pot_samp_str = f"T{str(hour).zfill(2)}:{str(minute).zfill(2)}"
            pot_time_list.append(pot_samp_str)
            continue



    # NN [P|A]: NN 타입

    samp_patterns = re.findall(r'\d+[P|A]:*\d+', search_samp_str)
    samp_short_patterns = re.findall(r'\d+[P|A]', search_samp_str)
    if (len(samp_patterns)>0) or (len(samp_short_patterns)>0):

        # NN [P|A] 만 있는 경우 -> 뒤에 minute을 00으로 붙여줌
        if ((len(samp_patterns)>0)==False) and (len(samp_short_patterns)>0):
            sp = samp_short_patterns[0] + '00'
        # raise ValueError
        else:
            sp = samp_patterns[0].replace(':','')
        if 'A' in sp:
            sp_split = sp.split('A')
            hour = int(sp_split[0])
            minute = int(sp_split[-1])

            if (hour > 12) or (minute >= 60):
                print('AM - (hour > 12) or (minute >= 60)')
                if hour <= 24:
                    pot_samp_str = f"T{str(hour).zfill(2)}:{str(minute).zfill(2)}"
                    pot_time_list.append(pot_samp_str)
                    continue

                str_hour = str(hour)
                if len(str_hour)==4:
                    hour = int(str(hour)[-2:])
                elif len(str_hour)==3:
                    sp_first = sp_split[0]
                    if ((int((sp_first[:2])) <= 31) and (int((sp_first[2:])) <= 12)) and ((int((sp_first[:1])) <= 31) and (int((sp_first[-2:])) <= 12)):
                        raise ValueError
                    else:
                        if (int((sp_first[:2])) <= 31) and (int((sp_first[2:])) <= 12):
                            hour = int((sp_first[2:]))
                        elif (int((sp_first[:1])) <= 31) and (int((sp_first[-2:])) <= 12):
                            hour = int((sp_first[-2:]))

        elif 'P' in sp:
            sp_split = sp.split('P')
            hour = int(sp_split[0]) + 12
            minute = int(sp_split[-1])

            if (hour >= 24) or (minute >= 60):
                print('PM - (hour >= 24) or (minute >= 60)')
                if hour <= 36:
                    pot_samp_str = 'T'+sp.replace('P',':')
                    pot_time_list.append(pot_samp_str)
                    continue

                str_hour = str(hour)
                if len(str_hour)==4:
                    hour = int(str(hour)[-2:])
                elif len(str_hour)==3:
                    sp_first = sp_split[0]
                    if ((int((sp_first[:2])) <= 31) and (int((sp_first[2:])) <= 12)) and ((int((sp_first[:1])) <= 31) and (int((sp_first[-2:])) <= 12)):
                        raise ValueError
                    else:
                        if (int((sp_first[:2])) <= 31) and (int((sp_first[2:])) <= 12):
                            hour = int((sp_first[2:])) + 12
                        elif (int((sp_first[:1])) <= 31) and (int((sp_first[-2:])) <= 12):
                            hour = int((sp_first[-2:])) + 12

                    # raise ValueError

        pot_samp_str = f"T{str(hour).zfill(2)}:{str(minute).zfill(2)}"
        pot_time_list.append(pot_samp_str)
        continue


    pot_time_list.append(pot_samp_str)
conc_result_df['POT_SAMP_TIME'] = pot_time_list
conc_result_df['POT_SAMP_MONTHDAY'] = pot_date_list


conc_result_df = conc_result_df[['ID', 'NAME', '보고일', '오더일', 'DRUG', 'CONC', 'POT_SAMP_MONTHDAY', 'POT_SAMP_TIME','POTENTIAL_SAMPLING_INFO']].copy()
conc_result_df.to_csv(f"{output_dir}/final_conc_df.csv", encoding='utf-8-sig', index=False)

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
pt_info = pd.read_csv(f"{resource_dir}/AMK_tdm.csv", encoding='euc-kr')
pt_info = pt_info.sort_values(['등록번호','검사일']).drop_duplicates(['등록번호'], ignore_index=True)
rename_dict = {'등록번호':'ID','환자명':'NAME','검사일':'TDM_REQ_DATE','작성일':'TDM_RES_DATE','성별':'SEX','몸무게':'WT','키':'HT','나이':'AGE','병동':'WARD','진료과':'DEP'}
pt_info = pt_info.rename(columns=rename_dict)
pt_info['ID'] = pt_info['ID'].astype(str)
# pt_info['TDM_REQ_DATE'] = pt_info['TDM_REQ_DATE'].map(lambda x:x.replace(' ','T'))
pt_info['TDM_REQ_DATE'] = pt_info['TDM_REQ_DATE'].map(lambda x:x.split(' ')[0])
pt_info = pt_info[list(rename_dict.values())]
pt_info.to_csv(f"{output_dir}/patient_info.csv", encoding='utf-8-sig', index=False)
# pt_info.columns

pt_files = glob.glob(f'{resource_dir}/lab/{prj_name}_lab(*).xlsx')

dose_result_df = pd.read_csv(f"{output_dir}/final_dose_df.csv")
dose_result_df['ID'] = dose_result_df['ID'].astype(str)
dose_result_df['DOSE_DT'] = dose_result_df['DATE']+'T'+dose_result_df['TIME']
dose_result_df['TIME'] = "T"+dose_result_df['TIME']
dose_result_df = dose_result_df.rename(columns={'DATE':'DOSE_DATE','TIME':'DOSE_TIME'})
dose_result_df = dose_result_df[['ID', 'NAME', 'DOSE_DATE','DOSE_TIME','DOSE_DT', 'DOSE']].copy()
uniq_dose_pids = list(dose_result_df.drop_duplicates(['ID'])['ID'].astype(str))


conc_result_df = pd.read_csv(f"{output_dir}/final_conc_df.csv")
conc_result_df = conc_result_df[['ID', 'NAME', '보고일', '오더일', 'DRUG', 'CONC', 'POT_SAMP_MONTHDAY', 'POT_SAMP_TIME']].copy()
conc_result_df['ID'] = conc_result_df['ID'].astype(str)
conc_result_df['POT_SAMP_MONTHDAY'] = conc_result_df['POT_SAMP_MONTHDAY'].replace(np.nan,'')
conc_result_df['POT_SAMP_TIME'] = conc_result_df['POT_SAMP_TIME'].replace(np.nan,'')
conc_result_df['POT채혈DT'] = conc_result_df['POT_SAMP_MONTHDAY'] + conc_result_df['POT_SAMP_TIME']
uniq_conc_pids = list(conc_result_df.drop_duplicates(['ID'])['ID'].astype(str))

sampling_result_df = pd.read_csv(f"{output_dir}/final_sampling_df.csv")
# sampling_result_df = sampling_result_df[['ID', 'NAME', '오더일','보고일','채혈DT']].copy()
sampling_result_df['ID'] = sampling_result_df['ID'].astype(str)
uniq_sampling_pids = list(sampling_result_df.drop_duplicates(['ID'])['ID'].astype(str))


final_df = list()
conc_samp_mismatch_pids = list()
no_dose_pids = list()
error_passing_pids = set()
for finx, fpath in enumerate(pt_files): #break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]
    # pt_df.append({'ID':pid,'NAME':pname})
    id_info_df = pt_info[pt_info['ID'] == pid].copy()
    min_date = (datetime.strptime(id_info_df['TDM_REQ_DATE'].iloc[0],'%Y-%m-%d')-timedelta(days=62)).strftime('%Y-%m-%d')
    max_date = (datetime.strptime(id_info_df['TDM_RES_DATE'].iloc[0],'%Y-%m-%d')+timedelta(days=92)).strftime('%Y-%m-%d')

    id_dose_df = dose_result_df[dose_result_df['ID'] == pid].copy()
    id_conc_df = conc_result_df[conc_result_df['ID'] == pid].copy()
    id_conc_df = id_conc_df[(id_conc_df['오더일']>=min_date)&(id_conc_df['오더일']<=max_date)&(id_conc_df['보고일']>=min_date)&(id_conc_df['보고일']<=max_date)]
    id_samp_df = sampling_result_df[sampling_result_df['ID'] == pid].copy()

    res_frag_df = id_conc_df.copy()
    for c in ['SAMP_DT','채혈DT', '라벨DT', '접수DT', '시행DT', '보고DT']:
        res_frag_df[c] = ''

    # if finx
    # raise ValueError
    # res_frag_df.columns

    # POT채혈DT 있는 경우
    pot_dt_rows = res_frag_df[(res_frag_df['POT채혈DT'] != '') & (~res_frag_df['POT채혈DT'].isna())]
    if len(pot_dt_rows)!=0:
        print(f"({finx}) {pname} / {pid} / POT채혈DT 결과 존재")
        for potdt_inx, potdt_row in pot_dt_rows.iterrows(): #break
            # pot_dt_rows['POT채혈DT']
            # POT채혈DT에 날짜도 써 있는 경우 -> 해당 날짜로 SAMP_DT 입력
            if bool(re.search(r'\d\d\d\d-\d\d-\d\dT\d\d:\d\d',potdt_row['POT채혈DT'])):
                if res_frag_df.at[potdt_inx, 'SAMP_DT'] == '':
                    res_frag_df.at[potdt_inx, 'SAMP_DT'] = potdt_row['POT채혈DT']
                else:
                    continue
                error_passing_pids.add(pid)
            # res_frag_df['POT채혈DT']
            # res_frag_df['SAMP_DT']
            # POT채혈DT에 시간만 써 있는 경우 -> SAMPLING DATA 있는 경우는 아래서 처리. 없는 경우는 여기서 처리
            elif bool(re.search(r'T\d\d:\d\d',potdt_row['POT채혈DT'])):
                if res_frag_df.at[potdt_inx, 'SAMP_DT'] == '':
                    if pid in uniq_sampling_pids:
                        print("POT채혈DT 결과 존재 / SAMPLING DATA 있어 아래서 처리")
                        error_passing_pids.add(pid)
                        pass
                    else:
                        print("POT채혈DT 결과 존재 하는데 / SAMPLING DATA 없는 경우 입니다")
                        # POT채혈DT, Dosing 기록으로 유추해야.

                        # 농도측정의 오더일과 보고일이 같은 날짜인 경우는 그 날짜로 사용
                        if potdt_row['오더일']==potdt_row['보고일']:
                            if res_frag_df.at[potdt_inx, 'SAMP_DT'] == '':
                                res_frag_df.at[potdt_inx, 'SAMP_DT'] = potdt_row['오더일'] + potdt_row['POT채혈DT']

                        # 농도측정의 오더일과 보고일이 다른 날짜들의 경우
                        ord_noteq_rep_rows = res_frag_df[(res_frag_df['오더일'] != res_frag_df['보고일'])&(res_frag_df['SAMP_DT']=='')].copy()
                        if len(ord_noteq_rep_rows)>0:
                            # 농도채혈 오더 난 날짜가 1개만 있을때
                            if len(ord_noteq_rep_rows['오더일'].drop_duplicates()) == 1:
                                for oneqr_inx, oneqr_row in ord_noteq_rep_rows.iterrows():
                                    res_frag_df.at[oneqr_inx, 'SAMP_DT'] = oneqr_row['오더일'] + oneqr_row['POT채혈DT']
                            # 농도채혈 오더 난 날짜가 1개 이상 존재할 때 있을때
                            else:
                                print("POT채혈DT 결과 존재 / SAMPLING DATA 없음 / 농도채혈 오더 난 날짜가 1개 이상 존재")
                                raise ValueError
                        error_passing_pids.add(pid)


                    # res_frag_df[['보고일', '오더일', 'CONC', ]]
                    # id_dose_df
                    # id_samp_df[['보고일', '오더일', '채혈DT', ]]
                    # res_frag_df.at[potdt_inx, 'SAMP_DT'] = potdt_row['POT채혈DT']
                else:
                    continue
            # POT채혈DT에 날짜만 써 있는 경우 -> SAMPLING DATA 있는 경우는 아래서 처리. 없는 경우는 여기서 처리
            elif bool(re.search(r'\d\d\d\d-\d\d-\d\d',potdt_row['POT채혈DT'])):
                if res_frag_df.at[potdt_inx, 'SAMP_DT']=='':
                    # id_info_df
                    mean_conc = res_frag_df['CONC'].mean()
                    date_dose_rows = id_dose_df[id_dose_df['DOSE_DATE']==potdt_row['POT채혈DT']].copy()
                    if len(date_dose_rows)<=2:
                        est_conc_dt_tups = ((datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')-timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                                            (datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')+timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
                                            )
                        # est_conc_dt_tups = (datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')-timedelta(minutes=30),date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')+timedelta(minutes=30))

                        if res_frag_df.at[potdt_inx, 'CONC'] < mean_conc:
                            res_frag_df.at[potdt_inx, 'SAMP_DT'] = potdt_row['POT채혈DT'] + est_conc_dt_tups[0]
                        else:
                            res_frag_df.at[potdt_inx, 'SAMP_DT'] = potdt_row['POT채혈DT'] + est_conc_dt_tups[1]
                    else:
                        print(f"({finx}) {pname} / {pid} / POT채혈DT에 날짜만 존재 / 해당 날짜 dose기록 3개 이상")
                        raise ValueError

            else:
                # 이 조건은 만족하는 게 없을 것이라 예상
                print(f"({finx}) {pname} / {pid} / 이 조건은 만족하는 게 없을 것이라 예상")
                raise ValueError
        # len(error_passing_pids)
        # if pid not in error_passing_pids:
        #     raise ValueError

#     final_df.append(res_frag_df)
# final_df = pd.concat(final_df, ignore_index=True)
# final_df.to_csv(f"{output_dir}/final_conc_df(with sampling).csv", encoding='utf-8-sig', index=False)

    # CONC 데이터만 존재 (SAMPLING TIME은 X, POT채혈DT는 X 예상)
    # 위에서 내려온 케이스 :
    #   POT채혈DT에 시간만 써 있는 경우 / SAMPLING DATA 존재시
    #   POT채혈DT에 날짜만 써 있는 경우 / SAMPLING DATA 존재시
    if (pid in uniq_conc_pids) and (pid not in uniq_sampling_pids):
        # POT채혈DT 존재하는 잘못된 케이스 있나 확인
        if len(res_frag_df[(res_frag_df['SAMP_DT']=='')&(res_frag_df['POT채혈DT']!='')]):
            print('POT채혈DT 존재하는 잘못된 케이스 있나 확인')
            raise ValueError
        # DOSING 데이터 없는 경우
        if len(id_dose_df) == 0:
            print(f"({finx}) {pname} / {pid} / Only Conc data / No Dose Data (Check!)")
            no_dose_pids.append(pid)
            continue
        else:
            print(f"({finx}) {pname} / {pid} / Only Conc data (Dosing: O, Sampling X, POT채혈DT: X")

        # res_frag_df.columns
        # id_info_df
        # res_frag_df[(res_frag_df['SAMP_DT']=='')]

        ord_eq_rep_rows = res_frag_df[(res_frag_df['SAMP_DT'] == '') & (res_frag_df['오더일'] == res_frag_df['보고일'])].copy()
        # 농도측정의 오더일과 보고일이 같은 날짜인 경우는 그 날짜로 사용

        for oer_inx, oer_row in ord_eq_rep_rows.iterrows():

            mean_conc = res_frag_df['CONC'].mean()
            date_dose_rows = id_dose_df[id_dose_df['DOSE_DATE'] == oer_row['POT채혈DT']].copy()
            if len(date_dose_rows) <= 2:
                est_conc_dt_tups = ((datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'], '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                                    (datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'], '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
                                    )
                # est_conc_dt_tups = (datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')-timedelta(minutes=30),date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')+timedelta(minutes=30))

                if res_frag_df.at[oer_inx, 'CONC'] < mean_conc:
                    res_frag_df.at[oer_inx, 'SAMP_DT'] = oer_row['오더일'] + est_conc_dt_tups[0]
                else:
                    res_frag_df.at[oer_inx, 'SAMP_DT'] = oer_row['오더일'] + est_conc_dt_tups[1]
            else:
                print(f"({finx}) {pname} / {pid} / Only Conc data (Dosing: O, Sampling X, POT채혈DT: X / 해당 날짜 dose기록 3개 이상")
                raise ValueError

        # 농도측정의 오더일과 보고일이 다른 날짜들의 경우
        ord_noteq_rep_rows = res_frag_df[(res_frag_df['SAMP_DT'] == '') & (res_frag_df['오더일'] != res_frag_df['보고일'])].copy()

        conc_order_date_list = list(res_frag_df['오더일'].unique())
        # 농도채혈 order가 난 날짜가 1개만 있을때 (dosing 타임과 농도데이터 개수 및 농도 값 고려해서 배열)
        if len(conc_order_date_list) == 1:
            id_dose_df = id_dose_df[id_dose_df['DOSE_DATE'] == order_date_list[0]].copy()

            # 측정된 농도 값이 2개 이하이면 DOSING 타임 및 CONC 값 고려하여 배열
            if len(res_frag_df) <= 2:
                est_dose_dt = datetime.strptime(id_dose_df.iloc[0]['DOSE_DT'], '%Y-%m-%dT%H:%M')
                est_conc_dt_tups = ((est_dose_dt - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                                    (est_dose_dt + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))
                res_frag_df = res_frag_df.sort_values(['CONC'], ignore_index=True)
                for resf_inx, resf_row in res_frag_df.iterrows():  # break
                    if resf_row['SAMP_DT'] == '':
                        res_frag_df.at[resf_inx, 'SAMP_DT'] = est_conc_dt_tups[resf_inx]
                    # res_frag_df['SAMP_DT']
            else:
                print(f"({finx}) {pname} / {pid} / 농도 측정 오더 날짜 1개, 측정된 농도 값이 2개 이상")
                raise ValueError

        # 농도 측정 오더 날짜 2개 이상
        else:
            print(f"({finx}) {pname} / {pid} / 농도 측정 오더 날짜 2개 이상")
            raise ValueError

#     final_df.append(res_frag_df)
# final_df = pd.concat(final_df, ignore_index=True)
# final_df.to_csv(f"{output_dir}/final_conc_df(with sampling).csv", encoding='utf-8-sig', index=False)




    # CONC, SAMP 데이터 둘 다 있는 경우
    elif (pid in uniq_conc_pids) and (pid in uniq_sampling_pids):
        continue
        # print(f"({finx}) {pname} / {pid} / CONC, SAMP 데이터 둘 다 존재")
        raise ValueError

        # id_info_df
        # res_frag_df[['보고일','오더일', 'CONC']]

        cdf = id_conc_df.sort_values(['보고일', 'CONC'])
        sdf = id_samp_df[['보고일', '채혈DT', '라벨DT', '접수DT','시행DT','보고DT']].copy()

        mdf = cdf.merge(sdf, on=['보고일'], how='outer')
        mdf = mdf.sort_values(['보고일', 'CONC', '채혈DT'])

        trough_mdf = mdf.drop_duplicates(['보고일'], keep='first')
        peak_mdf = mdf.drop_duplicates(['보고일'], keep='last')
        total_mdf = pd.concat([trough_mdf, peak_mdf]).drop_duplicates(['보고일', 'CONC', '채혈DT'], keep='last').sort_values(['보고일', 'CONC', '채혈DT'])
        # total_mdf.columns
        # cdf = id_conc_df[['보고일', 'POT채혈DT', 'CONC']].sort_values(['보고일', 'CONC'])
        # sdf = id_samp_df[['보고일', '채혈DT', '라벨DT', '접수DT']].copy()
        #
        # mdf = cdf.merge(sdf, on=['보고일'], how='outer')[['보고일', 'CONC', '채혈DT', 'POT채혈DT']]
        # mdf = mdf.sort_values(['보고일', 'CONC', '채혈DT'])
        #
        # trough_mdf = mdf.drop_duplicates(['보고일'], keep='first')
        # peak_mdf = mdf.drop_duplicates(['보고일'], keep='last')
        # total_mdf = pd.concat([trough_mdf, peak_mdf]).drop_duplicates(['보고일','CONC', '채혈DT'], keep='last').sort_values(['보고일', 'CONC', '채혈DT'])
        #


        # CONC 데이터는 오더일, 보고일 날짜 두개만 존재
        # SAMPLING 데이터는 여러개의 날짜 존재 두개만 존재

        if len(total_mdf)!=len(cdf):
            print(f"({finx}) {pname} / {pid} / No matched")
            conc_samp_mismatch_pids.append(pid)
            # len(conc_samp_mismatch_pids)
            # continue
        else:
            print(f"({finx}) {pname} / {pid} / matched")
            # if pid not in ['11116501', '10112328', '10143478', '10228470', '10533576', '10885385', '10914959', ]:
            #     raise ValueError
            """
            '10112328' : 04.17 의 샘플링 타임데이터는 있는데 농도데이터는 부재함.
            '10143478' : 2018-06-04 의 농도데이터가 3개 있는데, 샘플링 데이터는 2개라 max, min만 남기면 df 길이 달라짐
            '10228470' : 2005-04-22 샘플링 데이터는 존재, 농도 데이터는 그날 것 없음. (total_mdf 에 CONC가 NAN인 값 생김)
            '10533576' : 2004-09-14에 CONC가 0.5로 똑같은 데이터 2개 존재. 아마도 9-13일 채혈일듯. (보고일 기준으로만 하고 있는데, 중복된 데이터의 오더일은 다른 것으로 보아 하나는 2004-09-13 데이터인듯
            '10885385' : 샘플링 데이터 보고일기준 2003-12-30 는 1개, 2004-01-02 는 2개 인데, 농도데이터에서는 2, 1개로 되어 있음
            '10914959'
            '11116501'
            """
        final_df.append(total_mdf)

        # if finx==3:
        #     raise ValueError

        # cdf = id_conc_df[['오더일', '보고일', 'POT채혈DT', 'CONC']].copy()
        # sdf = id_samp_df[['오더일', '보고일', '채혈DT', '라벨DT', '접수DT']].copy()

        # mdf = cdf.merge(sdf, on=['오더일', '보고일'], how='outer')[['오더일', '보고일','CONC','채혈DT','POT채혈DT']]
        # mdf = mdf.sort_values(['오더일', '보고일', 'CONC', '채혈DT'])

        # trough_mdf = mdf.drop_duplicates(['오더일', '보고일'], keep='first')
        # peak_mdf = mdf.drop_duplicates(['오더일', '보고일'], keep='last')
        # total_mdf = pd.concat([trough_mdf, peak_mdf]).drop_duplicates(['오더일', '보고일', 'CONC', '채혈DT'], keep='last')


    else:
        print(f"({finx}) {pname} / {pid} / 그 어디도 X")
        raise ValueError


    # elif (pid not in uniq_conc_pids) or (pid not in uniq_sampling_pids):
    #     print(f"({finx}) {pname} / {pid} / No eighther data")
    #     continue
    # else:
    #     raise ValueError

final_df = pd.concat(final_df, ignore_index=True)
final_df.to_csv(f"{output_dir}/final_conc_df(with sampling).csv", encoding='utf-8-sig', index=False)


"""
## 새로운 로직 필요할듯

(1) 'POT채혈DT' 있는 경우 -> POT채혈DT로 시간 설정
# 날짜,시간 모두 존재: POT채혈DT로 채혈시간 설정 [CONFIRMED]
# 채혈시간만 존재: 
    - (보고일==오더일)인 경우: 오더일로 날짜 설정 후 채혈시간 적용 [CONFIRMED]
    - (보고일!=오더일)인 경우: [PASS]
        * '채혈DT' 있는 경우: [PASS]
        * '채혈DT' 없는 경우: [PASS]

(2) '채혈DT' 있는 경우


(보고일==오더일)인 경우는 채혈시간은 그날 NN:NN으로 일단 설정 
      # '채혈DT' 있는 경우 -> 그날의 Dosing 시점 및 '채혈DT' 
      # 채혈기록 없는 경우 -> 그날의 Dosing 시점 고려하여 Dosing 전후로 설정
# (2) (보고일=!오더일)인 경우
      # '채혈DT' 있는 경우 -> 그날의 
"""
