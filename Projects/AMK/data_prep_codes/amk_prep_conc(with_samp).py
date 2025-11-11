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


pt_info = pd.read_csv(f"{resource_dir}/AMK_tdm.csv", encoding='euc-kr')
pt_info = pt_info.sort_values(['등록번호','검사일'])
rename_dict = {'등록번호':'ID','환자명':'NAME','검사일':'TDM_REQ_DATE','작성일':'TDM_RES_DATE','성별':'SEX','몸무게':'WT','키':'HT','나이':'AGE','병동':'WARD','진료과':'DEP'}
pt_info = pt_info.rename(columns=rename_dict)
pt_info['ID'] = pt_info['ID'].astype(str)
# pt_info['TDM_REQ_DATE'] = pt_info['TDM_REQ_DATE'].map(lambda x:x.replace(' ','T'))
pt_info['TDM_REQ_DATE'] = pt_info['TDM_REQ_DATE'].map(lambda x:x.split(' ')[0])
pt_info = pt_info[list(rename_dict.values())]
pt_info_dup = pt_info.copy()
pt_info = pt_info.drop_duplicates(['ID'], ignore_index=True)
# pt_info['MAX_DT'] = pt_info[['TDM_REQ_DATE','TDM_RES_DATE']].max(axis=1).map(lambda x:(datetime.strptime(x,'%Y-%m-%d') + timedelta(days=92)).strftime('%Y-%m-%dT%H:%M'))
# pt_info['MIN_DT'] = pt_info[['TDM_REQ_DATE','TDM_RES_DATE']].min(axis=1).map(lambda x:(datetime.strptime(x,'%Y-%m-%d') - timedelta(days=61)).strftime('%Y-%m-%dT%H:%M'))
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
# uniq_sampling_pids = list(sampling_result_df.drop_duplicates(['ID'])['ID'].astype(str))
uniq_sampling_pids = list()

final_df = list()
conc_samp_mismatch_pids = list()
no_dose_pids = set()
ambiguous_dosing_time_pids = set()
partly_no_dose_pids = set()
re_dose_pids = set()
error_passing_pids = set()
more_conc_than_samp = set()
outpatient_list = ['14047969']
for finx, fpath in enumerate(pt_files): #break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    # if pid=='13431066':
    #     raise ValueError
    # else:
    #     continue
    # pt_df.append({'ID':pid,'NAME':pname})
    id_info_df = pt_info_dup[pt_info_dup['ID'] == pid].copy()
    min_date = (datetime.strptime(id_info_df['TDM_REQ_DATE'].iloc[0],'%Y-%m-%d')-timedelta(days=62)).strftime('%Y-%m-%d')
    max_date = (datetime.strptime(id_info_df['TDM_RES_DATE'].iloc[0],'%Y-%m-%d')+timedelta(days=92)).strftime('%Y-%m-%d')

    id_dose_df = dose_result_df[dose_result_df['ID'] == pid].copy()
    # raise ValueError
    if len(id_dose_df)==0:
        print('Dose 정보 부재한 사람 제외')
        no_dose_pids.add(pid)
        continue
        # id_dose_df = dose_df
    if len(id_dose_df[id_dose_df['DOSE_TIME']=='TNN:NN']) > 0:
        print('Dosing TIME이 NN:NN으로 불명확한 경우 제외')
        ambiguous_dosing_time_pids.add(pid)
        continue
    id_conc_df = conc_result_df[conc_result_df['ID'] == pid].copy()
    id_conc_df = id_conc_df[(id_conc_df['오더일']>=min_date)&(id_conc_df['오더일']<=max_date)&(id_conc_df['보고일']>=min_date)&(id_conc_df['보고일']<=max_date)]
    id_samp_df_ori = sampling_result_df[(sampling_result_df['ID'] == pid)].copy()
    id_samp_df = id_samp_df_ori[(id_samp_df_ori['채혈DT'] != '') & (~id_samp_df_ori['채혈DT'].isna()) & (id_samp_df_ori['시행DT']<=max_date)&(id_samp_df_ori['시행DT']>=min_date)]
    if len(id_samp_df)>0:
        uniq_sampling_pids.append(pid)

    samp_na = id_samp_df_ori[id_samp_df_ori['채혈DT'].isna()].copy()
    samp_na['채혈DT'] = samp_na['시행DT']
    samp_not_na = id_samp_df_ori[~id_samp_df_ori['채혈DT'].isna()].copy()
    id_samp_df_ori = pd.concat([samp_na, samp_not_na])
    # raise ValueError
    # id_samp_df_ori['채혈DT']



    # res_frag_df = id_conc_df.sort_values(['오더일'], ascending=False)
    res_frag_df = id_conc_df.copy()
    # res_frag_df = id_conc_df.sort_values(['오더일'], ascending=True)
    # pd.to_datetime(res_frag_df['오더일']).diff().dt.days
    for c in ['SAMP_DT','채혈DT', '라벨DT', '접수DT', '시행DT', '보고DT','REC_REASON']:
        res_frag_df[c] = ''
        # res_frag_df['CONC']


    ## 10059302 -> 오더비고 자료도 수집할 필요 있겠음.
    ## 10048706 -> 채혈시간의 peak/trough 위치가 있는데(2004-01-28T21:14 / 2004-01-29T20:39), 농도값은 서로 반대로 배치됨 -> 코드 수정필요
    ## 10053140 -> 채혈일시랑 라벨출력일시랑 같음. 채혈시간 기록을 완전히 믿기 어려울듯 -> trough 농도 레벨인데 dosing과 비교했을때 애매한 위치 (dosing 3시간 후) 측정됨


    samp_time_in_conc_etc_info_dict = {'16685088':['2008-10-29T09:05','2008-10-29T07:30','2008-10-21T09:50','2008-10-21T06:45','2008-10-18T09:30','2008-10-18T07:50','2008-10-14T10:35','2008-10-14T08:35','2008-10-04T11:30','2008-10-04T09:20'],
                                       '19868921':['2012-03-01T19:30','2012-03-01T17:30','2012-02-24T17:30','2012-02-24T19:30','2012-02-21T15:30','2012-02-21T14:00','2012-02-15T14:00','2012-02-15T12:30','2012-02-14T08:00','2012-02-13T11:45','2012-02-11T14:05','2012-02-11T12:20','2012-02-07T14:00','2012-02-07T12:30','2012-02-03T14:00','2012-02-03T12:30'],
                                       '21249383':['2013-10-02T04:40','2013-10-02T06:45','2013-08-19T11:10','2013-08-19T08:20','2013-08-12T09:59','2013-08-12T08:15','2013-08-05T12:00','2013-08-05T08:30','2013-07-29T10:45','2013-07-29T08:40','2013-07-25T11:50','2013-07-25T08:30','2013-07-18T21:00','2013-07-18T23:30','2013-07-18T10:30','2013-07-18T08:50','2013-07-12T09:30','2013-07-12T07:30','2013-07-08T06:15','2013-07-08T08:15','2013-07-04T04:25','2013-07-04T06:40'],
                                       '21393684':['2012-11-23T20:30','2012-11-23T22:00','2012-11-20T10:00','2012-11-20T09:00'],

                                       '10533576':['2004-09-14T14:00','2004-09-14T13:00','2004-09-13T14:30','2004-09-13T12:43','2004-09-08T01:33','2004-09-08T08:57'],
                                       '10610651':['2019-08-26T05:53','2019-08-26T04:53']
                                       # '10004009':[]
                                       # id_samp_df_ori['채혈DT']
                                       # res_frag_df[['CONC']]
                                       #  conc_ord_date_type
                                       }
    if pid in list(samp_time_in_conc_etc_info_dict.keys()):
        # res_frag_df[['보고일','오더일','CONC','SAMP_DT']]
        # res_frag_df['POT채혈DT'] = samp_time_in_conc_etc_info_dict[pid]
        if len(res_frag_df)==len(samp_time_in_conc_etc_info_dict[pid]):
            res_frag_df['POT채혈DT'] = samp_time_in_conc_etc_info_dict[pid]
        else:
            print('길이가 안 맞음. 다시 확인')
            raise ValueError

        # id_samp_df_ori
        # samp_date_rows
        # id_samp_df_ori

    samp_rearr_to_ori_list = ['19447177']
    if pid in samp_rearr_to_ori_list:
        id_samp_df = id_samp_df_ori.copy()

    # if finx
    # raise ValueError
    # res_frag_df.columns

    # POT채혈DT 있는 경우
    pot_dt_rows = res_frag_df[(res_frag_df['POT채혈DT']!='') & (~res_frag_df['POT채혈DT'].isna())].copy()
    if len(pot_dt_rows)!=0:
        print(f"({finx}) {pname} / {pid} / POT채혈DT 결과 존재")
        for potdt_inx, potdt_row in pot_dt_rows.iterrows(): #break
            # pot_dt_rows['CONC']
            # pot_dt_rows['POT채혈DT']
            # POT채혈DT에 날짜도 써 있는 경우 -> 해당 날짜로 SAMP_DT 입력
            if bool(re.search(r'\d\d\d\d-\d\d-\d\dT\d\d:\d\d',potdt_row['POT채혈DT'])):
                if res_frag_df.at[potdt_inx, 'SAMP_DT'] == '':
                    res_frag_df.at[potdt_inx, 'SAMP_DT'] = potdt_row['POT채혈DT']
                    res_frag_df.at[potdt_inx, 'REC_REASON'] += '랩결과옆(POT채혈날짜시간)/'
                else:
                    res_frag_df.at[potdt_inx, 'REC_REASON'] += '랩결과옆(POT채혈날짜시간)/ 이미 DATETIME 존재'
                    continue
                error_passing_pids.add(pid)
            # res_frag_df['POT채혈DT']
            # res_frag_df['SAMP_DT']
            # POT채혈DT에 시간만 써 있는 경우 -> SAMPLING DATA 있는 경우는 아래서 처리. 없는 경우는 여기서 처리
            elif bool(re.search(r'T\d\d:\d\d',potdt_row['POT채혈DT'])):
                if res_frag_df.at[potdt_inx, 'SAMP_DT'] == '':
                    if pid in uniq_sampling_pids:
                        print("POT채혈DT의 TIME만 결과 존재 / SAMPLING DATA 있어 아래서 처리")
                        error_passing_pids.add(pid)
                        res_frag_df.at[potdt_inx, 'REC_REASON'] += '랩결과옆(POT채혈시간만)/'
                    else:
                        print("POT채혈DT의 TIME만 결과 존재 하는데 / SAMPLING DATA 없는 경우 입니다")
                        # POT채혈DT, Dosing 기록으로 유추해야.

                        # 농도측정의 오더일과 보고일이 같은 날짜인 경우는 그 날짜로 사용
                        if potdt_row['오더일']==potdt_row['보고일']:
                            if res_frag_df.at[potdt_inx, 'SAMP_DT'] == '':
                                res_frag_df.at[potdt_inx, 'SAMP_DT'] = potdt_row['오더일'] + 'T' + potdt_row['POT채혈DT'].split('T')[-1]
                                res_frag_df.at[potdt_inx, 'REC_REASON'] += '랩결과옆(POT채혈시간만)/(NO SAMP)오더일==보고일로 추정/'

                        # 농도측정의 오더일과 보고일이 다른 날짜들의 경우
                        ord_noteq_rep_rows = res_frag_df[(res_frag_df['오더일'] != res_frag_df['보고일'])&(res_frag_df['SAMP_DT']=='')].copy()
                        if len(ord_noteq_rep_rows)>0:
                            # 농도채혈 오더 난 날짜가 1개만 있을때
                            if len(ord_noteq_rep_rows['오더일'].drop_duplicates()) == 1:
                                for oneqr_inx, oneqr_row in ord_noteq_rep_rows.iterrows():
                                    res_frag_df.at[oneqr_inx, 'SAMP_DT'] = oneqr_row['오더일'] + 'T' + potdt_row['POT채혈DT'].split('T')[-1]
                                    res_frag_df.at[oneqr_inx, 'REC_REASON'] += '랩결과옆(채혈시간만)/(NO SAMP)오더일!=보고일,오더일로 추정/'
                            # 농도채혈 오더 난 날짜가 1개 이상 존재할 때 있을때
                            else:
                                print("POT채혈DT 결과 존재 / SAMPLING DATA 없음 / 농도채혈 오더 난 날짜가 1개 이상 존재")
                                raise ValueError


                        # error_passing_pids.add(pid)


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
                    mean_conc = (res_frag_df['CONC'].max()+ res_frag_df['CONC'].min())/2
                    date_dose_rows = id_dose_df[id_dose_df['DOSE_DATE']==potdt_row['POT채혈DT']].copy()
                    if len(date_dose_rows)<=2:
                        est_conc_dt_tups = ((datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')-timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                                            (datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')+timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
                                            )
                        # est_conc_dt_tups = (datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')-timedelta(minutes=30),date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')+timedelta(minutes=30))

                        if res_frag_df.at[potdt_inx, 'CONC'] < mean_conc:
                            res_frag_df.at[potdt_inx, 'SAMP_DT'] = potdt_row['POT채혈DT'] + 'T' + est_conc_dt_tups[0].split('T')[-1]
                            res_frag_df.at[potdt_inx, 'REC_REASON'] += '랩결과옆(채혈날짜만)/(NO SAMP)평균농도이상,해당일Peak로배정/'
                        else:
                            res_frag_df.at[potdt_inx, 'SAMP_DT'] = potdt_row['POT채혈DT'] + 'T' + est_conc_dt_tups[1].split('T')[-1]
                            res_frag_df.at[potdt_inx, 'REC_REASON'] += '랩결과옆(채혈날짜만)/(NO SAMP)평균농도미만,해당일Trough로배정/'
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
        if len(res_frag_df[(res_frag_df['SAMP_DT']=='')&(res_frag_df['POT채혈DT']!='')])>0:
            print('POT채혈DT 존재하는 잘못된 케이스 있나 확인')
            raise ValueError
        # DOSING 데이터 없는 경우
        if len(id_dose_df) == 0:
            print(f"({finx}) {pname} / {pid} / Only Conc data / No Dose Data (Check!)")
            no_dose_pids.add(pid)
            continue
        else:
            print(f"({finx}) {pname} / {pid} / Only Conc data (Dosing: O, Sampling X, POT채혈DT: X")

        # res_frag_df.columns
        # id_info_df
        # res_frag_df[(res_frag_df['SAMP_DT']=='')]

        ord_eq_rep_rows = res_frag_df[(res_frag_df['SAMP_DT'] == '') & (res_frag_df['오더일'] == res_frag_df['보고일'])].copy()
        # 농도측정의 오더일과 보고일이 같은 날짜인 경우는 그 날짜로 사용
        for oer_inx, oer_row in ord_eq_rep_rows.iterrows(): #break

            mean_conc = (res_frag_df['CONC'].max()+res_frag_df['CONC'].min())/2
            date_dose_rows = id_dose_df[id_dose_df['DOSE_DATE'] == oer_row['오더일']].copy()
            date_dose_rows['TROUGH_DT'] = date_dose_rows['DOSE_DT'].map(lambda x:(datetime.strptime(x, '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))
            date_dose_rows['PEAK_DT'] = date_dose_rows['DOSE_DT'].map(lambda x:(datetime.strptime(x, '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))

            if len(date_dose_rows) == 0:
                pass
            elif (len(date_dose_rows) > 0) and (len(date_dose_rows) <= 5) and (len(date_dose_rows) >= 1):
            # elif (len(date_dose_rows) <= 2) and (len(date_dose_rows) >= 1):

                # est_conc_dt_tups = ((datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'], '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                #                     (datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'], '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
                #                     )
                # est_conc_dt_tups = (datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')-timedelta(minutes=30),date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')+timedelta(minutes=30))

                if res_frag_df.at[oer_inx, 'CONC'] < mean_conc:
                    dd_col = 'TROUGH_DT'
                else:
                    dd_col = 'PEAK_DT'
                for dd_inx, dd_row in date_dose_rows.iterrows():
                    if dd_row[dd_col] not in list(res_frag_df['SAMP_DT']):
                        res_frag_df.at[oer_inx, 'SAMP_DT'] = oer_row['오더일'] + 'T' + dd_row[dd_col].split('T')[-1]
                        res_frag_df.at[oer_inx, 'REC_REASON'] += f'(NO SAMP,YES DOSE)오더일==보고일로, 해당일{dd_col}로추정/'
                        break
                    else:
                        continue
                        # date_dose_rows[['DOSE','DOSE_DT']]
                    # res_frag_df['SAMP_DT']
                    # res_frag_df.at[oer_inx, 'SAMP_DT'] = oer_row['오더일'] + 'T' + est_conc_dt_tups[1].split('T')[-1]
            else:
                print(f"({finx}) {pname} / {pid} / Only Conc data (Dosing: O, Sampling X, POT채혈DT: X / 해당 날짜 dose기록 5개 이상")
                raise ValueError

        # 농도측정의 오더일과 보고일이 다른 날짜들의 경우
        ord_noteq_rep_rows = res_frag_df[(res_frag_df['SAMP_DT'] == '') & (res_frag_df['오더일'] != res_frag_df['보고일'])].copy()

        conc_order_date_list = list(ord_noteq_rep_rows['오더일'].unique())
        conc_report_date_list = list(ord_noteq_rep_rows['보고일'].unique())
        # if pid=='10042506':
        #     raise ValueError
        # 농도채혈 order가 난 날짜가 1개만 있을때 (dosing 타임과 농도데이터 개수 및 농도 값 고려해서 배열)
        if len(conc_order_date_list) == 0:
            pass
        elif len(conc_order_date_list) == 1:
            # 해당 Conc 에 매치되는 Dosing 오더 범위 설정 (가장 넓은 범위)
            temp_id_dose_df = id_dose_df[(id_dose_df['DOSE_DATE'] >= conc_order_date_list[0])&(id_dose_df['DOSE_DATE'] <= min(conc_report_date_list))].copy()

            # 오더일과 같은날 Dosing된 날짜 존재시 -> Conc 측정 날짜를 오더일에 되었다고 추정
            if len(temp_id_dose_df[(temp_id_dose_df['DOSE_DATE'] == conc_order_date_list[0])])>0:
                temp_id_dose_df = temp_id_dose_df[(temp_id_dose_df['DOSE_DATE'] == conc_order_date_list[0])].copy()
                ascending_orders = [True, True]
                rec_date_reason = '오더일'

            # 보고일과 같은날 Dosing된 날짜 존재시 -> Conc 측정 날짜를 보고일에 되었다고 추정
            elif len(temp_id_dose_df[(temp_id_dose_df['DOSE_DATE'] == conc_report_date_list[0])])>0:
                temp_id_dose_df = temp_id_dose_df[(temp_id_dose_df['DOSE_DATE'] == conc_report_date_list[0])].copy()
                ascending_orders = [True, True]
                rec_date_reason = '보고일'

            # 일치되는 날짜가 없을 경우 -> DOSE 날짜 범위 중 가장 Latest 에 Conc 측정되었다고 추정
            else:
                ascending_orders = [False, True]
                rec_date_reason = '범위내DOSING날짜없어LATEST투약일'

            mean_conc = (res_frag_df['CONC'].max()+res_frag_df['CONC'].min())/2
            # temp_id_dose_df = temp_id_dose_df.sort_values(['DOSE_DATE','DOSE_TIME'],ascending=ascending_orders)

            temp_id_dose_df = temp_id_dose_df.sort_values(['DOSE_DATE', 'DOSE_TIME'], ascending=ascending_orders)
            temp_id_dose_df['TROUGH_DT'] = temp_id_dose_df['DOSE_DT'].map(lambda x: (datetime.strptime(x, '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))
            temp_id_dose_df['PEAK_DT'] = temp_id_dose_df['DOSE_DT'].map(lambda x: (datetime.strptime(x, '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))

            if len(temp_id_dose_df)==0:
                print(f"({finx}) {pname} / {pid} / 농도 측정 오더일 1개, 일부 CONC 날짜의 DOSING DATA 부재함")
                partly_no_dose_pids.add(pid)
            else:
                # ord_noteq_rep_rows['CONC']
                # 측정된 농도 값이 2개 이하이면 DOSING 타임 및 CONC 값 고려하여 배열
                if len(ord_noteq_rep_rows) <= 2:
                    for resf_inx, resf_row in ord_noteq_rep_rows.iterrows():
                        if res_frag_df.at[resf_inx, 'CONC'] < mean_conc:
                            dd_col = 'TROUGH_DT'
                        else:
                            dd_col = 'PEAK_DT'
                        for dd_inx, dd_row in temp_id_dose_df.iterrows():
                            if dd_row[dd_col] not in list(res_frag_df['SAMP_DT']):
                                res_frag_df.at[resf_inx, 'SAMP_DT'] = dd_row[dd_col]
                                res_frag_df.at[resf_inx, 'REC_REASON'] += f'(NO SAMP,YES DOSE)오더일!=보고일로, {rec_date_reason}에{dd_col}로추정/'
                                break
                            else:
                                continue

                    # est_dose_dt = datetime.strptime(temp_id_dose_df.iloc[0]['DOSE_DT'], '%Y-%m-%dT%H:%M')
                    # est_conc_dt_tups = ((est_dose_dt - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                    #                     (est_dose_dt + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))
                    #
                    # if ((est_conc_dt_tups[0] in list(res_frag_df['SAMP_DT'])) or (est_conc_dt_tups[1] in list(res_frag_df['SAMP_DT']))) and (len(temp_id_dose_df)>1):
                    #     est_dose_dt = datetime.strptime(temp_id_dose_df.iloc[1]['DOSE_DT'], '%Y-%m-%dT%H:%M')
                    #     est_conc_dt_tups = ((est_dose_dt - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                    #                         (est_dose_dt + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))

                    # ord_noteq_rep_rows = ord_noteq_rep_rows.sort_values(['CONC'])
                    # for resf_inx, resf_row in ord_noteq_rep_rows.reset_index(drop=False).iterrows():  # break
                    #     if (resf_inx==0) and (est_conc_dt_tups[resf_inx] not in list(res_frag_df['SAMP_DT'])):
                    #         res_frag_df.at[resf_row['index'], 'SAMP_DT'] = est_conc_dt_tups[resf_inx]
                    #     elif (resf_inx==0) and (est_conc_dt_tups[resf_inx] in list(res_frag_df['SAMP_DT'])):
                    #         res_frag_df.at[resf_row['index'], 'SAMP_DT'] = est_conc_dt_tups[resf_inx+1]
                    #     elif (resf_inx == 1) and (est_conc_dt_tups[resf_inx] not in list(res_frag_df['SAMP_DT'])):
                    #         res_frag_df.at[resf_row['index'], 'SAMP_DT'] = est_conc_dt_tups[resf_inx]
                    #     else:
                    #         print(f'({finx}) {pname} / {pid} / 기록할 인덱스가 넘어감')
                    #         raise ValueError
                            
                        
                        # res_frag_df['SAMP_DT']
                else:
                    print(f"({finx}) {pname} / {pid} / 농도 측정 오더일 1개, 도스 3개 이상")
                    raise ValueError

        # 농도 측정 오더 날짜 2개 이상
        else:
            print(f"({finx}) {pname} / {pid} / 농도 측정 오더 날짜 2개 이상")
            mean_conc = (res_frag_df['CONC'].max()+res_frag_df['CONC'].min())/2
            tdm_relate_dates = set(id_info_df['TDM_RES_DATE']).union(set(id_info_df['TDM_REQ_DATE']))
            if len(tdm_relate_dates.intersection(set(ord_noteq_rep_rows['오더일']))) >= len(tdm_relate_dates.intersection(set(ord_noteq_rep_rows['보고일']))):
                conc_ord_date_type = '오더일'
            else:
                conc_ord_date_type = '보고일'
            for concord_date in ord_noteq_rep_rows[conc_ord_date_type].drop_duplicates(): #break
                concord_date_rows = ord_noteq_rep_rows[ord_noteq_rep_rows[conc_ord_date_type]==concord_date] #break
                # ord_noteq_rep_rows[ord_noteq_rep_rows]
                tdm_info_row = {'TDM_REQ_DATE':'0000-00-00','TDM_RES_DATE':'9999-99-99'}
                try: tdm_info_row = id_info_df[id_info_df['TDM_REQ_DATE'] == concord_date].iloc[0]
                except:
                    try: tdm_info_row = id_info_df[id_info_df['TDM_RES_DATE'] == concord_date].iloc[0]
                    except:
                        pass
                        # raise ValueError

                # len()
                tdm_date_tups = (tdm_info_row['TDM_REQ_DATE'],tdm_info_row['TDM_RES_DATE'])
                min_ord_date = concord_date_rows['오더일'].min()
                max_ord_date = concord_date_rows['보고일'].max()

                dose_ord_frag = id_dose_df[(id_dose_df['DOSE_DATE'] >= min_ord_date)&(id_dose_df['DOSE_DATE'] <= max_ord_date)].copy()
                dose_ord_frag = dose_ord_frag[(dose_ord_frag['DOSE_DATE'] >= tdm_date_tups[0])&(dose_ord_frag['DOSE_DATE'] <= tdm_date_tups[1])].copy()
                dose_ord_frag['TROUGH_DT'] = dose_ord_frag['DOSE_DT'].map(lambda x: (datetime.strptime(x, '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))
                dose_ord_frag['PEAK_DT'] = dose_ord_frag['DOSE_DT'].map(lambda x: (datetime.strptime(x, '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))

                # id_info_df
                # id_samp_df
                if len(dose_ord_frag)==0:
                    print(f"({finx}) {pname} / {pid} / 농도 측정된 날짜 범위 중 DOSE 기록이 없음")
                    re_dose_pids.add(pid)
                else:
                    # concord_date_rows
                    # dose_ord_frag

                    # 측정된 농도 값이 2개 이하이면 DOSING 타임 및 CONC 값 고려하여 배열
                    if len(concord_date_rows) <= 2:
                        # raise ValueError
                        print(f"({finx}) {pname} / {pid} / (부분 CONC 날짜 기준) 농도 측정 오더 날짜 2개 이하, DOSE ORD 존재")

                        # est_dose_dt = datetime.strptime(dose_ord_frag.iloc[0]['DOSE_DT'], '%Y-%m-%dT%H:%M')
                        # est_conc_dt_tups = ((est_dose_dt - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                        #                     (est_dose_dt + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))
                        # res_frag_df['SAMP_DT']
                        for resf_inx, resf_row in concord_date_rows.iterrows():
                            if res_frag_df.at[resf_inx, 'CONC'] < mean_conc:
                                dd_col = 'TROUGH_DT'
                            else:
                                dd_col = 'PEAK_DT'
                            for dd_inx, dd_row in dose_ord_frag.iterrows():
                                if dd_row[dd_col] not in list(res_frag_df['SAMP_DT']):
                                    res_frag_df.at[resf_inx, 'SAMP_DT'] = dd_row[dd_col]
                                    res_frag_df.at[resf_inx, 'REC_REASON'] += f'(NO SAMP,YES DOSE) 오더일!=보고일(농도측정오더일 2개 이상, {conc_ord_date_type}기준FITTING, {dd_col}로 추정/'
                                    break
                                else:
                                    # res_frag_df.at[resf_inx, 'REC_REASON'] += f'(NO SAMP,YES DOSE) 오더일!=보고일(농도측정오더일 2개 이상, {conc_ord_date_type}기준FITTING, {dd_col}로 추정/이미 DATETIME 존재'
                                    continue

                        # concord_date_rows = concord_date_rows.sort_values(['CONC']).reset_index(drop=False)
                        #
                        # for resf_inx, resf_row in concord_date_rows.iterrows():  # break
                        #     res_frag_df.at[resf_row['index'], 'SAMP_DT'] = est_conc_dt_tups[resf_inx]
                        #     # res_frag_df['SAMP_DT']
                    else:
                        print(f"({finx}) {pname} / {pid} / 농도 측정 오더 날짜 1개, 측정된 농도 값이 3개 이상")
                        raise ValueError


                    # raise ValueError

            # if pid not in ['10042506', '10150089']:
            #     print('정상작동중')
            #     raise ValueError

#     final_df.append(res_frag_df)
# final_df = pd.concat(final_df, ignore_index=True)
# final_df.to_csv(f"{output_dir}/final_conc_df(with sampling).csv", encoding='utf-8-sig', index=False)


    # print(f"{} {}")

    # CONC, SAMP 데이터 둘 다 있는 경우
    elif (pid in uniq_conc_pids) and (pid in uniq_sampling_pids):
        # continue
        # print(f"({finx}) {pname} / {pid} / CONC, SAMP 데이터 둘 다 존재")
        # raise ValueError

        cdf = res_frag_df[res_frag_df['SAMP_DT'] == ''].copy()
        pdf = cdf[cdf['POT채혈DT'] != ''].copy()
        ddf = id_dose_df.copy()
        # sdf = id_samp_df[['보고일', '채혈DT', '라벨DT', '접수DT', '시행DT', '보고DT']].copy()
        sdf = id_samp_df.copy()

    # sdf = id_samp_df.copy()
    #
    #     res_frag_df['CONC']
    #     id_samp_df['채혈DT']
    #     id_info_df

        print(f"({finx}) {pname} / {pid} / SAMP 데이터 존재")
        mean_conc = (res_frag_df['CONC'].max()+res_frag_df['CONC'].min())/2
        # tdm_relate_dates = set(id_info_df['TDM_RES_DATE']).union(set(id_info_df['TDM_REQ_DATE']))
        samp_relate_dates = set(sdf['보고일'])
        if len(samp_relate_dates.intersection(set(cdf['오더일']))) > len(samp_relate_dates.intersection(set(cdf['보고일']))):
            conc_ord_date_type = '오더일'
            print(f"선택된 날짜 타입이 '{conc_ord_date_type}'이 더 잘 맞음")
            raise ValueError
        elif len(samp_relate_dates.intersection(set(cdf['오더일']))) == len(samp_relate_dates.intersection(set(cdf['보고일']))):
            conc_ord_date_type = '오더일'
            other_type = '보고일'

            # conc_ord_date_type = '보고일'
            # other_type = '오더일'
        else:
            conc_ord_date_type = '보고일'
            other_type = '오더일'
        for concord_date in cdf[conc_ord_date_type].drop_duplicates(): # break
            concord_date_rows = cdf[cdf[conc_ord_date_type] == concord_date].copy()  # break
            samp_date_rows = sdf[sdf[conc_ord_date_type] == concord_date].copy()  # break

            # id_samp_df_ori['채혈DT']
            # concord_date_rows['CONC']
            # samp_date_rows['채혈DT']
            # SAMPLING TIME DATA가 부분적으로 없는 경우
            if len(samp_date_rows) == 0:
                # sdf = id_samp_df_ori.copy()
                print("SAMPLING TIME DATA가 부분적으로 없는 경우")
                samp_date_rows = id_samp_df_ori[id_samp_df_ori[conc_ord_date_type] == concord_date].copy()
                if len(samp_date_rows)>0:
                    pass
                else:

                    other_type_dates = list(set(concord_date_rows[other_type]).intersection(set(id_samp_df_ori[conc_ord_date_type]) - set(concord_date_rows[conc_ord_date_type])))

                    if len(other_type_dates)==1:
                        print("SAMPLING TIME DATA가 부분적으로 없는 경우 / '오더일' '보고일' 바꿔서 해볼 수 있는 경우")
                        samp_date_rows = id_samp_df_ori[id_samp_df_ori[other_type] == other_type_dates[0]].copy()
                        # samp_date_rows['채혈DT']
                        if len(samp_date_rows)==0:
                            print("SAMPLING TIME DATA가 부분적으로 없는 경우 / '오더일' '보고일' 바꿔서 해볼 수 있는 경우 / 바꿔도 날짜가 안나옴")
                            raise ValueError
                    elif len(other_type_dates)>1:
                        print("SAMPLING TIME DATA가 부분적으로 없는 경우 / '오더일' '보고일' 바꿔서 해볼 수 있는 경우 / 날짜가 여러 개 나와서 정하기 어려움")
                        raise ValueError
                    else:
                        print("SAMPLING TIME DATA가 부분적으로 없는 경우 / ori sample file에도 없음")
                        # id_dose_df
                        # id_samp_df['채혈DT']
                        # id_conc_df
                        # concord_date_rows['POT채혈DT']

                        if pid in outpatient_list:
                            print("외래 베이스로 투여된 환자라 정확한 dosing time 추정 어려움")
                            res_frag_df['REC_REASON']='외래 베이스 환자라 정확한 dosing time 추정 어려움'
                            continue
                        else:
                            raise ValueError
                # id_info_df
                # id_samp_df['채혈DT']
                # id_samp_df_ori['채혈DT']
                # id_conc_df

            # ord_noteq_rep_rows[ord_noteq_rep_rows]
            # tdm_info_row = {'TDM_REQ_DATE': '0000-00-00', 'TDM_RES_DATE': '9999-99-99'}
            # try: tdm_info_row = id_info_df[id_info_df['TDM_REQ_DATE'] == concord_date].iloc[0]
            # except:
            #     try:tdm_info_row = id_info_df[id_info_df['TDM_RES_DATE'] == concord_date].iloc[0]
            #     except:
            #         pass
                    # raise ValueError



            # len()
            # tdm_date_tups = (tdm_info_row['TDM_REQ_DATE'], tdm_info_row['TDM_RES_DATE'])
            # samp_date_tups = (samp_date_rows['TDM_REQ_DATE'], samp_date_rows['TDM_RES_DATE'])
            min_ord_date = samp_date_rows['오더일'].min()
            max_ord_date = samp_date_rows['보고일'].max()

            dose_ord_frag = id_dose_df[(id_dose_df['DOSE_DATE'] >= min_ord_date) & (id_dose_df['DOSE_DATE'] <= max_ord_date)].copy()
            # dose_ord_frag = dose_ord_frag[(dose_ord_frag['DOSE_DATE'] >= tdm_date_tups[0]) & (dose_ord_frag['DOSE_DATE'] <= tdm_date_tups[1])].copy()
            # dose_ord_frag['TROUGH_DT'] = dose_ord_frag['DOSE_DT'].map(lambda x: (datetime.strptime(x, '%Y-%m-%dT%H:%M') - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))
            # dose_ord_frag['PEAK_DT'] = dose_ord_frag['DOSE_DT'].map(lambda x: (datetime.strptime(x, '%Y-%m-%dT%H:%M') + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'))

            ## 농도 기록과 채혈 시간 기록이 일치할경우 -> Dose 투약 시점에 맞춰 Conc 크기에 따라 배열
            if (len(concord_date_rows)==len(samp_date_rows)) and (len(samp_date_rows) <= 2):
                if (len(samp_date_rows) == 1):
                    print('CONC 데이터 길이 == SAMP 데이터 길이 / samp_date_rows 기록이 1개만 존재')

                    # dose_before_samps = id_dose_df[(id_dose_df['DOSE_DT'] < samp_date_rows['채혈DT'].min())].copy()
                    # dose_after_samps = id_dose_df[(id_dose_df['DOSE_DT'] > samp_date_rows['채혈DT'].min())].copy()
                    # raise ValueError
                    concord_date_rows['채혈DT'] = samp_date_rows['채혈DT'].iloc[0]

                    for dd_inx, dd_row in concord_date_rows.iterrows():
                        if dd_row['채혈DT'] not in list(res_frag_df['SAMP_DT']):
                            # print('기록')
                            res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                            res_frag_df.at[dd_inx, 'REC_REASON'] += f'(YES SAMP,YES DOSE) 채혈 기록 1개, 농도 및 채혈 기록일치/'
                            # break
                        else:
                            # res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                            # res_frag_df['채혈DT']
                            # res_frag_df['SAMP_DT']
                            # res_frag_df['POT채혈DT']
                            res_frag_df.at[dd_inx, 'REC_REASON'] += f'(YES SAMP,YES DOSE) 채혈 기록 1개, 농도 및 채혈 기록일치/ 이미 DATETIME 존재'

                            print('이미 날짜 DT가 존재합니다')
                else:
                    # if pid=='10024058':
                    #     raise ValueError

                    # samp_date_rows['채혈DT']
                    dose_before_samps = id_dose_df[(id_dose_df['DOSE_DT'] < samp_date_rows['채혈DT'].min())].copy()
                    dose_between_samps = id_dose_df[(id_dose_df['DOSE_DT'] >= samp_date_rows['채혈DT'].min())&(id_dose_df['DOSE_DT'] <= samp_date_rows['채혈DT'].max())].copy()
                    dose_after_samps = id_dose_df[(id_dose_df['DOSE_DT'] > samp_date_rows['채혈DT'].min())].copy()
                    # 채혈 시간 기록 사이에 Dose 투약 기록이 있는 경우 -> 채혈 농도값으로 peak / trough에 해당하는 농도들 구분하고, '채혈DT'를 이에 맞춰 배정 후 SAMP_DT로 확정시킴 (DOSE_DT 앞 뒤에 있는게 맞는지 확인도하면 좋을듯)
                    if len(dose_between_samps)>=2:
                        print('CONC 데이터 길이 == SAMP 데이터 길이 / sampling 사이에 dose 기록이 2개 이상 존재')
                        # dose_between_samps['DOSE_DT']
                        temp_min_samp_dt_str = samp_date_rows['채혈DT'].min()
                        temp_max_samp_dt_str = samp_date_rows['채혈DT'].max()

                        temp_min_samp_dt = datetime.strptime(samp_date_rows['채혈DT'].min(),'%Y-%m-%dT%H:%M')
                        temp_max_samp_dt = datetime.strptime(samp_date_rows['채혈DT'].max(),'%Y-%m-%dT%H:%M')

                        before_samp_max_dt_str = dose_before_samps['DOSE_DT'].max()
                        if type(before_samp_max_dt_str)==str:
                            temp_dose_before_samps_dt = datetime.strptime(before_samp_max_dt_str,'%Y-%m-%dT%H:%M')
                        else:
                            temp_dose_before_samps_dt = datetime.strptime('0001-01-01T00:00','%Y-%m-%dT%H:%M')

                        temp_dose_between_samps_min_dt = datetime.strptime(dose_between_samps['DOSE_DT'].min(), '%Y-%m-%dT%H:%M')
                        temp_dose_between_samps_max_dt = datetime.strptime(dose_between_samps['DOSE_DT'].max(), '%Y-%m-%dT%H:%M')

                        after_samp_min_dt_str = dose_after_samps['DOSE_DT'].min()
                        if type(after_samp_min_dt_str)==str:
                            temp_dose_after_samps_dt = datetime.strptime(after_samp_min_dt_str,'%Y-%m-%dT%H:%M')
                        else:
                            temp_dose_after_samps_dt = datetime.strptime('9999-12-31T23:59','%Y-%m-%dT%H:%M')

                        # 직전까지의 DOSE_DT가 다음 DOSE_DT 보다 시간적으로 가까우면 -> 채혈DT중 min값은 peak conc
                        if ((temp_min_samp_dt - temp_dose_before_samps_dt).total_seconds()/60) < ((temp_dose_between_samps_min_dt - temp_min_samp_dt).total_seconds()/60):
                            min_samp_dt_conc_type = 'peak'
                        else:
                            min_samp_dt_conc_type = 'trough'

                        if ((temp_max_samp_dt - temp_dose_between_samps_max_dt).total_seconds()/60) < ((temp_dose_after_samps_dt - temp_max_samp_dt).total_seconds()/60):
                            max_samp_dt_conc_type = 'peak'
                        else:
                            max_samp_dt_conc_type = 'trough'
                        # id_conc_df
                        concval_saving_rows = concord_date_rows.copy()
                        concord_date_rows.at[concord_date_rows.index[0], 'SAMP_DT'] = temp_min_samp_dt_str
                        if min_samp_dt_conc_type=='peak':
                            concord_date_rows.at[concord_date_rows.index[0],'CONC'] = concval_saving_rows[concval_saving_rows['CONC'] >= mean_conc].iloc[-1]['CONC']
                        else:
                            concord_date_rows.at[concord_date_rows.index[0],'CONC'] = concval_saving_rows[concval_saving_rows['CONC'] < mean_conc].iloc[-1]['CONC']

                        concord_date_rows.at[concord_date_rows.index[-1], 'SAMP_DT'] = temp_max_samp_dt_str
                        if max_samp_dt_conc_type=='peak':
                            concord_date_rows.at[concord_date_rows.index[-1],'CONC'] = concval_saving_rows[concval_saving_rows['CONC'] >= mean_conc].iloc[0]['CONC']
                        else:
                            concord_date_rows.at[concord_date_rows.index[-1],'CONC'] = concval_saving_rows[concval_saving_rows['CONC'] < mean_conc].iloc[0]['CONC']

                        # concord_date_rows[['SAMP_DT','CONC']]
                        for dd_inx, dd_row in concord_date_rows.iterrows():
                            if dd_row['채혈DT'] not in list(res_frag_df['SAMP_DT']):
                                # print('기록')
                                res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                                res_frag_df.at[dd_inx, 'CONC'] = dd_row['CONC']
                                res_frag_df.at[dd_inx, 'REC_REASON'] += f'(YES SAMP,YES DOSE, 두 기록 일치) 채혈 기록 2개 이하, 채혈 사이 DOSING 기록 여러 개, 농도 값과 DOSING 고려해서 배열'
                                # break
                            else:
                                # res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                                # res_frag_df['채혈DT']
                                # res_frag_df['SAMP_DT']
                                # res_frag_df['POT채혈DT']

                                print('이미 날짜 DT가 존재합니다')

                        # res_frag_df['']
                        #
                        # res_frag_df[(res_frag_df[conc_ord_date_type] == concord_date) & (res_frag_df['SAMP_DT'] == '')]

                        # raise ValueError
                    else:

                        # if pid=='10004009':
                        #     raise ValueError

                        # if samp_date_rows:
                        # cond_check_conc_date_rows['CONC']
                        cond_check_conc_date_rows = concord_date_rows.sort_values(['CONC'], ascending=True)
                        COND1 = ((cond_check_conc_date_rows['CONC'].diff().map(np.abs) / (cond_check_conc_date_rows['CONC'].min())) >= 3).reset_index(drop=True) # 3배 이상의 급격한 CONC 변화
                        COND2 = ((cond_check_conc_date_rows['CONC'].diff().map(np.abs) / (cond_check_conc_date_rows['CONC'].min())) < 3).reset_index(drop=True)  # 3배 미만의 완만한 CONC 변화
                        COND3 = ((samp_date_rows['채혈DT'].map(lambda x:datetime.strptime(x,'%Y-%m-%dT%H:%M').timestamp()).diff())/60 <= 90).map(np.abs).reset_index(drop=True)  # 시간 차이가 90분 이내
                        COND4 = (len(dose_between_samps)==0) # Sampling time 사이에 DOSE 기록 없음
                        # # COND4 = (np.abs(res_frag_df['SAMP_DT'].map(lambda x:datetime.strptime(x,'%Y-%m-%dT%H:%M').timestamp()).diff())/60 <= 3) # 시간 차이가 3분 이내
                        #
                        if ((COND2 & COND3).sum() > 0) and COND4:  # 시간이 거의 일치하는데 상대적으로 완만한 변화 -> Toxic 후 감소하는 것으로 간주. 그러나 그 측정시간이 서로 너무 가까운 경우 예측이 어려우므로 하나는 날린다
                            arranged_conc_samp_rows = concord_date_rows.sort_values(['CONC'], ascending=False).iloc[0:1,:].copy()
                            arranged_conc_samp_rows['채혈DT'] = samp_date_rows.iloc[0]['채혈DT']
                            between_dosing_rec_count = '고려하며, 시간이 거의 일치하는데 상대적으로 완만한 변화 -> Toxic 후 감소하는 것으로 간주(시간 너무 가까우면 날림)'
                        elif ((COND1 & COND3).sum() > 0) and COND4:  # 시간이 거의 일치하는데 급격한 변화 -> 가장 가까운 Dosing time 기준 trough, peak로 간주
                            min_diff_dosing_dt = np.abs(id_dose_df['DOSE_DT'].map(lambda x:datetime.strptime(x,'%Y-%m-%dT%H:%M').timestamp()) - datetime.strptime(samp_date_rows.iloc[0]['채혈DT'],'%Y-%m-%dT%H:%M').timestamp())/3600
                            min_inx = min_diff_dosing_dt[min_diff_dosing_dt==min_diff_dosing_dt.min()].index[0]
                            min_diff_dt_str = id_dose_df[id_dose_df.index==min_inx].iloc[0]['DOSE_DT']
                            min_diff_dt = datetime.strptime(min_diff_dt_str,'%Y-%m-%dT%H:%M')
                            est_conc_dt_list = [(min_diff_dt - timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                                                (min_diff_dt + timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M'),
                                                ]
                            # est_conc_dt_tups = (datetime.strptime(date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')-timedelta(minutes=30),date_dose_rows.iloc[0]['DOSE_DT'],'%Y-%m-%dT%H:%M')+timedelta(minutes=30))

                            arranged_conc_samp_rows = cond_check_conc_date_rows.copy()
                            arranged_conc_samp_rows['채혈DT'] = est_conc_dt_list
                            between_dosing_rec_count = '고려하며, 시간이 거의 일치하는데 급격한 변화 -> 가장 가까운 Dosing time 기준 trough, peak로 간주'
                        else:
                            # CONC 데이터 길이 == SAMP 데이터 길이 / sampling 사이에 dose 기록이 한 개 존재
                            if len(dose_between_samps) == 1:
                                conc_ascending_sort = True
                                samp_ascending_sort = True
                                between_dosing_rec_count = '1개'
                                # print(f'CONC 데이터 길이 == SAMP 데이터 길이 / sampling 사이에 dose 기록이 한 개 존재')
                                # raise ValueError

                            # 채혈 시간 기록 사이에 Dose 투약 기록이 없는 경우 -> 채혈 농도값의 크기를 내림차순(큰->작은)으로 배열하고 '채혈DT'를 오름차순으로 맞춰 배정
                            else:
                                conc_ascending_sort = False
                                samp_ascending_sort = True
                                between_dosing_rec_count = '0개'
                            arranged_conc_samp_rows = concord_date_rows.sort_values(['CONC'], ascending=conc_ascending_sort)
                            arranged_conc_samp_rows['채혈DT'] = list(samp_date_rows.sort_values(['채혈DT'], ascending=samp_ascending_sort)['채혈DT'])

                        for dd_inx, dd_row in arranged_conc_samp_rows.iterrows():
                            if dd_row['채혈DT'] not in list(res_frag_df['SAMP_DT']):
                                print('기록')
                                res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                                res_frag_df.at[dd_inx, 'REC_REASON'] += f'(YES SAMP, YES DOSE, 두 기록 일치) 채혈 기록 2개 이하, 채혈 사이 DOSING 기록 {between_dosing_rec_count}, 농도값 고려하여 추정'
                                # break
                            else:
                                res_frag_df.at[dd_inx, 'REC_REASON'] += f'(YES SAMP, YES DOSE, 두 기록 일치) 채혈 기록 2개 이하, 채혈 사이 DOSING 기록 {between_dosing_rec_count}, 농도값 고려하여 추정 / 이미 DATETIME 존재'


                            #     res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                            # res_frag_df['CONC']
                                # res_frag_df['채혈DT']
                                # res_frag_df['SAMP_DT']
                                # res_frag_df['POT채혈DT']

                                print('이미 날짜 DT가 존재합니다')
                                # raise ValueError


            elif (len(concord_date_rows)==len(samp_date_rows)) and (len(samp_date_rows) > 2):
                print('샘플이 3 개 이상')
                dose_between_samps = id_dose_df[(id_dose_df['DOSE_DT'] >= samp_date_rows['채혈DT'].min()) & (id_dose_df['DOSE_DT'] <= samp_date_rows['채혈DT'].max())].copy().sort_values(['DOSE_DT'])
                # concord_date_rows
                # samp_date_rows
                # concord_date_rows['CONC']

                peak_concord_date_rows = concord_date_rows[concord_date_rows['CONC'] >= mean_conc].copy()
                trough_concord_date_rows = concord_date_rows[concord_date_rows['CONC'] < mean_conc].copy()

                dose_samp_date_rows = samp_date_rows.copy()
                dose_samp_date_rows['임시DT'] = dose_samp_date_rows['채혈DT']
                dose_samp_date_rows['임시CONC'] = np.nan

                temp_dose_between_samps = pd.DataFrame(columns=samp_date_rows.columns)
                temp_dose_between_samps['임시DT'] = dose_between_samps['DOSE_DT'].copy()
                temp_dose_between_samps['보고DT'] = ''
                temp_dose_between_samps['임시CONC'] = np.nan

                dose_samp_date_rows = pd.concat([dose_samp_date_rows, temp_dose_between_samps]).sort_values(['임시DT','보고DT'],ascending=[True,False])
                # dose_samp_date_rows['임시']

                prev_dosedt = '0000-00-00T00:00'
                for dosedt_inx, dosedt in enumerate(dose_between_samps['DOSE_DT']): #break
                    # dose_samp_date_rows['임시CONC'].iloc[0]
                    print('사이클')
                    try: next_dosedt = dose_between_samps[dose_between_samps['DOSE_DT'] > dosedt].iloc[0]['DOSE_DT']
                    except: next_dosedt = '9999-99-99T99:99'
                    # id_dose_df

                    prev_ds_frag_row = dose_samp_date_rows[(dose_samp_date_rows['임시DT'] > prev_dosedt)&(dose_samp_date_rows['임시DT'] <= dosedt)&(~dose_samp_date_rows['ID'].isna())].copy()
                    next_ds_frag_row = dose_samp_date_rows[(dose_samp_date_rows['임시DT'] <= next_dosedt) & (dose_samp_date_rows['임시DT'] > dosedt) & (~dose_samp_date_rows['ID'].isna())].copy()

                    if len(prev_ds_frag_row)==0:
                        if len(next_ds_frag_row)==0:
                            pass
                        # id_dose_df
                        elif len(next_ds_frag_row)==1:
                            print('prev_ds_frag_row는 0개 / next_ds_frag_row가 1 개')
                            raise ValueError
                            # temp_inx = next_ds_frag_row.index[0]
                            # dose_samp_date_rows.at[temp_inx, '임시CONC'] = peak_concord_date_rows.iloc[-1]['CONC']
                            # peak_concord_date_rows = peak_concord_date_rows.iloc[:-1, :]
                        elif len(next_ds_frag_row)==2:
                            print('prev_ds_frag_row는 0개 / next_ds_frag_row가 2 개')
                            raise ValueError
                        else:
                            print('prev_ds_frag_row는 0개 / next_ds_frag_row가 3 개 이상')
                            raise ValueError
                    elif len(prev_ds_frag_row)==1:
                        temp_inx = prev_ds_frag_row.index[0]
                        # dose_samp_date_rows['채혈DT']
                        # id_dose_df
                        dose_samp_date_rows.at[temp_inx, '임시CONC'] = trough_concord_date_rows.iloc[-1]['CONC']
                        trough_concord_date_rows = trough_concord_date_rows.iloc[:-1,:]
                        # peak_concord_date_rows
                    elif len(prev_ds_frag_row)==2:

                        # dose_samp_date_rows['채혈DT']
                        # id_dose_df
                        peak_inx = prev_ds_frag_row.index[0]
                        dose_samp_date_rows.at[peak_inx, '임시CONC'] = peak_concord_date_rows.iloc[-1]['CONC']
                        peak_concord_date_rows = peak_concord_date_rows.iloc[:-1, :]

                        trough_inx = prev_ds_frag_row.index[-1]
                        dose_samp_date_rows.at[trough_inx, '임시CONC'] = trough_concord_date_rows.iloc[-1]['CONC']
                        trough_concord_date_rows = trough_concord_date_rows.iloc[:-1,:]

                    else:
                        # concord_date_rows['CONC']

                        # dose_samp_date_rows
                        # peak_concord_date_rows['CONC']
                        # trough_concord_date_rows['CONC']


                        print('DOSE와 DOSE사이 TROUGH 샘플이 3 개 이상 - peak/trough 배열 조정 요망')
                        raise ValueError

                    # samps 사이 데이터상 마지막 dose인 경우
                    if dosedt_inx==len(dose_between_samps)-1:
                        if len(next_ds_frag_row)==0:
                            print('마지막 dose 처리 / next_ds_frag_row가 0 개')
                            pass
                        elif len(next_ds_frag_row)==1:
                            temp_inx = next_ds_frag_row.index[0]
                            # concord_date_rows['CONC']
                            dose_samp_date_rows.at[temp_inx, '임시CONC'] = peak_concord_date_rows.iloc[0]['CONC']
                            peak_concord_date_rows = peak_concord_date_rows.iloc[:-1,:]
                        elif len(next_ds_frag_row)==2:

                            # dose_samp_date_rows['채혈DT']
                            # id_dose_df
                            peak_inx = next_ds_frag_row.index[0]
                            dose_samp_date_rows.at[peak_inx, '임시CONC'] = peak_concord_date_rows.iloc[-1]['CONC']
                            peak_concord_date_rows = peak_concord_date_rows.iloc[:-1, :]

                            trough_inx = next_ds_frag_row.index[-1]
                            if len(trough_concord_date_rows)>0:
                                dose_samp_date_rows.at[trough_inx, '임시CONC'] = trough_concord_date_rows.iloc[-1]['CONC']
                            else:
                                dose_samp_date_rows.at[trough_inx, '임시CONC'] = peak_concord_date_rows.iloc[-1]['CONC']

                            trough_concord_date_rows = trough_concord_date_rows.iloc[:-1,:]
                        else:
                            for temp_inx in next_ds_frag_row.index:
                                dose_samp_date_rows.at[temp_inx, '임시CONC'] = peak_concord_date_rows.iloc[-1]['CONC']
                                peak_concord_date_rows = peak_concord_date_rows.iloc[:-1, :]

                                # trough_concord_date_rows['CONC']
                                # peak_concord_date_rows['CONC']
                                print('마지막 dose 처리 / next_ds_frag_row가 3 개 이상')
                            # raise ValueError
                        decided_dose_samp_date_rows = dose_samp_date_rows[~dose_samp_date_rows['ID'].isna()].copy()
                        decided_dose_samp_date_rows['CONC'] = decided_dose_samp_date_rows['임시CONC']
                        decided_dose_samp_date_rows['SAMP_DT'] = decided_dose_samp_date_rows['임시DT']

                        # decided_dose_samp_date_rows['CONC']
                        # concord_date_rows.columns
                        # concord_date_rows['CONC']

                        # rest_of_samp_dt_rows = res_frag_df[res_frag_df['SAMP_DT']==''].copy()
                        if len(concord_date_rows)==len(decided_dose_samp_date_rows):

                            for col in ['CONC','SAMP_DT',]:
                                concord_date_rows[col] = list(decided_dose_samp_date_rows[col])

                            for dd_inx, dd_row in concord_date_rows.iterrows():
                                if dd_row['채혈DT'] not in list(res_frag_df['SAMP_DT']):
                                    # print('기록')
                                    res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                                    res_frag_df.at[dd_inx, 'CONC'] = dd_row['CONC']
                                    res_frag_df.at[dd_inx, 'REC_REASON'] += f'(YES SAMP, YES DOSE, 두 기록 일치) 채혈 기록 3개 이상, 채혈 사이 DOSING 기록, 농도값 고려하여 추정'

                                    # break
                                else:
                                    res_frag_df.at[dd_inx, 'REC_REASON'] += f'(YES SAMP, YES DOSE, 두 기록 일치) 채혈 기록 3개 이상, 채혈 사이 DOSING 기록, 농도값 고려하여 추정 / 이미 DATETIME 존재'
                                    # res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                                    # res_frag_df['채혈DT']
                                    # res_frag_df['SAMP_DT']
                                    # res_frag_df['POT채혈DT']

                                    print('이미 날짜 DT가 존재합니다')
                        else:
                            # rest_of_samp_dt_rows['CONC']
                            # res_frag_df['SAMP_DT']
                            # id_samp_df
                            # conc_ord_date_type

                            print(f"len(res_frag_df)!=len(decided_dose_samp_date_rows)")
                            raise ValueError

                        # decided_dose_samp_date_rows['CONC']
                        # res_frag_df.columns
                        # res_frag_df['CONC']
                        # for
                        #     res_frag_df.at[temp_inx, '임시CONC'] = peak_concord_date_rows.iloc[0]['CONC']
                        # decided_dose_samp_date_rows['CONC']
                    prev_dosedt = dosedt
                # id_samp_df_ori
                # mean_conc


            #################################        아래는 날렸음        ################################################

            elif len(concord_date_rows) > len(samp_date_rows):   ################# 채혈정보랑 안 맞는 부정확한 정보라 날린다


                print(f'CONC 데이터 길이: {len(concord_date_rows)} > SAMP 데이터 길이: {len(samp_date_rows)}')
                more_conc_than_samp.add(pid)
                # len(more_conc_than_samp)
                continue


                # print("SAMPLING TIME DATA가 부분적으로 없는 경우")
                # check_samp_date_rows = id_samp_df_ori[id_samp_df_ori[conc_ord_date_type] == concord_date].copy()
                # if len(check_samp_date_rows) > 0:
                #     pass

                # if pid=='10533576':
                #
                #     res_frag_df[['SAMP_DT','CONC']]
                #     for concord_inx, concord_row in concord_date_rows:
                #
                #     concord_date_rows['오더일']==concord_date_rows['보고일']:




                dose_before_samps = id_dose_df[(id_dose_df['DOSE_DT'] < samp_date_rows['채혈DT'].min())].copy()
                dose_between_samps = id_dose_df[(id_dose_df['DOSE_DT'] >= samp_date_rows['채혈DT'].min()) & (id_dose_df['DOSE_DT'] <= samp_date_rows['채혈DT'].max())].copy()
                dose_after_samps = id_dose_df[(id_dose_df['DOSE_DT'] > samp_date_rows['채혈DT'].max())].copy()

                # 채혈 시간 기록 사이에 Dose 투약 기록이 있는 경우 -> 채혈 농도값으로 peak / trough에 해당하는 농도들 구분하고, '채혈DT'를 이에 맞춰 배정 후 SAMP_DT로 확정시킴 (DOSE_DT 앞 뒤에 있는게 맞는지 확인도하면 좋을듯)
                if len(dose_between_samps)>0:
                    # id_samp_df

                    # if len(dose_between_samps) >=2:
                    #     print(f'CONC 데이터 길이 ({len(concord_date_rows)}) > SAMP 데이터 길이 ({len(samp_date_rows)}) / sampling 사이에 dose 기록이 여러개 존재')
                    #     raise ValueError

                    # dose_between_samps
                    print(f'CONC 데이터 길이 ({len(concord_date_rows)}) > SAMP 데이터 길이 ({len(samp_date_rows)}) / sampling 사이에 dose 기록이 한 개 존재')


                    # raise ValueError

                    # samp_date_rows['채혈DT']
                    # concord_date_rows['CONC']

                    temp_samp_date_rows = samp_date_rows.copy()
                    temp_samp_date_rows['임시P_T'] = ''
                    # temp_samp_date_rows['채혈DT']
                    # 채혈시간 데이터가 peak일지 trough 일지 앞, 뒤 dosing time 고려하여 결정
                    for samp_inx, samp_row in temp_samp_date_rows.iterrows():  # break
                        samp_dt_str = samp_row['채혈DT']
                        samp_dt = datetime.strptime(samp_dt_str, '%Y-%m-%dT%H:%M')

                        dose_before_samps = id_dose_df[id_dose_df['DOSE_DT'] <= samp_dt_str].copy()
                        dose_after_samps = id_dose_df[id_dose_df['DOSE_DT'] > samp_dt_str].copy()

                        before_samp_max_dt_str = dose_before_samps['DOSE_DT'].max()
                        after_samp_min_dt_str = dose_after_samps['DOSE_DT'].min()

                        if type(before_samp_max_dt_str) == str:
                            temp_dose_before_samps_dt = datetime.strptime(before_samp_max_dt_str, '%Y-%m-%dT%H:%M')
                        else:
                            temp_dose_before_samps_dt = datetime.strptime('0001-01-01T00:00', '%Y-%m-%dT%H:%M')

                        if type(after_samp_min_dt_str) == str:
                            temp_dose_after_samps_dt = datetime.strptime(after_samp_min_dt_str, '%Y-%m-%dT%H:%M')
                        else:
                            temp_dose_after_samps_dt = datetime.strptime('9999-12-31T23:59', '%Y-%m-%dT%H:%M')

                        # 직전까지의 DOSE_DT가 다음 DOSE_DT 보다 시간적으로 가까우면 -> 채혈DT중 min값은 peak conc
                        if ((samp_dt - temp_dose_before_samps_dt).total_seconds() / 60) < ((temp_dose_after_samps_dt - samp_dt).total_seconds() / 60):
                            samp_dt_conc_type = 'peak'
                        else:
                            samp_dt_conc_type = 'trough'

                        temp_samp_date_rows.at[samp_inx, '임시P_T'] = samp_dt_conc_type


                    for dd_inx, dd_row in concord_date_rows.iterrows():  # break
                        # id_conc_df
                        # conc_ord_date_type
                        # concord_date_rows[['보고일','오더일','CONC']]
                        # res_frag_df
                        peak_samp_date_rows = temp_samp_date_rows[temp_samp_date_rows['임시P_T'] == 'peak'].copy()
                        trough_samp_date_rows = temp_samp_date_rows[temp_samp_date_rows['임시P_T'] == 'trough'].copy()
                        used_dt_list = list()
                        if dd_row['CONC'] >= mean_conc:
                            conc_level_type = 'peak'

                            if len(peak_samp_date_rows) != 0:
                                input_samp_date_rows = peak_samp_date_rows.copy()
                            else:
                                # input_samp_date_rows = trough_samp_date_rows.copy()
                                print('peak 부족')
                                raise ValueError
                        else:
                            conc_level_type = 'trough'

                            if len(trough_samp_date_rows) != 0:
                                input_samp_date_rows = trough_samp_date_rows.copy()
                            else:
                                # input_samp_date_rows = peak_samp_date_rows.copy()
                                print('trough 부족')
                                raise ValueError

                        # check_input_samp_date_rows = input_samp_date_rows.copy()
                        input_samp_dt = input_samp_date_rows['채혈DT'].iloc[-1]
                        res_frag_df.at[dd_inx, 'SAMP_DT'] = input_samp_dt
                        res_frag_df.at[dd_inx, 'REC_REASON'] = '날린 곳'
                        used_dt_list.append(input_samp_dt)

                        peak_samp_date_rows = peak_samp_date_rows[~peak_samp_date_rows['채혈DT'].isin(used_dt_list)]
                        trough_samp_date_rows = trough_samp_date_rows[~trough_samp_date_rows['채혈DT'].isin(used_dt_list)]
                        # if input_samp_dt not in list(res_frag_df['SAMP_DT']):
                        #     print('기록')
                        #     res_frag_df.at[dd_inx, 'SAMP_DT'] = input_samp_dt
                        #     used_dt_list.append(input_samp_dt)
                        #     # res_frag_df['SAMP_DT']
                        # else:
                        #     print('이미 날짜 DT가 존재합니다')
                        #     # raise ValueError
                    # print('뒤처리 필요?')
                    # id_samp_df['채혈DT']
                    # res_frag_df['SAMP_DT']
                    # res_frag_df['POT채혈DT']
                    raise ValueError

                # 채혈 시간 기록 사이에 Dose 투약 기록이 없는 경우 -> 채혈 농도값의 크기를 내림차순(큰->작은)으로 배열하고 '채혈DT'를 내림차순으로 맞춰 배정
                else:
                    # id_dose_df
                    # id_conc_df
                    # id_samp_df
                    arranged_conc_samp_rows = concord_date_rows.sort_values(['CONC'], ascending=True)
                    samp_date_rows = samp_date_rows.sort_values(['채혈DT'], ascending=True)

                    # arranged_conc_samp_rows['CONC']
                    # id_conc_df
                    # id_samp_df

                    # 첫 번째 row를 여러 번 복사
                    # samp_date_rows['채혈DT']
                    num_to_add=len(arranged_conc_samp_rows)-len(samp_date_rows)
                    first_row = samp_date_rows.iloc[-1]
                    rows_to_add = pd.DataFrame([first_row] * (num_to_add), columns=samp_date_rows.columns)

                    # 원래 DataFrame과 합치기
                    # samp_date_rows = samp_date_rows.iloc[:-1,:].copy()
                    modi_samp_date_rows = pd.concat([samp_date_rows, rows_to_add])

                    # arranged_conc_samp_rows[['채혈DT','CONC']]

                    # list(['채혈DT'])
                    # for dd_inx, dd_row in arranged_conc_samp_rows:

                    arranged_conc_samp_rows['채혈DT'] = list(modi_samp_date_rows['채혈DT'])
                    for dd_inx, dd_row in arranged_conc_samp_rows.iterrows():
                        if dd_row['채혈DT'] not in list(res_frag_df['SAMP_DT']):
                            print('기록')
                            res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                            res_frag_df.at[dd_inx, 'REC_REASON'] = '날린 곳'
                            # break
                        else:
                            # res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                            # res_frag_df['채혈DT']
                            # res_frag_df['SAMP_DT']
                            # res_frag_df['POT채혈DT']

                            print('이미 날짜 DT가 존재합니다')


    #################################        아래는 날렸음        ################################################
            elif len(concord_date_rows) < len(samp_date_rows):   ################# 채혈정보랑 안 맞는 부정확한 정보라 날린다
                print(f'CONC 데이터 길이: {len(concord_date_rows)} < SAMP 데이터 길이: {len(samp_date_rows)}')
                more_conc_than_samp.add(pid)
                continue
                # len(more_conc_than_samp)



                not_na_pot_dt_rows = concord_date_rows[~concord_date_rows['POT채혈DT'].isna()].copy()
                not_na_pot_dt_rows = not_na_pot_dt_rows[not_na_pot_dt_rows['POT채혈DT']!=''].copy()
                if len(not_na_pot_dt_rows)>0:
                    for dd_inx, dd_row in not_na_pot_dt_rows.iterrows(): # break

                        # dd_row['POT_SAMP_TIME']
                        if (dd_row['POT_SAMP_MONTHDAY']!='') and (dd_row['POT_SAMP_TIME']!=''):
                            samp_dt_input = dd_row['POT채혈DT']
                        elif (dd_row['POT_SAMP_MONTHDAY']=='') and (dd_row['POT_SAMP_TIME']!=''):
                            samp_dt_input = dd_row['보고일'] + dd_row['POT_SAMP_TIME']
                        elif (dd_row['POT_SAMP_MONTHDAY']!='') and (dd_row['POT_SAMP_TIME']==''):
                            samp_dt_input = dd_row['POT채혈DT'] + 'T00:00'
                            raise ValueError
                        else:
                            print('이런경우는 없을듯. 이미 앞전 조건이 pot dt가 있어야하는 조건')
                            raise ValueError



                        if samp_dt_input not in list(res_frag_df['SAMP_DT']):
                            # print('기록')
                            res_frag_df.at[dd_inx, 'SAMP_DT'] = samp_dt_input
                            res_frag_df.at[dd_inx, 'REC_REASON'] = '날린 곳'
                        else:

                            print('이미 날짜 DT가 존재합니다')
                else:
                    # concord_date_rows['CONC']
                    print(f'pot dt 도 없으면서 실제 / CONC 데이터 길이: {len(concord_date_rows)} < SAMP 데이터 길이: {len(samp_date_rows)}')
                    # res_frag_df
                    # samp_date_rows['채혈DT']
                    for conc_ord_row_inx in range(len(concord_date_rows)):#break
                    # if len(concord_date_rows)==1:
                    # conc_ord_row_inx=1
                        temp_concord_date_row = concord_date_rows.iloc[conc_ord_row_inx:conc_ord_row_inx+1, :].copy()
                        # concord_date_rows.iloc[conc_ord_row_inx:,:]
                        # concord_date_rows[conc_ord_date_type]
                        # concord_date_rows['CONC']
                        # mean_conc
                        # samp_date_rows
                        # id_dose_df
                        # id_dose_df['DOSE_DT']

                        temp_samp_date_rows = samp_date_rows.copy()
                        temp_samp_date_rows['임시P_T'] = ''

                        for samp_inx, samp_row in temp_samp_date_rows.iterrows(): #break
                            samp_dt_str = samp_row['채혈DT']
                            samp_dt = datetime.strptime(samp_dt_str, '%Y-%m-%dT%H:%M')

                            dose_before_samps = id_dose_df[id_dose_df['DOSE_DT'] <= samp_dt_str].copy()
                            dose_after_samps = id_dose_df[id_dose_df['DOSE_DT'] > samp_dt_str].copy()

                            before_samp_max_dt_str = dose_before_samps['DOSE_DT'].max()
                            after_samp_min_dt_str = dose_after_samps['DOSE_DT'].min()

                            if type(before_samp_max_dt_str) == str:
                                temp_dose_before_samps_dt = datetime.strptime(before_samp_max_dt_str, '%Y-%m-%dT%H:%M')
                            else:
                                temp_dose_before_samps_dt = datetime.strptime('0001-01-01T00:00', '%Y-%m-%dT%H:%M')

                            if type(after_samp_min_dt_str) == str:
                                temp_dose_after_samps_dt = datetime.strptime(after_samp_min_dt_str, '%Y-%m-%dT%H:%M')
                            else:
                                temp_dose_after_samps_dt = datetime.strptime('9999-12-31T23:59', '%Y-%m-%dT%H:%M')



                            # 직전까지의 DOSE_DT가 다음 DOSE_DT 보다 시간적으로 가까우면 -> 채혈DT중 min값은 peak conc
                            if ((samp_dt - temp_dose_before_samps_dt).total_seconds() / 60) < ((temp_dose_after_samps_dt - samp_dt).total_seconds() / 60):
                                samp_dt_conc_type = 'peak'
                            else:
                                samp_dt_conc_type = 'trough'

                            temp_samp_date_rows.at[samp_inx, '임시P_T'] = samp_dt_conc_type

                        # for concord_inx,
                        for dd_inx, dd_row in temp_concord_date_row.iterrows(): #break

                            peak_samp_date_rows = temp_samp_date_rows[temp_samp_date_rows['임시P_T']=='peak'].copy()
                            trough_samp_date_rows = temp_samp_date_rows[temp_samp_date_rows['임시P_T'] == 'trough'].copy()
                            used_dt_list = list()
                            if dd_row['CONC'] >= mean_conc:
                                conc_level_type = 'peak'

                                if len(peak_samp_date_rows)==0:
                                    input_samp_date_rows = trough_samp_date_rows.copy()
                                else:
                                    input_samp_date_rows = peak_samp_date_rows.copy()
                            else:
                                conc_level_type = 'trough'

                                if len(trough_samp_date_rows)==0:
                                    input_samp_date_rows = peak_samp_date_rows.copy()
                                else:
                                    input_samp_date_rows = trough_samp_date_rows.copy()

                            check_input_samp_date_rows = input_samp_date_rows[~input_samp_date_rows['채혈DT'].isin(used_dt_list)].copy()
                            if len(check_input_samp_date_rows) == 0:
                                input_samp_date_rows = input_samp_date_rows.iloc[-1:, :].copy()
                            input_samp_dt = input_samp_date_rows['채혈DT'].iloc[-1]
                            # concord_date_rows['CONC']
                            # res_frag_df['CONC']
                            if input_samp_dt not in list(res_frag_df['SAMP_DT']):
                                print('기록')
                                res_frag_df.at[dd_inx, 'SAMP_DT'] = input_samp_dt
                                res_frag_df.at[dd_inx, 'REC_REASON'] = '날린 곳'
                                used_dt_list.append(input_samp_dt)
                            else:
                                input_samp_dt = input_samp_date_rows['채혈DT'].iloc[-2]
                                print('이미 날짜 DT가 존재합니다 / 다음 농도로 대체하여 설정합니다')
                                raise ValueError


                        # print('여기 할 차례')
                        # raise ValueError
                    # else:
                    #     print('CONC 데이터 길이 2개 이상')
                    #     raise ValueError

            else:
                print('len(concord_date_rows), len(samp_date_rows), len(samp_date_rows) 그외 케이스')
                raise ValueError


    else:


    # print(f"({finx}) {pname} / {pid} / 그 어디도 X")
    # COND1 = (res_frag_df['CONC'].diff().map(np.abs) / res_frag_df['CONC'] >= 3) # 3배 이상의 급격한 CONC 변화
    # COND2 = (res_frag_df['CONC'].diff().map(np.abs) / res_frag_df['CONC'] < 3) # 3배 미만의 완만한 CONC 변화
    # COND3 = (np.abs(res_frag_df['SAMP_DT'].map(lambda x:datetime.strptime(x,'%Y-%m-%dT%H:%M').timestamp()).diff())/60 <= 15) # 시간 차이가 15분 이내
    # # COND4 = (np.abs(res_frag_df['SAMP_DT'].map(lambda x:datetime.strptime(x,'%Y-%m-%dT%H:%M').timestamp()).diff())/60 <= 3) # 시간 차이가 3분 이내
    #
    # if (COND2 & COND3).sum() > 0:  # 시간이 거의 일치하는데 상대적으로 완만한 변화 -> Toxic 후 감소하는 것으로 간주
    #     pass
    # elif (COND1 & COND3).sum() > 0:  # 시간이 거의 일치하는데 급격한 변화 -> 가장 가까운 Dosing time 기준 trough, peak로 간주
    #     pass
    # else:
    #     pass



    # if len(res_frag_df)>len(res_frag_df.drop_duplicates(['ID','SAMP_DT'])):
    #     res_frag_df[~res_frag_df['SAMP_DT'].isin(res_frag_df.drop_duplicates(['ID', 'SAMP_DT'],keep=False)['SAMP_DT'].unique())]

        """
        # 다시 확인해봐야할 사람들 
        10197507	이동주    2009-08-01T12:44 ~


        """
        # raise ValueError

    final_df.append(res_frag_df)
final_df = pd.concat(final_df, ignore_index=True)
final_df.to_csv(f"{output_dir}/final_conc_df(with sampling).csv", encoding='utf-8-sig', index=False)

print(f"Total patients in final_df: {len(final_df.drop_duplicates(['ID']))} (from {len(pt_info)}) / {len(final_df)} rows")
print(f"No dose data: {len(no_dose_pids)}")
print(f"Ambiguous dosing time data: {len(ambiguous_dosing_time_pids)}")
print(f"Partly no dose data: {len(partly_no_dose_pids)}")
print(f"Re-dose data: {len(re_dose_pids)}")
print(f"Error_passing data: {len(error_passing_pids)}")
# final_df.to_csv(f"{output_dir}/final_conc_df(with sampling).csv", encoding='utf-8-sig', index=False)

#
# final_df = pd.concat(final_df, ignore_index=True)
# final_df.to_csv(f"{output_dir}/final_conc_df(with sampling).csv", encoding='utf-8-sig', index=False)


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
