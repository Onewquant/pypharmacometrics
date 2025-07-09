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
partly_no_dose_pids = set()
re_dose_pids = set()
error_passing_pids = set()
for finx, fpath in enumerate(pt_files): #break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]
    # pt_df.append({'ID':pid,'NAME':pname})
    id_info_df = pt_info_dup[pt_info_dup['ID'] == pid].copy()
    min_date = (datetime.strptime(id_info_df['TDM_REQ_DATE'].iloc[0],'%Y-%m-%d')-timedelta(days=62)).strftime('%Y-%m-%d')
    max_date = (datetime.strptime(id_info_df['TDM_RES_DATE'].iloc[0],'%Y-%m-%d')+timedelta(days=92)).strftime('%Y-%m-%d')

    id_dose_df = dose_result_df[dose_result_df['ID'] == pid].copy()
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



    # res_frag_df = id_conc_df.sort_values(['오더일'], ascending=False)
    res_frag_df = id_conc_df.copy()
    # res_frag_df = id_conc_df.sort_values(['오더일'], ascending=True)
    # pd.to_datetime(res_frag_df['오더일']).diff().dt.days
    for c in ['SAMP_DT','채혈DT', '라벨DT', '접수DT', '시행DT', '보고DT']:
        res_frag_df[c] = ''

    # if finx
    # raise ValueError
    # res_frag_df.columns

    # POT채혈DT 있는 경우
    pot_dt_rows = res_frag_df[(res_frag_df['POT채혈DT'] != '') & (~res_frag_df['POT채혈DT'].isna())].copy()
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

            mean_conc = res_frag_df['CONC'].mean()
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
                        break
                    else:
                        continue
                    # res_frag_df['SAMP_DT']
                    # res_frag_df.at[oer_inx, 'SAMP_DT'] = oer_row['오더일'] + 'T' + est_conc_dt_tups[1].split('T')[-1]
            else:
                print(f"({finx}) {pname} / {pid} / Only Conc data (Dosing: O, Sampling X, POT채혈DT: X / 해당 날짜 dose기록 3개 이상")
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
                ascending_orders = [True,True]

            # 보고일과 같은날 Dosing된 날짜 존재시 -> Conc 측정 날짜를 보고일에 되었다고 추정
            elif len(temp_id_dose_df[(temp_id_dose_df['DOSE_DATE'] == conc_report_date_list[0])])>0:
                temp_id_dose_df = temp_id_dose_df[(temp_id_dose_df['DOSE_DATE'] == conc_report_date_list[0])].copy()
                ascending_orders = [True, True]

            # 일치되는 날짜가 없을 경우 -> DOSE 날짜 범위 중 가장 Latest 에 Conc 측정되었다고 추정
            else:
                ascending_orders = [False, True]

            mean_conc = res_frag_df['CONC'].mean()
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
            mean_conc = res_frag_df['CONC'].mean()
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
                                    break
                                else:
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

        print(f"({finx}) {pname} / {pid} / SAMP 데이터 존재")
        mean_conc = res_frag_df['CONC'].mean()
        # tdm_relate_dates = set(id_info_df['TDM_RES_DATE']).union(set(id_info_df['TDM_REQ_DATE']))
        samp_relate_dates = set(sdf['보고일'])
        if len(samp_relate_dates.intersection(set(cdf['오더일']))) > len(samp_relate_dates.intersection(set(cdf['보고일']))):
            conc_ord_date_type = '오더일'
            print(f"선택된 날짜 타입이 '{conc_ord_date_type}'이 더 잘 맞음")
            raise ValueError
        elif len(samp_relate_dates.intersection(set(cdf['오더일']))) == len(samp_relate_dates.intersection(set(cdf['보고일']))):
            conc_ord_date_type = '오더일'
            other_type = '보고일'
        else:
            conc_ord_date_type = '보고일'
            other_type = '오더일'
        for concord_date in cdf[conc_ord_date_type].drop_duplicates(): # break
            concord_date_rows = cdf[cdf[conc_ord_date_type] == concord_date].copy()  # break
            samp_date_rows = sdf[sdf[conc_ord_date_type] == concord_date].copy()  # break

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
                        if len(samp_date_rows)==0:
                            print("SAMPLING TIME DATA가 부분적으로 없는 경우 / '오더일' '보고일' 바꿔서 해볼 수 있는 경우 / 바꿔도 날짜가 안나옴")
                            raise ValueError
                    elif len(other_type_dates)>1:
                        print("SAMPLING TIME DATA가 부분적으로 없는 경우 / '오더일' '보고일' 바꿔서 해볼 수 있는 경우 / 날짜가 여러 개 나와서 정하기 어려움")
                        raise ValueError
                    else:
                        print("SAMPLING TIME DATA가 부분적으로 없는 경우 / ori sample file에도 없음")

                        raise ValueError
                # id_info_df
                # id_samp_df['채혈DT']
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
                            # break
                        else:
                            # res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                            # res_frag_df['채혈DT']
                            # res_frag_df['SAMP_DT']
                            # res_frag_df['POT채혈DT']

                            print('이미 날짜 DT가 존재합니다')
                else:
                    # samp_date_rows['채혈DT']
                    dose_before_samps = id_dose_df[(id_dose_df['DOSE_DT'] < samp_date_rows['채혈DT'].min())].copy()
                    dose_between_samps = dose_ord_frag[(dose_ord_frag['DOSE_DT'] >= samp_date_rows['채혈DT'].min())&(dose_ord_frag['DOSE_DT'] <= samp_date_rows['채혈DT'].max())].copy()
                    dose_after_samps = id_dose_df[(id_dose_df['DOSE_DT'] > samp_date_rows['채혈DT'].min())].copy()
                    # 채혈 시간 기록 사이에 Dose 투약 기록이 있는 경우 -> 채혈 농도값으로 peak / trough에 해당하는 농도들 구분하고, '채혈DT'를 이에 맞춰 배정 후 SAMP_DT로 확정시킴 (DOSE_DT 앞 뒤에 있는게 맞는지 확인도하면 좋을듯)
                    if len(dose_between_samps)>=2:
                        print('CONC 데이터 길이 == SAMP 데이터 길이 / sampling 사이에 dose 기록이 2개 이상 존재')
                        # dose_between_samps['DOSE_DT']
                        temp_min_samp_dt_str = samp_date_rows['채혈DT'].min()
                        temp_max_samp_dt_str = samp_date_rows['채혈DT'].max()

                        temp_min_samp_dt = datetime.strptime(samp_date_rows['채혈DT'].min(),'%Y-%m-%dT%H:%M')
                        temp_max_samp_dt = datetime.strptime(samp_date_rows['채혈DT'].max(),'%Y-%m-%dT%H:%M')

                        before_samp_min_dt_str = dose_before_samps['DOSE_DT'].max()
                        if type(before_samp_min_dt_str)==str:
                            temp_dose_before_samps_dt = datetime.strptime(before_samp_min_dt_str,'%Y-%m-%dT%H:%M')
                        else:
                            temp_dose_before_samps_dt = '0000-00-00T00:00'

                        temp_dose_between_samps_min_dt = datetime.strptime(dose_between_samps['DOSE_DT'].min(), '%Y-%m-%dT%H:%M')
                        temp_dose_between_samps_max_dt = datetime.strptime(dose_between_samps['DOSE_DT'].max(), '%Y-%m-%dT%H:%M')

                        after_samp_min_dt_str = dose_after_samps['DOSE_DT'].min()
                        if type(after_samp_min_dt_str)==str:
                            temp_dose_after_samps_dt = datetime.strptime(after_samp_min_dt_str,'%Y-%m-%dT%H:%M')
                        else:
                            temp_dose_after_samps_dt = '9999-99-99T99:99'

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

                    elif len(dose_between_samps)==1:

                        conc_ascending_sort = True
                        samp_ascending_sort = True

                        # print(f'CONC 데이터 길이 == SAMP 데이터 길이 / sampling 사이에 dose 기록이 한 개 존재')
                        # raise ValueError

                    # 채혈 시간 기록 사이에 Dose 투약 기록이 없는 경우 -> 채혈 농도값의 크기를 내림차순(큰->작은)으로 배열하고 '채혈DT'를 오름차순으로 맞춰 배정
                    else:
                        conc_ascending_sort = False
                        samp_ascending_sort = True

                    arranged_conc_samp_rows = concord_date_rows.sort_values(['CONC'], ascending=conc_ascending_sort)
                    arranged_conc_samp_rows['채혈DT'] = list(samp_date_rows.sort_values(['채혈DT'], ascending=samp_ascending_sort)['채혈DT'])
                    for dd_inx, dd_row in arranged_conc_samp_rows.iterrows():
                        if dd_row['채혈DT'] not in list(res_frag_df['SAMP_DT']):
                            print('기록')
                            res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                            # break
                        else:
                            res_frag_df.at[dd_inx, 'SAMP_DT'] = dd_row['채혈DT']
                            # res_frag_df['채혈DT']
                            # res_frag_df['SAMP_DT']
                            # res_frag_df['POT채혈DT']

                            print('이미 날짜 DT가 존재합니다')
                            # raise ValueError


            elif (len(concord_date_rows)==len(samp_date_rows)) and (len(samp_date_rows) > 2):
                print('샘플이 3 개 이상')
                dose_between_samps = dose_ord_frag[(dose_ord_frag['DOSE_DT'] >= samp_date_rows['채혈DT'].min()) & (dose_ord_frag['DOSE_DT'] <= samp_date_rows['채혈DT'].max())].copy().sort_values(['DOSE_DT'])
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

                    prev_ds_frag_row = dose_samp_date_rows[(dose_samp_date_rows['임시DT'] > prev_dosedt)&(dose_samp_date_rows['임시DT'] <= dosedt)&(~dose_samp_date_rows['ID'].isna())].copy()
                    next_ds_frag_row = dose_samp_date_rows[(dose_samp_date_rows['임시DT'] <= next_dosedt) & (dose_samp_date_rows['임시DT'] > dosedt) & (~dose_samp_date_rows['ID'].isna())].copy()

                    if len(prev_ds_frag_row)==0:
                        if len(next_ds_frag_row)==0:
                            pass
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
                        dose_samp_date_rows.at[temp_inx, '임시CONC'] = trough_concord_date_rows.iloc[-1]['CONC']
                        trough_concord_date_rows = trough_concord_date_rows.iloc[:-1,:]
                    else:
                        print('DOSE와 DOSE사이 TROUGH 샘플이 2 개 이상 - peak/trough 배열 조정 요망')
                        raise ValueError

                    if dosedt_inx==len(dose_between_samps)-1:
                        if len(next_ds_frag_row)==0:
                            print('마지막 dose 처리 / next_ds_frag_row가 0 개')
                            pass
                        elif len(next_ds_frag_row)==1:
                            temp_inx = next_ds_frag_row.index[0]
                            dose_samp_date_rows.at[temp_inx, '임시CONC'] = peak_concord_date_rows.iloc[0]['CONC']
                            peak_concord_date_rows = peak_concord_date_rows.iloc[:-1,:]
                        else:
                            print('마지막 dose 처리 / next_ds_frag_row가 2 개 이상')
                            raise ValueError
                        decided_dose_samp_date_rows = dose_samp_date_rows[~dose_samp_date_rows['ID'].isna()].copy()
                        decided_dose_samp_date_rows['CONC'] = decided_dose_samp_date_rows['임시CONC']
                        decided_dose_samp_date_rows['SAMP_DT'] = decided_dose_samp_date_rows['임시DT']
                        # concord_date_rows.columns
                        if len(res_frag_df)==len(decided_dose_samp_date_rows):
                            for col in ['CONC','SAMP_DT',]:
                                res_frag_df[col] = list(decided_dose_samp_date_rows[col])
                        else:
                            print(f"len(res_frag_df)!=len(decided_dose_samp_date_rows)")
                            raise ValueError

                        # decided_dose_samp_date_rows['CONC']
                        # res_frag_df.columns
                        # res_frag_df['CONC']
                        # for
                        #     res_frag_df.at[temp_inx, '임시CONC'] = peak_concord_date_rows.iloc[0]['CONC']
                        # decided_dose_samp_date_rows['CONC']
                    prev_dosedt = dosedt

                # mean_conc


                # raise ValueError

            elif len(concord_date_rows) > len(samp_date_rows):
                print(f'CONC 데이터 길이: {len(concord_date_rows)} > SAMP 데이터 길이: {len(samp_date_rows)}')
                dose_between_samps = dose_ord_frag[(dose_ord_frag['DOSE_DT'] >= samp_date_rows['채혈DT'].min()) & (dose_ord_frag['DOSE_DT'] <= samp_date_rows['채혈DT'].max())].copy()
                # 채혈 시간 기록 사이에 Dose 투약 기록이 있는 경우 -> 채혈 농도값으로 peak / trough에 해당하는 농도들 구분하고, '채혈DT'를 이에 맞춰 배정 후 SAMP_DT로 확정시킴 (DOSE_DT 앞 뒤에 있는게 맞는지 확인도하면 좋을듯)
                if len(dose_between_samps)>0:


                    if len(dose_between_samps) >=2:
                        print(f'CONC 데이터 길이 ({len(concord_date_rows)}) > SAMP 데이터 길이 ({len(samp_date_rows)}) / sampling 사이에 dose 기록이 여러개 존재')
                        raise ValueError

                    dose_between_samps
                    print(f'CONC 데이터 길이 ({len(concord_date_rows)}) > SAMP 데이터 길이 ({len(samp_date_rows)}) / sampling 사이에 dose 기록이 한 개 존재')
                    raise ValueError

                    samp_date_rows['채혈DT']
                    concord_date_rows['CONC']

                # 채혈 시간 기록 사이에 Dose 투약 기록이 없는 경우 -> 채혈 농도값의 크기를 내림차순(큰->작은)으로 배열하고 '채혈DT'를 내림차순으로 맞춰 배정
                else:
                    arranged_conc_samp_rows = concord_date_rows.sort_values(['CONC'], ascending=True)
                    samp_date_rows = samp_date_rows.sort_values(['채혈DT'], ascending=True)

                    # id_conc_df
                    # id_samp_df

                    # 첫 번째 row를 여러 번 복사
                    num_to_add=len(arranged_conc_samp_rows)
                    first_row = samp_date_rows.iloc[-1]
                    rows_to_add = pd.DataFrame([first_row] * (num_to_add-1), columns=samp_date_rows.columns)

                    # 원래 DataFrame과 합치기
                    # samp_date_rows = samp_date_rows.iloc[:-1,:].copy()
                    modi_samp_date_rows = pd.concat([samp_date_rows, rows_to_add])


                    # list(['채혈DT'])
                    # for dd_inx, dd_row in arranged_conc_samp_rows:
                    arranged_conc_samp_rows['채혈DT'] = list(modi_samp_date_rows['채혈DT'])
                # id_conc_df
                # res_frag_df['SAMP_DT']
                # raise ValueError
            elif len(concord_date_rows) < len(samp_date_rows):
                print(f'CONC 데이터 길이: {len(concord_date_rows)} < SAMP 데이터 길이: {len(samp_date_rows)}')
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
                        else:

                            print('이미 날짜 DT가 존재합니다')
                else:
                    print(f'pot dt 도 없으면서 실제 / CONC 데이터 길이: {len(concord_date_rows)} < SAMP 데이터 길이: {len(samp_date_rows)}')
                    if len(concord_date_rows)==1:
                        concord_date_rows[conc_ord_date_type]
                        samp_date_rows
                        id_dose_df
                        id_dose_df['DOSE_DT']
                        print('여기 할 차례')
                        raise ValueError
                    else:
                        print('CONC 데이터 길이 2개 이상')
                        raise ValueError

    ## 기타 추가 처리

    # (1) 같은 DT의 농도채혈시간이 존재하는경우 앞뒤 Dosing time 고려하여 외삽


            # id_info_df
            # id_samp_df


# """
#         # id_info_df
#         # res_frag_df[['보고일','오더일', 'CONC']]
#         # res_frag_df[res_frag_df['SAMP_DT']=='']
#         cdf = res_frag_df[res_frag_df['SAMP_DT']==''].copy()
#         pdf = cdf[cdf['POT채혈DT'] != ''].copy()
#         sdf = id_samp_df[['보고일', '채혈DT', '라벨DT', '접수DT','시행DT','보고DT']].copy()
#
#
#         mdf = cdf.merge(sdf, on=['보고일'], how='outer')
#         mdf = mdf.sort_values(['보고일', 'CONC', '채혈DT'])
#
#         trough_mdf = mdf.drop_duplicates(['보고일'], keep='first')
#         peak_mdf = mdf.drop_duplicates(['보고일'], keep='last')
#         total_mdf = pd.concat([trough_mdf, peak_mdf]).drop_duplicates(['보고일', 'CONC', '채혈DT'], keep='last').sort_values(['보고일', 'CONC', '채혈DT'])
#         # total_mdf.columns
#         # cdf = id_conc_df[['보고일', 'POT채혈DT', 'CONC']].sort_values(['보고일', 'CONC'])
#         # sdf = id_samp_df[['보고일', '채혈DT', '라벨DT', '접수DT']].copy()
#         #
#         # mdf = cdf.merge(sdf, on=['보고일'], how='outer')[['보고일', 'CONC', '채혈DT', 'POT채혈DT']]
#         # mdf = mdf.sort_values(['보고일', 'CONC', '채혈DT'])
#         #
#         # trough_mdf = mdf.drop_duplicates(['보고일'], keep='first')
#         # peak_mdf = mdf.drop_duplicates(['보고일'], keep='last')
#         # total_mdf = pd.concat([trough_mdf, peak_mdf]).drop_duplicates(['보고일','CONC', '채혈DT'], keep='last').sort_values(['보고일', 'CONC', '채혈DT'])
#         #
#
#
#         # CONC 데이터는 오더일, 보고일 날짜 두개만 존재
#         # SAMPLING 데이터는 여러개의 날짜 존재 두개만 존재
#
#         if len(total_mdf)!=len(cdf):
#             print(f"({finx}) {pname} / {pid} / No matched")
#             conc_samp_mismatch_pids.append(pid)
#             # len(conc_samp_mismatch_pids)
#             # continue
#         else:
#             print(f"({finx}) {pname} / {pid} / matched")
#             # if pid not in ['11116501', '10112328', '10143478', '10228470', '10533576', '10885385', '10914959', ]:
#             #     raise ValueError
#
#             # '10112328' : 04.17 의 샘플링 타임데이터는 있는데 농도데이터는 부재함.
#             # '10143478' : 2018-06-04 의 농도데이터가 3개 있는데, 샘플링 데이터는 2개라 max, min만 남기면 df 길이 달라짐
#             # '10228470' : 2005-04-22 샘플링 데이터는 존재, 농도 데이터는 그날 것 없음. (total_mdf 에 CONC가 NAN인 값 생김)
#             # '10533576' : 2004-09-14에 CONC가 0.5로 똑같은 데이터 2개 존재. 아마도 9-13일 채혈일듯. (보고일 기준으로만 하고 있는데, 중복된 데이터의 오더일은 다른 것으로 보아 하나는 2004-09-13 데이터인듯
#             # '10885385' : 샘플링 데이터 보고일기준 2003-12-30 는 1개, 2004-01-02 는 2개 인데, 농도데이터에서는 2, 1개로 되어 있음
#             # '10914959'
#             # '11116501'
#
#
#         final_df.append(total_mdf)
#
#         # if finx==3:
#         #     raise ValueError
#
#         # cdf = id_conc_df[['오더일', '보고일', 'POT채혈DT', 'CONC']].copy()
#         # sdf = id_samp_df[['오더일', '보고일', '채혈DT', '라벨DT', '접수DT']].copy()
#
#         # mdf = cdf.merge(sdf, on=['오더일', '보고일'], how='outer')[['오더일', '보고일','CONC','채혈DT','POT채혈DT']]
#         # mdf = mdf.sort_values(['오더일', '보고일', 'CONC', '채혈DT'])
#
#         # trough_mdf = mdf.drop_duplicates(['오더일', '보고일'], keep='first')
#         # peak_mdf = mdf.drop_duplicates(['오더일', '보고일'], keep='last')
#         # total_mdf = pd.concat([trough_mdf, peak_mdf]).drop_duplicates(['오더일', '보고일', 'CONC', '채혈DT'], keep='last')
#
# """



    else:
        print(f"({finx}) {pname} / {pid} / 그 어디도 X")
        # raise ValueError


    # elif (pid not in uniq_conc_pids) or (pid not in uniq_sampling_pids):
    #     print(f"({finx}) {pname} / {pid} / No eighther data")
    #     continue
    # else:
    #     raise ValueError

    final_df.append(res_frag_df)
final_df = pd.concat(final_df, ignore_index=True)
final_df.to_csv(f"{output_dir}/final_conc_df(with sampling_test).csv", encoding='utf-8-sig', index=False)
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
