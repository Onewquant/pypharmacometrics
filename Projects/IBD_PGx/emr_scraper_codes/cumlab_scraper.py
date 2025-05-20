from emr_scraper.tools import *

df = pd.read_csv("C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/anti_TNFa_patients.csv")

## Lab 누적조회 수집 (농도 및 PD marker (Calprotectin) 측정시간)

cumdata_list = ['Infliximab', 'Adalimumab', 'Calprotectin']

# try:
#     existing_df = pd.read_excel("C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/IBD_PGx_cumlab.xlsx")
#     existing_id_list = list(existing_df['EMR ID'])
# except:
#     existing_df = pd.DataFrame(columns=['EMR ID', 'name', 'IBD type', 'raw_CumLab'])
#     existing_id_list = list()

existing_files = glob.glob("C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/cumlab/IBD_PGx_cumlab(*).xlsx")
existing_id_list = [int(f.split('IBD_PGx_cumlab(')[-1].split('_')[0]) for f in existing_files]
# move_to_content_window(content='Hx')
# len(existing_id_list)

lab_df = list()

for inx, row in df.iterrows(): #break

    new_sleep_time = 20
    row['name'] = row['name'].split('(')[0]
    if row['name'] in ['이건수', '김규진', '김은숙', '김기백', '최수호']:
        new_sleep_time = 360
        # continue

    prev2_sp_lab, prev_sp_lab = '',''
#     print(1)
#     if inx==2:
#         break
    pid = row['EMR ID']
    if pid in existing_id_list:
        print(f"({inx}) {pid} / {row['name']} / {row['IBD type']} / 이미 자료가 존재합니다.")
        continue

    # pyautogui.position()

    # 환자ID 입력

    move_to_click_and_wait(x=20, y=245, sleep_time=3)   # 환자번호 입력
    move_to_click_and_wait(x=239, y=921, sleep_time=1)   # 검색버튼 누르기
    keyboard_input_prep(sleep_time=0)
    keyboard_input(input_str=pid)
    move_to_click_and_wait(x=90, y=968, sleep_time=1)   # '전체' 누르기
    move_to_click_and_wait(x=370, y=882, sleep_time=1)   # 조회 누르기

    move_to_click_and_wait(x=586, y=157, sleep_time=1, double_click=True)   # 환자선택 누르기

    move_to_click_and_wait(x=684, y=653, sleep_time=1)  # 환자선택사유 입력
    move_to_click_and_wait(x=1162, y=926, sleep_time=10)  # 환자선택사유 확인


    # 검사결과조회

    move_to_click_and_wait(x=1011, y=362, sleep_time=1)  # 기간을 전체로 변경
    move_to_click_and_wait(x=1041, y=361, sleep_time=3)  # 재조회

    # 필터로 적용 데이터 고르기

    delta = 26
    # delnum = 0
    cum_lab_df = list()
    inc_lab_df = pd.DataFrame(columns=['보고일', '오더일', '검사명', '검사결과'])
    # for delnum in range(22):

    for delnum in range(22):
    # for delnum in range(12,17):
        raw_sp_lab = ''

        y_pos = 455 + delnum * delta
        move_to_click_and_wait(x=857, y=y_pos, sleep_time=1)  # 조회 누르기
        time.sleep(3)

        # for i in range(30):
        #     move_to_click_and_wait(x=1474, y=979, sleep_time=2)     # 다음결과 N번 클릭


        while_count = 0
        while (raw_sp_lab == '') or ('검사명\t참고치\t' in raw_sp_lab.split('\r\n')[0]):

            move_to_click_and_wait(x=1230, y=451, sleep_time=1)  # 첫보고일란 선택
            pyautogui.dragTo(x=1231, y=474, duration=0.3, button='left')  # 드래그 시작위치
            pyautogui.keyDown('shift')
            time.sleep(1)
            move_to_click_and_wait(x=1296, y=451, sleep_time=1)  # 첫오더일 선택
            move_to_click_and_wait(x=1419, y=451, sleep_time=1)  # 첫검사명 선택
            move_to_click_and_wait(x=1562, y=451, sleep_time=1)  # 검사결과 선택
            move_to_click_and_wait(x=1664, y=451, sleep_time=1)  # 직전결과 선택
            move_to_click_and_wait(x=1762, y=451, sleep_time=1)  # 참고치 선택
            move_to_click_and_wait(x=1868, y=451, sleep_time=1)  # 결과비고 선택
            pyautogui.dragTo(x=1868, y=472, duration=1, button='left')  # 드래그 시작위치
            pyautogui.dragTo(x=1868, y=944, duration=1, button='left')  # 드래그 끝위치
            # move_to_click_and_wait(x=1869, y=944, sleep_time=1)  # 마지막보고일란 선택

            time.sleep(1)
            pyautogui.keyUp('shift')
            time.sleep(1)

            pyautogui.hotkey('ctrl', 'c')
            raw_sp_lab = pyperclip.paste()
            time.sleep(1)
            while_count+=1
            if while_count==3:break
            # print('copy and pasted')

        if while_count==3:
            move_to_click_and_wait(x=857, y=y_pos, sleep_time=2)  # 조회 누르기
            continue

        if (prev_sp_lab == prev2_sp_lab) and (prev2_sp_lab!='') and (prev_sp_lab!=''):
            print('다음사람으로 넘어갑니다')
            break

        # 랩항목 샘플에 원하는 키워드 포함되어 있는지 확인
        lab_tbl_to_df = pd.DataFrame([c.replace('\t\t','\t').split('\t') for c in raw_sp_lab.strip().split('\r\n')])
        lab_col_len = len(lab_tbl_to_df.columns)
        if lab_col_len==6:
            lab_tbl_to_df.columns = ['보고일', '오더일', '검사명', '검사결과', '직전결과', '참고치']
            raw_splab_df = lab_tbl_to_df.copy()
            raw_splab_df['결과비고'] = ''
        elif lab_col_len==7:
            lab_tbl_to_df.columns = ['보고일', '오더일', '검사명', '검사결과', '직전결과', '참고치','결과비고']
            raw_splab_df = lab_tbl_to_df.copy()
        else:
            move_to_click_and_wait(x=857, y=y_pos, sleep_time=2)  # 조회 누르기
            continue
        raw_splablist_ds = raw_splab_df['검사명'].drop_duplicates()
        raw_splablist_ds = raw_splablist_ds[~raw_splablist_ds.isna()]
        raw_splabstr = str(list(raw_splablist_ds)).upper()
        check_yn = False
        for clab in cumdata_list:
            if clab.upper() in raw_splabstr:
                check_yn = True

        # 필터 검사항목 모두 확인 완료시 다음사람으로 넘어감
        splab_check_df = raw_splab_df.copy()
        # splab_check_df = raw_splab_df[raw_splab_df['검사명'].map(lambda x:True if len(set(cumdata_list).intersection(set(x.split(' ')))) > 0 else False)].copy()
        splab_check_df['Existance_Check'] = True
        exist_check_df = inc_lab_df.merge(splab_check_df[['보고일', '오더일', '검사명', '검사결과','Existance_Check']], on=['보고일', '오더일', '검사명', '검사결과'], how='left')
        print(f"Existance_check: {((exist_check_df['Existance_Check']==True)*1).sum()}")
        if ((exist_check_df['Existance_Check']==True)*1).sum() > 0:
            print('다음사람으로 넘어갑니다 (필터 검사항목들 모두 확인 완료)')
            break
            # raise ValueError

        cumlab_frag_df = pd.DataFrame(columns=['DT','Lab','Value'])

        if check_yn:

            cum_delta = 26


            move_to_click_and_wait(x=1230, y=451, sleep_time=1, button='right')  # 첫보고일란 우클릭
            move_to_click_and_wait(x=1267, y=466, sleep_time=1)                  # 누적조회 누르기

            time.sleep(5)

            move_to_click_and_wait(x=855, y=467, sleep_time=new_sleep_time)  # 누적조회 시작항목 클릭
            pyautogui.keyDown('shift')
            move_to_click_and_wait(x=855, y=467 + (len(raw_splablist_ds)-1) * cum_delta , sleep_time=new_sleep_time)  # 누적조회 끝항목 클릭
            pyautogui.keyUp('shift')

            for i in range(3):
                pyautogui.hotkey('ctrl', 'c')
                raw_cum_lab = pyperclip.paste()
                time.sleep(1)

            # raw_cum_lab = input('직접 써 넣으세요.')

            raw_cum_lab_text = raw_cum_lab.strip().replace('\r\n*','*').replace('\r\nPositive','Positive').replace('\r\nEquivocal','Equivocal').replace('재검한 결과입니다.\r\n','재검한 결과입니다.')
            cumlab_frag_df = pd.DataFrame([c.split('\t') for c in raw_cum_lab_text.split('\r\n')]).copy()

            # cumlab_frag_df
            cumlab_frag_df = cumlab_frag_df.T.copy()
            # cumlab_frag_df.iloc[:,2]
            # cumlab_frag_df[['Value']].head(50)
            cumlab_frag_cols = ['DT']+[c.split('\t')[0] for c in cumlab_frag_df.iloc[0,1:]]
            cumlab_frag_df = cumlab_frag_df.iloc[2:,:].copy()
            cumlab_frag_df.columns = cumlab_frag_cols
            if '' in cumlab_frag_cols:
                cumlab_frag_df = cumlab_frag_df.drop('', axis=1)

            # if '차세대 의료정보시스템이(가) 응답하지 않습니다.' in cumlab_frag_df.columns:
            #     time.sleep(720)
            #     move_to_click_and_wait(x=855, y=467 + (len(raw_splablist_ds) - 1) * cum_delta,
            #                            sleep_time=new_sleep_time)  # 누적조회 끝항목 클릭
            #     pyautogui.keyUp('shift')
            #
            #     for i in range(3):
            #         pyautogui.hotkey('ctrl', 'c')
            #         raw_cum_lab = pyperclip.paste()
            #         time.sleep(1)
            #
            #     # raw_cum_lab = input('직접 써 넣으세요.')
            #
            #     raw_cum_lab_text = raw_cum_lab.strip().replace('\r\n*', '*').replace('\r\nPositive',
            #                                                                          'Positive').replace(
            #         '\r\nEquivocal', 'Equivocal').replace('재검한 결과입니다.\r\n', '재검한 결과입니다.')
            #     cumlab_frag_df = pd.DataFrame([c.split('\t') for c in raw_cum_lab_text.split('\r\n')]).copy()
            #
            #     # cumlab_frag_df
            #     cumlab_frag_df = cumlab_frag_df.T.copy()
            #     # cumlab_frag_df.iloc[:,2]
            #     # cumlab_frag_df[['Value']].head(50)
            #     cumlab_frag_cols = ['DT'] + [c.split('\t')[0] for c in cumlab_frag_df.iloc[0, 1:]]
            #     cumlab_frag_df = cumlab_frag_df.iloc[2:, :].copy()
            #     cumlab_frag_df.columns = cumlab_frag_cols
            #     if '' in cumlab_frag_cols:
            #         cumlab_frag_df = cumlab_frag_df.drop('', axis=1)

            cumlab_frag_df = pd.melt(cumlab_frag_df.dropna(axis=1, how='all').reset_index(drop=True), id_vars=['DT'], var_name='Lab', value_name='Value')

            # print(cumlab_frag_df.columns)
            # print(cumlab_frag_df.iloc[0])

            cum_lab_df.append(cumlab_frag_df)
            prev2_sp_lab = prev_sp_lab
            prev_sp_lab = raw_sp_lab
            move_to_click_and_wait(x=1890, y=401, sleep_time=3)  # 누적조회 끄기

            inc_lab_df = pd.concat([inc_lab_df, raw_splab_df])

        # print(f"{delnum} / {check_yn} / raw_sp_lab: {raw_sp_lab}")
        print(f"{delnum} / {check_yn} / 추가df: {len(cumlab_frag_df)} rows / 누적df: {len(cum_lab_df)}")

        move_to_click_and_wait(x=857, y=y_pos, sleep_time=2)  # 조회 누르기

        if len(cum_lab_df)>1:
            cdl_count=0
            cdl_comp_str = str(pd.concat(cum_lab_df, ignore_index=True)['Lab'].unique())
            for cdl_str in cumdata_list:
                if cdl_str in cdl_comp_str:
                    cdl_count+=1
            if (cdl_count>1) and ('Calprotectin' in cdl_comp_str):
                print('다음사람으로 넘어갑니다 (모든 필요 랩 누적수치 수집완료)')
                break

    cum_lab_count = len(cum_lab_df)
    cum_lab_df = pd.concat(cum_lab_df, ignore_index=True)
    cum_lab_df.to_excel(f"C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/cumlab/IBD_PGx_cumlab({pid}_{row['name']}).xlsx",index=False)
    print(f"({inx}) {pid} / {row['name']} / {row['IBD type']} / {cum_lab_count} dfs / {len(cum_lab_df)} rows")

    # raise ValueError

        # pyautogui.position()

            # move_to_click_and_wait(x=1230, y=451, sleep_time=2, double_click=True, button='right')         # Lab 데이터 첫보고일란 선택

    # move_to_click_and_wait(x=1230, y=451, sleep_time=1)     # 다음결과 N번 클릭
    #
    # pyautogui.moveTo(x=1891, y=454)                             # 드래그 시작위치
    # pyautogui.dragTo(x=1890, y=983, duration=1, button='left')  # 드래그 끝위치
    #
    # pyautogui.keyDown('shift')
    # time.sleep(1)
    # move_to_click_and_wait(x=1793, y=947, sleep_time=1)     # 다음결과 N번 클릭
    # pyautogui.keyUp('shift')

    # raw_lab=''
    # try: