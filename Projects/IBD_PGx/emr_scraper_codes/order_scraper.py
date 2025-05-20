from emr_scraper.tools import *

prj_dir = "C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/IBD-PGx"

df = pd.read_csv(f"{prj_dir}/anti_TNFa_patients.csv")

## ORDER 수집

existing_files = glob.glob(f"{prj_dir}/order/IBD_PGx_order(*).xlsx")
existing_id_list = [int(f.split('IBD_PGx_order(')[-1].split('_')[0]) for f in existing_files]
# move_to_content_window(content='Hx')
# len(existing_id_list)

# df['EMR ID']
prev_clipboard_text = ''
for inx, row in df.iterrows(): #break

    # if row['EMR ID'] not in [15322168, 19739357, 34835292, 37366865, 21618097, 36898756,36975211,37858047]:
    # continue if row['EMR ID'] not in ['14188505', '17677819', '21169146', '21201336', '24028105', '24106625',
    # '25269024','37858047']: continue if row['EMR ID'] =='10875838': break if row['name']!='임정진': continue

    # new_sleep_time = 20
    # row['name'] = row['name'].split('(')[0]
    # if row['name'] in ['이학준', ]:
    #     # new_sleep_time = 360
    #     print(f"({inx}) {pid} / {row['name']} / 코드에러로 넘어감. Fix 필요")
    #     continue

    prev2_sp_order, prev_sp_order = '',''
#     print(1)
#     if inx==2:
#         break
    row['name'] = row['name'].split('(')[0]
    if row['name']=='최덕영': row['name']='최수호'
    pid = row['EMR ID']
    if pid in existing_id_list:
        print(f"({inx}) {pid} / {row['name']} / {row['IBD type']} / 이미 자료가 존재합니다.")
        continue

    # 수집된 Lab 데이터의 초기날짜 가져오기

    pid_lab_df = pd.read_excel(f"{prj_dir}/lab/IBD_PGx_lab({pid}_{row['name']}).xlsx")


    # lab_report_dates = list()
    # for labdate in pid_lab_df['보고일'].unique():
    #     if (type(labdate) == str):
    #         if re.match(r'[\d]+-[\d]+-[\d]+',labdate):
    #             lab_report_dates.append(labdate)
    #
    # lab_order_dates = list()
    # for labdate in pid_lab_df['오더일'].unique():
    #     if (type(labdate) == str):
    #         if re.match(r'[\d]+-[\d]+-[\d]+', labdate):
    #             lab_order_dates.append(labdate)
    #
    # min_date = min(min(lab_report_dates),min(lab_order_dates))

    min_date = min(pid_lab_df['보고일'].min(), pid_lab_df['오더일'].min())
    # pid_lab_df['보고일'].unique()

    today_date = datetime.today().strftime('%Y-%m-%d')

    # 수집된 Order_frag 데이터의 마지막날짜 가져오기
    # raise ValueError
    existing_order_frag_files = glob.glob(f"{prj_dir}/order_frag/{pid}_{row['name']}/IBD_PGx_order({pid}_{row['name']}_*).xlsx")
    existing_max_order_frag_inx=np.nan
    # len(existing_order_frag_files)
    if len(existing_order_frag_files)>0:
        existing_max_order_frag_inx = max([int(ex_offile.split(row['name']+'_')[-1].split(')')[0]) for ex_offile in existing_order_frag_files])

    # raise ValueError
    # 7일 간격으로 나누기
    days_interval = 360
    # if row['name'] in ['한상희']:
    #     days_interval = 360
    # else:
    #     days_interval = 60

    min_dt = datetime.strptime(min_date, '%Y-%m-%d')
    today_dt = datetime.strptime(today_date, '%Y-%m-%d')

    intervals = list()
    current_start = min_dt

    while current_start < today_dt:
        current_end = min(current_start + timedelta(days=days_interval), today_dt)
        intervals.append({
            "start_date": current_start.strftime('%Y-%m-%d'),
            "end_date": current_end.strftime('%Y-%m-%d')
        })
        current_start = current_end + timedelta(days=1)

    intervals_df = pd.DataFrame(intervals)

    # pyautogui.position()

    # 환자ID

    move_to_click_and_wait(x=135, y=57, sleep_time=1)  # 입력란 활성화
    keyboard_input_prep(sleep_time=1)
    keyboard_input(input_str=pid, sleep_time=1)  # 환자번호 입력

    # 날짜 입력

    for date_inx, date_row in intervals_df.iterrows():

        if float(date_inx) <= existing_max_order_frag_inx:
            print(f"({inx}-{date_inx}) / {pid} / {row['name']} 이미 존재하는 order frag 입니다")
            continue

        # raise ValueError

        move_to_click_and_wait(x=559, y=56, sleep_time=1.5)  # 입력란 활성화
        keyboard_input_prep(sleep_time=1.5)
        # keyboard_input(input_str=today_date, sleep_time=1)  # 끝날짜 입력
        keyboard_input(input_str=date_row['end_date'], sleep_time=3)  # 끝날짜 입력


        move_to_click_and_wait(x=443, y=58, sleep_time=1.5)  # 입력란 활성화
        keyboard_input_prep(sleep_time=1.5)
        # keyboard_input(input_str=min_date, sleep_time=1)  # 시작날짜 입력
        keyboard_input(input_str=date_row['start_date'], sleep_time=3)  # 시작날짜 입력

        # move_to_click_and_wait(x=654, y=61, sleep_time=120)  # 조회 누르기
        # move_to_click_and_wait(x=654, y=61, sleep_time=90)  # 조회 누르기
        move_to_click_and_wait(x=654, y=61, sleep_time=60)  # 조회 누르기
        # move_to_click_and_wait(x=654, y=61, sleep_time=15)  # 조회 누르기
        # move_to_click_and_wait(x=654, y=61, sleep_time=30)  # 조회 누르기

        waiting_message1 = '[Window Title]\n차세대 의료정보시스템\n\n[Main Instruction]\n차세대 의료정보시스템이(가) 응답하지 않습니다.\n\n[Content]\n프로그램을 닫으면 정보를 잃을 수 있습니다.\n\n[프로그램 닫기] [프로그램 응답 대기]'
        waiting_message2 = '[Window Title]\r\n차세대 의료정보시스템\r\n\r\n[Main Instruction]\r\n차세대 의료정보시스템이(가) 응답하지 않습니다.\r\n\r\n[Content]\r\nWindows에서 프로그램 복원을 시도할 수 있습니다. 프로그램을 복원하거나 닫으면 정보를 잃을 수 있습니다.\r\n\r\n[프로그램 복원 시도] [프로그램 닫기] [프로그램 응답 대기]'
        waiting_message1_pos = {'x':863, 'y':279}
        waiting_message2_pos = {'x':870, 'y':253}

        ## 오더 뜨는 것 기다리기

        test_raw_order=''
        order_waiting_count = 0
        while test_raw_order in ['', waiting_message1, waiting_message2, min_date, today_date]:
            move_to_click_and_wait(x=350, y=217, sleep_time=5)  # 첫보고일란 선택
            test_raw_order=''
            pyperclip.copy("")
            for i in range(3):
                pyautogui.hotkey('ctrl', 'c')
                time.sleep(2)
            test_raw_order = pyperclip.paste()
            # print(f"test raw order : {test_raw_order}")
            if test_raw_order==waiting_message1:
                move_to_click_and_wait(x=863, y=279, sleep_time=2)  # waiting message 닫기
            elif test_raw_order==waiting_message1:
                move_to_click_and_wait(x=870, y=253, sleep_time=2)  # waiting message 닫기
            elif test_raw_order=='':
                # continue
                break
            elif len(test_raw_order.strip().split('\t')) > 1:
                # print('end')
                break
            time.sleep(30)
            print(f"({inx}-{date_inx}) / order waiting count: {order_waiting_count} / {pid} / {row['name']} / 오더 조회 기다리는 중")
            order_waiting_count+=1

        if test_raw_order=='':
            print(f"start: {date_row['start_date']} / end: {date_row['end_date']} / 다음 기간으로 넘어감")
            continue
        # else:
        #     pyperclip.copy("")
        #     move_to_click_and_wait(x=350, y=217, sleep_time=1)  # 첫보고일란 선택해제
        #     pyautogui.hotkey('ctrl', 'c')
        #     time.sleep(2)
        #     print(pyperclip.paste())
        # raise ValueError


        # 검사결과조회

        raw_order = ''
        raw_order_count = 0
        while raw_order in ['', waiting_message1, waiting_message2, min_date, today_date]:
            # pyautogui.position()
            # if raw_order_count == 0:
            #     pass
            # else:
            # raise ValueError

            move_to_click_and_wait(x=350, y=217, sleep_time=1)  # 첫보고일란 선택
            # pyautogui.scroll(100)
            # primer_delta = 25
            # for primer_order_inx in range(1,7):
            #     pyautogui.keyDown('shift')
            #     time.sleep(1)
            #     move_to_click_and_wait(x=696, y=217+primer_delta*primer_order_inx, sleep_time=1)  # 다음결과 N번 클릭
            #     pyautogui.scroll(100)
            #     pyautogui.keyUp('shift')
            pyautogui.keyDown('shift')
            pyautogui.dragTo(x=350, y=217, duration=1, button='left')  # 드래그 시작위치
            # pyautogui.scroll(100)
            pyautogui.dragTo(x=350, y=207, duration=1, button='left')  # 드래그 시작위치
            # pyautogui.dragTo(x=350, y=217, duration=1, button='left')  # 드래그 시작위치
            pyautogui.moveTo(x=353, y=217)  # 드래그 무빙
            pyautogui.dragTo(x=353, y=683, duration=1, button='left')  # 드래그 끝위치
            pyautogui.keyUp('shift')

            test_raw_order = ''                         # 여러 줄이 잘 드래그 되었는지 확인
            pyperclip.copy("")
            for i in range(3):
                pyautogui.hotkey('ctrl', 'c')
                time.sleep(2)
            test_raw_order = pyperclip.paste()

            if len(test_raw_order.split('\n'))<=2:
                if row['name'] in ['한상희']:
                    pass
                else:
                    continue
                # continue
            # else:
            # raise ValueError

            move_to_click_and_wait(x=1386, y=132, sleep_time=1)  # 첫보고일란 선택
            pyautogui.dragTo(x=1386, y=132, duration=1, button='left')  # 드래그 시작위치
            # pyautogui.scroll(100)
            pyautogui.moveTo(x=1386, y=132)  # 드래그 무빙
            pyautogui.dragTo(x=1387, y=683, duration=3, button='left')  # 드래그 끝위치

            pyautogui.keyDown('shift')
            time.sleep(1)
            move_to_click_and_wait(x=696, y=660, sleep_time=1)  # 다음결과 N번 클릭
            pyautogui.keyUp('shift')

            raw_order = ''
            pyperclip.copy("")
            for i in range(3):
                pyautogui.hotkey('ctrl', 'c')
                time.sleep(2)
                raw_order = pyperclip.paste()
                time.sleep(1)
            # print(f"raw order : {raw_order}")
            raw_order_str_list = raw_order.split('\n')
            # if len(raw_order_str_list[0].split('\t'))>1:
            #     break

            # move_to_click_and_wait(x=350, y=217, sleep_time=1)  # 첫보고일란 선택

            # print(raw_order)

            # move_to_click_and_wait(x=1383, y=658, sleep_time=1)  # 드래그 끝위치
            # pyautogui.dragTo(x=1383, y=658, duration=1, button='left')  # 드래그 끝위치
            # # pyautogui.scroll(-100)
            # pyautogui.moveTo(x=1386, y=132)  # 드래그 무빙
            # pyautogui.dragTo(x=1386, y=132, duration=3, button='left')  # 드래그 시작

            pyautogui.scroll(5000)

            # raise ValueError

            # move_to_click_and_wait(x=350, y=217, sleep_time=1)  # 첫보고일란 선택

        order_cols = ['VAC','처방지시','발행처','처방의','수납','약국_검사','주사시행처','Acting','변경의']
        order_frag_df = pd.DataFrame([order_row.split('\t') for order_row in raw_order_str_list], columns=order_cols)
        order_frag_df = order_frag_df.drop(columns=['VAC'])

        print(f"({inx}-{date_inx}) {pid} / {row['name']} / {row['IBD type']} / {len(order_frag_df)} rows / start: {date_row['start_date']} / end: {date_row['end_date']}")

        if not os.path.exists(f"{prj_dir}/order_frag/{pid}_{row['name']}"):
            os.makedirs(f"{prj_dir}/order_frag/{pid}_{row['name']}")
        order_frag_df.to_excel(f"{prj_dir}/order_frag/{pid}_{row['name']}/IBD_PGx_order({pid}_{row['name']}_{date_inx}).xlsx",index=False)
        time.sleep(1)
        # order_frag_df.columns
        # order_frag_df['처방의']

        # raise ValueError

        prev_clipboard_text = raw_order
    # pid=28275802
    # row={'name':'김민직'}
    order_frag_files = glob.glob(f"{prj_dir}/order_frag/{pid}_{row['name']}/IBD_PGx_order({pid}_{row['name']}_*).xlsx")
    order_df = list()
    for ordfile in order_frag_files:
        order_frag_df = pd.read_excel(ordfile)
        order_df.append(order_frag_df)
    # if len(order_df)=
    order_df = pd.concat(order_df, ignore_index=True)
    order_df.to_excel(f"{prj_dir}/order/IBD_PGx_order({pid}_{row['name']}).xlsx", index=False)


