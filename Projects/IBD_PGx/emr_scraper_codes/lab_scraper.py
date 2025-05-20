from emr_scraper.tools import *

df = pd.read_csv("C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/anti_TNFa_patients.csv")

## Lab 수집

existing_files = glob.glob("C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/IBD-PGx/lab/IBD_PGx_lab(*).xlsx")
existing_id_list = [int(f.split('IBD_PGx_lab(')[-1].split('_')[0]) for f in existing_files]
# move_to_content_window(content='Hx')
# len(existing_id_list)


prev_clipboard_text, raw_lab = '',''
for inx, row in df.iterrows(): #break
    # if row['EMR ID'] not in ['15322168', '19739357', '34835292', '37366865']:
    #     continue
    # if row['name']=='강혁':
    #     break
    row['name'] = row['name'].split('(')[0]

    # if row['name'] != '임성준':
    #     continue

    # new_sleep_time = 20
    # row['name'] = row['name'].split('(')[0]
    # if row['name'] in ['이건수', '김규진', '김은숙', '김기백', '최수호(최덕영 / 개명)']:
    #     new_sleep_time = 360
    #     # continue

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
    for i in range(50):
        move_to_click_and_wait(x=1474, y=979, sleep_time=2)     # 다음결과 N번 클릭
    # move_to_click_and_wait(x=1230, y=451, sleep_time=2, double_click=True, button='right')         # Lab 데이터 첫보고일란 선택
    # pyautogui.moveTo(x=1230, y=451)                             # Lab 데이터 첫보고일란 선택
    move_to_click_and_wait(x=1230, y=451, sleep_time=1)     # 첫보고일란 선택
    pyautogui.dragTo(x=1234, y=573, duration=1, button='left')  # 드래그 끝위치
    pyautogui.moveTo(x=1891, y=461)                             # 드래그 시작위치
    pyautogui.dragTo(x=1890, y=983, duration=15, button='left')  # 드래그 끝위치

    pyautogui.keyDown('shift')
    time.sleep(1)
    move_to_click_and_wait(x=1793, y=947, sleep_time=1)     # 다음결과 N번 클릭
    pyautogui.keyUp('shift')

    # raw_lab=''
    # try:
    while (prev_clipboard_text==raw_lab) or (raw_lab in ['소화기내과',pid,'',"''",'""','소화']) or (len(raw_lab)<10):
        pyautogui.hotkey('ctrl', 'c')

        for i in range(3):
            raw_lab = pyperclip.paste()
            time.sleep(1)

        time.sleep(5)
    # except KeyboardInterrupt:
    #     raise ValueError

    # pyautogui.position()

    print(f"({inx}) {pid} / {row['name']} / {row['IBD type']} / {len(raw_lab)} strings")
    # raise ValueError
    ldf = get_parsed_lab_df(value=raw_lab)
    ldf.to_excel(f"C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/lab/IBD_PGx_lab({pid}_{row['name']}).xlsx",index=False)
    time.sleep(1)

    prev_clipboard_text = raw_lab