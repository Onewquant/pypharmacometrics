from emr_scraper.tools import *

drug_name = 'AMK'
prj_dir = f"C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/{drug_name}"
if not os.path.exists(prj_dir):
    os.makedirs(prj_dir)

df = pd.read_csv(f"{prj_dir}/AMK_tdm.csv",encoding='euc-kr')
df = df.rename(columns={'등록번호':'EMR ID', '검사일':'TDM_DATE','환자명':'name', '담당부서':'DEP'})
df = df.drop_duplicates(['EMR ID'])

## Lab 수집
if not os.path.exists(f"{prj_dir}/lab"):
    os.makedirs(f"{prj_dir}/lab")
existing_files = glob.glob(f"{prj_dir}/lab/AMK_lab(*).xlsx")
existing_id_list = [int(f.split('AMK_lab(')[-1].split('_')[0]) for f in existing_files]
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
        print(f"({inx}) {pid} / {row['name']} / {row['DEP']} / 이미 자료가 존재합니다.")
        continue

    # from emr_scraper.tools import *
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
    # pyautogui.position()
    # 검사결과조회

    move_to_click_and_wait(x=968, y=346, sleep_time=1, double_click=True)  # 기간을 전체로 변경
    move_to_click_and_wait(x=994, y=347, sleep_time=3)  # 재조회
    for i in range(50):
        move_to_click_and_wait(x=1449, y=988, sleep_time=2)     # 다음결과 N번 클릭
    # move_to_click_and_wait(x=1230, y=451, sleep_time=2, double_click=True, button='right')         # Lab 데이터 첫보고일란 선택
    # pyautogui.moveTo(x=1230, y=451)                             # Lab 데이터 첫보고일란 선택
    move_to_click_and_wait(x=1147, y=417, sleep_time=1)  # 첫보고일란 선택
    pyautogui.dragTo(x=1148, y=437, duration=1, button='left')  # 드래그 끝위치
    pyautogui.moveTo(x=1891, y=461)                             # 드래그 시작위치
    pyautogui.dragTo(x=1890, y=983, duration=30, button='left')  # 드래그 끝위치

    pyautogui.keyDown('shift')
    time.sleep(1)
    move_to_click_and_wait(x=1704, y=947, sleep_time=1)     # 다음결과 N번 클릭
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

    print(f"({inx}) {pid} / {row['name']} / {row['DEP']} / {len(raw_lab)} strings")
    # raise ValueError
    ldf = get_parsed_lab_df(value=raw_lab)
    ldf.to_excel(f"{prj_dir}/lab/AMK_lab({pid}_{row['name']}).xlsx",index=False)
    time.sleep(1)

    prev_clipboard_text = raw_lab