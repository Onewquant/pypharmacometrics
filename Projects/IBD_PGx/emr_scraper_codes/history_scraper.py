from emr_scraper.tools import *

df = pd.read_csv("C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/anti_TNFa_patients.csv")

## History 수집

existing_files = glob.glob("C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/hx/IBD_PGx_hx(*).txt")
existing_id_list = [int(f.split('IBD_PGx_hx(')[-1].split('_')[0]) for f in existing_files]
# move_to_content_window(content='Hx')
# len(existing_id_list)



# move_to_content_window(content='Hx')

hx_df = list()
prev_clipboard_text, raw_hx = '', ''
for inx, row in df.iterrows():
    #     if inx==2:
    #         break
    row['name'] = row['name'].split('(')[0]
    pid = row['EMR ID']
    if pid in existing_id_list:
        print(f"({inx}) {pid} / {row['name']} / {row['IBD type']} / 이미 자료가 존재합니다.")
        continue

    # pyautogui.position()

    # 환자ID 입력

    move_to_click_and_wait(x=488, y=57, sleep_time=1)  # 환자번호 입력
    keyboard_input_prep(sleep_time=0)
    keyboard_input(input_str=pid)
    move_to_click_and_wait(x=528, y=57, sleep_time=3)  # 검색버튼 누르기
    move_to_click_and_wait(x=684, y=653, sleep_time=1)  # 환자선택사유 입력
    move_to_click_and_wait(x=1162, y=926, sleep_time=1)  # 환자선택사유 확인

    # 기록일자 '전체'로 입력(조건별 상세조회에서)
    move_to_click_and_wait(x=329, y=190, sleep_time=1)  # 기록일자 Toggle 활성화
    move_to_click_and_wait(x=329, y=209, sleep_time=1)  # '전체'로 선택

    # 진료과- 작성과(소화기내과) 기준으로 설정
    move_to_click_and_wait(x=211, y=262, sleep_time=0)  # 진료과를 '작성과' 기준으로
    move_to_click_and_wait(x=293, y=262, sleep_time=1)  # 입력란 활성화
    keyboard_input_prep(sleep_time=1)
    keyboard_input(input_str='소화기내과', sleep_time=7)  # 소화기내과 입력
    move_to_click_and_wait(x=307, y=341, sleep_time=1)  # Toggle에서 소화기내과 선택
    move_to_click_and_wait(x=363, y=369, sleep_time=3)  # 목록조회 버튼 누르기

    # pyautogui.position()

    # 기록목록 전체 조회
    move_to_click_and_wait(x=27, y=449, sleep_time=1)  # 기록목록 전체 체크
    move_to_click_and_wait(x=585, y=416, sleep_time=15)  # 조회 누르기
    #
    # raw_hx=''
    # try:
    while (prev_clipboard_text == raw_hx) or (raw_hx in ['소화기내과', pid, '', "''", '""', '소화']):
        raw_hx = get_mainpage_text(x=891, y=501)  # 내용 복사
        time.sleep(5)
    # except KeyboardInterrupt:
    #     raise ValueError

    print(f"({inx}) {pid} / {row['name']} / {row['IBD type']} / {len(raw_hx)} strings")

    with open(f"C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/hx/IBD_PGx_hx({pid}_{row['name']}).txt", "w", encoding='utf-8-sig') as f:
        f.write(raw_hx)

    # hx_df.append({'EMR ID': pid, 'name': row['name'], 'IBD type': row['IBD type'], 'raw_Hx': raw_hx})
    # pd.concat([existing_df, pd.DataFrame(hx_df)]).to_excel(
    #     "C:/Users/snubh/PycharmProjects/snubhprac/emr_scraper/IBD_PGx_hx.xlsx", index=False)
    time.sleep(1)

    prev_clipboard_text = raw_hx