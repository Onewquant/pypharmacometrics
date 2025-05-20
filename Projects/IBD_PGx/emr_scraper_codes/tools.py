import os.path
import time

import pandas as pd
import pyautogui
import pyperclip
import ctypes
import re
import numpy as np
import glob
from datetime import datetime, timedelta

def get_parsed_lab_df(value):

    raw_ldf_cols = ['보고일', '오더일', '검사명', '검사결과', '직전결과', '참고치', '결과비고']

    raw_ldf = pd.DataFrame([tbl_row.split('\t') for tbl_row in value.strip().split('\r\n') if tbl_row != ''])
    cur_rldf_cols = list(raw_ldf.columns)
    vld_rldf_cols = list()
    for i in range(len(cur_rldf_cols)):
        if i + 1 > len(cur_rldf_cols):
            vld_rldf_cols.append(str(i))
        else:
            vld_rldf_cols.append(raw_ldf_cols[i])

    raw_ldf.columns = vld_rldf_cols
    if (len(raw_ldf.columns) == 1) or (len(raw_ldf) <= 1):
        ldf = pd.DataFrame(columns=['date', 'dt'])
        return ldf

    for inx, rrow in raw_ldf.iterrows():
        if (rrow['검사명'] == 'WBC') and ('HPF' in rrow['참고치']):
            raw_ldf.at[inx, '검사명'] = 'u.WBC'
        elif (rrow['검사명'] == 'WBC') and ('mm³' in rrow['참고치']):
            raw_ldf.at[inx, '검사명'] = 'em.WBC'
        elif (rrow['검사명'] == 'RBC') and ('HPF' in rrow['참고치']):
            raw_ldf.at[inx, '검사명'] = 'u.RBC'

    raw_ldf['검사명'] = raw_ldf['검사명'].map(lambda x: x.strip() if type(x)==str else '')

    return raw_ldf

# current_position = pyautogui.position()
# pyautogui.doubleClick()
# pyautogui.rightClick()
# print(pyautogui.position())

# # 한글로 전환
# ctypes.windll.user32.PostMessageW(0xFFFF, 0x50, 0, 0x00000112)
#
# # 영어로 전환
# ctypes.windll.user32.PostMessageW(0xFFFF, 0x50, 0, 0x00000113)
# ', '.join(list(df['drug']))
# drug_str = ''
# for inx, drug in enumerate(df['drug']):
#     drug_str +

def get_mainpage_text(x, y, sleep_time=1):
    pyautogui.moveTo(x, y, duration=0)
    pyautogui.click()
    time.sleep(sleep_time)
    pyautogui.hotkey('ctrl', 'a')
    pyautogui.hotkey('ctrl', 'c')

    return pyperclip.paste()

def move_to_content_window(content):

    # pyautogui.position()

    pix_x_dict = {'Hx':906, 'Order':1006, 'Lab':711}
    working_list = (906, 1060)

    pyautogui.moveTo(x=working_list[0], y=working_list[1], duration=0)
    time.sleep(1)
    pyautogui.doubleClick()
    time.sleep(1)
    pyautogui.moveTo(x=pix_x_dict[content], y=975, duration=0)
    pyautogui.click()
    print(f"{content}로 이동함")

def keyboard_input_prep(sleep_time=0):
    pyautogui.hotkey('ctrl', 'a')
    pyautogui.press('delete')
    time.sleep(sleep_time)
def keyboard_input(input_str, sleep_time=0):
    pyperclip.copy(input_str)
    pyautogui.hotkey('ctrl', 'v')
    time.sleep(sleep_time)

def move_to_click_and_wait(x, y, sleep_time=1, double_click=False, button='left'):
    pyautogui.moveTo(x, y, duration=0)
    if double_click:
        pyautogui.doubleClick(button=button)
    else:
        pyautogui.click(button=button)
    time.sleep(sleep_time)
