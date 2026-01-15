import numpy as np
import pandas as pd
import msoffcrypto
import io, os, re

# result_type = 'R'

prj_name = 'Valoproic acid'
prj_dir = f'C:/Users/tldus/PycharmProjects/DILI/resource/DILI_{prj_name}'
resource_dir = f'{prj_dir}/CONC'
output_dir = f'C:/Users/tldus/PycharmProjects/DILI/resource/prep'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# # Dosing 정보 불러오기
#
# dose_result_df = pd.read_csv(f"{output_dir}/final_dose_df.csv")
# dose_result_df['ID'] = dose_result_df['ID'].astype(str)
# dose_result_df['DOSE_DT'] = dose_result_df['DATE']+'T'+dose_result_df['TIME']
# dose_result_df['TIME'] = "T"+dose_result_df['TIME']
# dose_result_df = dose_result_df.rename(columns={'DATE':'DOSE_DATE','TIME':'DOSE_TIME'})
# dose_result_df = dose_result_df[['ID', 'NAME', 'DOSE_DATE','DOSE_TIME','DOSE_DT', 'DOSE']].copy()
# uniq_dose_pids = list(dose_result_df.drop_duplicates(['ID'])['ID'].astype(str))

# 암호 걸린 파일 열기
file_path = f"{resource_dir}/202512-00178 BESTCare 자료제공 신청서 DILI 연구 데이터(진단검사)2.xlsx"
password = "snubhsnubh"

decrypted = io.BytesIO()
with open(file_path, "rb") as f:
    office_file = msoffcrypto.OfficeFile(f)
    office_file.load_key(password=password)   # 암호 입력
    office_file.decrypt(decrypted)

# pandas로 읽기
raw_df = pd.read_excel(decrypted, engine="openpyxl")
raw_df = raw_df[~raw_df['검사결과'].isna()].copy()


result_cols = ['ID', 'NAME', '보고일', '오더일', 'DRUG', 'CONC', 'SAMP_DT_ORI', 'SAMP_DT', 'REC_REASON']
# raw_df.columns

# raw_df = pd.read_excel(f"{resource_dir}/AMK_REQ_DATA/amk_req_drug_order_data.xlsx", engine="openpyxl")
# raw_df.columns
raw_df = raw_df.rename(columns={'환자번호':'ID','접수일자':'접수일','환자명':'NAME', '오더일자':'오더일','결과보고일자':'보고DT', '접수일시':'접수DT', '검사결과':'VALUE', '채혈일시':'SAMP_DT','라벨출력일시':'LABEL_DT'})
raw_df['ID'] = raw_df['ID'].astype(str)
raw_df['DATE'] = raw_df['접수일'].copy()
raw_df['접수DT'] = raw_df['접수일'] + 'T' + raw_df['접수DT']
# raw_df = raw_df.drop('SEQ',axis=1)
raw_df['LABEL_DT'] = raw_df['LABEL_DT'].dt.strftime("%Y-%m-%dT%H:%M:%S")
raw_df['SAMP_DT'] = raw_df['SAMP_DT'].dt.strftime("%Y-%m-%dT%H:%M:%S")


raw_df['NEW_SAMP_DT'] =raw_df.apply(lambda x:x['LABEL_DT'] if type(x['SAMP_DT'])==float else x['SAMP_DT'], axis=1)
# raw_df['NEW_SAMP_DT'] = raw_df['NEW_SAMP_DT'].map(lambda x:x[:-3])
raw_df['REC_REASON'] = ''

## LLOQ 고려
# raw_df['LLOQ'] = raw_df['VALUE'].map(lambda x: float(x.replace('<','').replace('(재검)','')) if '<' in x else np.nan).fillna(method='bfill').fillna(method='ffill')
# raw_df['BLQ'] = raw_df['VALUE'].map(lambda x: '<' in x)*1
# raw_df['BLQ'].sum()

## 농도결과 숫자로만 !
val_list = list()
for x in raw_df['VALUE']:
    if type(x)==str:
        # if '<' in x:
        #     val = 0.0
        # else:
        #     val = float(x.replace('<','').replace('>','').replace('(재검)',''))
        val = float(x.replace('<', '').replace('>', '').replace('(재검)', ''))
    else:
        val = x
    val_list.append(val)
raw_df['VALUE'] = val_list

## LLOQ 고려 안 할때
# raw_df['VALUE'] = raw_df['VALUE'].map(lambda x: '0' if '<' in x else x)
# raw_df['VALUE'] = raw_df['VALUE'].map(lambda x: float(x.replace('<','').replace('>','').replace('(재검)','')) if type(x)==str else x)


## 오더비고에 존재하는 채혈시각 반영
ordbigo_count = 0
ordbigo_patients = set()
for inx, row in raw_df.iterrows():# break
    if type(row['결과비고'])==float:
        continue
    else:
        # if '채혈시각' in row['오더비고']:
            # raise ValueError
        row['결과비고'] = row['결과비고'].replace('MD   시','12   시')
        ordbigo_samp_patterns = re.findall(r'채혈시각\s*\(\s*\d+\s*월\s*\d+\s*일\s*\d+\s*시\s*\d+\s*분\s*.*\)',row['결과비고'])
        # ordbigo_samp_patterns[0]
        if len(ordbigo_samp_patterns)==0:
            continue
        else:
            raise ValueError
            ordbigo_patients.add(row['ID'])
            ordbigo_count+=1
            year_str = row['NEW_SAMP_DT'][:4]
            month_str = ordbigo_samp_patterns[0].split('(')[-1].split('월')[0].strip().zfill(2)
            day_str = ordbigo_samp_patterns[0].split('월')[-1].split('일')[0].strip().zfill(2)
            hour_str = ordbigo_samp_patterns[0].split('일')[-1].split('시')[0].strip().zfill(2)
            minute_str = ordbigo_samp_patterns[0].split('시')[-1].split('분')[0].strip().zfill(2)
            ampm_str = ordbigo_samp_patterns[0].split('분')[-1].split(')')[0].strip().zfill(2)
            if 'PM' in ampm_str.upper():
                corrected_hour = int(hour_str) + 12
                if corrected_hour < 24:
                    hour_str=str(corrected_hour).zfill(2)
                else:
                    pass
                # raise ValueError
            elif 'AM' in ampm_str.upper():
                pass
                # raise ValueError
            elif 'MD' in ampm_str.upper():
                pass
                # raise ValueError
            elif ampm_str.upper()=='00':
                # raise ValueError
                pass
            else:
                print(ampm_str.upper())
                raise ValueError
            dt_str = f"{year_str}-{month_str}-{day_str}T{hour_str}:{minute_str}"
            # if (row['ID']=='11952372') and (dt_str == '2013-08-16T40:40'):
            #     dt_str = '2013-08-16T09:40'
            # elif (row['ID']=='13523484') and (dt_str == '2010-07-20T23:30'):
            #     dt_str = '2010-07-20T11:30'
            # elif (row['ID']=='13523484') and (dt_str == '2010-07-21T01:30'):
            #     dt_str = '2010-07-20T13:30'
            # elif (row['ID']=='11657204') and (dt_str == '2010-06-23T22:00'):
            #     # raise ValueError
            #     dt_str = '2010-06-23T10:00'
            raw_df.at[inx, 'NEW_SAMP_DT'] = dt_str
            raw_df.at[inx, 'REC_REASON'] = '결과비고반영'
            print(f'({inx}) / {row["ID"]} / {row["NAME"]} / {ordbigo_count} 번째 / {dt_str}')


## 결과비고에 존재하는 채혈시각 반영
resbigo_count = 0
resbigo_patients = set()
for inx, row in raw_df[(raw_df['REC_REASON']=='')].iterrows():# break
    if type(row['결과비고'])==float:
        continue
    else:
        # if (row['ID']=='11504533') and (row['VALUE']==37.08):
        #     raise ValueError

        resbigo_str = row['결과비고'].upper().replace('MN:','00:').replace('MD:','12:').replace('MD 시','12:').replace('MD ','12:').replace('12/27  MD','12/27 12:00').strip()

        ## 월/일 시간:분 A,P 패턴인 경우
        resbigo_samp_patterns = re.findall(r'\d+/\d+\s+\d+\:?\d*\s*[A|P]', resbigo_str)
        if len(resbigo_samp_patterns)==0:
            raw_df.at[inx, '결과비고'] = resbigo_str
        else:
            # raise ValueError
            resbigo_count += 1
            resbigo_patients.add(row['ID'])
            year_str = row['NEW_SAMP_DT'][:4]
            rs_pattern = resbigo_samp_patterns[0].replace('  ',' ')
            month_str = rs_pattern.split('/')[0].strip().zfill(2)
            day_str = rs_pattern.split('/')[-1].split(' ')[0].strip().zfill(2)
            if ':' not in rs_pattern:
                hour_str = rs_pattern.split(':')[0].split(' ')[1].replace('A','').replace('P','').strip().zfill(2)
                minute_str = '00'
                ampm_str = rs_pattern[-1]
            else:
                hour_str = rs_pattern.split(':')[0].split(' ')[-1].strip().zfill(2)
                minute_str = rs_pattern.split(':')[-1].replace('A','').replace('P','').strip().zfill(2)
                ampm_str = rs_pattern[-1]
            if 'P' in ampm_str:
                corrected_hour = int(hour_str) + 12
                if corrected_hour < 24:
                    hour_str=str(corrected_hour).zfill(2)
                else:
                    pass
                # print(ampm_str)
                # hour_str = str(int(hour_str) + 12).zfill(2)
                # raise ValueError
            elif 'A' in ampm_str:
                pass
                # print(ampm_str)
                # raise ValueError
            else:
                # print(ampm_str)
                raise ValueError
            dt_str = f"{year_str}-{month_str}-{day_str}T{hour_str}:{minute_str}"
            # if row['SAMP_DT']=='2003-11-30T21:04:59':
            #     raise ValueError
            raw_df.at[inx, 'NEW_SAMP_DT'] = dt_str
            raw_df.at[inx, 'REC_REASON'] = '결과비고반영(날짜_시간AP)'
            # raw_df.at[inx, '결과비고'] = resbigo_str
            print(f'({inx}) / {row["ID"]} / {row["NAME"]} / {resbigo_count} 번째 / {dt_str}')
            continue

        ## 월/일 시간:분 패턴인 경우
        resbigo_samp_patterns = re.findall(r'\d+/\d+\s+\d+\:?\d*\s*', resbigo_str)
        if len(resbigo_samp_patterns) == 0:
            pass
        else:
            # raise ValueError
            resbigo_count += 1
            resbigo_patients.add(row['ID'])
            year_str = row['NEW_SAMP_DT'][:4]
            rs_pattern = resbigo_samp_patterns[0].replace('  ', ' ')
            month_str = rs_pattern.split('/')[0].strip().zfill(2)
            day_str = rs_pattern.split('/')[-1].split(' ')[0].strip().zfill(2)
            if ':' not in rs_pattern:
                hour_str = rs_pattern.split(':')[0].split(' ')[1].replace('A', '').replace('P', '').strip().zfill(2)
                minute_str = '00'
            else:
                hour_str = rs_pattern.split(':')[0].split(' ')[1].strip().zfill(2)
                minute_str = rs_pattern.split(':')[-1].replace('A', '').replace('P', '').strip().zfill(2)
            dt_str = f"{year_str}-{month_str}-{day_str}T{hour_str}:{minute_str}"
            # if row['SAMP_DT']=='2003-11-30T21:04:59':
            #     raise ValueError
            raw_df.at[inx, 'NEW_SAMP_DT'] = dt_str
            raw_df.at[inx, 'REC_REASON'] = '결과비고반영(날짜_시간)'
            # raw_df.at[inx, '결과비고'] = resbigo_str
            print(f'({inx}) / {row["ID"]} / {row["NAME"]} / {resbigo_count} 번째 / {dt_str}')
            continue


        ## 시간:분 A,P 패턴인 경우
        # resbigo_str = resbigo_str.replace('AMK 450MG TID','')
        resbigo_samp_patterns = re.findall(r'\d+\:\d*\s*[A|P]', resbigo_str)
        if len(resbigo_samp_patterns)==0:
            pass
        else:
            # raise ValueError
            resbigo_count += 1
            resbigo_patients.add(row['ID'])
            date_str = row['NEW_SAMP_DT'].split('T')[0]
            rs_pattern = resbigo_samp_patterns[0].replace('  ',' ')
            # month_str = rs_pattern.split('/')[0].strip().zfill(2)
            # day_str = rs_pattern.split('/')[-1].split(' ')[0].strip().zfill(2)
            if ':' not in rs_pattern:
                hour_str = rs_pattern.split(':')[0].split(' ')[0].replace('A','').replace('P','').strip().zfill(2)
                minute_str = '00'
                ampm_str = rs_pattern[-1]
            else:
                hour_str = rs_pattern.split(':')[0].split(' ')[0].strip().zfill(2)
                minute_str = rs_pattern.split(':')[-1].replace('A','').replace('P','').strip().zfill(2)
                ampm_str = rs_pattern[-1]
            if 'P' in ampm_str:
                corrected_hour = int(hour_str) + 12
                if corrected_hour < 24:
                    hour_str=str(corrected_hour).zfill(2)
                else:
                    pass
                # raise ValueError
            elif 'A' in ampm_str:
                pass
                # print(ampm_str)
                # raise ValueError
            else:
                # print(ampm_str)
                raise ValueError
            dt_str = f"{date_str}T{hour_str}:{minute_str}"
            # if row['SAMP_DT']=='2003-11-30T21:04:59':
            #     raise ValueError
            raw_df.at[inx, 'NEW_SAMP_DT'] = dt_str
            raw_df.at[inx, 'REC_REASON'] = '결과비고반영(시간_분AP)'
            print(f'({inx}) / {row["ID"]} / {row["NAME"]} / {resbigo_count} 번째 / {dt_str}')
            continue

        ## 시간A,P분 패턴인 경우
        # resbigo_str = resbigo_str.replace('AMK 450MG TID','')
        resbigo_samp_patterns = re.findall(r'\d+[A|P]\d+\s*', resbigo_str)
        if len(resbigo_samp_patterns) == 0:
            pass
        else:
            # raise ValueError
            resbigo_count += 1
            resbigo_patients.add(row['ID'])
            date_str = row['NEW_SAMP_DT'].split('T')[0]
            rs_pattern = resbigo_samp_patterns[0].replace('  ', ' ')
            # month_str = rs_pattern.split('/')[0].strip().zfill(2)
            # day_str = rs_pattern.split('/')[-1].split(' ')[0].strip().zfill(2)
            if 'A' in rs_pattern:
                ap_str = 'A'
            elif 'P' in rs_pattern:
                ap_str = 'P'
            else:
                ap_str = ''
                raise ValueError
            hour_str = rs_pattern.split(ap_str)[0].strip().zfill(2)
            minute_str = rs_pattern.split(ap_str)[-1].replace(ap_str, '').strip().zfill(2)
            dt_str = f"{date_str}T{hour_str}:{minute_str}"
            # if row['SAMP_DT']=='2003-11-30T21:04:59':
            #     raise ValueError
            raw_df.at[inx, 'NEW_SAMP_DT'] = dt_str
            raw_df.at[inx, 'REC_REASON'] = '결과비고반영(시간AP분)'
            print(f'({inx}) / {row["ID"]} / {row["NAME"]} / {resbigo_count} 번째 / {dt_str}')
            continue


# 채혈시간이 중복되는 데이터 처리 -> DOSING TIME 고려하여 P/T 구분
dup_id_sampdate_df = (raw_df.dropna(subset=['SAMP_DT'])).groupby(['ID','NEW_SAMP_DT'], as_index=False).agg({'NAME':'count'}).copy()