from tools import *
from pynca.tools import *
import msoffcrypto
import io

# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Dosing 정보 불러오기

dose_result_df = pd.read_csv(f"{output_dir}/final_dose_df.csv")
dose_result_df['ID'] = dose_result_df['ID'].astype(str)
dose_result_df['DOSE_DT'] = dose_result_df['DATE']+'T'+dose_result_df['TIME']
dose_result_df['TIME'] = "T"+dose_result_df['TIME']
dose_result_df = dose_result_df.rename(columns={'DATE':'DOSE_DATE','TIME':'DOSE_TIME'})
dose_result_df = dose_result_df[['ID', 'NAME', 'DOSE_DATE','DOSE_TIME','DOSE_DT', 'DOSE']].copy()
uniq_dose_pids = list(dose_result_df.drop_duplicates(['ID'])['ID'].astype(str))

# 암호 걸린 파일 열기
file_path = f"{resource_dir}/AMK_REQ_DATA/amk_req_conc_data.xlsx"
password = "snubhsnubh"

decrypted = io.BytesIO()
with open(file_path, "rb") as f:
    office_file = msoffcrypto.OfficeFile(f)
    office_file.load_key(password=password)   # 암호 입력
    office_file.decrypt(decrypted)

# pandas로 읽기
raw_df = pd.read_excel(decrypted, engine="openpyxl")
result_cols = ['ID', 'NAME', '보고일', '오더일', 'DRUG', 'CONC', 'POT채혈DT', 'SAMP_DT', 'REC_REASON']


# raw_df = pd.read_excel(f"{resource_dir}/AMK_REQ_DATA/amk_req_drug_order_data.xlsx", engine="openpyxl")
# raw_df.columns
raw_df = raw_df.rename(columns={'환자번호':'ID','검사 접수일자':'DATE','환자명':'NAME', '오더일자':'오더일','검사 결과보고일자':'보고일', '검사 접수일시':'접수DT', '검사결과':'VALUE', '채혈일시':'SAMP_DT'})
raw_df['ID'] = raw_df['ID'].astype(str)
raw_df = raw_df.drop('SEQ',axis=1)
raw_df['SAMP_DT'] = raw_df['SAMP_DT'].dt.strftime("%Y-%m-%dT%H:%M:%S")
raw_df['NEW_SAMP_DT'] = raw_df['SAMP_DT'].map(lambda x:x[:-3])
raw_df['REC_REASON'] = ''
raw_df = raw_df[raw_df['VALUE'] != '중복 오더임'].copy()
raw_df['VALUE'] = raw_df['VALUE'].map(lambda x: float(x.replace('<','').replace('>','').replace('(재검)','')) if type(x)==str else x)

## 오더비고에 존재하는 채혈시각 반영
ordbigo_count = 0
for inx, row in raw_df.iterrows():# break
    if type(row['오더비고'])==float:
        continue
    else:
        # if '채혈시각' in row['오더비고']:
            # raise ValueError
        ordbigo_samp_patterns = re.findall(r'채혈시각\s*\(\s*\d+\s*월\s*\d+\s*일\s*\d+\s*시\s*\d+\s*분\s*.*\)',row['오더비고'])
        if len(ordbigo_samp_patterns)==0:
            continue
        else:
            ordbigo_count+=1
            year_str = row['SAMP_DT'][:4]
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
            if dt_str == '2013-08-16T40:40':
                dt_str = '2013-08-16T09:40'
            raw_df.at[inx, 'NEW_SAMP_DT'] = dt_str
            raw_df.at[inx, 'REC_REASON'] = '오더비고반영'
            print(f'({inx}) / {row["ID"]} / {row["NAME"]} / {ordbigo_count} 번째 / {dt_str}')


## 결과비고에 존재하는 채혈시각 반영
resbigo_count = 0
for inx, row in raw_df[(raw_df['REC_REASON']=='')].iterrows():# break
    if type(row['결과비고'])==float:
        continue
    else:
        resbigo_str = row['결과비고'].upper().replace('MN:','00:').replace('MD:','12:').replace('MD 시','12:').replace('MD ','12:').replace('12/27  MD','12/27 12:00').strip()

        ## 월/일 시간:분 A,P 패턴인 경우
        resbigo_samp_patterns = re.findall(r'\d+/\d+\s+\d+\:?\d*\s*[A|P]', resbigo_str)
        if len(resbigo_samp_patterns)==0:
            raw_df.at[inx, '결과비고'] = resbigo_str
        else:
            # raise ValueError
            resbigo_count += 1
            year_str = row['SAMP_DT'][:4]
            rs_pattern = resbigo_samp_patterns[0].replace('  ',' ')
            month_str = rs_pattern.split('/')[0].strip().zfill(2)
            day_str = rs_pattern.split('/')[-1].split(' ')[0].strip().zfill(2)
            if ':' not in rs_pattern:
                hour_str = rs_pattern.split(':')[0].split(' ')[1].replace('A','').replace('P','').strip().zfill(2)
                minute_str = '00'
                ampm_str = rs_pattern[-1]
            else:
                hour_str = rs_pattern.split(':')[0].split(' ')[1].strip().zfill(2)
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
            year_str = row['SAMP_DT'][:4]
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
            date_str = row['SAMP_DT'].split('T')[0]
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
            date_str = row['SAMP_DT'].split('T')[0]
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
dup_id_sampdate_df = raw_df.groupby(['ID','NEW_SAMP_DT'], as_index=False).agg({'NAME':'count'}).copy()
# test_df = dup_id_sampdate_df.groupby('ID', as_index=False).agg({'NAME':'sum'}).copy()
# test_df[test_df['NAME']]
# dup_id_sampdate_df[dup_id_sampdate_df['NAME']>2]
used_id_dosedt = list()
dupsampdt_count = 0

for inx, row in dup_id_sampdate_df[dup_id_sampdate_df['NAME']>1].copy().iterrows(): #break
    dupsampdt_count+=1
    dup_dt = row['NEW_SAMP_DT']
    dup_date = dup_dt.split('T')[0]
    changing_df = raw_df[(raw_df['ID']==row['ID'])&(raw_df['NEW_SAMP_DT']==dup_dt)].copy()
    mean_value = changing_df['VALUE'].mean()
    changing_df['TEMP_PT'] = (changing_df['VALUE'] > mean_value)*1

    id_dose_df = dose_result_df[dose_result_df['ID']==row['ID']].copy()
    if len(id_dose_df)==0:
        # No dosing record
        continue

    id_dupdate_dose_df = id_dose_df[id_dose_df['DOSE_DATE']==dup_date].copy()

    if len(id_dupdate_dose_df)==0:

        if row['ID'] == '10042359':
            id_dupdate_dose_df = id_dose_df[id_dose_df['DOSE_DATE'] == (datetime.strptime(dup_dt.split('T')[0], '%Y-%m-%d') - timedelta(days=1)).strftime('%Y-%m-%d')].copy()
            # continue
        # if row['ID']=='10014664':
        #     id_dupdate_dose_df = id_dose_df[id_dose_df['DOSE_DATE']==(datetime.strptime(dup_dt.split('T')[0],'%Y-%m-%d') + timedelta(days=1)).strftime('%Y-%m-%d')].copy()
        #     continue
        min_dose_date = id_dose_df['DOSE_DATE'].min()
        max_dose_date = id_dose_df['DOSE_DATE'].max()

        # 근처 투약력 있는 경우 -> 다음 투약 일로 추정
        if (min_dose_date <= dup_date) and (max_dose_date >= dup_date):
            id_dupdate_dose_df = id_dose_df[id_dose_df['DOSE_DATE'] == (datetime.strptime(dup_dt.split('T')[0], '%Y-%m-%d') + timedelta(days=1)).strftime('%Y-%m-%d')].copy()
            if len(id_dupdate_dose_df)==0:
                for cg_inx, cg_row in changing_df.iterrows():  # break
                    raw_df.at[cg_inx, 'NEW_SAMP_DT'] = np.nan
                    raw_df.at[cg_inx, 'REC_REASON'] = '투약력부재로유추불가'
                continue
        # 근처 투약력 투약력 없는 경우
        else:
            for cg_inx, cg_row in changing_df.iterrows():  # break
                raw_df.at[cg_inx, 'NEW_SAMP_DT'] = np.nan
                raw_df.at[cg_inx, 'REC_REASON'] = '투약력부재로유추불가'
            # ('10137505', '10125999', '10216802')
            # No dosing record
            continue

    ## 중복된 DOSE 기록이 없게 DOSE DT 선택

    for ddt_inx, dose_dt in enumerate(id_dupdate_dose_df['DOSE_DT']):
        id_dosedt_str = f"{row['ID']}DT{dose_dt}"
        if dose_dt not in used_id_dosedt:
            used_id_dosedt.append(id_dosedt_str)
            break
        else:
            ddt_inx += 1
            continue
    if (ddt_inx==len(id_dupdate_dose_df)):
        raise ValueError

    pre_dt = (datetime.strptime(dose_dt,'%Y-%m-%dT%H:%M')-timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')
    post_dt = (datetime.strptime(dose_dt,'%Y-%m-%dT%H:%M')+timedelta(minutes=30)).strftime('%Y-%m-%dT%H:%M')

    for cg_inx, cg_row in changing_df.iterrows(): #break
        if cg_row['TEMP_PT']==0:
            raw_df.at[cg_inx, 'NEW_SAMP_DT'] = pre_dt
            raw_df.at[cg_inx, 'REC_REASON'] = '채혈시간중복(TROUGH)'
            print(f'({cg_inx}) / {cg_row["ID"]} / {cg_row["NAME"]} / {dupsampdt_count} 번째 / {pre_dt}')
        else:
            raw_df.at[cg_inx, 'NEW_SAMP_DT'] = post_dt
            raw_df.at[cg_inx, 'REC_REASON'] = '채혈시간중복(PEAK)'
            print(f'({cg_inx}) / {cg_row["ID"]} / {cg_row["NAME"]} / {dupsampdt_count} 번째 / {post_dt}')


# dup_id_sampdate_df = raw_df.groupby(['ID','NEW_SAMP_DT'], as_index=False).agg({'NAME':'count'}).copy()
final_conc_df = raw_df.copy()
final_conc_df['SAMP_DT'] = final_conc_df['NEW_SAMP_DT'].copy()
final_conc_df['DRUG'] = 'amikacin'
final_conc_df['CONC'] = final_conc_df['VALUE'].copy()
final_conc_df = final_conc_df[['ID', 'NAME', '보고일', '오더일', 'DRUG', 'CONC', 'SAMP_DT', 'REC_REASON']].copy()
final_conc_df = final_conc_df[~final_conc_df['SAMP_DT'].isna()].copy()
final_conc_df.to_csv(f"{output_dir}/final_conc_df(with sampling).csv", encoding='utf-8-sig', index=False)