from tools import *
from pynca.tools import *

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Orders
drug_order_set = set()

order_files = glob.glob(f'{resource_dir}/order/{prj_name}_order(*).xlsx')
dose_result_df = list()
wierd_result_df = list()
result_cols = ['ID','NAME','DATETIME','DRUG','DOSE','ROUTE','ACTING','PERIOD','ETC_INFO','PLACE']
no_dup_cols = [c for c in result_cols if c!='NAME']
for finx, fpath in enumerate(order_files): #break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    print(f"({finx}) {pname} / {pid}")

    # if pid in (,):       # lab, order 파일 다시 수집 필요
    #     continue

    fdf = pd.read_excel(fpath)
    raise ValueError

    # fdf.columns
    # fdf.to_csv(f"{outcome_dir}/error_dose_df.csv", encoding='utf-8-sig', index=False)

    dose_df = fdf[fdf['처방지시'].map(lambda x: (('amikacin' in x.lower()) or ('infliximab' in x.lower()) or ('ustekinumab' in x.lower())) and ('quantification' not in x.lower()))].copy()
    # dose_df.iloc[0]['처방지시']
    # dose_df.iloc[1]['처방지시']
    # dose_df.iloc[2]

    dose_df = dose_df[(~dose_df['Acting'].isna())].copy()

    dose_df['처방지시비고'] = dose_df['처방지시'].map(lambda x:x.split(' : ')[-1] if len(x.split(' : '))>1 else '')

    dose_df['ID'] = pid
    dose_df['NAME'] = pname

    # if pname=='김옥순': raise ValueError

    #### 약국_검사가 NA인 것은 제외하고 진행함 (추후 필요시 추가)

    # dose_df[dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test.xlsx",index=False)
    # dose_df[~dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test2.xlsx",index=False)
    dose_df = dose_df[(~dose_df['약국_검사'].isna())].copy()

    #### 성분명이 IBD Biologics 인 경우만으로 필터링 (Infliximab, Adalimumab)

    regex_pattern = r'\(infliximab|\(adalimumab|\(ustekinumab'
    dose_df = dose_df[dose_df['처방지시'].map(lambda x: bool(re.search(regex_pattern, x, flags=re.IGNORECASE)))].copy()

    dose_df['DT1'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[0].split('[')[-1].replace(' ','T'))
    dose_df['DT2'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[-1].split('[')[-1].replace(' ','T')[:-3])
    dose_df['DATETIME'] = dose_df[['DT1','DT2']].min(axis=1)
    dose_df['ETC_INFO'] = dose_df['처방지시비고'].copy()
    dose_df['DRUG'] = dose_df['처방지시'].map(lambda x: re.search(regex_pattern, x, flags=re.IGNORECASE).group().lower().replace('(','').replace(')',''))
    dose_df['ROUTE'] = dose_df['처방지시'].map(lambda x: 'IV' if " [SC] " not in x else 'SC')
    dose_df['PLACE'] = dose_df['주사시행처']
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+',x.split('▣')[-1].split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split(' ')[-1].strip())[0] if " [SC] " not in x else x.split('(Infliximab)')[-1].split('(Adalimumab)')[-1].split('▣')[-1].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+', x.split('mg')[1].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0] if (" [SC] " not in x) else x.split('(Infliximab')[-1].split('(Adalimumab')[-1].split('(Ustekinumab')[-1].split('▣')[-1].split('srg')[0].split('via')[0].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # mg_inx_dict = {'ustekinumab':0,'infliximab':1,'adalimumab':1}
    dose_series = list()
    for dose_inx, dose_row in dose_df.iterrows():
        x = dose_row['처방지시']
        if (" [SC] " not in x):
            # raise ValueError
            try: dose_val = int(re.findall(r'\d+mg',x)[1].replace('mg','').strip())
            except:
                dose_val = int(re.findall(r'\d+mg', x)[0].replace('mg', '').strip()) * int(re.findall(r'\d+ via', x)[0].replace('via', '').strip())
            # dose_val = re.findall(r'\d+',x.split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0]
            # raise ValueError
            dose_series.append(dose_val)
        else:
            try: dose_val = int(re.findall(r'\d+mg',x)[1].replace('mg','').strip())
            except:
                if '(Ustekinumab' in x:
                    # raise ValueError
                    dose_val = int(re.findall(r'\d+mg',x)[0].replace('mg','').strip()) * int(re.findall(r'\d+ srg',x)[0].replace('srg','').strip())
                elif '(Adalimumab' in x:
                    # raise ValueError
                    # if 'ml pen' in x:
                    #     x = x.replace('ml pen','ml 1 pen')
                    try: pen_count = int(re.findall(r'\d+ pen',x)[0].replace('pen','').strip())
                    except: pen_count = int(re.findall(r'\d+ srg',x)[0].replace('srg','').strip())
                    dose_val = int(re.findall(r'\d+mg',x)[0].replace('mg','').strip()) * pen_count
                    # raise ValueError
                else:
                    raise ValueError
                    dose_val = x.split('(Infliximab')[-1].split('(Adalimumab')[-1].split('(Ustekinumab')[-1].split('▣')[-1].split('srg')[0].split('via')[0].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip()
            dose_series.append(dose_val)
    # raise ValueError
    dose_df['DOSE'] = dose_series
    # dose_df['DOSE'] = dose_df['DOSE'].map({'1 pen': 40, '2 pen': 80, '2 pen': 160})
    # (1) [원내] Remsima 100mg inj (Infliximab Korea) ...
    dose_df['ACTING'] = dose_df['Acting']

    dose_df['PERIOD'] = dose_df['처방지시'].map(lambda x: 'x1' if " [SC] " not in x else x.split(' [SC] ')[-1].split(':')[0].strip())
    # dose_df.to_csv(f"{output_dir}/dose_df_lhj.csv", encoding='utf-8-sig', index=False)
    # dose_df.loc[3203,'처방지시']

    dose_result_df.append(dose_df[result_cols].drop_duplicates(no_dup_cols))

    # drug_order_set = drug_order_set.union(set(dose_df['처방지시'].map(lambda x:''.join(x.split(':')[0].replace('  ',' ').split(') ')[1:]).replace('[원내]','').replace('[D/C]','').replace('[보류]','').replace('[반납]','').replace('[Em] ','').strip()).drop_duplicates()))


dose_result_df = pd.concat(dose_result_df, ignore_index=True).sort_values(['ID','DATETIME'], ascending=[True,False])


## DATETIME 추가 조정
dt_series = list()
for inx, row in dose_result_df.iterrows():
    if 'Y' in row['ACTING']:
        ymdt_patterns = re.findall(r"[\d][\d][\d][\d]-[\d][\d]-[\d][\d] [\d][\d]:[\d][\d]",row['ACTING'])
        if len(ymdt_patterns)>0:
            new_dt = ymdt_patterns[0].replace(' ','T')
        else:
            time_patterns = re.findall(r"[\d][\d]:[\d][\d]",row['ACTING'])
            new_dt = row['DATETIME'].split('T')[0] + 'T' + time_patterns[0]
        dt_series.append(new_dt)
    else:
        dt_series.append(row['DATETIME'])

dose_result_df['DATETIME'] = dt_series

## PERIOD 추가 조정

dose_result_df['PERIOD'] = dose_result_df['PERIOD'].map(lambda x:x.replace('2주 간격 유지용량 ','').replace('#1 ','').replace('3/20 ','').replace('매주 ','').replace(' X1 Day','').replace('1회 ','').replace('/wk ','/1wks').replace('x1 ','1'))


## H, C, Z는 제거

dose_result_df = dose_result_df[dose_result_df['ACTING'].map(lambda x: False if (('H' in x) or ('Z' in x) or ('C' in x)) else True)].reset_index(drop=True)
dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)
dose_result_df['ETC_INFO_TREATED'] = ''

## 처방비고 반영 및 기타 정리