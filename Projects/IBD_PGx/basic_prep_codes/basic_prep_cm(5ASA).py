from tools import *
from pynca.tools import *

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Orders
drug_order_set = set()

order_files = glob.glob(f'{resource_dir}/order/IBD_PGx_order(*).xlsx')
# len(order_files)
dose_result_df = list()
wierd_result_df = list()
result_cols = ['ID','NAME','DATETIME','DRUG','DOSE','ROUTE','ACTING','PERIOD','PLACE']
no_dup_cols = [c for c in result_cols if c!='NAME']
for finx, fpath in enumerate(order_files): #break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    # if pname=='이학준':
    #     raise ValueError

    # if pid in ('32411481',):
    #     pass
    # else:
    #     continue

    # if pid in ("15322168", "19739357", "34835292", "37366865", "21618097", "36898756", "36975211", "37858047"):       # lab, order 파일 다시 수집 필요
    #     continue

    fdf = pd.read_excel(fpath)

    # fdf.columns
    # fdf.to_csv(f"{outcome_dir}/error_dose_df.csv", encoding='utf-8-sig', index=False)

    # 5-ASA: mesalamine, basalazide, sulfasalazine
    dose_df = fdf[fdf['처방지시'].map(lambda x: (('mesalazine' in x.lower()) or ('mesalamine' in x.lower()) or ('basalazide' in x.lower()) or ('sulfasalazine' in x.lower())) and ('quantification' not in x.lower()))].copy()
    # dose_df.iloc[0]['처방지시']
    # dose_df.iloc[1]['처방지시']
    # dose_df.iloc[2]
    # dose_df['처방지시'].unique()

    # dose_df['처방지시'].map(lambda x: '[D/C]' in x)
    # dose_df['처방지시'].map(lambda x: '[D/C]' in x)

    # if len(dose_df)>0:
    #     raise ValueError
    # else:
    #     continue


    dose_df = dose_df[(~dose_df['Acting'].isna())].copy()

    # dose_df['처방지시'] = dose_df['처방지시'].map(lambda x:x.replace('임상 SUITE_램시마 주 100mg','Remsima 100mg inj (infliximab) ▣'))
    # dose_df['처방지시비고'] = dose_df['처방지시'].map(lambda x:x.split(' : ')[-1] if len(x.split(' : '))>1 else '')

    dose_df['ID'] = pid
    dose_df['NAME'] = pname

    # if pname=='김옥순': raise ValueError

    #### 약국_검사가 NA인 것은 제외하고 진행함 (추후 필요시 추가)

    # dose_df[dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test.xlsx",index=False)
    # dose_df[~dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test2.xlsx",index=False)
    dose_df = dose_df[(~dose_df['약국_검사'].isna())].copy()
    if len(dose_df)==0:
        print(f"({finx}) {pname} / {pid} / No administration records (5-ASA)")
        continue
    else:
        print(f"({finx}) {pname} / {pid}")

    #### 성분명이 IBD Biologics 인 경우만으로 필터링 (Infliximab, Adalimumab)

    regex_pattern = r'\(mesalamine|\(mesalazine|\(basalazide|\(sulfasalazine'
    dose_df = dose_df[dose_df['처방지시'].map(lambda x: bool(re.search(regex_pattern, x, flags=re.IGNORECASE)))].copy()

    dose_df['DT1'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[0].split('[')[-1].replace(' ','T'))
    dose_df['DT2'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[-1].split('[')[-1].replace(' ','T')[:-3])
    dose_df['DATETIME'] = dose_df[['DT1','DT2']].min(axis=1)
    # dose_df['DATE'] = dose_df['처방지시비고'].copy()
    dose_df['DRUG'] = dose_df['처방지시'].map(lambda x: re.search(regex_pattern, x, flags=re.IGNORECASE).group().lower().replace('(','').replace(')',''))
    dose_df['ROUTE'] = dose_df['처방지시'].map(lambda x: x.split('[')[-1].split(']')[0] if (" [P.O] " in x) or ("[per Rectal]" in x) or ("[Apply]" in x) else 'ETC')
    # dose_df['처방지시'].iloc[0]
    dose_df['PLACE'] = dose_df['주사시행처']
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+',x.split('▣')[-1].split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split(' ')[-1].strip())[0] if " [SC] " not in x else x.split('(Infliximab)')[-1].split('(Adalimumab)')[-1].split('▣')[-1].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+', x.split('mg')[1].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0] if (" [SC] " not in x) else x.split('(Infliximab')[-1].split('(Adalimumab')[-1].split('(Ustekinumab')[-1].split('▣')[-1].split('srg')[0].split('via')[0].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # mg_inx_dict = {'ustekinumab':0,'infliximab':1,'adalimumab':1}
    dose_series = list()
    for dose_inx, dose_row in dose_df.iterrows():
        x = dose_row['처방지시'].replace('1g(','1000mg(').replace('4g/100ml','4000mg(').replace('2g','2000mg')
        if (" [P.O] " in x) or ("[per Rectal]" in x) or ("[Apply]" in x):
            # raise ValueError
            try: dose_val = float(re.findall(r'\d+[\.]?\d*mg',x)[0].replace('mg','').strip())
            except:
                raise ValueError
                dose_val = float(re.findall(r'\d+[\.]?\d*mg', x)[0].replace('mg', '').strip()) * float(re.findall(r'\d+[\.]?\d* via', x)[0].replace('via', '').strip())
            # dose_val = re.findall(r'\d+',x.split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split('Stelara')[-1].split(' ')[-1].strip())[0]
            # raise ValueError
            dose_series.append(dose_val)
        else:
            print('PO 아님')
            raise ValueError
        #     if ('infliximab' in x.lower()) and (not ('Remsima 120mg pen (Infliximab) ▣ 120mg [SC]' in x)) and (not ('Remsima 120mg pen (Infliximab) 120mg [SC]' in x)):
        #         raise ValueError
        #     else:continue
        # continue
        # pid
        #     try: dose_val = float(re.findall(r'\d+[\.]?\d*mg',x)[1].replace('mg','').strip())
        #     except:
        #         if '(Ustekinumab' in x:
        #             # raise ValueError
        #             dose_val = float(re.findall(r'\d+[\.]?\d*mg',x)[0].replace('mg','').strip()) * float(re.findall(r'\d+[\.]?\d* srg',x)[0].replace('srg','').strip())
        #         elif '(Adalimumab' in x:
        #             # raise ValueError
        #             # if 'ml pen' in x:
        #             #     x = x.replace('ml pen','ml 1 pen')
        #             try: pen_count = float(re.findall(r'\d+[\.]?\d* pen',x)[0].replace('pen','').strip())
        #             except: pen_count = float(re.findall(r'\d+[\.]?\d* srg',x)[0].replace('srg','').strip())
        #             dose_val = float(re.findall(r'\d+[\.]?\d*mg',x)[0].replace('mg','').strip()) * pen_count
        #             # raise ValueError
        #         else:
        #             # raise ValueError
        #             dose_val = x.split('(Infliximab')[-1].split('(Adalimumab')[-1].split('(Ustekinumab')[-1].split('▣')[-1].split('srg')[0].split('via')[0].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip()
        #     dose_series.append(dose_val)
    # raise ValueError
    dose_df['DOSE'] = dose_series
    # dose_df['DOSE'] = dose_df['DOSE'].map({'1 pen': 40, '2 pen': 80, '2 pen': 160})
    # (1) [원내] Remsima 100mg inj (Infliximab Korea) ...
    dose_df['ACTING'] = dose_df['Acting']

    dose_df['PERIOD'] = dose_df['처방지시'].map(lambda x: re.findall(r'X\d+ Days',x)[0] if bool(re.search(r'X\d+ Days', x)) else 'X1 Days')
    # dose_df['처방지시'].unique()
    # dose_df.to_csv(f"{output_dir}/dose_df_lhj.csv", encoding='utf-8-sig', index=False)
    # dose_df.loc[3203,'처방지시']

    ## ADALIMUMAB DOSE 데이터 Correction
    # if pid=='37590846':
    #     add_rows = pd.DataFrame([dose_df.sort_values(['ID','DATETIME']).iloc[0],dose_df.sort_values(['ID','DATETIME']).iloc[0]])
    #     add_rows['DATETIME'] = ['2023-11-23T16:22','2023-12-06T16:22']
    #     dose_df = pd.concat([dose_df, add_rows])
    # if pid=='26675590':
    #     add_rows = pd.DataFrame([dose_df.sort_values(['ID','DATETIME']).iloc[0],dose_df.sort_values(['ID','DATETIME']).iloc[0]])
    #     add_rows['DATETIME'] = ['2023-12-14T16:03','2023-12-28T16:03']
    #     dose_df = pd.concat([dose_df, add_rows])
    # if pid=='27364521':
    #     dose_df = dose_df.sort_values(['ID', 'DATETIME'],ignore_index=True)
    #     dose_df.at[0,'DOSE']=160
    #     dose_df = dose_df[dose_df.index != 1].copy()
    # if pid=='32411481':
    #     raise ValueError
        # dose_df.iloc[20]
        # dose_df.iloc[21]
        # dose_df = dose_df.sort_values(['ID', 'DATETIME'],ignore_index=True)
        # dose_df.at[0,'DOSE']=160
        # dose_df = dose_df[dose_df.index != 1].copy()



        # dose_df[['DOSE','DATETIME']]
        # dose_df = pd.concat([dose_df, add_rows])


        # dose_df[['DATETIME','ETC_INFO']]
        # raise ValueError


    dose_result_df.append(dose_df[result_cols].drop_duplicates(no_dup_cols))

    # drug_order_set = drug_order_set.union(set(dose_df['처방지시'].map(lambda x:''.join(x.split(':')[0].replace('  ',' ').split(') ')[1:]).replace('[원내]','').replace('[D/C]','').replace('[보류]','').replace('[반납]','').replace('[Em] ','').strip()).drop_duplicates()))

# dose_result_df = dose_result_df[0]
dose_result_df = pd.concat(dose_result_df, ignore_index=True).sort_values(['ID','DATETIME'], ascending=[True,False])
# dose_result_df['DATE'] = dose_result_df['DATETIME']

# dose_result_df['DRUG'].unique()
# dose_result_df.iloc[0]
# dose_result_df[dose_result_df['ID'].isin(('37590846',))]

## DATETIME 추가 조정
# dt_series = list()
# for inx, row in dose_result_df.iterrows():
#     if 'Y' in row['ACTING']:
#         ymdt_patterns = re.findall(r"[\d][\d][\d][\d]-[\d][\d]-[\d][\d] [\d][\d]:[\d][\d]",row['ACTING'])
#         if len(ymdt_patterns)>0:
#             new_dt = ymdt_patterns[0].replace(' ','T')
#         else:
#             time_patterns = re.findall(r"[\d][\d]:[\d][\d]",row['ACTING'])
#             new_dt = row['DATETIME'].split('T')[0] + 'T' + time_patterns[0]
#         dt_series.append(new_dt)
#     else:
#         dt_series.append(row['DATETIME'])
#
# dose_result_df['DATETIME'] = dt_series

## PERIOD 추가 조정

dose_result_df['PERIOD'] = dose_result_df['PERIOD'].map(lambda x:x.replace('2주 간격 유지용량 ','').replace('#1 ','').replace('3/20 ','').replace('매주 ','').replace(' X1 Day','').replace('1회 ','').replace('/wk ','/1wks').replace('x1 ','1'))

## H, C, Z는 제거
# dose_result_df['ACTING'].unique()
dose_result_df = dose_result_df[dose_result_df['ACTING'].map(lambda x: False if (('H' in x) or ('Z' in x) or ('C' in x) or ('N' in x.replace('PEN',''))) else True)].reset_index(drop=True)
dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,False], ignore_index=True)
dose_result_df['ETC_INFO_TREATED'] = ''

# ADDL, II 추가 및 다음 DOSE 보다 넘어가는 ADDL 확인

dose_result_df['ADDL'] = dose_result_df['PERIOD'].map(lambda x:int(x.split('X')[-1].split(' ')[0].strip()))
dose_result_df['II'] = 24
dose_result_df['DATETIME'] = dose_result_df['DATETIME'].replace('2016-09-07T17:','2016-09-07T17:00')

# ADDL을 하나씩 다 넣어보자
addl_dose_result_df = dose_result_df.sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True).copy()
addl_added_df = list()
for inx, row in addl_dose_result_df[addl_dose_result_df['ADDL']>0].iterrows():
    # dose_result_df[dose_result_df.index==inx]
    # if row['DATETIME']=='2021-07-19T12:05':
        # raise ValueError

    max_addl = row['ADDL']
    init_datetime = row['DATETIME']
    addl_frag_df = list()
    for addl_num in range(0,max_addl):
        row['ETC_INFO_TREATED'] = f'ADDL반영_{addl_num}'
        row['PERIOD'] = 'x1'
        row['ADDL']=0
        row['DATETIME'] = (datetime.strptime(init_datetime,'%Y-%m-%dT%H:%M')+timedelta(addl_num*row['II']/24)).strftime('%Y-%m-%dT%H:%M')
        if addl_num==0:
            continue
        else:
            addl_frag_df.append(pd.DataFrame([row]))

    if len(addl_frag_df)!=0:
        print(f'{row["ID"]} / {row["DATETIME"]} / ADDL ({max_addl})')
        addl_added_df.append(pd.concat(addl_frag_df))
    else:
        print(f'{row["ID"]} / {row["DATETIME"]} / No ADDL')

    addl_dose_result_df.at[inx,'ADDL']=0
    addl_dose_result_df.at[inx,'ETC_INFO_TREATED'] = f"ADDL존재_{max_addl}"

addl_added_df = pd.concat(addl_added_df)


# 정리 및 저장

# addl_dose_result_df['ROUTE'].unique()
# addl_dose_result_df[addl_dose_result_df['ROUTE']=='per Rectal']

addl_dose_result_df = pd.concat([addl_dose_result_df, addl_added_df]).sort_values(['ID','DATETIME'], ascending=[True,True], ignore_index=True)
addl_dose_result_df['DATETIME'] = addl_dose_result_df['DATETIME'].map(lambda x:x.split('T')[0])
addl_dose_result_df = addl_dose_result_df.drop(['ADDL','ETC_INFO_TREATED','II'],axis=1).rename(columns={'DATETIME':'DATE'})
addl_dose_result_df.to_csv(f"{output_dir}/cm_df(5ASA).csv", encoding='utf-8-sig', index=False)

