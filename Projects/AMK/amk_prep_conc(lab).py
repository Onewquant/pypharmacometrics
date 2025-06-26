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
lab_set = set()

conc_files = glob.glob(f'{resource_dir}/lab/{prj_name}_lab(*).xlsx')
conc_result_df = list()
for finx, fpath in enumerate(conc_files): #break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    print(f"({finx}) {pname} / {pid}")

    # if pid in ("15322168", "19739357", "34835292", "37366865", "21618097", "36898756", "36975211", "37858047"):       # lab, order 파일 다시 수집 필요
    #     continue

    fdf = pd.read_excel(fpath)
    # lab_set = lab_set.union(set(fdf[fdf['검사명'].map(lambda x:'amikacin' in x.lower())]['검사명'].unique()))

    # fdf['검사명'].unique()
    # break

    # fdf.columns

    # cumlab_set = cumlab_set.union(set(fdf['Lab'].unique()))
    # [c for c in cumlab_set if ('infliximab' in c.lower()) or ('adalimumab' in c.lower())]
    # lab_set = lab_set.union(set(fdf['검사명'].unique()))
    # [c for c in lab_set if ('infliximab' in c.lower()) or ('adalimumab' in c.lower())]

    # 'Infliximab 정량: 재검한 결과입니다.'
    # 'Infliximab Quantification'
    # 'Adalimumab Quantification'

    # conc_df.columns
    conc_df = fdf[fdf['검사명'].isin(['Amikacin'])].copy()     # Amikacin 농도 측정 항목만 추출
    """
    try:
        conc_df['검사결과'].astype(float)
    except:
        print(f"({finx}) {pname} / {pid} / float 아님 추가")
        conc_result_df.append(conc_df[['검사결과']])
        continue        
    print(f"({finx}) {pname} / {pid}")
    continue
    conc_result_df = pd.concat(conc_result_df, ignore_index=True)
    
    strconc_inx_list = list()
    for inx, res in enumerate(conc_result_df['검사결과']): #break
        try: float(res.replace('<','').replace('>',''))
        except:strconc_inx_list.append(inx)

    strconc_df = conc_result_df.iloc[strconc_inx_list].copy()
    for inx, row in strconc_df.reset_index(drop=True).iterrows():
        print(f"({inx}) {row['검사결과']}")

    strconc_df.to_csv(f"{output_dir}/str_conc_df(lab).csv", encoding='utf-8-sig', index=False)
    """


    # fdf = fdf[~fdf['Acting'].isna()].copy()
    # dose_df = fdf[fdf['처방지시'].map(lambda x: (('adalimumab' in x.lower()) or ('infliximab' in x.lower())) and ('quantification' not in x.lower()))].copy()
    # dose_df['처방지시비고'] = dose_df['처방지시'].map(lambda x:x.split(' : ')[-1] if len(x.split(' : '))>1 else '')
    #
    #
    conc_df['ID'] = pid
    conc_df['NAME'] = pname
    conc_df['DRUG'] = conc_df['검사명'].map(lambda x:x.split(' ')[0].lower())
    conc_df['검사결과'] = conc_df['검사결과'].replace('중복 오더임',np.nan)
    conc_df['CONC'] = conc_df['검사결과'].map(lambda x:float(x.split('(')[0].replace('<','').replace('>','').strip()) if str(x).lower()!='nan' else np.nan)
    conc_df['ETC_INFO'] = conc_df['검사결과'].map(lambda x:x.split('(')[-1].strip() if str(x).lower() != 'nan' else np.nan)


    conc_df = conc_df[~conc_df['CONC'].isna()].copy()
    # conc_df[] = conc_df[['보고일','오더일']].max(axis=1)

    conc_result_df.append(conc_df[['ID','NAME','보고일','오더일','DRUG','CONC', 'ETC_INFO']])
    # conc_result_df['ID'].drop_duplicates()
    # # if pname=='김옥순': raise ValueError
    #
    # #### 약국_검사가 NA인 것은 제외하고 진행함 (추후 필요시 추가)
    #
    # # dose_df[dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test.xlsx",index=False)
    # # dose_df[~dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test2.xlsx",index=False)
    # dose_df = dose_df[~dose_df['약국_검사'].isna()].copy()
    #
    # #### 성분명이 IBD Biologics 인 경우만으로 필터링 (Infliximab, Adalimumab)
    #
    # regex_pattern = r'\(infliximab|adalimumab\)'
    # dose_df = dose_df[dose_df['처방지시'].map(lambda x: bool(re.search(regex_pattern, x, flags=re.IGNORECASE)))].copy()
    #
    # dose_df['DT1'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[0].split('[')[-1].replace(' ','T'))
    # dose_df['DT2'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[-1].split('[')[-1].replace(' ','T')[:-3])
    # dose_df['DATETIME'] = dose_df[['DT1','DT2']].min(axis=1)
    # dose_df['ETC_INFO'] = dose_df['처방지시비고'].copy()
    # dose_df['DRUG'] = dose_df['처방지시'].map(lambda x: re.search(regex_pattern, x, flags=re.IGNORECASE).group().lower().replace('(','').replace(')',''))
    # dose_df['ROUTE'] = dose_df['처방지시'].map(lambda x: 'IV' if " [SC] " not in x else 'SC')
    # dose_df['DOSE'] = dose_df['처방지시'].map(lambda x: re.findall(r'\d+',x.split('▣')[-1].split('mg')[0].split('Remsima')[-1].split('Humira')[-1].split(' ')[-1].strip())[0] if " [SC] " not in x else x.split('(Infliximab)')[-1].split('(Adalimumab)')[-1].split('▣')[-1].split('mg')[0].split(':')[0].split(' [SC] ')[0].strip())
    # dose_df['PERIOD'] = dose_df['처방지시'].map(lambda x: 'x1' if " [SC] " not in x else x.split(' [SC] ')[-1].split(':')[0].strip())
    #
    # # dose_df.loc[3203,'처방지시']
    # dose_result_df.append(dose_df[['ID','NAME','DATETIME','DRUG','DOSE','ROUTE','PERIOD','ETC_INFO']].copy())

    # drug_order_set = drug_order_set.union(set(dose_df['처방지시'].map(lambda x:''.join(x.split(':')[0].replace('  ',' ').split(') ')[1:]).replace('[원내]','').replace('[D/C]','').replace('[보류]','').replace('[반납]','').replace('[Em] ','').strip()).drop_duplicates()))

conc_result_df = pd.concat(conc_result_df, ignore_index=True).drop_duplicates(['ID','보고일','오더일','DRUG','CONC'])
conc_result_df.to_csv(f"{output_dir}/conc_df(lab).csv", encoding='utf-8-sig', index=False)

# ot_list = list()
# for inx_ot, order_text in enumerate(drug_order_set):
#     if '[SC]' not in order_text:
#         ot_list.append(order_text)
#         continue
#
# for inx_ot, order_text in enumerate(ot_list):
#     print(f"({inx_ot}) {order_text}")
# len(order_text)

# ['비고','처방지시', '발행처', '발행의', '수납', '약국/검사', '주사시행처', 'Acting', '변경의']
#
# len(order_files)
#
#
# df = pd.read_csv(f'{resource_dir}/glpharma_CONC.csv')
# seq_df = pd.read_csv(f'{prj_dir}/glpharma_SEQUENCE.csv')