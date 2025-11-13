from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


# pid_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/CCM.csv")
# pid_df = pid_df.rename(columns={'Deidentification_ID':'UID','환자번호':'REQ_ID'}).drop(['Column1', 'Column2'],axis=1)

pid_decode_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/재식별 파일.csv")
pid_decode_df = pid_decode_df.rename(columns={'환자번호':'UID','Deidentification_ID':'PID'})
pid_decode_df['UID'] = pid_decode_df['UID'].map(lambda x: x.split('-')[0])

mic_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/MICROBIOLOGY.csv", encoding='euc-kr')
# mic_df.columns
# mic_df['NAME'].unique()

mic_df = mic_df.rename(columns={'ID':'UID'})
mic_df = mic_df[~mic_df['RESULT'].isna()].copy()
# mic_df['UID'] = mic_df['UID'].map(lambda x: x.split('-')[0])
# for c in mic_df['UID']:

mic_df = mic_df[~(mic_df['CODE'].isin(['L41071','L4107','L7347','L2585','L26002','L731002','L2584']))].copy() # 결핵균 약제내성 검사 제외
mic_df = mic_df[~(mic_df['CODE'].isin(['L2588','L25137','L2537','L25371','L25361','L25373','L4067','L43242']))].copy() # pharmacogenomics 및 기타 약제검사 제외

afb_df = mic_df[(mic_df['CODE'].isin(['L4106']))].copy() # 항산균 검사
tb_df = mic_df[(mic_df['CODE'].isin(['L2510']))].copy() # tb 검사
ntm_df = mic_df[(mic_df['CODE'].isin(['L2518','L2508']))].copy() # ntm 검사

respvirus_df = mic_df[(mic_df['CODE'].isin(['L25159','L25158','L25198','L25199','L25104','L2513','L25140']))].copy() # 호흡기바이러스/폐렴 검사
parainfluenza_df = mic_df[(mic_df['CODE'].isin(['L25132']))].copy() # parainfluenza 검사
diar_df = mic_df[(mic_df['CODE'].isin(['L25116','L25117','L25118','L25153','L25154']))].copy() # Acute diarrhea 검사
legionella_df = mic_df[(mic_df['CODE'].isin(['L2507']))].copy() # legionella pcr 검사
fungi_df = mic_df[(mic_df['CODE'].isin(['L25125','L2566']))].copy() # fungi 검사
vre_df = mic_df[(mic_df['CODE'].isin(['L1704','L1703','L1702','L1701','L1705','L1706']))].copy() # vre 검사
mening_df = mic_df[(mic_df['CODE'].isin(['L25152']))].copy() # vre 검사

# AFB 결과
afb_df['RESULT'] = afb_df['RESULT'].map(lambda x: 'Mycobacterium tuberculosis' if ('동정결과]\nM. tuberculosis complex isolated' in x) else x)
afb_df['RESULT'] = afb_df['RESULT'].map(lambda x: 'NonTuberculous Mycobacteria (NTM)' if ('동정결과]\nNonTuberculous Mycobacteria (NTM) isolated' in x) else x)
afb_df['RESULT'] = afb_df['RESULT'].map(lambda x: 'Mycobacterium tuberculosis' if ('동정결과 : Mycobacterium tuberculosis isolated' in x) else x)
afb_df['RESULT'] = afb_df['RESULT'].map(lambda x: 'NonTuberculous Mycobacteria (NTM)' if ('동정결과 : NonTuberculous Mycobacteria (NTM) isolated' in x) else x)

afb_df['RESULT'] = afb_df['RESULT'].map(lambda x:x.split('분리 균명')[-1].split('동정결과]')[-1].split('동정결과 :')[-1].split('결핵결과 :')[-1].split('* 두 가지 종류의 배지(Solid and Liquid)')[0].strip()) # .split('\n\n')[0]
afb_df['RESULT'] = afb_df['RESULT'].map(lambda x: 'Negative' if ('QNS(검체량 부족)' in x) or ('QNS(검체량부족)' in x) or ('QNS (검체량부족)' in x) or ('QNS (검체량 부족)' in x) or ('Q.N.S.(검체량 부족)' in x) or ('No Acid-Fast bacilli are observed' in x) or ('No acid fast bacilli isolated' in x) or ('Contamination' in x) else x)

afb_df['RESULT'] = afb_df['RESULT'].map(lambda x: 'Mycobacterium tuberculosis' if ('Mycobacterium tuberculosis complex isolated' in x) or ('M. tuberculosis complex isolated' in x) or ('Mycobacterium tuberculosis isolated' in x) else x)
afb_df['RESULT'] = afb_df['RESULT'].map(lambda x: 'NonTuberculous Mycobacteria (NTM)' if ('NonTuberculous Mycobacteria (NTM) isolated' in x) else x)
afb_df['RESULT'] = afb_df['RESULT'].map(lambda x: 'MOTT (Mycobacterium other than tuberculosis)' if ('Mycobacterium other than tuberculosis isolated' in x) else x)
afb_df['RESULT'] = afb_df['RESULT'].map({'Negative':'Negative','42 Days':'Negative','57 Days':'Negative','most closely related to Nocardia 43 Days':'Nocardia','NonTuberculous Mycobacteria (NTM)':'NonTuberculous Mycobacteria (NTM)','Mycobacterium tuberculosis':'Mycobacterium tuberculosis','Mycobacterium other than tuberculosis':'Mycobacterium other than tuberculosis'})

# TB 검사
tb_df['RESULT'] = tb_df['RESULT'].map(lambda x:x.split('판독결과:')[-1].split('[검체]')[0].split('------------------------------------------------------------')[0].strip()) # .split('\n\n')[0]
tb_df['RESULT'] = tb_df['RESULT'].map(lambda x: 'Negative' if ('Negative for M.tuberculosis complex' in x) or ('Negative for M.tuberculosis' in x) else x)
tb_df['RESULT'] = tb_df['RESULT'].map(lambda x: 'Mycobacterium tuberculosis' if ('Positive for M.tuberculosis complex' in x) or ('Positive for M.tuberculosis' in x) else x)

# NTM 검사
ntm_df['RESULT'] = ntm_df['RESULT'].map(lambda x:x.split('판독결과 :')[-1].split('판독결과:')[-1].split('[검체]')[0].split('------------------------------------------------------------')[0].strip()) # .split('\n\n')[0]
ntm_df['RESULT'] = ntm_df['RESULT'].map(lambda x: 'Negative' if ('Negative for ' in x) or ('No amplification' in x) else x.replace('Positive for ','').split(' + '))
ntm_df = ntm_df.explode('RESULT')
ntm_df['RESULT'] = ntm_df['RESULT'].map(lambda x:x.split('(formally')[0].strip())
for ntm_str in ['Mycobacterium avium','MAC[Mycobacterium avium complex]','MAC(Mycobacterium avium complex)']:
    ntm_df['RESULT'] = ntm_df['RESULT'].replace(ntm_str,'Mycobacterium avium complex')
for ntm_str in ['Mycobacterium abscessus sensu stricto','Mycobacterium abscessus complex','Mycobacterium abscessus']:
    ntm_df['RESULT'] = ntm_df['RESULT'].replace(ntm_str,'Mycobacterium abscessus complex')
for ntm_str in ['NonTuberculous Mycobacteria, Unidentified','NTM (NonTuberculous Mycobacteria)']:
    ntm_df['RESULT'] = ntm_df['RESULT'].replace(ntm_str,'NonTuberculous Mycobacteria (NTM)')
ntm_df['RESULT'] = ntm_df['RESULT'].replace('most closely related to Nocardia','Nocardia')

# Respiratory virus
respvirus_df['RESULT'] = respvirus_df['RESULT'].map(lambda x:x.split('[판독결과]')[-1].split('판독결과:')[-1].split('[검사방법]')[0].strip().split('\n')) # .split('\n\n')[0]
respvirus_df = respvirus_df.explode('RESULT')
respvirus_df['RESULT'] = respvirus_df['RESULT'].map(lambda x: 'Negative' if ('Negative for ' in x) else x.replace('Positive for ',''))


# Parainfluenza 검사
parainfluenza_df['RESULT'] = parainfluenza_df['RESULT'].map(lambda x:x.split('[판독결과]')[-1].split('판독결과:')[-1].split('[검사방법]')[0].strip().split('\n')) # .split('\n\n')[0]
parainfluenza_df = parainfluenza_df.explode('RESULT')
parainfluenza_df['RESULT'] = parainfluenza_df['RESULT'].map(lambda x: 'Negative' if ('Negative for ' in x) else x.replace('Positive for ',''))

# Diarrhea 검사
diar_df['RESULT'] = diar_df['RESULT'].map(lambda x:x.split('[판독결과]')[-1].split('판독결과:')[-1].split('[검사방법]')[0].strip().split('\n')) # .split('\n\n')[0]
diar_df = diar_df.explode('RESULT')
diar_df['RESULT'] = diar_df['RESULT'].map(lambda x: 'Negative' if ('Negative for ' in x) else x.replace('Positive for ',''))

# Legionella 검사
legionella_df['RESULT'] = legionella_df['RESULT'].map(lambda x:x.split('[판독결과]')[-1].split('판독결과:')[-1].split('[검사방법]')[0].strip().split('\n')) # .split('\n\n')[0]
legionella_df = legionella_df.explode('RESULT')
legionella_df['RESULT'] = legionella_df['RESULT'].map(lambda x: 'Negative' if ('Negative for ' in x) else x.replace('Positive for ',''))

# Fungi 검사
fungi_df['RESULT'] = fungi_df['RESULT'].map(lambda x:x.split('[판독결과]')[-1].split('판독결과:')[-1].split('[검사개요]')[0].split('[검체]')[0].strip().split('\n')) # .split('\n\n')[0]
fungi_df = fungi_df.explode('RESULT')
fungi_df['RESULT'] = fungi_df['RESULT'].map(lambda x: 'Negative' if ('Negative for ' in x) or ('Unable to identify due to probably mixed condition'==x) or ('No amplification'==x) or ('(검사진행을 하였으나 염기서열분석이 불가능함)'==x) or ('(검사진행을 하였으나 유전자 증폭이 되지 않아 염기서열분석이 불가능함)'==x) else x.replace('Positive for ',''))
fungi_df['RESULT'] = fungi_df['RESULT'].replace('Pseudallescheria spp','Pseudallescheria spp.')
# fungi_df
# fungi_df['RESULT'].iloc[0]

# VRE 검사
vre_df['RESULT'] = vre_df['RESULT'].map(lambda x:x.split('<수정 후 결과>')[-1].split('동정결과:')[-1].split('COMMENT')[0].split('항생제 감수성결과')[0].split('본 검사실은 대한진단검사의학회/진단검사의학재단의')[0].split('--------------------------------------------------')[0].strip()) # .split('\n\n')[0]
for vre_str in ['No Vancomycin-resistance Enterococci(VRE) isolated','No Vancomycin-resistant Enterococci (VRE) isolated']:
    vre_df['RESULT'] = vre_df['RESULT'].replace(vre_str,'Negative')

# Meningitis 검사
# mening_df['RESULT'] = mening_df['RESULT'].map(lambda x:x.split('<수정 후 결과>')[-1].split('동정결과:')[-1].split('COMMENT')[0].split('항생제 감수성결과')[0].split('본 검사실은 대한진단검사의학회/진단검사의학재단의')[0].split('--------------------------------------------------')[0].strip()) # .split('\n\n')[0]
mening_df['RESULT'] = mening_df['RESULT'].map(lambda x:x.split('[판독결과]')[-1].split('판독결과 :')[-1].split('판독결과:')[-1].split('[검체]')[0].strip().split('\n')) # .split('\n\n')[0]
mening_df = mening_df.explode('RESULT')
mening_df['RESULT'] = mening_df['RESULT'].map(lambda x: 'Negative' if ('Negative for ' in x) else x.replace('Positive for ',''))
mening_df['RESULT'] = mening_df['RESULT'].map(lambda x: 'Negative' if (': Negative' in x) else x.replace(': Positive',''))

for inx, x in enumerate(mening_df['RESULT'].unique()):
    print(f"({inx}) {x}")
    # break


other_mic_df = mic_df[~(mic_df['CODE'].isin(['L4106','L2518','L2508','L25125','L2566','L25159','L25158','L25198','L25199','L25116','L25117','L2507','L25104','L2513','L1704','L25118','L2513','L1702','L1703','L25153','L25152','L1701','L1705','L25154','L25132','L2510','L25140','L1706']))].copy()
other_mic_df['RESULT'] = other_mic_df['RESULT'].map(lambda x:x.split('분리 균명')[-1].split('항생제 감수성결과')[0].split('동정결과')[-1].split('동정 결과')[-1].split('(최종보고')[0].split('판독결과')[-1].split('\n\n[검')[0].replace(':','').replace('균주명  ','').replace(']','').split('중간보고')[0].split('[검체')[0].split('약제명        절대농도(㎍/ml)      판정결과')[0].split('[검사방법')[0].split('Comment 본 검사는 CRE surveillance를 목적으로 시행하였음')[0].split('항생제 감수성 검사결과')[0].split('* 항생제 감수성 결과가 있는 경우')[0].split('카바페넴내성 장내세균 확인검사')[0].split('Comment')[0].split('COMMENT')[0].split('----------------------------------------------')[0].strip()) # .split('\n\n')[0]
other_mic_df['RESULT'] = other_mic_df.apply(lambda x: x['NAME'].split('[분자진단]')[0] if 'IU/mL' in x['RESULT'] else x['RESULT'], axis=1)
other_mic_df['RESULT'] = other_mic_df['RESULT'].map(lambda x: 'Negative' if ('Negative for' in x) or ('Unable to identify due to probably mixed condition' in x) or ('No amplification' in x) or ('possible Normal flora' in x)else x)

## Negative 인 문구 define
neg_factor_list = list()
# for c in ['Rare(<1)','A few(1~5)']:
#     for g in ['G(-) rod','G(+) rod','Yeast like cell']:
for c in ['Rare(<1)','Afew(1~5)']:
    for g in ['G(-)rod','G(+)rod','G(-)cocci','G(+)cocci','G(-)cocci,chains','G(+)cocci,chains','Yeastlikecell']:
        cg_str = f"{c}{g}"
        neg_factor_list.append(cg_str)

neg_list = set()
for a in neg_factor_list:
    for b in neg_factor_list:
        if a==b:
            continue
        else:
            # neg_list.add(f"{a}, {b}")
            neg_list.add(f"{a},{b}")

row_inx_list = list()
for inx, row in other_mic_df.iterrows():
    row_bool = False
    for neg_str in neg_list:
        # if neg_str in row['RESULT']:
        if neg_str in row['RESULT'].replace(' ',''):
            row_bool=True
            row_inx_list.append(inx)
            print(f"({len(row_inx_list)}) Negative // {other_mic_df.at[inx,'RESULT']}")
            other_mic_df.at[inx,'RESULT'] = 'Negative'
            break
    # if row_bool:

for inx, x in enumerate(other_mic_df['RESULT'].unique()):
    print(f"({inx}) {x}")


micbio_df = pd.concat([afb_df, tb_df, ntm_df, respvirus_df, parainfluenza_df, diar_df, legionella_df,fungi_df,vre_df, mening_df, other_mic_df])
if not os.path.exists(f'{output_dir}'):
    os.mkdir(f'{output_dir}')
micbio_df.to_csv(f"{output_dir}/final_microbiology_df.csv", encoding='utf-8-sig', index=False)

#################

gumun = '2월 11일 Blood 검체 접종  BAP plate에서 Acid fast bacilli 의심되는 col'
other_mic_df[other_mic_df['RESULT'].map(lambda x:gumun in x)]['CODE'].unique()
other_mic_df.columns

surgery_df['PROC_CATNUM'] = surgery_df['CODE'].astype(float).map(lambda x: 1 if x in (35.2, 35.22, 35.24, 35.28, 36.1, 36.11, 36.12, 37.12, 37.31, 37.33, 37.34, 37.5, 37.62, 37.66, 37.8) else 2)
surgery_df['PROC'] = surgery_df['NAME'].copy()
surgery_df = surgery_df[['UID','DATE','PROC_CATNUM','PROC']].copy()

# surgery_df.columns
icu_mortal_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/PROCEDURE_AND_OUTCOME/LOS,MORTALITY,ICU.csv")
icu_mortal_df = icu_mortal_df.rename(columns={'ID':'UID'})
icu_mortal_df['UID'] = icu_mortal_df['UID'].map(lambda x: x.split('-')[0])
# icu_mortal_df.columns
# icu_mortal_df[['UID','LOS','ICU','ADD_DATE','DIS_DATE']]
icu_mortal_df = icu_mortal_df[~icu_mortal_df['ICU'].isna()][['UID','LOS','ICU','ADD_DATE','DIS_DATE']].reset_index(drop=True)
icu_mortal_df['PROC_CATNUM'] = 3
icu_mortal_df['PROC'] = 'ICU'
# DIS_DATE 비어있는 데이터 채우기 (ADD_DATE + ICU 이용 days)
for inx, row in icu_mortal_df.iterrows():
    if type(row['DIS_DATE'])!=str:
        icu_mortal_df.at[inx, 'DIS_DATE'] = (datetime.strptime(row['ADD_DATE'],'%Y-%m-%d')+timedelta(days=int(row['LOS']))).strftime('%Y-%m-%d')

icu_mortal_df['DATE'] = icu_mortal_df.apply(lambda x:pd.date_range(x['ADD_DATE'],x['DIS_DATE']).strftime('%Y-%m-%d').tolist(), axis=1)
icu_mortal_df = icu_mortal_df[['UID','DATE','PROC_CATNUM', 'PROC']].explode('DATE').drop_duplicates(['UID','DATE']).sort_values(['UID','DATE'], ignore_index=True)

# icu_mortal_df[icu_mortal_df['ADD_DATE'].isna()]

# x = {'ADD_DATE':'2008-11-11','DIS_DATE':'2009-05-26'}

mv_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/PROCEDURE_AND_OUTCOME/MV.csv", encoding='euc-kr')
mv_df = mv_df.rename(columns={'ID':'UID'})
mv_df['UID'] = mv_df['UID'].map(lambda x: x.split('-')[0])
mv_df = mv_df[mv_df['MLCODE'] > 0][['UID','DATE']].copy()
mv_df['PROC_CATNUM'] = 4
mv_df['PROC'] = 'Mechanical Ventilation'
mv_df = mv_df[['UID','DATE','PROC_CATNUM', 'PROC']].copy()


proc_df = pd.concat([surgery_df, icu_mortal_df, mv_df])
proc_df = proc_df.drop_duplicates(['UID','DATE','PROC_CATNUM'])

proc_df['UID'] = proc_df['UID'].map(str)
proc_df = proc_df.merge(pid_decode_df, on=['UID'],how='left')
proc_df['UID'] = proc_df['PID'].copy()
proc_df = proc_df.drop(['PID'], axis=1)

proc_index_dict = {1:(365,1),2:(365,1),3:(365,1),4:(14,1)}
proc_df['PROC_PERIOD_FROM'] = proc_df['PROC_CATNUM'].map(lambda x:proc_index_dict[x][0])
proc_df['PROC_PERIOD_UNTIL'] = proc_df['PROC_CATNUM'].map(lambda x:proc_index_dict[x][1])


if not os.path.exists(f'{output_dir}'):
    os.mkdir(f'{output_dir}')
proc_df.to_csv(f"{output_dir}/final_procedure_df.csv", encoding='utf-8-sig', index=False)
