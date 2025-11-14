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
neg_degree_list = ['Rare(<1)','Afew(1~5)']
bact_type_list = ['G(-)rod','G(+)rod','G(-)cocci','G(+)cocci','G(-)cocci,chains','G(+)cocci,chains','Yeastlikecell']

for c in neg_degree_list:
    for g in bact_type_list:
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
    # Negative 구문이 2개 들어있는 경우
    for neg_str in neg_list:
        if neg_str in row['RESULT'].replace(' ',''):
            row_bool=True
            row_inx_list.append(inx)
            print(f"({len(row_inx_list)}) Negative // {other_mic_df.at[inx,'RESULT']}")
            other_mic_df.at[inx,'RESULT'] = 'Negative'
            break

    # Negative 구문이 1개 들어있고 나머지 균에 대한 구문은 없는 경우
    if not row_bool:
        for neg_factor_str in neg_factor_list:
            if (neg_factor_str in row['RESULT'].replace(' ', '')):
                no_other_bact = True
                for bactype in bact_type_list:
                    if (bactype not in neg_factor_str) and (bactype in row['RESULT']):
                        no_other_bact = False
                    # pass
                if no_other_bact:
                    # if other_mic_df.at[inx, 'RESULT'] not in ['Many(>10) WBC,A few(1~5) Yeast like cell']:
                    #     raise ValueError
                    row_bool = True
                    row_inx_list.append(inx)
                    print(f"({len(row_inx_list)}) Negative2 // {other_mic_df.at[inx, 'RESULT']}")
                    other_mic_df.at[inx, 'RESULT'] = 'Negative'
                    break
    # if row_bool:




micbio_df = pd.concat([afb_df, tb_df, ntm_df, respvirus_df, parainfluenza_df, diar_df, legionella_df,fungi_df,vre_df, mening_df, other_mic_df])
# if not os.path.exists(f'{output_dir}'):
#     os.mkdir(f'{output_dir}')
# micbio_df.to_csv(f"{output_dir}/final_microbiology_df.csv", encoding='utf-8-sig', index=False)



##############################

micbio_index_dict = {1:"Pseudomonas aeruginosa",
                     2:"Acinetobacter baumannii",
                     3:"Klebsiella pneumoniae",
                     4:"Escherichia coli",
                     5:"Enterobacter species",
                     6:"Serratia marcescens",
                     7:"Mycobacterium tuberculosis",
                     8:"NonTuberculous Mycobacterium",
                     9:"Others",
                     10:"Empirical"
                     }

# micbio_df = pd.read_csv(f"{output_dir}/final_microbiology_df.csv")
# for vre_str in ['No Vancomycin-resistance Enterococci(VRE) isolated','No Vancomycin-resistant Enterococci (VRE) isolated']:
#     vre_df['RESULT'] = vre_df['RESULT'].replace(vre_str,'Negative')

micbio_df = micbio_df[~micbio_df['RESULT'].isna()].copy()
micbio_df['RESULT'] = micbio_df['RESULT'].replace('No aerobic or anaerobic bacteria isolated','Negative')
micbio_df['RESULT'] = micbio_df['RESULT'].map(lambda x: 'Negative' if 'No Bacteria' in x else x)
catnum_list = list()
for inx, row in micbio_df.iterrows():
    shrink_res_str = row['RESULT'].replace(' ','').lower()
    if 'pseudomonasaeruginosa' in shrink_res_str:
        catnum_list.append(1)
        continue
    elif 'acinetobacterbaumannii' in shrink_res_str:
        catnum_list.append(2)
        continue
    elif 'klebsiellapneumoniae' in shrink_res_str:
        catnum_list.append(3)
        continue
    elif 'escherichiacoli' in shrink_res_str:
        catnum_list.append(4)
        continue
    elif 'enterobacter' in shrink_res_str:
        # raise ValueError
        catnum_list.append(5)
        continue
    elif 'klebsiellaaerogenes' in shrink_res_str:
        # raise ValueError
        catnum_list.append(5)
        continue
    elif 'serratiamarcescens' in shrink_res_str:
        # raise ValueError
        catnum_list.append(6)
        continue
    elif ('tuberculosis' in shrink_res_str) and ('mott' not in shrink_res_str):
        # raise ValueError
        catnum_list.append(7)
        continue
    elif ('nontuberculous' in shrink_res_str):
        # raise ValueError
        catnum_list.append(8)
        continue
    elif (('mycobacterium' in shrink_res_str) and ('tuberculosis' not in shrink_res_str)) or ('nontuberculous' in shrink_res_str) or ('mott' in shrink_res_str):
        # raise ValueError
        catnum_list.append(8)
        continue
    elif ('negative'==shrink_res_str) or ('nobacteria' in shrink_res_str) or ('throatnormalflora' in shrink_res_str) or ('nogrowth' in shrink_res_str) or ('noaerobicoranaerobicbacteriaisolated' in shrink_res_str) or ('nomicroorganism' in shrink_res_str) or ('nowbc' in shrink_res_str) or (shrink_res_str in ['rare(<1)wbc','afew(1~5)wbc','moderate(5~10)wbc','many(>10)wbc','stain불가(q.n.s.검체량부족)','[','', 'nowbcnowbc','nowbc, noorgani','nowbc,nonowbc,noorganismorganism']):
        # raise ValueError
        catnum_list.append(10)
        continue
    else:
        # raise ValueError
        catnum_list.append(9)
        continue

micbio_df['MICBIO_CATNUM'] = catnum_list
micbio_df['MICBIO'] = micbio_df['MICBIO_CATNUM'].map(micbio_index_dict)
# micbio_df.columns
micbio_df['DATE'] = micbio_df['DATETIME'].map(lambda x:x.split(' ')[0])
if not os.path.exists(f'{output_dir}'):
    os.mkdir(f'{output_dir}')
micbio_df[['UID','DATE','MICBIO_CATNUM','MICBIO','CODE','NAME']].to_csv(f"{output_dir}/final_microbiology_df.csv", encoding='utf-8-sig', index=False)


# micbio_df['MICBIO_CATNUM'].describe()
# sns.displot(micbio_df['MICBIO_CATNUM'])
micbio_stats = micbio_df.groupby('MICBIO_CATNUM', as_index=False)['UID'].agg('count')
micbio_stats['MICBIO'] = micbio_stats['MICBIO_CATNUM'].map(micbio_index_dict)
print(micbio_stats[['MICBIO_CATNUM','MICBIO','UID']])
# = micbio_df['RESULT'].map(lambda x: 'Negative' if 'No Bacteria' in x else x)

#

# for inx, x in enumerate(micbio_df['RESULT'].unique()):
#     print(f"({inx}) {x}")
