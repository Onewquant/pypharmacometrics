from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## DEMO Covariates Loading

# demo_df = pd.read_csv(f"{output_dir}/patient_info.csv",encoding='utf-8-sig')
# demo_df['AGE'] = demo_df['AGE'].map(lambda x: float(x.replace('개월',''))/12 if '개월' in x else float(x.replace('세','')))
# demo_df = demo_df.rename(columns={'ID':'UID'})
# demo_df = demo_df[['UID', 'AGE', 'SEX', 'HT', 'WT']].copy()
# demo_df['UID'] = demo_df['UID'].astype(str)
# demo_df['HT'] = demo_df['HT'].replace('44.5/45',np.nan).replace('R)600',np.nan).astype(float)
# demo_df['WT'] = demo_df['WT'].replace('71.4(+0.8)',72.2).astype(float)

"""
# 10036912 - 키가 44.5/45 이렇게 되어 있음.
# 10070501 - 키가 비어있음
"""

## LAB Covariates Loading
# lab_covar_rawcols = ['UID', 'DATETIME', 'Albumin', 'AST', 'AST(GOT)', 'ALT', 'ALT(GPT)', 'CRP', 'T. Bil.', 'γ-GT', 'Glucose', 'iCa','Cr (S)', 'Creatinine']
# lab_covar_modelcols = ['UID', 'DATETIME', 'ALB', 'AST', 'ALT','TBIL','GGT', 'CRP', 'CREATININE','ICAL']
#
# uid_fulldt_flist = glob.glob(f"{output_dir}/uid_lab_df/uid_fulldt_df(*).csv")
# totlab_df = list()
# for inx, fpath in enumerate(uid_fulldt_flist):
#     print(f'({inx}) {fpath.split("uid_lab_df")[-1][1:]}')
#     uid_lab_df = pd.read_csv(fpath)[lab_covar_rawcols].copy()
#     totlab_df.append(uid_lab_df)
# totlab_df = pd.concat(totlab_df, ignore_index=True)
# # totlab_df = totlab_df[lab_covar_rawcols].copy()
# for c in list(totlab_df.columns)[2:]:  # break
#     totlab_df[c] = totlab_df[c].map(lambda x: x if type(x) == float else float(re.findall(r'[\d]+.*[\d]*', str(x))[0]))
# totlab_df['UID'] = totlab_df['UID'].astype(str)
# totlab_df['AST'] = totlab_df[['AST', 'AST(GOT)']].max(axis=1)
# totlab_df['ALT'] = totlab_df[['ALT', 'ALT(GPT)']].max(axis=1)
# totlab_df['CREATININE'] = totlab_df[['Cr (S)', 'Creatinine']].max(axis=1)
# # set(totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'].unique())
# totlab_df = totlab_df.rename(columns={'Albumin': 'ALB','T. Bil.':'TBIL','γ-GT':'GGT','Glucose':'GLU','iCa':'ICAL'})
# totlab_df = totlab_df[lab_covar_modelcols].copy()
# # totlab_df = totlab_df[['UID', 'DATETIME', 'ALB', 'AST', 'ALT', 'CRP', 'CALPRTSTL', 'CREATININE']].copy()
# for c in list(totlab_df.columns)[2:]:  # break
#     totlab_df[c] = totlab_df[c].map(lambda x: x if type(x) == float else float(re.findall(r'[\d]+.*[\d]*', str(x))[0]))
#
# totlab_df = totlab_df.drop_duplicates(['UID','DATETIME'])
# totlab_df.to_csv(f"{output_dir}/totlab_df.csv", encoding='utf-8-sig', index=False)
#
#



totlab_df = pd.read_csv(f"{output_dir}/totlab_df.csv")

# [c for c in totlab_df.columns.unique() if 'ada' in c.lower()]

# totlab_df[~totlab_df['ALT'].isna()]['ALT']

## Modeling Data Loading
drug='amk'
modeling_datacheck_dir = f"{output_dir}/amk_modeling_datacheck"
covar_modeling_df = pd.read_csv(f'{modeling_datacheck_dir}/{drug}_modeling_datacheck.csv')

# covar_modeling_df = pd.read_csv(f"{output_dir}/{drug}_modeling_df_filt.csv")
covar_modeling_df['UID']= covar_modeling_df['UID'].astype(str)
covar_modeling_df['DATETIME_ORI'] = covar_modeling_df['DATETIME'].copy()
covar_modeling_df['DATETIME'] = covar_modeling_df['DATETIME'].map(lambda x:x.split('T')[0])

# modeling_df['UID'].drop_duplicates()
# raise ValueError

## DEMO Covariates Loading2

# age_df = pd.read_csv(f"{output_dir}/patient_info.csv",encoding='utf-8-sig')
# age_df = age_df.rename(columns={'ID':'UID'})
# age_df = age_df[['UID', 'SEX', ]].copy()
# age_df['UID'] = age_df['UID'].astype(str)

# ['AGE','HT', 'WT']

demo_df = pd.read_csv(f"{output_dir}/final_req_ptinfo_data.csv",encoding='utf-8-sig')
demo_df['UID'] = demo_df['UID'].astype(str)
demo_df['SEX'] = demo_df['SEX'].map({'남':'M','여':'F'})

# pt_info['AGE'] = pt_info['AGE'].map(lambda x: float(x.replace('개월',''))/12 if '개월' in x else float(x.replace('세','')))
# adult_pids = pt_info[pt_info['AGE'] >= 19].copy()['ID']
mindt_df = covar_modeling_df.groupby('UID', as_index=False)['DATETIME'].min()
mindt_df['MIN_DATE'] = mindt_df['DATETIME'].map(lambda x:x.split('T')[0])
mindt_df = mindt_df[['UID','MIN_DATE']].copy()
# mindt_df[mindt_df['UID']=='11356291']

# vs_df_ori = pd.read_csv(f"{output_dir}/final_req_vs_data.csv",encoding='utf-8-sig')
vs_df = pd.read_csv(f"{output_dir}/final_req_bodysize_data.csv",encoding='utf-8-sig')
vs_df['UID'] = vs_df['UID'].astype(str)
vs_df = vs_df.merge(mindt_df, on=['UID'], how='left')
vs_df = vs_df.dropna(subset=['MIN_DATE'])
vs_df = vs_df.drop_duplicates(['UID','DATETIME','VS'])
# pd.pivot_table(vs_df, values='VALUE', index=['UID','DATETIME','MIN_DATE'], aggfunc=["mean", "max", "min"])
vs_df['VALUE'] = vs_df['VALUE'].map(lambda x:float(re.findall(r'\d+[\.]?\d*',str(x))[0]) if len(re.findall(r'\d+[\.]?\d*',str(x)))!=0 else np.nan)
# vs_df[vs_df['UID']=='11356291']
# vs_df_ori[(vs_df_ori['UID']==29009921)&(vs_df_ori['VS']=='WT')]

## BMI 너무 크거나 작은 값 -> 제거
max_bmi = 50
min_bmi = 10
vs_df = vs_df[~((vs_df['VS']=='BMI')&((vs_df['VALUE']>=max_bmi)| (vs_df['VALUE']<=min_bmi)))].reset_index(drop=True)


## WT 너무 작은 값 -> 제거
vs_df = vs_df[~((vs_df['VALUE'] < 20) & (vs_df['VS']=='WT'))].reset_index(drop=True)

## WT 너무 큰 값 -> 소수점 조정
for inx, row in vs_df[(vs_df['VALUE'] > 200) & (vs_df['VS']=='WT')].iterrows():
    new_val = 60
    if (row['VALUE'] >= 10**3) and (row['VALUE'] < 10**4):
        new_val = row['VALUE']/100
    elif (row['VALUE'] >= 10**2) and (row['VALUE'] < 10**3):
        new_val = row['VALUE']/10
    else:
        raise ValueError
    vs_df.at[inx, 'VALUE'] = new_val

## WT 값에 HT가 들어가 있는 경우 조정 -> VS 항목을 HT로 바꿈
for inx in vs_df[(vs_df['VALUE'] > 130) & (vs_df['VS']=='WT')].index:
    vs_df.at[inx,'VS'] = 'HT'

## HT 너무 작은 값 -> 제거
vs_df = vs_df[~((vs_df['VS'] == 'HT') & (vs_df['VALUE'] < 30))].reset_index(drop=True)

## HT 너무 큰 값 -> 소수점 조정
for inx, row in vs_df[(vs_df['VALUE']>200)&(vs_df['VS']=='HT')].iterrows():
    new_val = 160
    if (row['VALUE'] >= 10**3) and (row['VALUE'] < 10**4):
        new_val = row['VALUE']/10
    else:
        raise ValueError
    vs_df.at[inx, 'VALUE'] = new_val

## HT 값에 WT가 들어가 있는 경우 조정 -> VS 항목을 WT로 바꿈
for inx in vs_df[(vs_df['VALUE'] >= 30) & (vs_df['VALUE'] < 100) & (vs_df['VS']=='HT')].index:
    vs_df.at[inx,'VS'] = 'WT'

vs_df = vs_df.drop_duplicates(['UID','DATETIME','VS'], ignore_index=True)
# vs_df_ori = pd.read_csv(f"{output_dir}/final_req_vs_data.csv",encoding='utf-8-sig')
# vs_df_ori[(vs_df_ori['UID']==11437323)&(vs_df_ori['VS']=='WT')]

vs_df = vs_df.pivot(index=['UID','DATETIME','MIN_DATE'], columns='VS', values='VALUE').reset_index(drop=False)
vs_df.columns.name = None
vs_df = vs_df.sort_values(['UID','DATETIME'])

for fill_col in ['HT','BMI','WT']: # 순서 바꾸면 안 됨! (키가 가장 보존적 var)

    # fill_col = 'WT'

    new_vs_df = list()
    for uid, uid_vs_df in vs_df.groupby('UID'):
        uid_vs_df[fill_col] = uid_vs_df[fill_col].fillna(method='ffill').fillna(method='bfill')
        new_vs_df.append(uid_vs_df)
    vs_df = pd.concat(new_vs_df)

    for inx, row in vs_df[(~vs_df[['BMI','HT','WT']].isna()*1).sum(axis=1)==2].iterrows():
        na_colname = row[row.isna()].index[0]
        if na_colname=='BMI':
            row[na_colname] = row['WT']/((row['HT']/100)**2)
        elif na_colname=='WT':
            row[na_colname] = row['BMI']*((row['HT']/100)**2)
        elif na_colname=='HT':
            row[na_colname] = 100*(row['WT']/row['BMI'])**0.5
        else:
            raise ValueError
        vs_df.at[inx, na_colname] = row[na_colname]

    # if fill_col=='WT':
    #     raise ValueError

# vs_df[vs_df['UID']=='11356291']
vs_df = vs_df[(~vs_df[['BMI','HT','WT']].isna()*1).sum(axis=1)>0].copy()
# vs_df = vs_df[vs_df['DATETIME'] <= vs_df['MIN_DATE']].copy()
recent_vs_df = vs_df[vs_df['DATETIME'] <= vs_df['MIN_DATE']].drop_duplicates(['UID'], keep='last', ignore_index=True)
recent_vs_df = recent_vs_df.drop(['DATETIME','BMI'],axis=1)[['UID','MIN_DATE','HT','WT']].copy()
# recent_vs_df[recent_vs_df['UID']=='11356291']
# recent_vs_df = vs_df[vs_df['DATETIME'] <= vs_df['MIN_DATE']].drop_duplicates(['UID'], keep='last', ignore_index=True)

after_vs_df = vs_df[vs_df['DATETIME'] > vs_df['MIN_DATE']].drop_duplicates(['UID'], keep='first', ignore_index=True)
after_vs_df = after_vs_df.drop(['DATETIME','BMI'],axis=1)[['UID','MIN_DATE','HT','WT']].copy()
# demo_df.merge()

nearest_vs_df = pd.concat([recent_vs_df, after_vs_df]).drop_duplicates(['UID'], keep='first')
# nearest_vs_df[nearest_vs_df['UID']=='11356291']
# nearest_vs_df[(nearest_vs_df['WT']>100)]

# ptinfo_df[ptinfo_df['UID']=='11356291']
ptinfo_df = pd.read_csv(f"{output_dir}/patient_info.csv",encoding='utf-8-sig')
# demo_df['AGE'] = demo_df['AGE'].map(lambda x: float(x.replace('개월',''))/12 if '개월' in x else float(x.replace('세','')))
ptinfo_df = ptinfo_df.rename(columns={'ID':'UID'})
ptinfo_df = ptinfo_df[['UID', 'SEX', 'HT', 'WT']].copy()
ptinfo_df['UID'] = ptinfo_df['UID'].astype(str)
ptinfo_df['HT'] = ptinfo_df['HT'].replace('44.5/45',np.nan).replace('R)600',np.nan).astype(float)
ptinfo_df['WT'] = ptinfo_df['WT'].replace('71.4(+0.8)',72.2).astype(float)
ptinfo_df = ptinfo_df.merge(mindt_df, on='UID',how='left')
ptinfo_df = ptinfo_df[['UID','MIN_DATE','HT','WT']].copy()
supplement_ptinfo_df = ptinfo_df[ptinfo_df['UID'].isin(list(set(ptinfo_df['UID']).difference(set(nearest_vs_df['UID']))))].copy()
# supplement_ptinfo_df[supplement_ptinfo_df['UID']=='11356291']
# nearest_vs_df.drop_duplicates(['UID'], keep='first')
# supplement_ptinfo_df['HT']

integ_vs_df = pd.concat([nearest_vs_df, supplement_ptinfo_df]).dropna(subset=['MIN_DATE'])
integ_vs_df = integ_vs_df.fillna(integ_vs_df.median(numeric_only=True))
# integ_vs_df[integ_vs_df['UID']=='11356291']

demo_df = integ_vs_df.merge(demo_df, on='UID',how='left')
demo_df['AGE'] = ((pd.to_datetime(demo_df['MIN_DATE'])-pd.to_datetime(demo_df['BIRTH_DATE'])).dt.total_seconds()/(365.25*86400)).map(int)
demo_df = demo_df[['UID', 'AGE', 'SEX', 'HT', 'WT']].copy()

# new_vs_df = list()
# for uid, uid_vs_df in vs_df.groupby('UID'):
#     uid_vs_df = uid_vs_df.fillna(method='ffill').fillna(method='bfill')
#     new_vs_df.append(uid_vs_df)
# vs_df = pd.concat(new_vs_df)

# vs_df.isna().sum()

#     # raise ValueError
# new_vs_df = list()
# for uid, uid_vs_df in vs_df.groupby('UID'):
#     uid_vs_df = uid_vs_df.fillna(method='ffill').fillna(method='bfill')
#     new_vs_df.append(uid_vs_df)
# vs_df = pd.concat(new_vs_df)

# vs_df[(~vs_df[['BMI', 'HT', 'WT']].isna() * 1).sum(axis=1) == 1]
#
#
# for uid, uid_vs_df in vs_df.groupby('UID'):
#     raise ValueError
# vs_df = vs_df[vs_df['DATETIME'] <= vs_df['MIN_DATE']].copy()


# demo_df = demo_df.merge(mindt_df, on=['UID'], how='left')
# demo_df = demo_df.dropna(subset=['MIN_DATE'])
# demo_df['AGE'] = ((pd.to_datetime(demo_df['MIN_DATE'])-pd.to_datetime(demo_df['BIRTH_DATE'])).dt.total_seconds()/(365.25*86400)).map(int)


"""
# list(totlab_df.columns)
# ['UID', 'DATETIME', '1,25-VitD3', '25-OH VIT D TOTAL', '25-OH Vit. D3(초진용)', '25-OH Vit.D (D3/D2)', '25-OH Vit.D (Total)', '25-OH Vit.D2', '25-OH Vit.D3', '3-methoxytyramine (Plasma)', '5-HIAA (24h urine)', 'ABO', 'ACL, IgG', 'ACL, IgM', 'ACTH', 'ACTHSTbase', 'ADA', 'AFP', 'AFP-L3(%)', 'ALT', 'ALT(GPT)', 'AMA', 'ANC', 'ANCA', 'ARR', 'ASCA IgA & IgG', 'ASO', 'AST', 'AST(GOT)', 'AT III', 'Ab scr', 'Acetoacetate', 'Adalimumab Quantification', 'Albumin', 'Aldo(B)', 'Alk. phos', 'Alk. phos.', 'Ammo', 'Amylase(Flu)', 'Amylase(S)', 'Anion gap', 'Anti - ccp', 'Anti Mullerian Hormone (AMH)', 'Anti RNP Ab', 'Anti Sm Antibody', 'Anti ds DNA', 'Anti-HBs', 'Anti-HCV', 'Anti-Infliximab Ab [정밀면역검사] (정량)', 'Anti-LKM', 'Anti-TSH receptor', 'ApoA1', 'ApoB', 'AtypicalLc', 'B2-MG(S)', 'B2-MG(U)', 'BE', 'BIL', 'BLD', 'BNP', 'BST', 'BUN', 'Bact.', 'Bacteria', 'Band.neut.', 'Basophil', 'Beta hydroxybutyric acid', 'Blast', 'C-Peptide(S)', 'C. difficile GDH', 'C. difficile toxin', 'C.diffcile toxin', 'C3', 'C4', 'CA 125', 'CA 15-3', 'CA 19-9', 'CBC (em) (diff), RDW제외', 'CEA', 'CK', 'CK(CPK)', 'CK-MB (em)', 'CMV Ag', 'CMV IgM', 'CORT30', 'CORT60', 'CORTbasal', 'CPEPSTbase', 'CRP', 'CTx (C-telopeptide)', 'Ca', 'Ca, total', 'Calcium', 'Calprotectin (Serum)', 'Calprotectin (Stool)', 'Cast', 'Casts', 'Cerulo', 'Chol.', 'Chromogranin A', 'Cl', 'Cl(U) (em)', 'Codfish IgE (F3)', 'Collagen 4', 'Color', 'Cortisol', 'Cortisol(S)', 'Cotinine', "Cow's milk IgE (F2)", 'Cr (S)', 'Cr (U)', 'Cr(u)', 'Creatinine', 'Cryoglobu.', 'Crystal', 'Crystals', 'Cu (24h Urine)', 'Cu (Serum)', 'Cystatin C', 'D-dimer', 'D-dimer(em)', 'D. Bil.', 'D. farinae IgE (D2)', 'D. farinae IgG4 (D2)', 'D. pteronyssinus IgE (D1)', 'D. pteronyssinus IgG4 (D1)', 'DHEA-S', 'DRBC', 'Delta ratio', 'Dopamine', 'ECP', 'ESR', 'Egg white IgE (F1)', 'Entero C', 'Eos.count', 'Eosinophil', 'Epinephrine', 'Erythropoiet', 'Erythropoietin', 'Estradiol', 'F. acid', 'FANA', 'FANA Titer', 'FBS', 'FDP정량', 'FSH', 'FT3', 'Factor 10', 'Factor 11', 'Factor 12', 'Factor 2', 'Factor 5', 'Factor 8', 'Factor 9', 'Ferritin', 'Fibrinogen', 'Folate', 'Free Fatty Acid', 'Free PSA', 'Free PSA%', 'Free T3', 'Free T4', 'GLU', 'Gastrin', 'Glu', 'Glucagon', 'Glucose', 'H. pylori IgG Ab', 'H.pylo IgG', 'HA', 'HAV Ab IgG', 'HAV Ab(IgG)', 'HAV Ab(IgM)', 'HAVAb IgM', 'HBcAb IgM', 'HBcAb Total (IgM+IgG) (진단검사의학)', 'HBcAb(IgG)', 'HBeAb', 'HBeAg', 'HBs Ag', 'HBsAb', 'HBsAg', 'HCO3-', 'HCV Ab', 'HDL Chol.', 'HE4 (Human Epididymis Protein 4)', 'HFR', 'HGH', 'HIV Ag/Ab', 'HIV Ag/Ab (건강검진용)', 'Hapto', 'Hb', 'HbA1c-IFCC', 'HbA1c-IFCC (건강증진센타용)', 'HbA1c-NGSP', 'HbA1c-NGSP (건강증진센타용)', 'HbA1c-eAG', 'HbA1c-eAG (건강증진센타용)', 'Hct', 'Helm,ova', 'Human Hb', 'I/Ⅱratio', 'IGF 1', 'IGFBP 3', 'IMA (Ischemia Modified Albumin Test) (em)', 'IRF', 'Ig A', 'Ig G', 'IgE(Total)', 'IgG sub 4', 'IgG sub1', 'IgG sub2', 'IgG sub3', 'IgG sub4', 'IgM', 'Imm.cell', 'Imm.lympho', 'Imm.mono', 'Infliximab Quantification', 'Influenza A&B Ag', 'Insulin', 'Insulin Ab', 'Interleukin-6', 'Iron', 'K', 'K (U) (em)', 'KET', 'Ketone (Beta-hydroxybutyrate)', 'LA', 'LD', 'LD(LDH)', 'LD-R', 'LD-R(B)', 'LDL Chol.', 'LH', 'LUC', 'Lactate', 'Lactic', 'Li', 'Lipase', 'Lp(a)', 'Lymph', 'Lympho', 'Lymphocyte', 'M/C ratio', 'MCH', 'MCHC', 'MCV', 'MFR', 'MN(Mononuclear cell)', 'MPV', 'MTX', 'Metamyelo', 'Metanephrine (Plasma)', 'Mg', 'Mic-Ab', 'MicroAlb', 'Mixing test (PT, aPTT 제외)', 'Monocyte', 'Mycopl. Ab', 'Mycoplasma Ab IgG', 'Mycoplasma Ab IgM', 'Myelocyte', 'Myoglobin', 'Myoglobin (em)', 'N. fat', 'NCV2019 응급용 선별검사 (교수용)', 'NCV2019 입원 선별 1단계 [분자진단]', 'NCV2019 입원 선별 2단계 [분자진단]', 'NIT', 'NMP 22', 'NSE', 'Na', 'Na(U) (em)', 'Norepinephrine', 'Normetanephrine (Plasma)', 'Normoblast', 'O2 CT', 'O2 SAT', 'Osmo-S', 'Osmo-U', 'Osteocalcin', 'Other', 'Others', 'O₂SAT', 'P', 'P. cyst', 'P. troph', 'P1NP', 'PCT', 'PDW', 'PIVKA II', 'PLT', 'PMN(Polymorphonuclear cell)', 'PP2', 'PRO', 'PSA', 'PSA (건증용)', 'PT %', 'PT % (MIX)', 'PT INR', 'PT INR (MIX)', 'PT sec', 'PT sec (MIX)', 'PTH', 'PTH(intact)', 'Pb(Blood)', 'Peanut IgE (F13)', 'PepsinogenⅠ', 'PepsinogenⅡ', 'Phosphorus', 'Pl. Hb', 'Plas.cell', 'Platelet', 'Poly', 'PreAlb', 'Prealbumin', 'Procalcitonin', 'Prolactin', 'Promyelo', 'Prostate Health Index', 'Protein', 'Protein C activity', 'Protein S activity', 'Protein/Creatinine ratio', 'Protozo', 'Pyruvate', 'RBC', 'RDW(CV)', 'RDW(SD)', 'RF', 'ROMA (postmenopausal)', 'ROMA (premenopausal)', 'RPR (VDRL, Auto) (serum)', 'RPR 정량 (serum)', 'RT', 'RTE', 'Renin(B)', 'Reticulocyte', 'Rh D', 'S.G', 'SAA (Serum Amyloid A) (em)', 'SARS-CoV-2 Ag(외래, 응급실, 중환자실 전용)', 'SCC(TA-4)', 'SG', 'SHBG', 'SQE', 'SS-A/Ro(52) Ab', 'SS-A/Ro(60) Ab', 'SS-B/La Ab', 'Salt intake', 'Seg.neut.', 'Selenium', 'Smooth muscle Ab', 'Sodium (random urine)', 'Soybean IgE (F14)', 'Sperm', 'T CO2', 'T. Bil.', 'T. Protein', 'T.B', 'T.P Ab Total (IgM+IgG) 정밀면역', 'T3', 'TCO2', 'TG', 'TIBC', 'TRE', 'TSH', 'Testosterone', 'Total IgE (진단검사의학과)', 'Total ketone', 'Toxo IgM', 'Toxocariasis (ELISA, IgG antibody)', 'Toxopla.Ab', 'Transferrin', 'Troponin I (em) [ng/mL]', 'Troponin I (em) [pg/mL]', 'Troponin T [ng/mL]', 'Troponin T [pg/mL]', 'Tryptase', 'Tur', 'Turbidity', 'UA', 'UIBC', 'URO', 'Uric acid', 'Urine S. pneumoniae Ag', 'VWF rel Ag', 'VWFrist co', 'VZV Ab IgG', 'Valproic', 'Vanco', 'Vitamin B12', 'Vitamin B₁ (Thiamin)', 'WBC', 'WBC(s)', 'Wheat IgE (F4)', 'X-mat', 'Yeast like organism', 'Zn (serum)', 'aPTT', 'aPTT (MIX)', 'alpha1-Antitrypsin (stool)', 'd1', 'd2', 'eGFR (CKD-EPI Cr-Cys)', 'eGFR (CKD-EPI Cys)', 'eGFR-CKD-EPI', 'eGFR-Cockcroft-Gault', 'eGFR-MDRD', 'eGFR-Schwartz(소아)', 'em.WBC', 'epine(U)', 'hCG', 'hsCRP', 'i6', 'iCa', 'iMg', 'metanephrine(U)', 'mito Ab', 'norepine(U)', 'p2PSA', 'pCO₂', 'pH', 'pO₂', 'proBNP (em)', 'tHcyst', 'u.RBC', 'u.WBC', 'β2 GP1-IgG', 'β2 GP1-IgM', 'γ-GT', '검체처리 1단계(채혈TUBE1~5개)', '마약선별검사', '모발 중금속 및 미네랄 40종 검사', '연구검체보관', '연구용채혈+검체처리1단계(채혈TUBE1~5개)', '연구용채혈+검체처리2단계', '절대단구수', '절대림프구수']
"""
# len(modeling_df)
totlab_df['UID']= totlab_df['UID'].astype(str)
covar_modeling_df = covar_modeling_df.merge(totlab_df, on=['UID','DATETIME'], how='left')
covar_modeling_df = covar_modeling_df.merge(demo_df, on=['UID'], how='left')
# modeling_df['UID'].iloc[0]
# demo_df['UID'].iloc[0]

## Covariates의 NA value 처리 (ffill 먼저 시도, 없으면 bfill, 그것도 없으면 전체의 median 값)
keep_ori_covar_cols = ['DATETIME_ORI','REC_REASON']
keep_ori_covar_df = covar_modeling_df[keep_ori_covar_cols].copy()

md_df_list = list()
for md_inx,md_df in covar_modeling_df.drop(keep_ori_covar_cols, axis=1).groupby(['UID']):
    md_df = md_df.sort_values(['TIME']).fillna(method='ffill').fillna(method='bfill')
    md_df_list.append(md_df)
covar_modeling_df = pd.concat(md_df_list).reset_index(drop=True)
covar_modeling_df = covar_modeling_df.fillna(covar_modeling_df.median(numeric_only=True))

covar_modeling_df = pd.concat([covar_modeling_df, keep_ori_covar_df], axis=1)


## Covariates의 NA value 처리 (ffill 먼저 시도, 없으면 bfill, 그것도 없으면 전체의 median 값)

# md_df_list = list()
# for md_inx, md_df in modeling_df.groupby(['UID']):
#     # md_df = md_df.sort_values(['DATETIME']).fillna(method='ffill')
#     md_df = md_df.sort_values(['DATETIME'])
#     md_df_list.append(md_df)
# modeling_df = pd.concat(md_df_list).reset_index(drop=True)

# modeling_df.fillna('.', inplace=True)
# modeling_df['UID'].drop_duplicates()
# raise ValueError

## Modeling Data Saving
# data_check_cols = ['ID','UID','NAME','DATETIME','TIME','DV','MDV','AMT','DUR','CMT','IBD_TYPE','ADDED_ADDL'] + list(modeling_df.loc[:,'DRUG':].iloc[:,1:].columns)
# modeling_df[data_check_cols].to_csv(f'{output_dir}/{drug}_{mode_str}_datacheck_covar.csv', index=False, encoding='utf-8-sig')
# raise ValueError

modeling_covar_dir = f"{output_dir}/amk_modeling_covar"
if not os.path.exists(modeling_covar_dir):
    os.mkdir(modeling_covar_dir)
# covar_modeling_df.columns
# covar_modeling_df['DATETIME_ORI']
right_covar_col = 'TDM_REQ_DATE'
datacheck_cols = ['ID',	'UID', 'NAME', 'DATETIME', 'TIME', 'DV', 'MDV', 'AMT', 'RATE', 'CMT','LLOQ','BLQ','REC_REASON'] + list(covar_modeling_df.loc[:,right_covar_col:].iloc[:,1:].columns)
covar_modeling_df[datacheck_cols].to_csv(f'{modeling_covar_dir}/{drug}_modeling_datacheck_covar.csv', index=False, encoding='utf-8-sig')
# covar_modeling_df = covar_modeling_df.drop(keep_ori_covar_cols, axis=1)
modeling_cols = ['ID','NAME','TIME','TAD','DV','MDV','CMT','AMT','RATE','UID','LLOQ','BLQ'] + list(covar_modeling_df.loc[:,right_covar_col:].iloc[:,1:].columns)

# modeling_df['AGE'] = modeling_df.apply(lambda x: int((datetime.strptime(x['DATETIME'],'%Y-%m-%d') - datetime.strptime(x['AGE'],'%Y-%m-%d')).days/365.25), axis=1)
covar_modeling_df['SEX'] = covar_modeling_df['SEX'].map({'M':0,'F':1})

# ## 최근 투약 기준 시간으로 TAD 설정
# covar_modeling_df['TAD'] = np.nan
#
# # 농도값 바로 전에 투약기록 있는 경우
# tad_cond0 = (covar_modeling_df['MDV']==0)
# tad_cond1 = (covar_modeling_df['MDV'].shift(1).fillna(1.0)==1)
# tad_cond2 = (covar_modeling_df['ID']==covar_modeling_df['ID'].shift(1))
# tad_cond3 = (covar_modeling_df['TAD'].isna())
# for inx in covar_modeling_df[tad_cond0&tad_cond1&tad_cond2&tad_cond3].index:
#     covar_modeling_df.at[inx,'TAD'] = covar_modeling_df.at[inx,'TIME']-covar_modeling_df.at[inx-1,'TIME']
#
# # 농도값 바로 전에 농도기록 있는 경우
# tad_cond0 = (covar_modeling_df['MDV']==0)
# tad_cond1 = (covar_modeling_df['MDV'].shift(1).fillna(1.0)==1)
# tad_cond2 = (covar_modeling_df['ID']==covar_modeling_df['ID'].shift(1))
# tad_cond3 = (covar_modeling_df['TAD'].isna())
# tad_cond4 = (~(covar_modeling_df['TAD'].shift(1).isna()))
# while len(covar_modeling_df[tad_cond0&(~tad_cond1)&tad_cond2&tad_cond3&tad_cond4])!=0:
#
# # covar_modeling_df[covar_modeling_df['TAD']<0]
# # covar_modeling_df.to_csv(f'{output_dir}/{drug}_checkcheck.csv',index=False, encoding='utf-8-sig')
#     for inx in covar_modeling_df[tad_cond0&(~tad_cond1)&tad_cond2&tad_cond3&tad_cond4].index:
#         covar_modeling_df.at[inx, 'TAD'] = covar_modeling_df.at[inx-1, 'TAD'] + covar_modeling_df.at[inx, 'TIME'] - covar_modeling_df.at[inx - 1, 'TIME']
#
#     tad_cond0 = (covar_modeling_df['MDV'] == 0)
#     tad_cond1 = (covar_modeling_df['MDV'].shift(1).fillna(1.0) == 1)
#     tad_cond2 = (covar_modeling_df['ID'] == covar_modeling_df['ID'].shift(1))
#     tad_cond3 = (covar_modeling_df['TAD'].isna())
#     tad_cond4 = (~(covar_modeling_df['TAD'].shift(1).isna()))
#
# covar_modeling_df.columns
# covar_modeling_df['TAD'] = covar_modeling_df['TAD'].replace(np.nan,0)
covar_modeling_df['TAD'] = add_time_after_dosing_column(df=covar_modeling_df)
covar_modeling_df['TAD'] = covar_modeling_df['TAD'].map(lambda x:round(x,7))
covar_modeling_df['TIME'] = covar_modeling_df['TIME'].map(lambda x:round(x,7))

to_change_rate = covar_modeling_df[(covar_modeling_df['TAD'] <= 0.5)&(covar_modeling_df['TAD'] != 0)&(covar_modeling_df['MDV'] == 0)].copy()
for tad_inx in (to_change_rate.index): #break
    rate_change_inx = tad_inx-1
    covar_modeling_df.at[rate_change_inx,'RATE'] = round(float(covar_modeling_df.at[rate_change_inx,'AMT'])/(covar_modeling_df.at[tad_inx,'TAD']*0.9),7)
# raise ValueError


covar_modeling_df = covar_modeling_df[modeling_cols].sort_values(['ID','TIME'], ignore_index=True)
covar_modeling_df = covar_modeling_df.drop(['NAME'],axis=1)
modeling_input_line = str(list(covar_modeling_df.columns)).replace("', '"," ")

print(f"Mode: {modeling_input_line}")

## 데이터 추가 조정
##############
# raise ValueError
# GFR < 15 대상자 삭제
covar_modeling_df = covar_modeling_df[~(covar_modeling_df['ID'].isin([156, 159, 255, 263, 449, 816, 1081, 1851, 2093, 2377, 840]))].copy()

# Wrong Data row 수정
to_be_adj = [
             # T/P 두번 기록, DV 움직임 비정상적
             # {'ID':377,'TIME':190.92,'DV':35.0, 'NEW_TIME':179.033},
             {'ID':367,'TIME':60.2,'DV':17.36, 'NEW_TIME':73},
             {'ID':353,'TIME':103.17,'DV':17.11, 'NEW_TIME':120},
]
# covar_modeling_df[(covar_modeling_df['ID']==row['ID'])&(covar_modeling_df['TIME'].map(lambda x:round(x,1)) == round(row['TIME'],1))&(covar_modeling_df['DV']==str(row['DV']))].copy()
# covar_modeling_df[(covar_modeling_df['ID']==row['ID'])]['TAD']
# to_be_adj_inx = list()
for row in to_be_adj:
    # raise ValueError
    # covar_modeling_df[(covar_modeling_df['ID'] == row['ID'])]
    # covar_modeling_df[(covar_modeling_df['ID'] == row['ID']) & (covar_modeling_df['TIME'].map(lambda x:round(x,1)) == round(row['TIME'],1))]
    adj_row = covar_modeling_df[(covar_modeling_df['ID']==row['ID'])&(covar_modeling_df['TIME'].map(lambda x:round(x,1)) == round(row['TIME'],1))&(covar_modeling_df['DV']==str(row['DV']))].copy()
    adj_inx = adj_row.index[0]
    covar_modeling_df.at[adj_inx, 'TIME'] = row['NEW_TIME']
    # raise ValueError
    # to_be_adj_inx.append()

# Wrong Data row 수정2
to_be_adj2 = [
             {'ID':777,'SEARCH_VARI':'TIME','SEARCH_VAL':498.21,'CHANGE_VAR':'TAD','PREV_VAL':0.0666667,'NEW_VAL':1},
             {'ID':777,'SEARCH_VARI':'TIME','SEARCH_VAL':498.21,'CHANGE_VAR':'TIME','PREV_VAL':498.21,'NEW_VAL':499.15},
             # {'ID':377,'SEARCH_VARI':'TIME','SEARCH_VAL':179.033,'CHANGE_VAR':'TAD','PREV_VAL':12.8833333,'NEW_VAL':1},

             {'ID':367,'SEARCH_VARI':'TIME','SEARCH_VAL':73,'CHANGE_VAR':'TAD','PREV_VAL':12.2,'NEW_VAL':1},
             # {'ID': 367, 'SEARCH_VARI': 'TIME', 'SEARCH_VAL': 73, 'CHANGE_VAR': 'TAD', 'PREV_VAL': 12.2, 'NEW_VAL': 1},

             {'ID':763,'SEARCH_VARI':'TIME','SEARCH_VAL':0.0,'CHANGE_VAR':'AMT','PREV_VAL':70,'NEW_VAL':700},
             {'ID':763,'SEARCH_VARI':'TIME','SEARCH_VAL':22.25,'CHANGE_VAR':'AMT','PREV_VAL':70,'NEW_VAL':700},
             {'ID':763,'SEARCH_VARI':'TIME','SEARCH_VAL':48.01667,'CHANGE_VAR':'AMT','PREV_VAL':70,'NEW_VAL':700},
             {'ID':763,'SEARCH_VARI':'TIME','SEARCH_VAL':74.41667,'CHANGE_VAR':'AMT','PREV_VAL':70,'NEW_VAL':700},

             {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':0.0,'CHANGE_VAR':'AMT','PREV_VAL':100,'NEW_VAL':1000.0},
             {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':24,'CHANGE_VAR':'AMT','PREV_VAL':100,'NEW_VAL':1000.0},
             {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':48,'CHANGE_VAR':'AMT','PREV_VAL':100,'NEW_VAL':1000.0},
             {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':72,'CHANGE_VAR':'AMT','PREV_VAL':100,'NEW_VAL':1000.0},
             {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':96,'CHANGE_VAR':'AMT','PREV_VAL':100,'NEW_VAL':1000.0},

             {'ID':656,'SEARCH_VARI':'TIME','SEARCH_VAL':57.11667,'CHANGE_VAR':'AMT','PREV_VAL':100,'NEW_VAL':1000.0},
             {'ID':656,'SEARCH_VARI':'TIME','SEARCH_VAL':68.45,'CHANGE_VAR':'AMT','PREV_VAL':100,'NEW_VAL':1000.0},
             {'ID':656,'SEARCH_VARI':'TIME','SEARCH_VAL':81.2,'CHANGE_VAR':'AMT','PREV_VAL':100,'NEW_VAL':1000.0},

             # 추후 아래서 삭제 작업도 필요한 목록들
             {'ID':935,'SEARCH_VARI':'TIME','SEARCH_VAL':126.9667,'CHANGE_VAR':'AMT','PREV_VAL':600,'NEW_VAL':700.0},
             {'ID':1278,'SEARCH_VARI':'TIME','SEARCH_VAL':71.866667,'CHANGE_VAR':'AMT','PREV_VAL':500,'NEW_VAL':600.0},
             {'ID':81,'SEARCH_VARI':'TIME','SEARCH_VAL':0,'CHANGE_VAR':'AMT','PREV_VAL':500,'NEW_VAL':525.0},
             {'ID':1679,'SEARCH_VARI':'TIME','SEARCH_VAL':204.4667,'CHANGE_VAR':'AMT','PREV_VAL':500,'NEW_VAL':530.0},

             # # 위 사람들 Rate 조정
             #
             # {'ID':763,'SEARCH_VARI':'TIME','SEARCH_VAL':0.0,'CHANGE_VAR':'RATE','PREV_VAL':70,'NEW_VAL':700},
             # {'ID':763,'SEARCH_VARI':'TIME','SEARCH_VAL':22.25,'CHANGE_VAR':'RATE','PREV_VAL':70,'NEW_VAL':700},
             # {'ID':763,'SEARCH_VARI':'TIME','SEARCH_VAL':48.01667,'CHANGE_VAR':'RATE','PREV_VAL':70,'NEW_VAL':700},
             # {'ID':763,'SEARCH_VARI':'TIME','SEARCH_VAL':74.41667,'CHANGE_VAR':'RATE','PREV_VAL':70,'NEW_VAL':700},
             #
             # {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':0.0,'CHANGE_VAR':'RATE','PREV_VAL':100,'NEW_VAL':1000.0},
             # {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':24,'CHANGE_VAR':'RATE','PREV_VAL':100,'NEW_VAL':1000.0},
             # {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':48,'CHANGE_VAR':'RATE','PREV_VAL':100,'NEW_VAL':1000.0},
             # {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':72,'CHANGE_VAR':'RATE','PREV_VAL':100,'NEW_VAL':1000.0},
             # {'ID':151,'SEARCH_VARI':'TIME','SEARCH_VAL':96,'CHANGE_VAR':'RATE','PREV_VAL':100,'NEW_VAL':1000.0},
             #
             # {'ID':656,'SEARCH_VARI':'TIME','SEARCH_VAL':57.11667,'CHANGE_VAR':'RATE','PREV_VAL':100,'NEW_VAL':1000.0},
             # {'ID':656,'SEARCH_VARI':'TIME','SEARCH_VAL':68.45,'CHANGE_VAR':'RATE','PREV_VAL':100,'NEW_VAL':1000.0},
             # {'ID':656,'SEARCH_VARI':'TIME','SEARCH_VAL':81.2,'CHANGE_VAR':'RATE','PREV_VAL':100,'NEW_VAL':1000.0},
             #
             # # 추후 아래서 삭제 작업도 필요한 목록들
             # {'ID':935,'SEARCH_VARI':'TIME','SEARCH_VAL':126.9667,'CHANGE_VAR':'RATE','PREV_VAL':600,'NEW_VAL':700.0},
             # {'ID':1278,'SEARCH_VARI':'TIME','SEARCH_VAL':71.866667,'CHANGE_VAR':'RATE','PREV_VAL':500,'NEW_VAL':600.0},
             # {'ID':81,'SEARCH_VARI':'TIME','SEARCH_VAL':0,'CHANGE_VAR':'RATE','PREV_VAL':500,'NEW_VAL':525.0},
             # {'ID':1679,'SEARCH_VARI':'TIME','SEARCH_VAL':204.4667,'CHANGE_VAR':'RATE','PREV_VAL':500,'NEW_VAL':530.0},
]

# to_be_adj_inx = list()
for row in to_be_adj2:
    # raise ValueError
    # if row['ID']==935:
    # covar_modeling_df[(covar_modeling_df['ID'] == 1679)&(covar_modeling_df['TIME'] > 192)][['ID', 'TIME', 'AMT']]
    # if row['ID'] == 935:
    #     covar_modeling_df[(covar_modeling_df['ID'] == 935)][['ID','TIME','AMT']]
    #     covar_modeling_df[(covar_modeling_df['ID'] == 81)][['ID','TIME','AMT']]
    #     covar_modeling_df[(covar_modeling_df['ID'] == 1278)][['ID','TIME','AMT']]
    #     covar_modeling_df[(covar_modeling_df['ID'] == 1679)][['ID','TIME','AMT']]
    #     covar_modeling_df[(covar_modeling_df['ID'] == row['ID'])][['ID','TIME','AMT']]
    #     raise ValueError
    if row['CHANGE_VAR'] in ['TAD', 'TIME']:
    # covar_modeling_df[(covar_modeling_df['ID'] == row['ID']) & (covar_modeling_df[row['SEARCH_VARI']].map(lambda x:round(x,1)) == round(row['SEARCH_VAL'],1))]
        adj2_row = covar_modeling_df[(covar_modeling_df['ID']==row['ID'])&(covar_modeling_df[row['SEARCH_VARI']].map(lambda x:round(x,1)) == round(row['SEARCH_VAL'],1))&(covar_modeling_df['DV']!='.')].copy()
        adj2_inx = adj2_row.index[0]
        covar_modeling_df.at[adj2_inx, row['CHANGE_VAR']] = row['NEW_VAL']
    elif row['CHANGE_VAR'] in ['AMT']:
        adj2_row = covar_modeling_df[(covar_modeling_df['ID']==row['ID'])&(covar_modeling_df[row['SEARCH_VARI']].map(lambda x:round(x,1)) == round(row['SEARCH_VAL'],1))&(covar_modeling_df['AMT']!='.')].copy()
        adj2_inx = adj2_row.index[0]
        covar_modeling_df.at[adj2_inx, row['CHANGE_VAR']] = row['NEW_VAL']
        covar_modeling_df.at[adj2_inx, 'RATE'] = row['NEW_VAL']/0.5
        # covar_modeling_df.at[adj2_inx+1, 'TIME'] - covar_modeling_df.at[adj2_inx, 'TIME']
        # covar_modeling_df.loc[adj2_inx]
        # covar_modeling_df.loc[adj2_inx+1]
        # raise ValueError
    else:
        continue
    # adj2_inx = adj2_row.index[0]
    # covar_modeling_df.at[adj2_inx, row['CHANGE_VAR']] = row['NEW_VAL']



# Wrong Data row 삭제
to_be_del = [
             # T/P 두번 기록, DV 움직임 비정상적
             {'ID':237,'TIME':417.15,'DV':27.2},
             {'ID':2638,'TIME':24.417,'DV':18.1},
             {'ID':371,'TIME':58,'DV':16.87},
             {'ID':1214,'TIME':663.283,'DV':3.8},

             # TAD>100 값
             {'ID':1034,'TIME':671.5,'DV':1.9},
             {'ID':1214,'TIME':759.1166667,'DV':5.3},
             {'ID':1630,'TIME':140.5833333,'DV':0.2},
             {'ID':1630,'TIME':212.7833333,'DV':0.3},

             # CWRE > 15
             {'ID': 353, 'TIME': 120, 'DV': 17.11},

             # 추가 중복데이터
             {'ID': 935, 'TIME': 126.85, 'AMT': 100.0},
             {'ID': 1278, 'TIME': 71.866667, 'AMT': 100.0},
             {'ID': 81, 'TIME': 3, 'AMT': 25.0},
             {'ID': 1679, 'TIME': 197.6333, 'AMT': 30.0},
             ]

to_be_del_inx = list()
for row in to_be_del:
    # raise ValueError
    # covar_modeling_df[(covar_modeling_df['ID'] == row['ID'])]
    # covar_modeling_df[(covar_modeling_df['ID'] == row['ID']) & (covar_modeling_df['TIME'].map(lambda x:round(x,1)) == round(row['TIME'],1))]
    if 'DV' in row:
        del_row = covar_modeling_df[(covar_modeling_df['ID']==row['ID'])&(covar_modeling_df['TIME'].map(lambda x:round(x,1)) == round(row['TIME'],1))&(covar_modeling_df['DV']==str(row['DV']))].copy()
    elif 'AMT' in row:
        # raise ValueError
        del_row = covar_modeling_df[(covar_modeling_df['ID']==row['ID'])&(covar_modeling_df['TIME'].map(lambda x:round(x,1)) == round(row['TIME'],1))&(covar_modeling_df['AMT']==str(row['AMT']))].copy()
    else:
        continue
    to_be_del_inx.append(del_row.index[0])
covar_modeling_df = covar_modeling_df[~(covar_modeling_df.index.isin(to_be_del_inx))].copy()

# DV==0, TIME==0 인 데이터는 모두 MDV==1 로 설정
for inx, row in covar_modeling_df.iterrows(): #break
    if (row['DV']=='0.0') and (row['TIME']==0):
        covar_modeling_df.at[inx,'MDV']=1

    # 한 대상자의 컬럼 하나 전체 수정
    if row['ID'] == 572:
        covar_modeling_df.at[inx, 'HT'] = 146.9

covar_modeling_df = covar_modeling_df.sort_values(['ID','TIME'], ignore_index=True)

##############
# raise ValueError
covar_modeling_df['CRCL'] = covar_modeling_df.apply(lambda x: ((140 - x['AGE']) * x['WT'] * 1) / (72 * x['CREATININE']) if x['SEX']==0 else  ((140 - x['AGE']) * x['WT'] * 0.85) / (72 * x['CREATININE']) ,axis=1)
no_ml_lab_pids = ['10042359', '10646256', '10844269', '11018784', '11505839', '11566803', '11666727', '11747950', '14258992', '15391683', '15957379', '18345874', '19380320', '19447177', '21971589', '26893691', '29819692', '32719154', '33246217', '33291185', '35679534', '39336264']
# inner_one_day_pids = ['25524226', '24961411', '10617861', '15499525', '22006666', '25389067', '19551122', '18347926', '14009765', '24311845', '30514726', '10190892', '10451629', '11895215', '16574645', '34728899', '11584452', '11845188', '24913861', '28650698', '23809356', '10006221', '13115597', '21210190', '26599247', '13158741', '14783830', '18146774', '26367192', '28080857', '18991455', '10474848', '24925408', '26668770', '35441638', '21249383', '10963177', '25351536', '11675891', '24785011', '32592760', '26809209', '23084924', '26948351']
inner_one_day_pids = ['25524226', '24961411', '10617861', '15499525', '22006666', '25389067', '19551122', '18347926', '24311845', '30514726', '14009765', '10190892', '10451629', '11895215', '16574645', '34728899', '11845188', '11584452', '24913861', '28650698', '23809356', '10006221', '21210190', '13115597', '26599247', '13158741', '14783830', '18146774', '26367192', '28080857', '18991455', '24925408', '10474848', '26668770', '35441638', '21249383', '10963177', '25351536', '24785011', '11675891', '32592760', '26809209', '23084924', '26948351', '10059364', '11667977', '14437586', '34023291', '10042359', '10646256', '10844269', '11505839', '11566803', '11666727', '11747950', '14258992', '15391683', '15957379', '19380320', '19447177', '21971589', '26893691', '29819692', '33246217', '33291185', '35679534', '39336264', '11018784', '18345874', '32719154', '10016994', '10228937', '10492062', '10535224', '10752283', '10820542', '11053826', '11428893', '11433600', '11497541', '11505459', '12019052', '12319497', '12687044', '13421634', '14035656', '14276767', '16021831', '17137032', '17938211', '18128891', '18395963', '19436230', '20198660', '20340959', '21435711', '21641345', '22344306', '22507400', '22945800', '24047733', '24425674', '24470991', '24608811', '24671198', '24902304', '25348422', '25351598', '25529726', '25578627', '25693300', '25736904', '26116332', '26347903', '26349031', '26903310', '26996381', '27243491', '27546657', '27682894', '27966097', '28099332', '28340566', '29030545', '29780662', '29819692', '29822232', '30193725', '30531163', '30790227', '31049427', '31090944', '31324515', '31521842', '31609665', '32111729', '32322536', '32369744', '32647057', '33018582', '33391294', '33880561', '34664320', '35536170', '35708146', '35791984', '35849131', '36075902', '36291937', '36463251', '36631775', '37560977', '37955638', '38489374', '39847005', '10819333']

covar_modeling_df = covar_modeling_df[~((covar_modeling_df['UID'].isin(inner_one_day_pids)) | (covar_modeling_df['UID'].isin(no_ml_lab_pids)))].copy()

for col in ['CRCL','WT']:
    covar_modeling_df[col]=covar_modeling_df[col].map(lambda x:round(x,6)*(10**6)/(10**6))
# raise ValueError
covar_modeling_df.to_csv(f'{modeling_covar_dir}/{drug}_modeling_df_covar.csv',index=False, encoding='utf-8')
covar_modeling_df.to_csv(f'{nonmem_dir}/{drug}_modeling_df_covar.csv',index=False, encoding='utf-8')

# (covar_modeling_df['MDV']==0).sum()

# covar_modeling_df = pd.read_csv(f'{output_dir}/amk_modeling_datacheck_covar.csv')
# covar_modeling_df['ID'].drop_duplicates()
# raise ValueError
####### NONMEM SDTAB

# nmsdtab_df = pd.read_csv(f"{nonmem_dir}/run/sdtab011",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
# modeling_datacheck_df = pd.read_csv(f"{output_dir}/amk_modeling_covar/amk_modeling_datacheck_covar.csv")
# nmsdtab_df = nmsdtab_df.merge(modeling_datacheck_df[['ID','UID']].drop_duplicates(['ID']), on=['ID'], how='left')
# nmsdtab_df['ID'] = nmsdtab_df['ID'].astype(int)
# # nmsdtab_df['TDM_YEAR'] = nmsdtab_df['TDM_YEAR'].astype(int)
# under_pred_df = nmsdtab_df[(nmsdtab_df['DV'] >= 10)&(nmsdtab_df['IPRED'] < 10)].copy()
# over_pred_df = nmsdtab_df[(nmsdtab_df['DV'] < 10)&(nmsdtab_df['IPRED'] >= 10)].copy()
#
#
# under_pred_df.to_csv(f"{modeling_covar_dir}/amk_modeling_underpred.csv",index=False, encoding='utf-8-sig')
# over_pred_df.to_csv(f"{modeling_covar_dir}/amk_modeling_overpred.csv",index=False, encoding='utf-8-sig')
# # under_pred_df['ID'].drop_duplicates()
#
# mis_pred_df = pd.concat([under_pred_df, over_pred_df])
# mis_pred_df['ID'].drop_duplicates()
# #
# # covar_modeling_df = covar_modeling_df[~(covar_modeling_df['ID'].isin(mis_pred_df['ID'].drop_duplicates()))]
# # # covar_modeling_df['SEX'] = covar_modeling_df['SEX'].map({'M':1,'F':2})
# # covar_modeling_df = covar_modeling_df.drop(['NAME'], axis=1)[['ID','TIME','TAD','DV','MDV','CMT','AMT','RATE','UID'] + list(covar_modeling_df.loc[:,'ALB':].iloc[:,:].columns)]
# # covar_modeling_df.to_csv(f"{nonmem_dir}/amk_modeling_df_covar_filt.csv",index=False, encoding='utf-8-sig')
#
#
# ####### NONMEM SDTAB
# # nmsdtab_df.columns
# modeling_covar_dir = f"{output_dir}/amk_modeling_datacheck"
# modeling_datacheck_df = pd.read_csv(f"{modeling_covar_dir}/amk_modeling_datacheck.csv")
# nmsdtab_df = pd.read_csv(f"{nonmem_dir}/run/sdtab011",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
# nmsdtab_df['ID'] = nmsdtab_df['ID'].astype(int)
#
# prep_conc_df = modeling_datacheck_df[modeling_datacheck_df['MDV']==0].copy()
# prep_conc_df['REC_REASON'] = prep_conc_df['REC_REASON'].replace(np.nan,'Vacant').map(lambda x:x.split('(')[0])
# nm_conc_df = nmsdtab_df[nmsdtab_df['MDV']==0]
# merge_conc_df = nm_conc_df.merge(prep_conc_df[['ID','TIME','REC_REASON']], on=['ID','TIME'], how='left')
#
# # merge_conc_df['REC_REASON'].unique()
# filt_gdf = merge_conc_df.copy()
# # filt_gdf.columns
# # filt_gdf = merge_conc_df[~merge_conc_df['REC_REASON'].isin(['오더비고반영', '결과비고반영(시간_분AP)', '결과비고반영(날짜_시간AP)', '결과비고반영(날짜_시간)'])].copy()
#
# # merge_conc_df['ID'].drop_duplicates()
# # # filt_gdf['ID'].drop_duplicates()
# # merge_conc_df[(merge_conc_df['REC_REASON']!='Vacant')]['ID'].drop_duplicates()
# # merge_conc_df[(merge_conc_df['REC_REASON']=='Vacant')]['ID'].drop_duplicates()
#
# import matplotlib.font_manager as fm
# plt.rc('font', family='Malgun Gothic')
# plt.rcParams['axes.unicode_minus'] = False
# sns.scatterplot(data=filt_gdf, x='IPRED',y='DV', hue='REC_REASON')
#
#
#
# # covar_modeling_df.drop(['NAME','DATETIME'], axis=1).applymap(lambda x:x if (type(x)==str) else np.nan).dropna(axis=1)
# # covar_modeling_df['DV'].unique()
# # covar_modeling_df.columns
# # covar_modeling_df['ID'].drop_duplicates()
# # str(covar_modeling_df.columns).replace("', '",' ').replace("',\n       '",' ')
#
# """
# # 결과에서 Covariate 일부분이 비어있는 이유
# # - basic_prep_lab_covar에서 각 사람마다의 lab 수치가 존재하는 날짜 기준으로 date range 를 생성했는데,
# # - 이게 order data (dosing date range)와 안 맞을 수 있다
# """
#
# # modeling_df[modeling_df['AST'].isna()]['UID'].unique()
#
