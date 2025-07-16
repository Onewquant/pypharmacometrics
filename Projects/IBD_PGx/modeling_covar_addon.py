from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

# rdf = pd.read_csv(f'{output_dir}/infliximab_induction_modeling_df.csv')
# rdf['ALB'].median()
# rdf['CRP'].median()
# rdf['CREATININE'].median()

## DEMO Covariates Loading

demo_df = pd.read_csv(f"{resource_dir}/demo/IBD_PGx_demo.csv")
demo_df = demo_df.rename(columns={'EMR ID':'UID','birthdate':'AGE','sex':'SEX','name':'NAME'})
demo_df = demo_df[['UID','AGE', 'SEX']].copy()
demo_df['UID'] = demo_df['UID'].astype(str)

demo2_df = pd.read_csv(f"{resource_dir}/demo2/IBD_PGx_demo2.csv")
demo2_df = demo2_df.rename(columns={'EMR ID':'UID','height':'HT','weight':'WT','bmi':'BMI','bsa':'BSA'})
demo2_df = demo2_df[['UID','HT', 'WT', 'BMI','BSA']].copy()
demo2_df['UID'] = demo2_df['UID'].astype(str)

## LAB Covariates Loading

totlab_df = pd.read_csv(f"{output_dir}/lab_df.csv")
totlab_df = totlab_df[['UID', 'DATETIME', 'Albumin', 'AST', 'AST(GOT)', 'ALT', 'ALT(GPT)', 'CRP', 'Calprotectin (Serum)', 'Calprotectin (Stool)', 'eGFR-CKD-EPI', 'Cr (S)', 'Creatinine','Anti-Infliximab Ab [정밀면역검사] (정량)']].copy()
totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'] = totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'].map(lambda x: float(re.findall(r'\d*\.\d+|\d+',str(x))[0]) if str(x)!='nan' else np.nan)
for c in list(totlab_df.columns)[2:]:  # break
    totlab_df[c] = totlab_df[c].map(lambda x: x if type(x) == float else float(re.findall(r'[\d]+.*[\d]*', str(x))[0]))
totlab_df['UID'] = totlab_df['UID'].astype(str)
totlab_df['AST'] = totlab_df[['AST', 'AST(GOT)']].max(axis=1)
totlab_df['ALT'] = totlab_df[['ALT', 'ALT(GPT)']].max(axis=1)
totlab_df['CREATININE'] = totlab_df[['Cr (S)', 'Creatinine']].max(axis=1)
# set(totlab_df['Anti-Infliximab Ab [정밀면역검사] (정량)'].unique())
totlab_df = totlab_df.rename(columns={'Albumin': 'ALB', 'Calprotectin (Stool)': 'CALPRTSTL', 'Calprotectin (Serum)': 'CALPRTSER', 'Anti-Infliximab Ab [정밀면역검사] (정량)':'INFATI'})
totlab_df = totlab_df[['UID', 'DATETIME', 'ALB', 'AST', 'ALT', 'CRP', 'CALPRTSTL', 'CREATININE','INFATI']].copy()
# totlab_df = totlab_df[['UID', 'DATETIME', 'ALB', 'AST', 'ALT', 'CRP', 'CALPRTSTL', 'CREATININE']].copy()
for c in list(totlab_df.columns)[2:]:  # break
    totlab_df[c] = totlab_df[c].map(lambda x: x if type(x) == float else float(re.findall(r'[\d]+.*[\d]*', str(x))[0]))

totlab_df = totlab_df.drop_duplicates(['UID','DATETIME'])

# [c for c in totlab_df.columns.unique() if 'ada' in c.lower()]

# totlab_df[~totlab_df['ALT'].isna()]['ALT']

## Modeling Data Loading
for drug in ['infliximab',]:
# for drug in ['infliximab','adalimumab']:
#     for mode_str in ['integrated','induction']:
    for mode_str in ['integrated',]:
    # for mode_str in ['integrated']:
        modeling_df = pd.read_csv(f'{output_dir}/{drug}_{mode_str}_datacheck.csv')
        modeling_df['UID']= modeling_df['UID'].astype(str)
        modeling_df['DATETIME'] = modeling_df['DATETIME'].map(lambda x:x.split('T')[0])


        # raise ValueError
        """
        # list(totlab_df.columns)
        # ['UID', 'DATETIME', '1,25-VitD3', '25-OH VIT D TOTAL', '25-OH Vit. D3(초진용)', '25-OH Vit.D (D3/D2)', '25-OH Vit.D (Total)', '25-OH Vit.D2', '25-OH Vit.D3', '3-methoxytyramine (Plasma)', '5-HIAA (24h urine)', 'ABO', 'ACL, IgG', 'ACL, IgM', 'ACTH', 'ACTHSTbase', 'ADA', 'AFP', 'AFP-L3(%)', 'ALT', 'ALT(GPT)', 'AMA', 'ANC', 'ANCA', 'ARR', 'ASCA IgA & IgG', 'ASO', 'AST', 'AST(GOT)', 'AT III', 'Ab scr', 'Acetoacetate', 'Adalimumab Quantification', 'Albumin', 'Aldo(B)', 'Alk. phos', 'Alk. phos.', 'Ammo', 'Amylase(Flu)', 'Amylase(S)', 'Anion gap', 'Anti - ccp', 'Anti Mullerian Hormone (AMH)', 'Anti RNP Ab', 'Anti Sm Antibody', 'Anti ds DNA', 'Anti-HBs', 'Anti-HCV', 'Anti-Infliximab Ab [정밀면역검사] (정량)', 'Anti-LKM', 'Anti-TSH receptor', 'ApoA1', 'ApoB', 'AtypicalLc', 'B2-MG(S)', 'B2-MG(U)', 'BE', 'BIL', 'BLD', 'BNP', 'BST', 'BUN', 'Bact.', 'Bacteria', 'Band.neut.', 'Basophil', 'Beta hydroxybutyric acid', 'Blast', 'C-Peptide(S)', 'C. difficile GDH', 'C. difficile toxin', 'C.diffcile toxin', 'C3', 'C4', 'CA 125', 'CA 15-3', 'CA 19-9', 'CBC (em) (diff), RDW제외', 'CEA', 'CK', 'CK(CPK)', 'CK-MB (em)', 'CMV Ag', 'CMV IgM', 'CORT30', 'CORT60', 'CORTbasal', 'CPEPSTbase', 'CRP', 'CTx (C-telopeptide)', 'Ca', 'Ca, total', 'Calcium', 'Calprotectin (Serum)', 'Calprotectin (Stool)', 'Cast', 'Casts', 'Cerulo', 'Chol.', 'Chromogranin A', 'Cl', 'Cl(U) (em)', 'Codfish IgE (F3)', 'Collagen 4', 'Color', 'Cortisol', 'Cortisol(S)', 'Cotinine', "Cow's milk IgE (F2)", 'Cr (S)', 'Cr (U)', 'Cr(u)', 'Creatinine', 'Cryoglobu.', 'Crystal', 'Crystals', 'Cu (24h Urine)', 'Cu (Serum)', 'Cystatin C', 'D-dimer', 'D-dimer(em)', 'D. Bil.', 'D. farinae IgE (D2)', 'D. farinae IgG4 (D2)', 'D. pteronyssinus IgE (D1)', 'D. pteronyssinus IgG4 (D1)', 'DHEA-S', 'DRBC', 'Delta ratio', 'Dopamine', 'ECP', 'ESR', 'Egg white IgE (F1)', 'Entero C', 'Eos.count', 'Eosinophil', 'Epinephrine', 'Erythropoiet', 'Erythropoietin', 'Estradiol', 'F. acid', 'FANA', 'FANA Titer', 'FBS', 'FDP정량', 'FSH', 'FT3', 'Factor 10', 'Factor 11', 'Factor 12', 'Factor 2', 'Factor 5', 'Factor 8', 'Factor 9', 'Ferritin', 'Fibrinogen', 'Folate', 'Free Fatty Acid', 'Free PSA', 'Free PSA%', 'Free T3', 'Free T4', 'GLU', 'Gastrin', 'Glu', 'Glucagon', 'Glucose', 'H. pylori IgG Ab', 'H.pylo IgG', 'HA', 'HAV Ab IgG', 'HAV Ab(IgG)', 'HAV Ab(IgM)', 'HAVAb IgM', 'HBcAb IgM', 'HBcAb Total (IgM+IgG) (진단검사의학)', 'HBcAb(IgG)', 'HBeAb', 'HBeAg', 'HBs Ag', 'HBsAb', 'HBsAg', 'HCO3-', 'HCV Ab', 'HDL Chol.', 'HE4 (Human Epididymis Protein 4)', 'HFR', 'HGH', 'HIV Ag/Ab', 'HIV Ag/Ab (건강검진용)', 'Hapto', 'Hb', 'HbA1c-IFCC', 'HbA1c-IFCC (건강증진센타용)', 'HbA1c-NGSP', 'HbA1c-NGSP (건강증진센타용)', 'HbA1c-eAG', 'HbA1c-eAG (건강증진센타용)', 'Hct', 'Helm,ova', 'Human Hb', 'I/Ⅱratio', 'IGF 1', 'IGFBP 3', 'IMA (Ischemia Modified Albumin Test) (em)', 'IRF', 'Ig A', 'Ig G', 'IgE(Total)', 'IgG sub 4', 'IgG sub1', 'IgG sub2', 'IgG sub3', 'IgG sub4', 'IgM', 'Imm.cell', 'Imm.lympho', 'Imm.mono', 'Infliximab Quantification', 'Influenza A&B Ag', 'Insulin', 'Insulin Ab', 'Interleukin-6', 'Iron', 'K', 'K (U) (em)', 'KET', 'Ketone (Beta-hydroxybutyrate)', 'LA', 'LD', 'LD(LDH)', 'LD-R', 'LD-R(B)', 'LDL Chol.', 'LH', 'LUC', 'Lactate', 'Lactic', 'Li', 'Lipase', 'Lp(a)', 'Lymph', 'Lympho', 'Lymphocyte', 'M/C ratio', 'MCH', 'MCHC', 'MCV', 'MFR', 'MN(Mononuclear cell)', 'MPV', 'MTX', 'Metamyelo', 'Metanephrine (Plasma)', 'Mg', 'Mic-Ab', 'MicroAlb', 'Mixing test (PT, aPTT 제외)', 'Monocyte', 'Mycopl. Ab', 'Mycoplasma Ab IgG', 'Mycoplasma Ab IgM', 'Myelocyte', 'Myoglobin', 'Myoglobin (em)', 'N. fat', 'NCV2019 응급용 선별검사 (교수용)', 'NCV2019 입원 선별 1단계 [분자진단]', 'NCV2019 입원 선별 2단계 [분자진단]', 'NIT', 'NMP 22', 'NSE', 'Na', 'Na(U) (em)', 'Norepinephrine', 'Normetanephrine (Plasma)', 'Normoblast', 'O2 CT', 'O2 SAT', 'Osmo-S', 'Osmo-U', 'Osteocalcin', 'Other', 'Others', 'O₂SAT', 'P', 'P. cyst', 'P. troph', 'P1NP', 'PCT', 'PDW', 'PIVKA II', 'PLT', 'PMN(Polymorphonuclear cell)', 'PP2', 'PRO', 'PSA', 'PSA (건증용)', 'PT %', 'PT % (MIX)', 'PT INR', 'PT INR (MIX)', 'PT sec', 'PT sec (MIX)', 'PTH', 'PTH(intact)', 'Pb(Blood)', 'Peanut IgE (F13)', 'PepsinogenⅠ', 'PepsinogenⅡ', 'Phosphorus', 'Pl. Hb', 'Plas.cell', 'Platelet', 'Poly', 'PreAlb', 'Prealbumin', 'Procalcitonin', 'Prolactin', 'Promyelo', 'Prostate Health Index', 'Protein', 'Protein C activity', 'Protein S activity', 'Protein/Creatinine ratio', 'Protozo', 'Pyruvate', 'RBC', 'RDW(CV)', 'RDW(SD)', 'RF', 'ROMA (postmenopausal)', 'ROMA (premenopausal)', 'RPR (VDRL, Auto) (serum)', 'RPR 정량 (serum)', 'RT', 'RTE', 'Renin(B)', 'Reticulocyte', 'Rh D', 'S.G', 'SAA (Serum Amyloid A) (em)', 'SARS-CoV-2 Ag(외래, 응급실, 중환자실 전용)', 'SCC(TA-4)', 'SG', 'SHBG', 'SQE', 'SS-A/Ro(52) Ab', 'SS-A/Ro(60) Ab', 'SS-B/La Ab', 'Salt intake', 'Seg.neut.', 'Selenium', 'Smooth muscle Ab', 'Sodium (random urine)', 'Soybean IgE (F14)', 'Sperm', 'T CO2', 'T. Bil.', 'T. Protein', 'T.B', 'T.P Ab Total (IgM+IgG) 정밀면역', 'T3', 'TCO2', 'TG', 'TIBC', 'TRE', 'TSH', 'Testosterone', 'Total IgE (진단검사의학과)', 'Total ketone', 'Toxo IgM', 'Toxocariasis (ELISA, IgG antibody)', 'Toxopla.Ab', 'Transferrin', 'Troponin I (em) [ng/mL]', 'Troponin I (em) [pg/mL]', 'Troponin T [ng/mL]', 'Troponin T [pg/mL]', 'Tryptase', 'Tur', 'Turbidity', 'UA', 'UIBC', 'URO', 'Uric acid', 'Urine S. pneumoniae Ag', 'VWF rel Ag', 'VWFrist co', 'VZV Ab IgG', 'Valproic', 'Vanco', 'Vitamin B12', 'Vitamin B₁ (Thiamin)', 'WBC', 'WBC(s)', 'Wheat IgE (F4)', 'X-mat', 'Yeast like organism', 'Zn (serum)', 'aPTT', 'aPTT (MIX)', 'alpha1-Antitrypsin (stool)', 'd1', 'd2', 'eGFR (CKD-EPI Cr-Cys)', 'eGFR (CKD-EPI Cys)', 'eGFR-CKD-EPI', 'eGFR-Cockcroft-Gault', 'eGFR-MDRD', 'eGFR-Schwartz(소아)', 'em.WBC', 'epine(U)', 'hCG', 'hsCRP', 'i6', 'iCa', 'iMg', 'metanephrine(U)', 'mito Ab', 'norepine(U)', 'p2PSA', 'pCO₂', 'pH', 'pO₂', 'proBNP (em)', 'tHcyst', 'u.RBC', 'u.WBC', 'β2 GP1-IgG', 'β2 GP1-IgM', 'γ-GT', '검체처리 1단계(채혈TUBE1~5개)', '마약선별검사', '모발 중금속 및 미네랄 40종 검사', '연구검체보관', '연구용채혈+검체처리1단계(채혈TUBE1~5개)', '연구용채혈+검체처리2단계', '절대단구수', '절대림프구수']
        """
        # len(modeling_df)

        modeling_df = modeling_df.merge(totlab_df, on=['UID','DATETIME'], how='left')
        modeling_df = modeling_df.merge(demo_df, on=['UID'], how='left')
        modeling_df = modeling_df.merge(demo2_df, on=['UID'], how='left')
        # modeling_df['UID'].iloc[0]
        # demo_df['UID'].iloc[0]

        ## Covariates의 NA value 처리 (ffill 먼저 시도, 없으면 bfill, 그것도 없으면 전체의 median 값)

        md_df_list = list()
        for md_inx,md_df in modeling_df.groupby(['UID']):
            md_df = md_df.sort_values(['DATETIME']).fillna(method='ffill').fillna(method='bfill')
            md_df_list.append(md_df)
        modeling_df = pd.concat(md_df_list).reset_index(drop=True)

        modeling_df.fillna(modeling_df.median(numeric_only=True), inplace=True)



        ## Covariates의 NA value 처리 (ffill 먼저 시도, 없으면 bfill, 그것도 없으면 전체의 median 값)

        # md_df_list = list()
        # for md_inx, md_df in modeling_df.groupby(['UID']):
        #     # md_df = md_df.sort_values(['DATETIME']).fillna(method='ffill')
        #     md_df = md_df.sort_values(['DATETIME'])
        #     md_df_list.append(md_df)
        # modeling_df = pd.concat(md_df_list).reset_index(drop=True)

        # modeling_df.fillna('.', inplace=True)

        # raise ValueError

        ## Modeling Data Saving
        # data_check_cols = ['ID','UID','NAME','DATETIME','TIME','DV','MDV','AMT','DUR','CMT','IBD_TYPE','ADDED_ADDL'] + list(modeling_df.loc[:,'DRUG':].iloc[:,1:].columns)
        # modeling_df[data_check_cols].to_csv(f'{output_dir}/{drug}_{mode_str}_datacheck_covar.csv', index=False, encoding='utf-8-sig')
        right_covar_col = 'TIME(WEEK)' if 'TIME(WEEK)' in modeling_df.columns else 'DRUG'
        datacheck_cols = ['ID',	'UID', 'NAME', 'DATETIME','TIME(WEEK)','TIME(DAY)','TIME', 'DV', 'MDV', 'AMT', 'DUR', 'CMT', 'IBD_TYPE', 'ROUTE','DRUG', 'ADDED_ADDL'] + list(modeling_df.loc[:,right_covar_col:].iloc[:,1:].columns)
        modeling_df[datacheck_cols].to_csv(f'{output_dir}/{drug}_{mode_str}_datacheck(covar).csv', index=False, encoding='utf-8-sig')

        modeling_cols = ['ID','TIME','DV','MDV','AMT','DUR','CMT','IBD_TYPE'] + list(modeling_df.loc[:,right_covar_col:].iloc[:,1:].columns)

        modeling_df['IBD_TYPE'] = modeling_df['IBD_TYPE'].map({'CD':1,'UC':2})
        modeling_df['AGE'] = modeling_df.apply(lambda x: int((datetime.strptime(x['DATETIME'],'%Y-%m-%d') - datetime.strptime(x['AGE'],'%Y-%m-%d')).days/365.25), axis=1)
        modeling_df['SEX'] = modeling_df['SEX'].map({'남':1,'여':2})

        modeling_df = modeling_df[modeling_cols].sort_values(['ID','TIME'], ignore_index=True)

        modeling_input_line = str(list(modeling_df.columns)).replace("', '"," ")

        print(f"Mode: {mode_str} / {modeling_input_line}")

        # if mode_str=='maintenance':  # Time decrease가 생김... 확인 !
        #     modeling_df["prev_time"] = modeling_df.groupby("ID")["TIME"].shift(1)
        #     modeling_df["time_decrease"] = (modeling_df["TIME"] < modeling_df["prev_time"])
        #     problem_rows = modeling_df[modeling_df["time_decrease"]]
        #     print(problem_rows)

        # modeling_df['CALPRTSTL'].median()
        # modeling_df['CREATININE'].median()
        # modeling_df['ALB'].median()
        # modeling_df['CRP'].median()
        # raise ValueError
        # modeling_df['A_0FLG'] = (modeling_df['ID'].shift(1)!=modeling_df['ID'])*1


        modeling_df.to_csv(f'{output_dir}/{drug}_{mode_str}_modeling_df.csv',index=False, encoding='utf-8-sig')
        modeling_df[~((modeling_df['TIME'] == 0) & (modeling_df['DV'] == '0.0'))].copy().to_csv(f'{output_dir}/{drug}_{mode_str}_modeling_df_without_zero_dv.csv',index=False, encoding='utf-8-sig')
        modeling_df['TIME']= modeling_df['TIME']/24
        modeling_df['DUR'] = modeling_df['DUR'].map(lambda x: float(x)/24 if x!='.' else x)
        modeling_df.to_csv(f'{output_dir}/{drug}_{mode_str}_modeling_df_dayscale.csv',index=False, encoding='utf-8')



"""
# 결과에서 Covariate 일부분이 비어있는 이유
# - basic_prep_lab_covar에서 각 사람마다의 lab 수치가 존재하는 날짜 기준으로 date range 를 생성했는데, 
# - 이게 order data (dosing date range)와 안 맞을 수 있다 
"""

# modeling_df[modeling_df['AST'].isna()]['UID'].unique()

