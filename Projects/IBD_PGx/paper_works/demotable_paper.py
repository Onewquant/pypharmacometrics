from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
paperwork_dir = f"{prj_dir}/paper_works"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

demo_res_table = list()
total_samples = np.nan
total_patients_with_samples = np.nan
total_patients = np.nan
rsid_df = pd.read_csv(f'{paperwork_dir}/rsid_dosage_matrix_with_alleles.csv')
# rsid_df = rsid_df.drop_duplicates(subset=['s']).copy()
rsid_df = rsid_df.dropna(subset=['UID'])
rsid_df['UID'] = rsid_df['UID'].map(lambda x:str(int(x)))
# set(md_df['UID'])-set(rsid_df['UID'])

for drug in ['all','infliximab', 'adalimumab']:
# for drug in ['all',]:
    # drug = 'infliximab'
    added_filename_str = '(for pda)'
    vacant_df = pd.DataFrame()
    # md_dict = {'integrated':vacant_df,'induction':vacant_df, 'maintenance':vacant_df}
    # md_demo_dict = {'integrated':vacant_df,'induction':vacant_df, 'maintenance':vacant_df}
    md_dict = {'integrated':vacant_df,}

    for mode_str in md_dict.keys():
        # md_df = pd.read_csv(f'{output_dir}/{drug}_{mode_str}_modeling_df.csv')
        if drug=='all':
            md_df1 = pd.read_csv(f'{output_dir}/modeling_df_covar/infliximab_{mode_str}_datacheck(covar){added_filename_str}.csv')
            md_df2 = pd.read_csv(f'{output_dir}/modeling_df_covar/adalimumab_{mode_str}_datacheck(covar){added_filename_str}.csv')
            common_columns = set(md_df1.columns).intersection(set(md_df2.columns))
            # md_df1['UID']
            # len(md_df1.columns)
            # raise ValueError
            # md_df.columns
            md_df = pd.concat([md_df1,md_df2])[list(md_df1.columns)].sort_values(['ID','DATETIME'])
            md_df['IBD_TYPE'] = md_df['IBD_TYPE'].map({'CD': 0, 'UC': 1})
            md_df['AGE'] = md_df.apply(lambda x: int((datetime.strptime(x['DATETIME'], '%Y-%m-%d') - datetime.strptime(x['AGE'], '%Y-%m-%d')).days / 365.25),axis=1)
            md_df['SEX'] = md_df['SEX'].map({'남': 0, '여': 1})
            md_df['ROUTE'] = md_df['ROUTE'].map({'IV': 1, 'SC': 2, '.': '.'})
        else:
            md_df = pd.read_csv(f'{nonmem_dir}/{drug}_{mode_str}_modeling_df_dayscale{added_filename_str}.csv')
        # raise ValueError
        # md_df = pd.read_csv(f'{output_dir}/{drug}_{mode_str}_datacheck.csv')

        md_df['Pediatric'] = (md_df['AGE'] < 19).map({False: 0, True: 1})
        md_df['BMI'] = md_df['WT']/((md_df['HT']/100)**2)
        # raise ValueError

        md_dict[mode_str] = md_df.copy()

        # raise ValueError

        mddemo_dict = dict()
        # mddemo_dict['Drug'] = drug
        # mddemo_dict['Mode'] = mode_str

        # mddemo_dict['Data spec, n (%)'] = ""
        if drug=='all':
            id_col = 'UID'
            md_df[id_col]=md_df[id_col].astype(str)
        else:
            id_col = 'ID'
        # raise ValueError

        mddemo_dict['Demographic categorical variables, n (%)'] = ""

        if drug=='all':
            total_patients = len(md_df[id_col].unique())
        mddemo_dict['Subtotal patients'] = f"{len(md_df[id_col].unique())} ({round(100 * len(md_df[id_col].unique())/total_patients,2)})"

        sex_series = md_df.drop_duplicates([id_col])['SEX'].copy()
        mddemo_dict['Female'] = f"{(sex_series == 1).sum()} ({round(((sex_series == 1).sum()) * 100 / len(sex_series), 2)})"

        pediatric_series = md_df.drop_duplicates([id_col])['Pediatric'].copy()
        mddemo_dict['Pediatric (Age < 19)'] = f"{(pediatric_series == 1).sum()} ({round(((pediatric_series == 1).sum()) * 100 / len(pediatric_series), 2)})"

        mddemo_dict['Demographic continuous variables, Mean (SD)'] = ""

        age_series = md_df[md_df['MDV']==1].drop_duplicates([id_col])['AGE'].copy()
        mddemo_dict['Age at the 1st Dose'] = f"{round(np.mean(age_series), 2)} ({round(np.std(age_series), 2)})"

        height_series = md_df.drop_duplicates([id_col])['HT'].copy()
        mddemo_dict['Height'] = f"{round(np.mean(height_series), 2)} ({round(np.std(height_series), 2)})"

        weight_series = md_df.drop_duplicates([id_col])['WT'].copy()
        mddemo_dict['Weight'] = f"{round(np.mean(weight_series), 2)} ({round(np.std(weight_series), 2)})"

        bmi_series = md_df.drop_duplicates([id_col])['BMI'].copy()
        mddemo_dict['BMI'] = f"{round(np.mean(bmi_series), 2)} ({round(np.std(bmi_series), 2)})"

        ## LAB findings

        mddemo_dict['Laboratory test, mean (SD)'] = ""
        labcol_dict = {'ALB':'Albumin (g/dL)', 'AST':'AST (IU/L)', 'ALT':'ALT (IU/L)', 'CRP':'CRP (mg/dL)', 'FCAL':'Fecal calprotectin (mg/kg)', 'CREATININE':'Serum creatinine (mg/dL)'}
        cols = ['Albumin (g/dL)', 'AST (IU/L)', 'ALT (IU/L)', 'CRP (mg/dL)', 'Fecal calprotectin (mg/kg)', 'Serum creatinine (mg/dL)']
        first_md_df = md_df.drop_duplicates(subset=[id_col])
        for col, res_col in labcol_dict.items():
            mean_val = first_md_df[col].mean()
            sd_val = first_md_df[col].std()

            # 소수점 자리수는 필요에 따라 조절 (여기선 2자리)
            value_str = f"{mean_val:.2f} ({sd_val:.2f})"
            mddemo_dict[res_col] = value_str

        ## Diagnosis

        ibdtype_series = md_df.drop_duplicates([id_col])['IBD_TYPE'].copy()
        mddemo_dict['Diagnosis, n (%)'] = ""
        mddemo_dict['CD'] = f"{(ibdtype_series == 0).sum()} ({round(((ibdtype_series == 0).sum()) * 100 / len(ibdtype_series), 2)})"
        mddemo_dict['UC'] = f"{(ibdtype_series == 1).sum()} ({round(((ibdtype_series == 1).sum()) * 100 / len(ibdtype_series), 2)})"

        ## Treatment phase


        mddemo_dict['Treatment phase, n (%)'] = ""
        maintenance_only_ids = md_df[(md_df['MDV'] == 0) & (md_df['TIME'] == 0) & (~md_df['DV'].isin(['0.0', '.']))][id_col]
        induction_df = md_df[~md_df[id_col].isin(maintenance_only_ids)].copy()
        subtotal_n = len(md_df[id_col].unique())
        induction_n = len(induction_df[id_col].unique())
        mddemo_dict['Whole phases'] = f"{induction_n} ({round(100 * induction_n / subtotal_n, 2)})"
        mddemo_dict['Maintenance only'] = f"{(subtotal_n - induction_n)} ({round(100 * (subtotal_n - induction_n) / subtotal_n, 2)})"

        ## Sampling
        mddemo_dict['Blood Sampling, n (%)'] = ""
        sampling_df = md_df[(md_df['MDV']!=1)&(~((md_df['DV']== '0.0') & (md_df['TIME'] == 0)))].copy()
        sampling_desc = sampling_df.groupby(id_col).agg(DV_COUNT=('DV', 'count'), DV_MIN=('DV', 'min')).reset_index(drop=False)
        if drug=='all':
            total_samples = len(sampling_df)
            total_patients_with_samples = len(sampling_desc)

        mddemo_dict['Total samples'] = f"{len(sampling_df)} ({round(100 * len(sampling_df) / total_samples, 2)})"
        mddemo_dict['Patients with samples'] = f"{len(sampling_desc)} ({round(100 * len(sampling_desc) / total_patients_with_samples, 2)})"
        # TL_cutoff = 5 if drug=='infliximab' else 8
        TL_cutoff = 5
        mddemo_dict[f'TL < {TL_cutoff}, n (%)'] = f"{np.sum(sampling_df['DV'].astype(float) < TL_cutoff)} ({round((np.sum(sampling_df['DV'].astype(float) < TL_cutoff)*100/len(sampling_df['DV'])),2)})"
        mddemo_dict[f'TL ≥ {TL_cutoff}, n (%)'] = f"{np.sum(sampling_df['DV'].astype(float) >= TL_cutoff)} ({round((np.sum(sampling_df['DV'].astype(float) >= TL_cutoff)*100/len(sampling_df['DV'])),2)})"

        # mddemo_dict['ATI'] = ""
        # not_na_ati_df = md_df[~md_df['INFATI'].isna()].copy()
        not_na_ati_df = md_df[md_df['ADA'] != 0].copy()
        not_na_ati_ids = not_na_ati_df.drop_duplicates([id_col])[id_col].reset_index(drop=True)
        ati_series = md_df.sort_values([id_col, 'ADA'], ascending=[True, False]).drop_duplicates([id_col])['ADA'].copy()
        mddemo_dict['Anti-drug antibody'] = f"{len(not_na_ati_ids)} ({round(len(not_na_ati_ids) * 100 / subtotal_n, 2)})"
        mddemo_dict['Samples/person'] = f"{round(np.mean(sampling_desc['DV_COUNT']), 2)} ({round(np.std(sampling_desc['DV_COUNT']), 2)})"


        # mddemo_dict['Patients with ATI positive, n (%)'] = f"{(ati_series >= 10).sum()} ({round(((ati_series >= 10).sum()) * 100 / len(not_na_ati_ids), 2)})"
        #
        # mddemo_dict['total ATI samples, n (%)'] = f"{len(not_na_ati_df)} ({round(100, 2)})"
        # mddemo_dict['high ATI samples, n (%)'] = f"{(not_na_ati_df['INFATI'] >= 10).sum()} ({round(((not_na_ati_df['INFATI'] >= 10).sum()) * 100 / len(not_na_ati_df), 2)})"
        # mddemo_dict['intermediate ATI samples, n (%)'] = f"{((not_na_ati_df['INFATI'] < 10)&(not_na_ati_df['INFATI'] > 2.5)).sum()} ({round(((not_na_ati_df['INFATI'] < 10)&(not_na_ati_df['INFATI'] > 2.5)).sum() * 100 / len(not_na_ati_df), 2)})"
        # mddemo_dict['LLOQ ATI samples, n (%)'] = f"{(not_na_ati_df['INFATI'] <= 2.5).sum()} ({round(((not_na_ati_df['INFATI'] <= 2.5).sum()) * 100 / len(not_na_ati_df), 2)})"

        # not_na_ati_df[(not_na_ati_df['INFATI'] < 10)&(not_na_ati_df['INFATI'] > 2.5)][['ID','INFATI']]
        # raise ValueError
        # not_na_ati_df = md_df[~md_df['INFATI'].isna()].copy()
        # not_na_ati_df[not_na_ati_df['INFATI'] > 2.5].drop_duplicates(['ID'])['INFATI']
        # not_na_ati_df[not_na_ati_df['INFATI'] >= 10]

        if (drug=='adalimumab'):
            mddemo_dict['Anti-drug antibody'] = f"0 (0.0)"
            # mddemo_dict['ATI high (value >= 10), n (%)'] = f"0 (0.0)"
            # mddemo_dict['ATI intermediate (10 > value >= 2.5), n (%)'] = f"0 (0.0)"
            # mddemo_dict['ATI negative (value < 2.5), n (%)'] = f"0 (0.0)"


        ## Dosing
        mddemo_dict['Dosing route, n (%)'] = ""

        dosing_df = md_df[md_df['MDV']==1].copy()
        pop_form_type_df = dosing_df['CMT'].value_counts()
        ind_form_type_df = dosing_df.groupby(id_col)['CMT'].value_counts().unstack(fill_value=0).reset_index(drop=False)
        ind_form_type_df.columns.name=None

        # if drug=='all':
        mddemo_dict['Total routes'] = f"{pop_form_type_df.sum()} ({round(pop_form_type_df.sum() * 100 / pop_form_type_df.sum(), 2)})"
        mddemo_dict['Subcutaneous'] = f"{pop_form_type_df[1]} ({round(pop_form_type_df[1] * 100 / pop_form_type_df.sum(), 2)})"
        try: mddemo_dict['Intravenous'] = f"{pop_form_type_df[2]} ({round(pop_form_type_df[2] * 100 / pop_form_type_df.sum(), 2)})"
        except: mddemo_dict['Intravenous'] = f"0 ({round(0 * 100 / pop_form_type_df.sum(), 2)})"

        res_df = pd.DataFrame([mddemo_dict]).T
        res_df.columns = [drug]

        # sampling_desc = sampling_df.groupby('ID').agg(DV_COUNT=('DV', 'count'), DV_MIN=('DV', 'min')).reset_index(drop=False)
        #
        # mddemo_dict['Samples/person, mean (SD)'] = f"{round(np.mean(sampling_desc['DV_COUNT']), 2)} ({round(np.std(sampling_desc['DV_COUNT']), 2)})"
        # mddemo_dict['TL < 3, n (%)'] = f"{np.sum(sampling_df['DV'].astype(float) < 3)} ({round((np.sum(sampling_df['DV'].astype(float) < 3)*100/len(sampling_df['DV'])),2)})"
        # mddemo_dict['TL ≥ 3, n (%)'] = f"{np.sum(sampling_df['DV'].astype(float) >= 3)} ({round((np.sum(sampling_df['DV'].astype(float) >= 3)*100/len(sampling_df['DV'])),2)})"

        res_df = pd.DataFrame([mddemo_dict]).T.reset_index(drop=False)
        res_df.columns = ['Characteristics', f"{drug}_{mode_str}"]
        demo_res_table.append(res_df)
# demo_res_table.columns
# for len(demo_res_table)
demo_res_table = (demo_res_table[0].merge(demo_res_table[1], on=['Characteristics'], how='left')).merge(demo_res_table[2], on=['Characteristics'], how='left')
demo_res_table.to_csv(f'{paperwork_dir}/demographic_table(paper).csv',index=False, encoding='utf-8-sig')