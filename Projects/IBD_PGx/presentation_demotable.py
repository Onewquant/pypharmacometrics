from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

# demo_res_table = list()
# for drug in ['infliximab', 'adalimumab']:
drug = 'infliximab'
vacant_df = pd.DataFrame()
# md_dict = {'integrated':vacant_df,'induction':vacant_df, 'maintenance':vacant_df}
# md_demo_dict = {'integrated':vacant_df,'induction':vacant_df, 'maintenance':vacant_df}
md_dict = {'integrated':vacant_df,}
md_demo_dict = {'integrated':vacant_df,}
total_patients = 140

for mode_str in md_dict.keys():
    md_df = pd.read_csv(f'{output_dir}/{drug}_{mode_str}_modeling_df.csv')
    # md_df = pd.read_csv(f'{output_dir}/{drug}_{mode_str}_datacheck.csv')
    md_dict[mode_str] = md_df.copy()

    mddemo_dict = dict()
    mddemo_dict['Drug'] = drug
    # mddemo_dict['Mode'] = mode_str

    mddemo_dict['Data spec'] = ""
    mddemo_dict['Subtotal(dose,TL exists), n (%)'] = f"{len(md_df['ID'].unique())} ({round(100 * len(md_df['ID'].unique())/total_patients,2)})"
    # if mode_str=='integrated':
    induction_df = pd.read_csv(f'{output_dir}/{drug}_induction_modeling_df.csv')
    subtotal_n = len(md_df['ID'].unique())
    induction_n = len(induction_df['ID'].unique())
    mddemo_dict['Whole phases'] = f"{induction_n} ({round(100 * induction_n/subtotal_n,2)})"
    mddemo_dict['Maintenance only'] = f"{(subtotal_n-induction_n)} ({round(100 * (subtotal_n-induction_n)/subtotal_n,2)})"

    mddemo_dict['Demographics'] = ""
    age_series = md_df[md_df['MDV']=='1'].drop_duplicates(['ID'])['AGE'].copy()
    mddemo_dict['AGE at the 1st Dose, mean (SD)'] = f"{round(np.mean(age_series), 2)} ({round(np.std(age_series), 2)})"

    sex_series = md_df.drop_duplicates(['ID'])['SEX'].copy()
    mddemo_dict['Female, n (%)'] = f"{(sex_series==1).sum()} ({round(((sex_series==1).sum()) * 100 /len(sex_series), 2)})"

    height_series = md_df.drop_duplicates(['ID'])['HT'].copy()
    mddemo_dict['Height(recent), mean (SD)'] = f"{round(np.mean(height_series), 2)} ({round(np.std(height_series), 2)})"

    weight_series = md_df.drop_duplicates(['ID'])['WT'].copy()
    mddemo_dict['Weight(recent), mean (SD)'] = f"{round(np.mean(weight_series), 2)} ({round(np.std(weight_series), 2)})"

    ibdtype_series = md_df.drop_duplicates(['ID'])['IBD_TYPE'].copy()
    mddemo_dict['Diagnosis, n (%)'] = ""
    mddemo_dict['CD'] = f"{(ibdtype_series == 1).sum()} ({round(((ibdtype_series == 1).sum()) * 100 / len(ibdtype_series), 2)})"
    mddemo_dict['UC'] = f"{(ibdtype_series == 2).sum()} ({round(((ibdtype_series == 2).sum()) * 100 / len(ibdtype_series), 2)})"

    ## Sampling
    mddemo_dict['Blood Sampling'] = ""
    sampling_df = md_df[(md_df['MDV']!='1')&(~((md_df['DV']== '0.0') & (md_df['TIME'] == 0)))].copy()
    sampling_desc = sampling_df.groupby('ID').agg(DV_COUNT=('DV', 'count'), DV_MIN=('DV', 'min')).reset_index(drop=False)

    mddemo_dict['Total samples, n'] = f"{len(sampling_df)}"
    mddemo_dict['Samples/person, mean (SD)'] = f"{round(np.mean(sampling_desc['DV_COUNT']), 2)} ({round(np.std(sampling_desc['DV_COUNT']), 2)})"
    mddemo_dict['TL < 3, n (%)'] = f"{np.sum(sampling_df['DV'].astype(float) < 3)} ({round((np.sum(sampling_df['DV'].astype(float) < 3)*100/len(sampling_df['DV'])),2)})"
    mddemo_dict['TL ≥ 3, n (%)'] = f"{np.sum(sampling_df['DV'].astype(float) >= 3)} ({round((np.sum(sampling_df['DV'].astype(float) >= 3)*100/len(sampling_df['DV'])),2)})"


    mddemo_dict['ATI'] = ""
    not_na_ati_ids = md_df[~md_df['INFATI'].isna()].drop_duplicates(['ID'])['ID'].reset_index(drop=True)
    ati_series = md_df.sort_values(['ID', 'INFATI'], ascending=[True, False]).drop_duplicates(['ID'])['INFATI'].copy()
    mddemo_dict['Patients with measured ATI, n (%)'] = f"{len(not_na_ati_ids)} ({round(len(not_na_ati_ids) * 100 / subtotal_n, 2)})"
    mddemo_dict['ATI positive (>= 10), n (%)'] = f"{(ati_series >= 10).sum()} ({round(((ati_series >= 10).sum()) * 100 / len(not_na_ati_ids), 2)})"

    ## Dosing
    mddemo_dict['Dosing'] = ""

    dosing_df = md_df[md_df['MDV']=='1'].copy()
    pop_form_type_df = dosing_df['CMT'].value_counts()
    ind_form_type_df = dosing_df.groupby('ID')['CMT'].value_counts().unstack(fill_value=0).reset_index(drop=False)
    ind_form_type_df.columns.name=None

    mddemo_dict['SC, n (%)'] = f"{pop_form_type_df[1]} ({round(pop_form_type_df[1] * 100 / pop_form_type_df.sum(), 2)})"
    mddemo_dict['IV, n (%)'] = f"{pop_form_type_df[2]} ({round(pop_form_type_df[2] * 100 / pop_form_type_df.sum(), 2)})"

    res_df = pd.DataFrame([mddemo_dict]).T
    res_df.columns = [drug]

    # sampling_desc = sampling_df.groupby('ID').agg(DV_COUNT=('DV', 'count'), DV_MIN=('DV', 'min')).reset_index(drop=False)
    #
    # mddemo_dict['Samples/person, mean (SD)'] = f"{round(np.mean(sampling_desc['DV_COUNT']), 2)} ({round(np.std(sampling_desc['DV_COUNT']), 2)})"
    # mddemo_dict['TL < 3, n (%)'] = f"{np.sum(sampling_df['DV'].astype(float) < 3)} ({round((np.sum(sampling_df['DV'].astype(float) < 3)*100/len(sampling_df['DV'])),2)})"
    # mddemo_dict['TL ≥ 3, n (%)'] = f"{np.sum(sampling_df['DV'].astype(float) >= 3)} ({round((np.sum(sampling_df['DV'].astype(float) >= 3)*100/len(sampling_df['DV'])),2)})"



    res_df = pd.DataFrame([mddemo_dict]).T.reset_index(drop=False)
    res_df.columns = ['Characteristics', 'Value']


res_df