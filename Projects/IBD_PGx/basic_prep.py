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
for finx, fpath in enumerate(order_files):

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    fdf = pd.read_excel(fpath)

    # fdf.columns


    fdf = fdf[~fdf['Acting'].isna()].copy()
    dose_df = fdf[fdf['처방지시'].map(lambda x: (('adalimumab' in x.lower()) or ('infliximab' in x.lower())) and ('quantification' not in x.lower()))].copy()
    dose_df['처방지시비고'] = dose_df['처방지시'].map(lambda x:x.split(' : ')[-1] if len(x.split(' : '))>1 else '')


    dose_df['ID'] = pid
    dose_df['NAME'] = pname


    #### 약국_검사가 NA인 것은 제외하고 진행함 (추후 필요시 추가)

    # dose_df[dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test.xlsx",index=False)
    # dose_df[~dose_df['약국_검사'].isna() & dose_df['처방지시'].map(lambda x: '(infliximab)' in x.lower())].to_excel(f"{output_dir}/test2.xlsx",index=False)
    dose_df = dose_df[~dose_df['약국_검사'].isna()].copy()

    #### 성분명이 IBD Biologics 인 경우만으로 필터링 (Infliximab, Adalimumab)

    regex_pattern = r'\(infliximab|adalimumab\)'
    dose_df = dose_df[dose_df['처방지시'].map(lambda x: bool(re.search(regex_pattern, x, flags=re.IGNORECASE)))].copy()

    dose_df['DT1'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[0].split('[')[-1].replace(' ','T'))
    dose_df['DT2'] = dose_df['약국_검사'].map(lambda x:x.split(']   ')[-1].split('[')[-1].replace(' ','T')[:-3])
    dose_df['DATETIME'] = dose_df[['DT1','DT2']].min(axis=1)
    dose_df['ETC_INFO'] = dose_df['처방지시비고'].copy()
    dose_df['DRUG'] = dose_df['처방지시'].map(lambda x: re.search(regex_pattern, x, flags=re.IGNORECASE).group().lower().replace('(','').replace(')',''))

    # dose_df[['ID','NAME','DATETIME','DRUG','ETC_INFO']]

    # dose_df['DRUG'].map(lambda x:x.split(' : ')[0]).unique()
    drug_order_set = drug_order_set.union(set(dose_df['DRUG'].drop_duplicates()))
    # dose_df = dose_df[dose_df]



# ['비고','처방지시', '발행처', '발행의', '수납', '약국/검사', '주사시행처', 'Acting', '변경의']
#
# len(order_files)
#
#
# df = pd.read_csv(f'{resource_dir}/glpharma_CONC.csv')
# seq_df = pd.read_csv(f'{prj_dir}/glpharma_SEQUENCE.csv')