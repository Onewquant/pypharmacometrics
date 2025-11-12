from tools import *
from pynca.tools import *

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
# prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Orders
lab_set = set()

sampling_files = glob.glob(f'{resource_dir}/sampling/{prj_name}_sampling(*).xlsx')
sampling_result_df = list()
no_sampling_list = list()
for finx, fpath in enumerate(sampling_files):  # break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]


    # if pid in ("19447177",):
    #     raise ValueError
    # else:
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
    conc_df = fdf[fdf['검사명'].isin(['Amikacin 농도'])].copy()  # Amikacin 농도 측정 항목만 추출
    conc_df['검사보고일시'] = conc_df['검사보고일시'].replace('_x000D_', np.nan)
    # conc_df = conc_df[(~conc_df['검사채혈일시'].isna()) & (~conc_df['검사시행일시'].isna()) & (~conc_df['검사보고일시'].isna()) & (~conc_df['처방발행일시'].isna())].copy()
    conc_df = conc_df[(~conc_df['검사시행일시'].isna()) & (~conc_df['검사보고일시'].isna()) & (~conc_df['처방발행일시'].isna())].copy()

    # samp_na = conc_df[conc_df['검사채혈일시'].isna()].copy()
    # samp_na['검사채혈일시'] = samp_na['검사접수일시']
    # samp_not_na = conc_df[~conc_df['검사채혈일시'].isna()].copy()
    # conc_df = pd.concat([samp_na, samp_not_na])
    # raise ValueError

    if len(conc_df)==0:
        print(f"({finx}) {pname} / {pid} / No sampling time data")
        no_sampling_list.append(pid)
        continue
    else:
        print(f"({finx}) {pname} / {pid}")


    conc_df['ID'] = pid
    conc_df['NAME'] = pname
    conc_df['DRUG'] = conc_df['검사명'].map(lambda x: x.split(' ')[0].lower())
    for k,v in {'오더DT':'처방발행일시', '라벨DT':'라벨출력일시','채혈DT':'검사채혈일시','접수DT':'검사접수일시','시행DT':'검사시행일시','보고DT':'검사보고일시'}.items():
        conc_df[k] = conc_df[v].map(lambda x: x.replace(' ','T').replace('_x000D_','') if type(x)==str else np.nan)


    conc_df['오더일'] = conc_df['접수DT'].map(lambda x:x.split('T')[0] if type(x)==str else np.nan)
    # conc_df['오더일'] = conc_df['접수DT'].map(lambda x: x.split('T')[0] if type(x) == str else np.nan)
    conc_df['보고일'] = conc_df['보고DT'].map(lambda x: x.split('T')[0] if type(x) == str else np.nan)

    # conc_df = conc_df[~conc_df['CONC'].isna()].copy()
    sampling_result_df.append(conc_df[['ID', 'NAME', '오더일','보고일', 'DRUG', '오더DT','라벨DT', '채혈DT','접수DT', '시행DT', '보고DT']])

sampling_result_df = pd.concat(sampling_result_df, ignore_index=True).drop_duplicates(['ID', '보고일', '오더일', '채혈DT', '접수DT'])
# sampling_result_df = pd.concat(sampling_result_df, ignore_index=True).drop_duplicates(['ID', '보고일', '오더일', '채혈DT'])
sampling_result_df = sampling_result_df.drop_duplicates(['ID', 'NAME', '오더일','보고일', 'DRUG', '오더DT','라벨DT', '채혈DT','접수DT', '시행DT', '보고DT'])
sampling_result_df.to_csv(f"{output_dir}/final_sampling_df.csv", encoding='utf-8-sig', index=False)
# sampling_result_df[['ID', 'NAME', '오더일','보고일','채혈DT']]