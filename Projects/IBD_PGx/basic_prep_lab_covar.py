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
total_lab_cols = set()

lab_files = glob.glob(f'{resource_dir}/lab/IBD_PGx_lab(*).xlsx')
lab_result_df = list()
for finx, fpath in enumerate(lab_files):# break

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    print(f"({finx}) {pname} / {pid}")

    # if pid in ("15322168", "19739357", "34835292", "37366865", "21618097", "36898756", "36975211", "37858047"):       # lab, order 파일 다시 수집 필요
    #     continue
    # fdf['검사명'].drop_duplicates().sort_values()
    fdf = pd.read_excel(fpath)
    total_lab_cols = total_lab_cols.union(set(fdf['검사명'].drop_duplicates().sort_values()))

total_lab_cols = list(total_lab_cols)
total_lab_cols.sort()
total_lab_cols = ['UID', 'DATETIME'] + total_lab_cols
print(f"# 총 랩 검사명 파악 완료 / total_lab_cols: {len(total_lab_cols)} 개\n")
# total_lab_cols


print(f"# 각 환자별 날짜별 랩수치 파악 시작\n")

lab_files = glob.glob(f'{resource_dir}/lab/IBD_PGx_lab(*).xlsx')
lab_result_df = list()
for finx, fpath in enumerate(lab_files):

    pid = fpath.split('(')[-1].split('_')[0]
    pname = fpath.split('_')[-1].split(')')[0]

    print(f"({finx}) {pname} / {pid}")

    # fdf.columns
    fdf = pd.read_excel(fpath)
    fdf['DATETIME'] = fdf[['보고일', '오더일']].max(axis=1)
    fdf['LAB'] = fdf['검사명']
    fdf['VALUE'] = fdf['검사결과']
    fdf = fdf.drop_duplicates(['DATETIME','LAB'], keep='last', ignore_index=True)
    fpv_df = fdf.pivot_table(index='DATETIME', columns='LAB', values='VALUE', aggfunc='min').reset_index(drop=False).fillna(method='ffill')
    # fpv_df = fdf.pivot_table(index='DATETIME', columns='LAB', values='VALUE', aggfunc='min').reset_index(drop=False)
    fpv_df.columns.name = None
    fpv_df['UID'] = pid
    ind_lab_cols = list(fpv_df.columns)
    for c in set(total_lab_cols).difference(set(ind_lab_cols)):
        fpv_df[c] = np.nan
    if finx==0:
        lab_result_df = fpv_df[total_lab_cols].copy()
    else:
        lab_result_df = pd.concat([lab_result_df, fpv_df[total_lab_cols].copy()])


print(f"# 날짜 전체로 매칭 시작\n")

full_result_df = list()
count = 0
for uid, uid_df in lab_result_df.groupby('UID',as_index=False): #break


    min_lab_date = uid_df['DATETIME'].min()
    max_lab_date = uid_df['DATETIME'].max()

    print(f"({count}) {uid} / Date Range : ({min_lab_date},  {max_lab_date})")

    uid_fulldt_df = pd.DataFrame(columns=['UID','DATETIME'])
    uid_fulldt_df['DATETIME'] = pd.date_range(start=min_lab_date,end=max_lab_date).astype(str)
    uid_fulldt_df['UID'] = uid

    uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left').fillna(method='ffill')
    # uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')

    if count==0:
        full_result_df = uid_fulldt_df.copy()
    else:
        full_result_df = pd.concat([full_result_df, uid_fulldt_df.copy()])

    count+=1

full_result_df = full_result_df.reset_index(drop=True)
# full_result_df.drop_duplicates(['UID','DATETIME'])
full_result_df.to_csv(f"{output_dir}/lab_df.csv", encoding='utf-8-sig', index=False)

