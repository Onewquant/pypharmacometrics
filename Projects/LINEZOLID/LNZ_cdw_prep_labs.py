from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr


prj_name = 'LINEZOLID'
prj_dir = './Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

lab_file_list = glob.glob(f"{resource_dir}/lab_frag_data/L_*.csv")
lab_data_list = list()
for finx, f in enumerate(lab_file_list): #break

    lab_frag_df = pd.read_csv(f, encoding='utf-8-sig')
    lab_frag_df = lab_frag_df.rename(columns={'환자번호':'UID','검사 접수일자':'DATE','검사 접수일시':'DATETIME', '검사명':'LAB', '항목별 검사결과':'VALUE', '채혈일시':'SAMPLING_DT', '검사 오더일자':'ORD_DATE','검사 시행일자':'PRAC_DATE','검사 결과보고일자':'REP_DATE','검사 코드':'LAB_CODE'})
    lab_frag_df = lab_frag_df.dropna(subset=['VALUE'])
    lab_frag_df['UID'] = lab_frag_df['UID'].map(lambda x:x.split('-')[0])
    lab_frag_df['DATE'] = lab_frag_df.apply(lambda x: x['DATE'] if str(x['SAMPLING_DT'])=='nan' else x['SAMPLING_DT'].split(' ')[0], axis=1)
    # raise ValueError
    lab_name = f.split('L_')[-1].replace('.csv','')
    print(f"({finx}) {lab_name} / {len(lab_frag_df)} rows")
    lab_frag_df.to_csv(f"{resource_dir}/lnz_cdw_lab_data.csv", index=False, mode='a', header=True if finx==0 else False)
# lab_frag_df.columns
#
df = pd.read_csv(f"{resource_dir}/lnz_cdw_lab_data.csv")

# Orders
total_lab_cols = list(df['LAB'].drop_duplicates().sort_values())
total_lab_cols = ['UID', 'DATE'] + total_lab_cols
print(f"# 총 랩 검사명 파악 완료 / total_lab_cols: {len(total_lab_cols[2:])} 개\n")

# lab_list_df = pd.DataFrame(columns=['LAB'])
# lab_list_df['LAB'] = total_lab_cols
# lab_list_df.to_csv(f"{output_dir}/lnz_lablist_df.csv", encoding='utf-8-sig', index=False)
# lab_list_df = pd.read_csv(f"{output_dir}/lnz_lablist_df.csv")
# total_lab_cols = lab_list_df['LAB'].to_list()

print(f"# 각 환자별 날짜별 랩수치 파악 시작\n")

# lab_files = glob.glob(f'{resource_dir}/lab/{prj_name}_lab(*).xlsx')
lab_result_df = list()
uid_count = 0
for pid, fdf in df.groupby('UID'): #break

    print(f"({uid_count}) {pid}")

    fdf['VALUE'] = pd.to_numeric(fdf['VALUE'], errors='coerce') # 숫자 아닌 랩은 NAN으로 변환
    fdf = fdf.drop_duplicates(['DATE','LAB'], keep='last', ignore_index=True)
    fpv_df = fdf.pivot_table(index='DATE', columns='LAB', values='VALUE', aggfunc='min').reset_index(drop=False).fillna(method='ffill')
    # fpv_df = fdf.pivot(index='DATETIME', columns='LAB', values='VALUE').reset_index(drop=False).fillna(method='ffill')

    # fpv_df = fdf.pivot_table(index='DATETIME', columns='LAB', values='VALUE', aggfunc='min').reset_index(drop=False)
    fpv_df.columns.name = None
    fpv_df['UID'] = [pid]*len(fpv_df)
    ind_lab_cols = list(fpv_df.columns)
    for c in set(total_lab_cols).difference(set(ind_lab_cols)):
        fpv_df[c] = [np.nan]*len(fpv_df)
    if uid_count==0:
        lab_result_df = fpv_df[total_lab_cols].copy()
    else:
        lab_result_df = pd.concat([lab_result_df, fpv_df[total_lab_cols].copy()])

    uid_count+=1

print(f"# 날짜 전체로 매칭 시작\n")

full_result_df = list()
count = 0
for uid, uid_df in lab_result_df.groupby('UID',as_index=False): #break

    min_lab_date = uid_df['DATE'].min()
    max_lab_date = uid_df['DATE'].max()

    print(f"({count}) {uid} / Date Range : ({min_lab_date},  {max_lab_date})")

    uid_fulldt_df = pd.DataFrame(columns=['UID','DATE'])
    uid_fulldt_df['DATE'] = pd.date_range(start=min_lab_date,end=max_lab_date).astype(str)
    uid_fulldt_df['UID'] = uid

    uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATE'], how='left').fillna(method='ffill')

    if count == 0:
        header = True
    else:
        header = False
    uid_fulldt_df.to_csv(f"{output_dir}/lnz_final_lab_df.csv", encoding='utf-8-sig', index=False, header=header, mode='a')

    count += 1
