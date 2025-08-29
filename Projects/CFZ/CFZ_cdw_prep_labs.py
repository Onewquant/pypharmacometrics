from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr



prj_name = 'CFZ'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/CFZ'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

df1 = pd.read_csv(f"{resource_dir}/cfz_cdw_lab_data.csv")
df2 = pd.read_csv(f"{resource_dir}/cfz_cdw_lab_data2.csv")
# df2 = df2.rename(columns={'환자번호':'UID','검사 접수일자':'DATE','검사 접수일시':'DATETIME', '검사명':'LAB','검사 세부 항목명':'LAB_DETAIL', '항목별 검사결과':'VALUE', '채혈일시':'SAMPLING_DT'})


df=pd.concat([df1,df2])
df = df.rename(columns={'환자번호':'UID','검사 접수일자':'DATE','검사 접수일시':'DATETIME', '검사명':'LAB','검사 세부 항목명':'LAB_DETAIL', '항목별 검사결과':'VALUE', '채혈일시':'SAMPLING_DT'})
df = df.dropna(subset=['VALUE'])
df['UID'] = df['UID'].map(lambda x:x.split('-')[0])
# df['SAMPLING_DT'].unique()
# df[df['SAMPLING_DT'].isna()]
# df['SAMPLING_DT'].map(str)
df['DATE'] = df.apply(lambda x: x['DATE'] if str(x['SAMPLING_DT'])=='nan' else x['SAMPLING_DT'].split(' ')[0], axis=1)


df_urbc = df[df['LAB']=='RBC'].copy()
df_urbc = df_urbc[~df_urbc['VALUE'].isna()].copy()
df_urbc = df_urbc[~df_urbc['VALUE'].isin(['.','q.n.s','+'])].copy()
df_urbc_new = list()
# df_urbc_cnt = list()
# df_urbc_cnt_ul = list()
# len(df_urbc)
for inx, row in df_urbc.iterrows(): #break

    value_split = re.findall(r'\d+[\.]?\d*',row['VALUE'])
    if len(value_split)==1:
        urbc_mean_cnt = pd.Series(value_split).astype(float).mean()
        urbc_cnt_ul = np.nan
    elif len(value_split)==2:
        if len(re.findall(r'\(\d+[\.]?\d*/㎕\)',row['VALUE']))==1:
            urbc_mean_cnt = float(value_split[0])
            urbc_cnt_ul = float(value_split[-1])
        else:
            urbc_mean_cnt = pd.Series(value_split).astype(float).mean()
            urbc_cnt_ul = np.nan
    elif len(value_split)==3:
        urbc_mean_cnt = pd.Series(value_split[:2]).astype(float).mean()
        urbc_cnt_ul = float(value_split[-1])
    elif len(value_split)==4:
        urbc_mean_cnt = pd.Series(value_split[:2]).astype(float).mean()
        urbc_cnt_ul = float(value_split[-1])
    else:
        raise ValueError

    print(f"({inx}) {row['VALUE']} / num count: {len(value_split)} / {urbc_mean_cnt} / {urbc_cnt_ul}")
    # value_split = row['VALUE'].split('(')
    # if len(value_split)==1:
    #     urbc_mean_cnt = pd.Series(re.findall(r'\d+[\.]?\d*',row['VALUE'])).mean()
    #     urbc_cnt_ul = np.nan
    # elif len(value_split)==2:
    #     pre_split = value_split[0]
    #     urbc_mean_cnt = float(re.findall(r'\d+[\.]?\d*', pre_split)[0])
    #     post_split = value_split[1]
    #     urbc_cnt_ul = float(re.findall(r'\d+[\.]?\d*', post_split)[0])
    # else:
    #     raise ValueError
    new_row1 = row.copy()
    new_row1['LAB'] = "RBC_count"
    new_row1['VALUE'] = urbc_mean_cnt

    new_row2 = row.copy()
    new_row2['LAB'] = "RBC_ul"
    new_row2['VALUE'] = urbc_cnt_ul

    if not np.isnan(new_row1['VALUE']):
        df_urbc_new.append(new_row1)
    if not np.isnan(new_row2['VALUE']):
        df_urbc_new.append(new_row2)

# df_urbc['VALUE'] = df_urbc['VALUE'].map(lambda x:x.split('(')[-1].split('/')[0])
df_urbc_new = pd.DataFrame(df_urbc_new)
df_blood = df[df['LAB']=='Blood'].copy()
df_blood['VALUE'] = df_blood['VALUE'].map({'-':0.0, '3+':3.0, '+/-':0.5, '1+':1.0, '2+':2.0, '?':np.nan, '1+(재검)':1.0, '2+(재검)':2.0,  '3+재검':3.0, '/':np.nan, '+':1.0, '.':np.nan, '1+(2번째검체로 재검)':1.0, '3+(재검)':3.0, ' 3+(재검)':3.0, '-+/-':0.0, '2+재검':2.0})
df = pd.concat([df_urbc_new, df_blood])


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
total_uid_count = len(df['UID'].drop_duplicates())
for pid, fdf in df.groupby('UID'): #break

    print(f"({uid_count} / {total_uid_count}) {pid}")

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

lab_result_df.sort_values(['UID','DATE']).to_csv(f"{output_dir}/cfz_semifinal_lab_df.csv", encoding='utf-8-sig', index=False, mode='w')

print(f"# 날짜 전체로 매칭 시작\n")

full_result_df = list()
count = 0
total_uid_count = len(lab_result_df['UID'].drop_duplicates())
for uid, uid_df in lab_result_df.groupby('UID',as_index=False): #break

    min_lab_date = uid_df['DATE'].min()
    max_lab_date = uid_df['DATE'].max()

    print(f"({count} / {total_uid_count}) {uid} / Date Range : ({min_lab_date},  {max_lab_date})")

    uid_fulldt_df = pd.DataFrame(columns=['UID','DATE'])
    uid_fulldt_df['DATE'] = pd.date_range(start=min_lab_date,end=max_lab_date).astype(str)
    uid_fulldt_df['UID'] = uid

    uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATE'], how='left').fillna(method='ffill')
    # uid_fulldt_df = uid_fulldt_df.merge(uid_df, on=['UID','DATETIME'], how='left')

    if count == 0:
        header = True
    else:
        header = False
    uid_fulldt_df.to_csv(f"{output_dir}/cfz_final_lab_df.csv", encoding='utf-8-sig', index=False, header=header, mode='a')

    # if not os.path.exists(f"{output_dir}/uid_lab_df"):
    #     os.mkdir(f"{output_dir}/uid_lab_df")
    # uid_fulldt_df.to_csv(f"{output_dir}/uid_lab_df/uid_fulldt_df({uid}).csv", encoding='utf-8-sig', index=False)

    count += 1

# full_result_df = full_result_df.reset_index(drop=True)
# full_result_df.to_csv(f"{output_dir}/cfz_final_lab_df.csv", encoding='utf-8-sig', index=False)

