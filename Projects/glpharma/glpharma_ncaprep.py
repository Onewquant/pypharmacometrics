from tools import *
from pynca.tools import *

result_type = 'Phoenix'
result_type = 'R'

prj_name = 'GLPHARMA'
ip_name = 'W2406'
prj_dir = 'Projects/glpharma/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

df = pd.read_csv(f'{prj_dir}/glpharma_CONC.csv')
seq_df = pd.read_csv(f'{prj_dir}/glpharma_SEQUENCE.csv')

df = df.drop(columns=['No.']).melt(id_vars=['PERIOD','NTIME'], var_name='ID', value_name='CONC')
seq_df = seq_df.T.iloc[1:].reset_index(drop=False).rename(columns={'index':'ID',0:'SEQUENCE'})

df = df.merge(seq_df, on=['ID'], how='left')
df['DOSE'] = 150
df['ATIME'] = df['NTIME']
df['DRUG'] = ip_name
df['TRT'] = ''
df['TRT'] = np.where((((df['PERIOD']==1)&(df['SEQUENCE']==1))|((df['PERIOD']==2)&(df['SEQUENCE']==2))), 'R', 'T')

# 채혈 되지 않은 대상자 분석에서 제외
df = df[df['CONC']!='-'].sort_values(['TRT','ID','NTIME'],ascending=[False,True,True],ignore_index=True)


result_cols = ['ID','DRUG','TRT','DOSE','ATIME','NTIME','CONC','PERIOD','SEQUENCE']
unit_row_dict = {'DOSE': 'mg', 'NTIME': 'h', 'ATIME': 'h', 'CONC': 'ug/mL'}

df = df[result_cols].reset_index(drop=True)





first_float_index_dict = dict()
for inx, frag_df in df.groupby(['ID','TRT']):
    first_index = frag_df.index[0]
    first_float_index = -1
    for finx, row in frag_df.iterrows():
        try:
            float(row['CONC'])
            first_float_index = finx
            break
        except:
            continue

    first_float_index_dict.update({f"{inx[0]}_{inx[1]}" : (first_index ,first_float_index)})

for inx, row in df.iterrows():
    key = f"{row['ID']}_{row['TRT']}"

    if (first_float_index_dict[key][0]<=inx) and (first_float_index_dict[key][1] > inx) and (row['CONC'] in ['ND']):
        df.at[inx,'CONC'] = 0.0
    elif (first_float_index_dict[key][0]<=inx) and (first_float_index_dict[key][1] > inx) and (row['CONC'] in ['BQL']):
        df.at[inx,'CONC'] = np.nan
    elif (row['NTIME'] != 0) and (row['CONC'] in ['BQL','ND']):
        df.at[inx, 'CONC'] = np.nan
    elif (row['NTIME'] != 0) and (row['CONC'] not in ['BQL','ND']):
        df.at[inx, 'CONC'] = float(df.at[inx, 'CONC'])
    else:
        raise ValueError

# df['CONC'] = df['CONC'].replace('ND',np.nan).replace('BQL',np.nan).replace('-',np.nan)

# raise ValueError
if result_type == 'Phoenix':
    df = df[~df['CONC'].isna()].reset_index(drop=True)
    df['CONC'] = df['CONC'].map(lambda x: str(x) if not np.isnan(x) else '.')
    prep_df = df[result_cols]


    additional_row = dict()
    for c in list(prep_df.columns):
        try:
            additional_row[c] = unit_row_dict[c]
        except:
            additional_row[c] = ''
    prep_df = pd.concat([pd.DataFrame([additional_row], index=['', ]), prep_df])
elif result_type == 'R':
    prep_df = df.dropna()

prep_df = prep_df[result_cols].sort_values(['TRT','ID','PERIOD','NTIME'])

result_file_name = f"{prj_name}_ConcPrep_{ip_name}_{result_type}.csv"
result_file_path = f"{output_dir}/{result_file_name}"
prep_df.to_csv(result_file_path, header=True, index=False)
