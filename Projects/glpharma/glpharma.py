from tools import *
from pynca.tools import *

result_type = 'Phoenix'
result_type = 'R'

prj_dir = 'Projects/glpharma/resource'

df = pd.read_csv(f'{prj_dir}/glpharma_CONC.csv')
seq_df = pd.read_csv(f'{prj_dir}/glpharma_SEQUENCE.csv')

df = df.drop(columns=['No.']).melt(id_vars=['DRUG','PERIOD','NTIME'], var_name='ID', value_name='CONC')
seq_df = seq_df.T.iloc[1:].reset_index(drop=False).rename(columns={'index':'ID',0:'SEQUENCE'})


df = df.merge(seq_df, on=['ID'], how='left')

# 채혈 되지 않은 대상자 분석에서 제외
df = df[df['CONC']!='-'].sort_values(['DRUG','ID','NTIME'],ascending=[False,True,True],ignore_index=True)


df = df[['ID','DRUG','NTIME','CONC','PERIOD','SEQUENCE']].reset_index(drop=True)




first_float_index_dict = dict()
for inx, frag_df in df.groupby(['ID','DRUG']):
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
    key = f"{row['ID']}_{row['DRUG']}"

    if (first_float_index_dict[key][0]<=inx) and (first_float_index_dict[key][1] > inx) and (row['CONC'] in ['BQL','ND']):
        df.at[inx,'CONC'] = 0.0
    elif (row['NTIME'] != 0) and (row['CONC'] in ['BQL','ND']):
        df.at[inx, 'CONC'] = np.nan
    elif (row['NTIME'] != 0) and (row['CONC'] not in ['BQL','ND']):
        df.at[inx, 'CONC'] = float(df.at[inx, 'CONC'])
    else:
        raise ValueError

# df['CONC'] = df['CONC'].replace('ND',np.nan).replace('BQL',np.nan).replace('-',np.nan)

# raise ValueError
if result_type == 'Phoenix':
    df['CONC'] = df['CONC'].map(lambda x: str(x) if not np.isnan(x) else '.')
    prep_df = df[['ID', 'DOSE', 'NTIME', 'ATIME', 'CONC', 'PERIOD', 'FEEDING', 'DRUG']]

    unit_row_dict = {'DOSE': 'mg', 'NTIME': 'h', 'ATIME': 'h', 'CONC': 'ng/mL'}
    additional_row = dict()
    for c in list(prep_df.columns):
        try:
            additional_row[c] = unit_row_dict[c]
        except:
            additional_row[c] = ''
    prep_df = pd.concat([pd.DataFrame([additional_row], index=['', ]), prep_df])
elif result_type == 'R':
    prep_df = prep_df.dropna()

prep_df = prep_df[['ID', 'DOSE', 'NTIME', 'ATIME', 'CONC', 'PERIOD', 'FEEDING', 'DRUG']].sort_values(['DRUG','ID','PERIOD','NTIME'])

result_file_name = f"CKD383_ConcPrep_({result_type}).csv"
result_file_path = f"{output_dir}/{result_file_name}"
prep_df.to_csv(result_file_path, header=True, index=False)
