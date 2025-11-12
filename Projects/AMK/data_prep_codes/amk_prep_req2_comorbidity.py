from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr

# result_type = 'Phoenix'
# result_type = 'R'

prj_name = 'AMK'
prj_dir = f'./Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


# pid_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/CCM.csv")
# pid_df = pid_df.rename(columns={'Deidentification_ID':'UID','환자번호':'REQ_ID'}).drop(['Column1', 'Column2'],axis=1)

pid_decode_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/재식별 파일.csv")
pid_decode_df = pid_decode_df.rename(columns={'환자번호':'UID','Deidentification_ID':'PID'})
pid_decode_df['UID'] = pid_decode_df['UID'].map(lambda x: x.split('-')[0])

df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/CM.csv")
df = df.rename(columns={'ID':'UID'})
df['UID'] = df['UID'].map(lambda x: x.split('-')[0])
# df.columns
dmhtn_df = pd.read_csv(f"{resource_dir}/[AMK_AKI_ML_DATA]/CM_DM_HTN.csv", encoding='euc-kr')
dmhtn_df = dmhtn_df.rename(columns={'ID':'UID','DATETIME':'DATE'})
dmhtn_df['UID'] = dmhtn_df['UID'].map(lambda x: x.split('-')[0])
dmhtn_df['DATE'] = dmhtn_df['DATE'].map(lambda x:x.split(' ')[0])

cm_index_df = pd.read_excel(f"{resource_dir}/[AMK_AKI_ML_DATA]/CM_index_table.xlsx")
cm_index_df['CM_PERIOD_FROM'] = cm_index_df['CM_PERIOD'].map(lambda x:int(x.split('d')[0].split('-')[-1]))
cm_index_df['CM_PERIOD_UNTIL'] = cm_index_df['CM_PERIOD'].map(lambda x:int(x.split('~-')[-1].split('d')[0]))
cm_index_df = cm_index_df.rename(columns={'CAT_NUM':'CM_CATNUM'})
cm_index_df = cm_index_df[['CM_CATNUM','CM_CAT','CM_PERIOD_FROM','CM_PERIOD_UNTIL']].copy()
# cm_index_df.columns
#
# df.columns
# df['NAME2']

df = df.melt(id_vars=['UID','DATE','FIRSTDATE', 'NAME', 'NAME2', 'CODE'], value_vars=['MLCODE1', 'MLCODE2', 'MLCODE3', 'MLCODE4', 'MLCODE5', 'MLCODE6', 'MLCODE7', 'MLCODE8', 'MLCODE9', 'MLCODE10', 'MLCODE11', 'MLCODE12', 'MLCODE13', 'MLCODE14', 'MLCODE15', 'MLCODE16', 'MLCODE17', 'MLCODE18', 'MLCODE19', 'MLCODE20', 'MLCODE21'], var_name='CM_CATNUM', value_name="VALUE")
df = df[df['VALUE']>0][['UID','DATE','FIRSTDATE','CM_CATNUM','VALUE']].copy()
df['DATE'] = df['DATE'].map(lambda x:f"{x[-4:]}-{x[3:5]}-{x[0:2]}")
df['FIRSTDATE'] = df['FIRSTDATE'].map(lambda x:f"{x[-4:]}-{x[3:5]}-{x[0:2]}")
df['CM_CATNUM'] = df['CM_CATNUM'].map(lambda x:int(x.replace('MLCODE','')))
df = df.drop(['VALUE'],axis=1)

dmhtn_df = dmhtn_df.drop(['NAME', 'CODE'],axis=1).melt(id_vars=['UID','DATE'],
                                            value_vars=['MLCODE1','MLCODE2'],
                                            var_name='CM_CATNUM', value_name="VALUE"
                                            )
dmhtn_df['VALUE'] = dmhtn_df['VALUE'].replace('.','0').map(int)
dmhtn_df = dmhtn_df[dmhtn_df['VALUE'] > 0].drop(['VALUE'],axis=1).copy()
dmhtn_df['CM_CATNUM'] = dmhtn_df['CM_CATNUM'].map(lambda x:int(x.replace('MLCODE','')))
dmhtn_df['FIRSTDATE'] = dmhtn_df['DATE'].copy()
dmhtn_df = dmhtn_df[list(df.columns)].copy()
# dmhtn_df['UID'].iloc[0]
cm_df = pd.concat([df, dmhtn_df])
cm_df = cm_df.drop_duplicates(['UID','DATE','CM_CATNUM']).sort_values(['UID','DATE','CM_CATNUM'])


cm_df = cm_df.merge(cm_index_df, on=['CM_CATNUM'], how='left')

cm_df['UID'] = cm_df['UID'].map(str)
cm_df = cm_df.merge(pid_decode_df, on=['UID'],how='left')
cm_df['UID'] = cm_df['PID'].copy()
cm_df = cm_df.drop(['PID'], axis=1)

if not os.path.exists(f'{output_dir}'):
    os.mkdir(f'{output_dir}')
cm_df.to_csv(f"{output_dir}/final_comorbidity_df.csv", encoding='utf-8-sig', index=False)
# cm_df.columns

# Comorbidity가 존재하면 데이터가 있는 것으로 생각하면 됨