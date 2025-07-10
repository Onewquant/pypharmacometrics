from tools import *
from pynca.tools import *

df = pd.read_csv("C:/Users/ilma0/NONMEMProjects/VancoObesity/1003.csv")
# df = df.iloc[:,:-1].copy()
df = df.drop(['PH','RATE.1','UA'], axis=1)
df.columns

core_df = df.loc[:,:'RATE'].copy()
covar_df = df.loc[:,'SEX':].copy()

for col in ['AST','CRP','ANC']:

    covar_df[col] = covar_df['AST'].map(lambda x:float(str(x).replace('>','').replace('<','')) if type(x)==str else x)

covar_df = covar_df.replace('nan',np.nan)
covar_df = covar_df.fillna(df.mean(numeric_only=True))

uid_dict = {val:inx+1001 for inx, val in enumerate(df['ID'].unique())}

df = pd.concat([core_df, covar_df],axis=1)
df['ID'] = df['ID'].map(uid_dict)
# df['CRP'].unique()
# df.isna().sum()





df.to_csv("C:/Users/ilma0/NONMEMProjects/VancoObesity/1004.csv", index=False)


# set(df['ANC'].unique())