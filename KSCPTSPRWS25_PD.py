from tools import *
from pynca.tools import *
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor


project_dict = {'KSCPTSPRWS25MULTI':['SGLT2INH'], 'KSCPTSPRWS25SINGLE':['SGLT2INH']}
modeling_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25"
modeling_prepconc_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/modeling_prep_data"
resource_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/resource"
results_dir_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/results"

## AUC 구하기

mdprep_conc_files = glob.glob(modeling_prepconc_dir_path+'/MDP_*.csv')
mdprep_conc_df = pd.read_csv(mdprep_conc_files[0]).copy()
mdprep_conc_df = mdprep_conc_df[(mdprep_conc_df['MDV']!=1)&(~(mdprep_conc_df['TAD'].isin([96, 120])))].copy()
mdprep_conc_df = mdprep_conc_df.drop(columns=['MDV','AMT','RATE','DUR','CMT'])
mdprep_conc_df['DAY'] = mdprep_conc_df['TAD'].map(lambda x:1 if x < 144 else 7)

adj_col_df = list()
for inx, df_frag in mdprep_conc_df.groupby(['ID','DAY','UID','GRP']):

    adj_tad_df = df_frag[['ID','UID','GRP', 'TAD', 'TIME']].copy()
    if inx[1]==1:
        adj_tad_df['ADJTAD'] = adj_tad_df['TAD'].copy()
        adj_tad_df['ADJTIME'] = adj_tad_df['TIME'].copy()
        adj_col_df.append(adj_tad_df)
    else:
        min_tad_rows = adj_tad_df[adj_tad_df['TAD']==144]
        if len(min_tad_rows)==0:
            raise ValueError
        else:
            min_tad = min_tad_rows.iloc[0]['TAD']
            min_time = min_tad_rows.iloc[0]['TIME']
        adj_tad_df['ADJTAD'] = adj_tad_df['TAD']-min_tad
        adj_tad_df['ADJTIME'] = adj_tad_df['TIME']-min_time
        adj_col_df.append(adj_tad_df)

adj_col_df = pd.concat(adj_col_df, ignore_index=True)
mdprep_conc_df = mdprep_conc_df.merge(adj_col_df, on=['ID','UID','GRP','TAD','TIME'], how='left')
mdprep_conc_df = mdprep_conc_df[mdprep_conc_df['ADJTAD'] <= 24].copy()
mdprep_conc_df['Short_UID'] = mdprep_conc_df['UID'].map(lambda x:int(x[-2:]))
mdprep_conc_df['UID'] = mdprep_conc_df['UID'].map(lambda x:x[1:])
mdprep_conc_df['DV'] = mdprep_conc_df['DV'].astype(float)

# mdprep_conc_df.to_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/modeling_prep_data/calculate_iAUC.csv", index=False, encoding='utf-8-sig')

# mdprep_conc_df['TADADJ'] =

iAUC_df=pd.DataFrame({'Name':['T[0-6h]','T[6-12h]','T[12-24h]','T[0-24h]'], 'Start':[0,6,12,0], 'End':[6,12,24,24]})
# result = tblNCA(concData=mdprep_conc_df,key=['ID','DAY'],colTime="A_TIME",colConc='CONC',down = "Linear",iAUC=iAUC_df,dose=0.5,slopeMode='BEST',colStyle='pw')
iauc_result = tblNCA(concData=mdprep_conc_df,key=['ID','UID','DAY'],colTime="ADJTIME",colConc='DV',down = "Log",iAUC=iAUC_df,dose=0.5,slopeMode='BEST',colStyle='pw')
iauc_res_df = iauc_result.iloc[1:,:-1].melt(id_vars=list(iauc_result.iloc[:,:-1].columns[:-4]), var_name='AUCTIMEINT', value_name='IAUC')

## PD와 연결

# iauc_result = pd.read_csv(modeling_prepconc_dir_path+'/sglt2i_iauc.csv')
# iauc_result = iauc_result.rename(columns={'Short_UID':'UID'})
# iauc_result['UID'] = iauc_result['UID'].map(lambda x:'S'+str(x).zfill(3))
# iauc_res_df = iauc_result.melt(id_vars=list(iauc_result.columns[:-4]), var_name='AUCTIMEINT', value_name='IAUC')
# iauc_result['AUClast'] = iauc_result['AUCLST']
# if 'AUCLST' in iauc_res_df.columns: iauc_res_df['AUClast'] = iauc_result['AUCLST']

iauc_res_df['AUCTIMEINT'] = iauc_res_df['AUCTIMEINT'].map(lambda x:x.split('[')[-1].split('h]')[0])
iauc_res_df['N_Day'] = iauc_res_df['DAY'].copy()

## Blood PD dataset 정리

bpd_df = pd.read_excel(f"{resource_dir_path}/PD prep_Final_0707_BASE group.xlsx")
bpd_df = bpd_df.iloc[1:]
bpd_df = bpd_df.replace('.',np.nan)
bpd_df['DAY'] = bpd_df['N_DAY'].astype(float)
bpd_df['N_TIME'] = bpd_df['N_TIME'].astype(float)
bpd_df['C_Serum GLU'] = bpd_df['C_Serum GLU'].astype(float)
bpd_df = bpd_df[(bpd_df['DAY'].isin([1,7]))&(bpd_df['N_TIME'] <= 24)].copy()

pg_zero_df = bpd_df[bpd_df['N_TIME']==0][['ID','DAY','C_Serum GLU']].rename(columns={'C_Serum GLU':'PG_ZERO'})

glu_nca_result = tblNCA(concData=bpd_df,key=['ID','DAY'],colTime="N_TIME",colConc='C_Serum GLU',down = "Log",dose=0.1,slopeMode='BEST',colStyle='pw')
# glu_nca_result.columns
s_glu_df = glu_nca_result[['ID','DAY','AUClast']].iloc[1:]
s_glu_df['PG_AVG'] = s_glu_df['AUClast'].astype(float)/24
s_glu_df = s_glu_df.merge(pg_zero_df, on=['ID','DAY'], how='left')
s_glu_df = s_glu_df.rename(columns={'AUClast':'AUC_GLU','DAY':'N_Day','ID':'UID'})

# bpd_df[bpd_df['ID']=='S031']

## Urine PD dataset 정리

upd_df = pd.read_excel(f"{resource_dir_path}/PD urine_prep_Final_0707_prep_BASE group.xlsx")
# upd_df = upd_df.iloc[1:].copy()
# upd_df.columns
upd_df = upd_df[upd_df['N_Time'].isin(['0-24','0-6','6-12','12-24'])]
upd_df = upd_df.rename(columns={'ID':'UID'})
upd_df['DAY'] = upd_df['N_Day']
upd_df['GRP'] = upd_df['GRP'].map(int)
upd_df['AUCTIMEINT'] = upd_df['N_Time']
upd_df['Amount_GLU'] = upd_df['Amount_GLU'].astype(float)
upd_df['Volume'] = upd_df['Volume'].astype(float)
upd_addcol_df = upd_df.groupby(['UID','DAY']).agg({'Amount_GLU':'sum','Volume':'sum'}).reset_index(drop=False).rename(columns={'Amount_GLU':'UGE24','Volume':'VOLUME24'})
upd_baseline_df = upd_addcol_df[upd_addcol_df['DAY']==-1].copy().rename(columns={'UGE24':'UGEbase'}).drop(['DAY','VOLUME24'],axis=1)
upd_addcol_df = upd_addcol_df.merge(upd_baseline_df, on=['UID'],how='left')
upd_df = upd_df.merge(upd_addcol_df, on=['UID','DAY'],how='left')

# daily auc만으로 정리
upd_daily_df = upd_df.copy()
upd_daily_df = upd_daily_df[upd_daily_df['Accumulated Amount_GLU']==upd_daily_df['UGE24']].copy()
upd_daily_df['N_Time'] = '0-24'
upd_daily_df['Volume'] = upd_daily_df['VOLUME24']
upd_daily_df['Amount_GLU'] = upd_daily_df['UGE24']
upd_daily_df['AUCTIMEINT'] = upd_daily_df['N_Time']
upd_daily_df = upd_daily_df.drop_duplicates(['UID','N_Day','N_Time'])


# 공변량 붙이기
covar_df = mdprep_conc_df[['ID','UID']+list(mdprep_conc_df.columns)[7:15]].drop_duplicates(['ID'])
upd_daily_df = upd_daily_df.merge(covar_df, on=['UID'],how='left')
upd_daily_df = upd_daily_df.merge(s_glu_df, on=['UID','N_Day'],how='left')
upd_daily_df = upd_daily_df.merge(iauc_res_df[['UID','N_Day','AUCTIMEINT','AUClast','Cmax','Tmax']], on=['UID','N_Day','AUCTIMEINT'],how='left')
# iauc_res_df.columns
upd_daily_df['EFFECT0'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])/upd_daily_df['PG_ZERO']
upd_daily_df['EFFECT1'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])/upd_daily_df['PG_AVG']
upd_daily_df['EFFECTgfr'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])/(upd_daily_df['PG_AVG']*upd_daily_df['GFR'])
upd_daily_df['EFFECTdelta'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])


# upd_daily_df['GFR_WEIGHTED_AUC'] = upd_daily_df['GFR']*upd_daily_df['AUClast']

# upd_daily_df.columns


# EDA
pre_df = upd_daily_df[upd_daily_df['DAY']==-1].sort_values(['GRP','UID'])
nss_df = upd_daily_df[upd_daily_df['DAY']==1].sort_values(['GRP','UID'])
ss_df = upd_daily_df[upd_daily_df['DAY']==7].sort_values(['GRP','UID'])
integ_df = pd.concat([nss_df, ss_df], ignore_index=True)

pre_df.to_csv(f"{results_dir_path}/KSCPTSPRWS25_PD_BASE.csv", index=False)
nss_df.to_csv(f"{results_dir_path}/KSCPTSPRWS25_PD_NSS.csv", index=False)
ss_df.to_csv(f"{results_dir_path}/KSCPTSPRWS25_PD_SS.csv", index=False)
integ_df.to_csv(f"{results_dir_path}/KSCPTSPRWS25_PD_INTEG(NSS_SS).csv", index=False)
#
# nss_df['TIME'] = 24
# nss_zero_df = nss_df.copy()
# nss_zero_df['TIME'] = 0
# nss_zero_df['EFFECT1'] = 0
# pd_adj_df = pd.concat([nss_df, nss_zero_df]).sort_values(['UID','TIME'])
# pd_adj_df['DV'] = pd_adj_df['EFFECT1']
# pd_adj_df.to_csv(f"{results_dir_path}/KSCPTSPRWS25_PD_NSS_ADJ.csv", index=False)

pd_df = nss_df.copy()
# pd_df = integ_df.copy()

## Multivariable linear regression

# 상관계수 행렬 계산
corr_matrix = pd_df[['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'GFR', 'TBIL', 'AUC_GLU', 'PG_AVG', 'PG_ZERO', 'AUClast']].corr().abs()

# # 히트맵 그리기
# plt.figure(figsize=(15, 10))
# sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap='coolwarm', vmin=-1, vmax=1)
# plt.title("Correlation Matrix Heatmap")
# plt.show()


# 상삼각행렬 추출 (자기 자신과 중복 제거)
lower = corr_matrix.where(np.tril(np.ones(corr_matrix.shape), k=-1).astype(bool))

# 상관계수 0.9 이상인 변수 중 하나 제거
to_drop = [column for column in lower.columns if any(lower[column] > 0.7)]

# 독립변수(X), 종속변수(y)
total_cand = ['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'GFR', 'TBIL', 'AUC_GLU', 'PG_AVG', 'PG_ZERO', 'AUClast']

X = pd_df[list(set(total_cand)-set(to_drop))]
# X = nss_df[list(set(total_cand))]
for c in X.columns:
    X[c]=X[c].astype(float)

y = pd_df['UGE24']

# VIF 계산
vif_data = pd.DataFrame(columns=['Variable','VIF'])
vif_data['Variable'] = list(X.columns)
vif_data['VIF'] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
no_vif_cols = list(vif_data[vif_data['VIF'] < 10]['Variable'])

# 회귀모형 적합
total_model = sm.OLS(y, X[list(vif_data['Variable'])]).fit()
sel_model = sm.OLS(y, X[no_vif_cols]).fit()

# 결과 출력
print(total_model.summary())
print(sel_model.summary())


# nss_df['EFFECT'] = (nss_df['UGE24']-nss_df['BASELINE'])/nss_df['PG_AVG']

# nss_df.columns

# scatter plot 그리기
plt.figure(figsize=(15, 10))
# x = 'AUClast'
# x = 'Cmax'
x = 'GFR'
# x = 'GFR_WEIGHTED_AUC'
# y = 'EFFECTgfr'
# y = 'EFFECT0'
y = 'EFFECT1'
# y = 'EFFECTgfr'
# y = 'EFFECTdelta'
hue = 'GRP'

# sns.scatterplot(data=ss_df, x='AUClast', y='EFFECT0', hue='GRP')
# sns.scatterplot(data=ss_df, x='AUClast', y='EFFECT1', hue='GRP')
# sns.scatterplot(data=nss_df, x='AUClast', y='EFFECT0', hue='GRP')
# sns.scatterplot(data=nss_df, x='AUClast', y='EFFECT1', hue='GRP')
sns.scatterplot(data=integ_df, x=x, y=y, hue=hue)
plt.title(f'{x} vs {y} by {hue}')
plt.xlabel(x)
plt.ylabel(y)
plt.grid(True)
plt.tight_layout()
plt.show()


X = nss_df[[x]]
X_const = sm.add_constant(X)
y_vals = nss_df[y]

model = sm.OLS(y_vals, X_const).fit()

intercept, slope = model.params
r_squared = model.rsquared
p_value = model.pvalues[x]


## Fitting

from scipy.optimize import curve_fit

# Sigmoid Emax 모델 함수 정의
def sigmoid_emax(conc, Emax, EC50, H):
    return Emax * conc**H / (EC50**H + conc**H)

gfr = nss_df['GFR']
effect = nss_df['EFFECT1']

# curve fitting
popt, pcov = curve_fit(sigmoid_emax, gfr, effect, p0=[0.5, 50, 1])  # 초기값 설정

# 피팅된 파라미터 출력
Emax_fit, EC50_fit, H_fit = popt
print(f"Emax: {Emax_fit:.2f}, EC50: {EC50_fit:.2f}, Hill coefficient: {H_fit:.2f}")

# 예측값 계산
gfr_fit = np.linspace(0.1, 120, 100)
effect_fit = sigmoid_emax(gfr_fit, *popt)

# 시각화
plt.figure(figsize=(15, 10))
plt.scatter(gfr, effect, label='Observed data', color='black')
plt.plot(gfr_fit, effect_fit, label='Sigmoid Emax Fit', color='blue')
plt.xlabel('GFR')
plt.ylabel('Effect1')
plt.title('Sigmoid Emax Model Fit (Python)')
plt.legend()
plt.grid(True)
plt.show()

