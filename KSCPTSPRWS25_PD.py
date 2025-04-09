from tools import *
from pynca.tools import *
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
# iauc_res_df.columns
# iauc_res_df[iauc_res_df['DAY']==1][['ID', 'N_Samples', 'Tmax', 'Cmax', 'AUClast', 'Vz_F_obs', 'Cl_F_obs']]

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
mmhba1c_df = bpd_df.copy()
mmhba1c_df = mmhba1c_df[(bpd_df['DAY'].isin([7]))].copy()[['ID','GRP','N_DAY','HbA1c']].drop_duplicates(['ID','HbA1c']).dropna().reset_index(drop=True)
mmhba1c_df = mmhba1c_df[['ID','N_DAY','HbA1c']].copy()
mmhba1c_df = mmhba1c_df.rename(columns={'ID':'UID', 'HbA1c':'HbA1c_7d','N_DAY':'N_Day'})
bpd_df = bpd_df[(bpd_df['DAY'].isin([1,7]))&(bpd_df['N_TIME'] <= 24)].copy()
# bpd_df[(bpd_df['DAY'].isin([7]))]

pg_zero_df = bpd_df[bpd_df['N_TIME']==0][['ID','DAY','C_Serum GLU','HbA1c']].rename(columns={'C_Serum GLU':'PG_ZERO', 'HbA1c':'HbA1c_base'})

glu_nca_result = tblNCA(concData=bpd_df,key=['ID','DAY'],colTime="N_TIME",colConc='C_Serum GLU',down = "Log",dose=0.1,slopeMode='BEST',colStyle='pw')
# glu_nca_result.columns
s_glu_df = glu_nca_result[['ID','DAY','AUClast']].iloc[1:]
s_glu_df['PG_AVG'] = s_glu_df['AUClast'].astype(float)/24
s_glu_df = s_glu_df.merge(pg_zero_df, on=['ID','DAY'], how='left')
s_glu_df = s_glu_df.rename(columns={'AUClast':'AUC_GLU','DAY':'N_Day','ID':'UID'})
# s_glu_df['PG_AVG'].mean()
# s_glu_df['PG_ZERO'].mean()

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
covar_df = mdprep_conc_df[['ID','UID']+list(mdprep_conc_df.columns)[7:18]].drop_duplicates(['ID'])
# covar_df.drop(columns=['ID']).rename(columns={'UID':'ID'}).to_csv(f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/sglt2i_dataset/KSCPTSPRWS25_SGLT2i_COVAR.csv',index=False)

upd_daily_df = upd_daily_df.merge(covar_df, on=['UID'],how='left')
upd_daily_df = upd_daily_df.merge(s_glu_df, on=['UID','N_Day'],how='left')
upd_daily_df = upd_daily_df.merge(iauc_res_df[['UID','N_Day','AUCTIMEINT','AUClast','Cmax','Tmax']], on=['UID','N_Day','AUCTIMEINT'],how='left')
upd_daily_df = upd_daily_df.merge(mmhba1c_df, on=['UID'],how='left')

# iauc_res_df.columns
# upd_daily_df['EFFECT0'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])/upd_daily_df['PG_ZERO']
upd_daily_df['eGFR'] = upd_daily_df['GFR']
upd_daily_df['eGFRxTBIL'] = upd_daily_df['GFR']*upd_daily_df['TBIL']
upd_daily_df['eGFRxTBILxAUC'] = upd_daily_df['GFR']*upd_daily_df['TBIL']*(upd_daily_df['AUClast'].map(np.log).map(float))
upd_daily_df['dUGE'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])
upd_daily_df['EFFECT0'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])/upd_daily_df['PG_ZERO']
upd_daily_df['EFFECT1'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])/upd_daily_df['PG_AVG']
upd_daily_df['EFFECTgfr'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])/(upd_daily_df['PG_AVG']*upd_daily_df['GFR'])

# upd_daily_df['EFFECTgfr'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])/(upd_daily_df['PG_AVG']*upd_daily_df['GFR'])
# upd_daily_df['EFFECTdelta'] = (upd_daily_df['UGE24']-upd_daily_df['UGEbase'])


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
# pd_df[['UID', 'GRP', 'UGE24', 'UGEbase', 'dUGE','PG_AVG', 'PG_ZERO', 'EFFECT0', 'EFFECT1']]
# pd_df['dUGE'].mean()
# pd_df.columns
# pd_df = integ_df.copy()



## Multivariable linear regression

# 상관계수 행렬 계산
# corr_matrix = pd_df[['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'GFR', 'TBIL', 'AUClast']].corr().abs()
# cov_cand = ['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'GFR', 'TBIL', 'AUClast']
cov_cand = ['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'eGFR', 'AST', 'ALT', 'SODIUM', 'TBIL', 'AUClast']
# cov_cand = ['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'GFR', 'AST','ALT','SODIUM','AUClast']
# cov_cand = ['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'AST', 'ALT', 'SODIUM', 'TBIL', 'AUClast']
corr_matrix = pd_df[cov_cand].corr().abs()

# 히트맵 그리기
plt.figure(figsize=(18, 15))
ax = sns.heatmap(
    corr_matrix,
    annot=True,
    fmt=".2f",
    cmap='coolwarm',
    vmin=-1,
    vmax=1,
    annot_kws={"size": 12}
)
fig_title = "[WSCT] Correlation matrix heatmap of the covariates"
plt.title(fig_title, fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

# colorbar 글씨 크기 조정
ax.collections[0].colorbar.ax.tick_params(labelsize=15)

plt.savefig(f"{results_dir_path}/{fig_title}.png")
plt.cla()
plt.clf()
plt.close()



# 상삼각행렬 추출 (자기 자신과 중복 제거)
lower = corr_matrix.where(np.tril(np.ones(corr_matrix.shape), k=-1).astype(bool))
# upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))

# 상관계수 0.9 이상인 변수 중 하나 제거
to_drop = [column for column in lower.columns if any(lower[column] >= 0.7)]
# to_drop = [column for column in upper.columns if any(upper[column] >= 0.7)]

# 독립변수(X), 종속변수(y)
# total_cand = ['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'GFR', 'TBIL', 'AUC_GLU', 'PG_AVG', 'PG_ZERO', 'AUClast']
# total_cand = ['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'GFR', 'TBIL', 'AUClast']
total_cand = cov_cand

X = pd_df[list(set(total_cand)-set(to_drop))].copy()
# X = nss_df[list(set(total_cand))]
for c in X.columns:
    X[c]=X[c].astype(float)

# VIF 계산
vif_data = pd.DataFrame(columns=['Covariates', 'VIF'])
vif_data['Covariates'] = list(X.columns)
vif_data['VIF'] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]

# VIF 값 오름차순 정렬
df_vif_sorted = vif_data.sort_values(by='VIF', ascending=False)
print(df_vif_sorted)

# VIF 수평 막대그래프 그리기
plt.figure(figsize=(17, 15))
plt.barh(df_vif_sorted['Covariates'], df_vif_sorted['VIF'], color='royalblue')
plt.xlabel('VIF Value', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
fig_title = '[WSCT] Variance Inflation Factor (VIF) of covariates'
plt.title(fig_title, fontsize=15)
plt.grid(axis='x', linestyle='--', alpha=0.5)

# VIF 값 표시
for index, value in enumerate(df_vif_sorted['VIF']):
    plt.text(value + 0.2, index, f"{value:.2f}", va='center', fontsize=13, color='black')

plt.tight_layout()

# 저장
plt.savefig(f"{results_dir_path}/{fig_title}.png", dpi=300)

plt.cla()
plt.clf()
plt.close()



effect_col_list = ['EFFECT1', 'EFFECTgfr']

for effect_col in effect_col_list:
    y = pd_df[effect_col].copy()

    no_vif_cols = list(vif_data[vif_data['VIF'] < 10]['Covariates'])


    # 다중 회귀모형 적합
    total_model = sm.OLS(y, X[list(vif_data['Covariates'])]).fit()
    sel_model = sm.OLS(y, X[no_vif_cols]).fit()

    # 결과 출력
    print(total_model.summary())
    print(sel_model.summary())

    # OLS 모델
    result_df = ols_result_df(total_model, X, y)
    create_pdf_report(total_model, result_df,filepath=f"{results_dir_path}/[WSCT] Multivariable LR results ({effect_col}).pdf")
    result_df = ols_result_df(sel_model, X[no_vif_cols], y)
    create_pdf_report(sel_model, result_df,filepath=f"{results_dir_path}/[WSCT] Univariable LR results ({effect_col}).pdf")


# nss_df['EFFECT'] = (nss_df['UGE24']-nss_df['BASELINE'])/nss_df['PG_AVG']

# nss_df.columns

# scatter plot 그리기
plt.figure(figsize=(15, 12))

# for x in ['AUClast','Cmax','GFR','TBIL','ALT','eGFRxTBIL','eGFRxTBILxAUC']:
for x in ['AUClast', 'Cmax', 'eGFR', 'TBIL', 'ALT']:
    # x='ALT'
    for y in effect_col_list:
        hue = 'GRP'

        X = nss_df[[x]].copy()
        X_const = sm.add_constant(X).applymap(float)
        y_vals = nss_df[y].map(float)

        model = sm.OLS(y_vals, X_const).fit()

        intercept, slope = model.params
        r_squared = model.rsquared
        p_value = model.pvalues[x]


        sns.scatterplot(data=nss_df, x=x, y=y, hue=hue, marker='o')
        fig_title = f'[WSCT] {x} vs {y} by {hue}\nR-squared:{r_squared:.4f}, p-value: {p_value:.4f}\nbeta: {slope:.4f}, intercept: {intercept:.4f} '
        plt.title(fig_title, fontsize=14)
        plt.xlabel(x, fontsize=14)
        plt.ylabel(y, fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=14)
        plt.grid(True)
        plt.tight_layout()
        # plt.show()

        plt.savefig(f"{results_dir_path}/{fig_title.split('R-squared')[0].strip()}.png")  # PNG 파일로 저장

        plt.cla()
        plt.clf()
        plt.close()


x_col = 'eGFR'
# x_col = 'eGFRxTBIL'
effect_col = 'EFFECT1'

X = nss_df[[x_col]]
X_const = sm.add_constant(X)
y_vals = nss_df[effect_col]

model = sm.OLS(y_vals, X_const).fit()

intercept, slope = model.params
r_squared = model.rsquared
p_value = model.pvalues[x_col]


## Fitting

from scipy.optimize import curve_fit

# Sigmoid Emax 모델 함수 정의
def sigmoid_emax(conc, E0, Emax, EC50, H):
    return E0+Emax * conc**H / (EC50**H + conc**H)


x_sigemax = nss_df[x_col]
effect = nss_df[effect_col]

# curve fitting
popt, pcov = curve_fit(sigmoid_emax, x_sigemax, effect, p0=[0.02, 0.5, 50, 1])  # 초기값 설정

# 피팅된 파라미터 출력
popt = np.array([0.0251, 0.435, 59.6, 7.36])
E0_fit, Emax_fit, EC50_fit, H_fit = popt
print(f"E0: {E0_fit:.2f}, Emax: {Emax_fit:.2f}, EC50: {EC50_fit:.2f}, Hill coefficient: {H_fit:.2f}")

# 적합도 판별

# 예측값
y_pred = sigmoid_emax(x_sigemax, * popt)

# R-squared 계산
ss_res = np.sum((y_vals - y_pred) ** 2)
ss_tot = np.sum((y_vals - np.mean(y_vals)) ** 2)
r_squared = 1 - (ss_res / ss_tot)

# RMSE
rmse = np.sqrt(np.mean((y_vals - y_pred) ** 2))

# AIC (Optional)
n = len(y_vals)
k = len(popt)
aic = n * np.log(ss_res / n) + 2 * k

print(f"Effect fitting / R-squared: {r_squared:.4f}")
print(f"Effect fitting / RMSE: {rmse:.4f}")
print(f"Effect fitting / AIC: {aic:.2f}")

# 예측값 계산
gfr_fit = np.linspace(0.1, int(np.max(x_sigemax))+1, 100)
effect_fit = sigmoid_emax(gfr_fit, * popt)

# 시각화

y = effect_col

plt.figure(figsize=(15, 10))
plt.scatter(x_sigemax, effect, label='Ovserved', color='black', marker='o')
plt.plot(gfr_fit, effect_fit, label='Model-fitting', color='royalblue')
plt.xlabel(x_col, fontsize=14)
plt.ylabel(y, fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# fig_title = f'[WSCT] Sigmoid Emax Model Fit ({x_col} vs {y})\nE0: {E0_fit:.3f}, Emax: {Emax_fit:.3f}, EC50: {EC50_fit:.2f}, Hill coefficient: {H_fit:.2f}\nR-squared: {r_squared:.4f}, RMSE: {rmse:.4f}, AIC: {aic:.2f}'
fig_title = f'[WSCT] Sigmoid Emax Model Fit ({x_col} vs {y})\nE0: {E0_fit:.3f}, Emax: {Emax_fit:.3f}, EC50: {EC50_fit:.2f}, Hill coefficient: {H_fit:.2f}, OFV: -126.849'
plt.title(fig_title, fontsize=14)
plt.legend(fontsize=14)
plt.grid(True)
plt.show()

plt.savefig(f"{results_dir_path}/{fig_title.split('E0:')[0].strip()}.png")  # PNG 파일로 저장

plt.cla()
plt.clf()
plt.close()

# ## PyMC (NLME fitting)
#
# import pymc as pm
# import numpy as np
#
# # 데이터 준비
# # N: 총 데이터 포인트 수
# # J: 총 개인 수
# # subj: 각 데이터 포인트의 개인 인덱스 (0부터 J-1까지의 정수)
# # conc: 농도 데이터 (길이 N의 배열)
# # effect: 효과 데이터 (길이 N의 배열)
#
# with pm.Model() as model:
#     # Hyperpriors for group means
#     E0 = pm.Normal('E0', mu=0, sigma=1)
#     Emax = pm.Normal('Emax', mu=0, sigma=1)
#     EC50 = pm.Normal('EC50', mu=0, sigma=1)
#     H = pm.Normal('H', mu=1, sigma=0.5)
#     sigma = pm.HalfNormal('sigma', sigma=1)
#
#     # Between-subject variability
#     omega_Emax = pm.HalfNormal('omega_Emax', sigma=1)
#     omega_EC50 = pm.HalfNormal('omega_EC50', sigma=1)
#
#     # Individual parameters
#     eta_Emax = pm.Normal('eta_Emax', mu=0, sigma=1, shape=J)
#     eta_EC50 = pm.Normal('eta_EC50', mu=0, sigma=1, shape=J)
#
#     Emax_ind = pm.Deterministic('Emax_ind', Emax * pm.math.exp(omega_Emax * eta_Emax))
#     EC50_ind = pm.Deterministic('EC50_ind', EC50 * pm.math.exp(omega_EC50 * eta_EC50))
#
#     # Expected value
#     mu = E0 + Emax_ind[subj] * conc ** H / (EC50_ind[subj] ** H + conc ** H)
#
#     # Likelihood
#     effect_obs = pm.Normal('effect_obs', mu=mu, sigma=sigma, observed=effect)
#
#     # Sampling
#     trace = pm.sample(2000, chains=4)
#
# # 결과 출력
# pm.summary(trace)

## Data For PD-Endpoint Simulation

nss_df = nss_df.reset_index(drop=True)
nss_df['dUGEc'] = nss_df[effect_col]
nss_df['ID'] = nss_df.index
nss_df[['ID','UID','GRP','eGFR','TBIL','ALT','AUClast','eGFRxTBIL','dUGEc','PG_AVG','PG_ZERO','HbA1c_base']].to_csv(f"{results_dir_path}/KSCPTSPR25_PD_Endpoint_Sim_data.csv", index=False)


ss_df = ss_df.reset_index(drop=True)
ss_df['dUGEc'] = ss_df[effect_col]
ss_df['ID'] = ss_df.index
ss_df[['ID','UID','GRP','eGFR','TBIL','ALT','AUClast','eGFRxTBIL','dUGEc','PG_AVG','PG_ZERO','HbA1c_base']].to_csv(f"{results_dir_path}/KSCPTSPR25_PD_Endpoint_Sim_data_MM.csv", index=False)


# nss_df['EFFECT0'].median()
# nss_df['EFFECT1'].median()

hue = 'GRP'

for reltup in [('PG_ZERO', 'PG_AVG'),('eGFR', 'TBIL'),('eGFR', 'ALT')]:
    x = reltup[0]
    y = reltup[1]

    X = nss_df[[x]].copy()
    X_const = sm.add_constant(X).applymap(float)
    y_vals = nss_df[y].map(float)

    model = sm.OLS(y_vals, X_const).fit()

    intercept, slope = model.params
    r_squared = model.rsquared
    p_value = model.pvalues[x]


    sns.scatterplot(data=nss_df, x=x, y=y, hue=hue)
    fig_title = f'[WSCT] {x} vs {y} by {hue}\nR-squared:{r_squared:.4f}, p-value: {p_value:.4f}\nbeta: {slope:.4f}, intercept: {intercept:.4f} '
    plt.title(fig_title)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.grid(True)
    plt.tight_layout()
    # plt.show()

    plt.savefig(f"{results_dir_path}/{fig_title.split('R-squared')[0].strip()}.png")  # PNG 파일로 저장

    plt.cla()
    plt.clf()
    plt.close()
