## 시간에 따른 모든 사람의 PD 경향성 그려보기

from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr, mannwhitneyu
import statsmodels.api as sm
# from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## 모델링 데이터셋 로딩

df = pd.read_csv(f'{output_dir}/modeling_df_covar/infliximab_integrated_modeling_df(for pda).csv')
# df['ADA']
count_df = df.groupby('ID',as_index=False).agg({'SEX':'mean','IBD_TYPE':'mean'})
count_df['IBD_TYPE'] = count_df['IBD_TYPE'].map({0.0:'CD',1.0:'UC'})
count_df['SEX'] = count_df['SEX'].map({0.0:'Male',1.0:'Female'})
count_df = count_df.groupby(['SEX','IBD_TYPE'],as_index=False).agg({'ID':'count'})
f_ibd_tot = count_df[count_df['SEX']=='Female']['ID'].sum()
m_ibd_tot = count_df[count_df['SEX']=='Male']['ID'].sum()
count_df = count_df.append(pd.DataFrame([{'SEX':'Female','IBD_TYPE':'ALL_IBD','ID':f_ibd_tot},
                              {'SEX':'Male','IBD_TYPE':'ALL_IBD','ID':m_ibd_tot}
                              ]))


count_df = count_df.rename(columns={'ID':'N'})
count_df['N_TOT'] = count_df['SEX'].map({'Female':f_ibd_tot,'Male':m_ibd_tot})
count_df = count_df.reset_index(drop=True)
count_df.loc[count_df[(count_df['SEX']=='Female')&(count_df['IBD_TYPE']=='ALL_IBD')].index[0], 'N_TOT'] = f_ibd_tot+m_ibd_tot
count_df.loc[count_df[(count_df['SEX']=='Male')&(count_df['IBD_TYPE']=='ALL_IBD')].index[0], 'N_TOT'] = f_ibd_tot+m_ibd_tot

count_df['pct'] = 100*count_df['N']/count_df['N_TOT']
count_df['N'] = count_df.apply(lambda x:f"{x['N']} ({round(x['pct'],1)})", axis=1)
count_df = count_df.rename(columns={'IBD_TYPE':'ITEM'})
count_df = count_df[['SEX','ITEM','N']].pivot(index='ITEM', columns='SEX', values='N').copy()
# count_df
# count_df = count_df[['N']].T.copy()
count_df.columns.name = None
count_df.index.name = 'ITEM'
count_df = count_df.reset_index(drop=False)

count_pval_df = pd.DataFrame([{'ITEM':'ALL', 'p_value':np.nan},
                              {'ITEM':'CD', 'p_value':np.nan},
                              {'ITEM':'UC', 'p_value':np.nan}])



## Demographics

demo_df = df.groupby('ID',as_index=False).agg({'SEX':'mean','AGE':'first','BMI':'median','WT':'median','HT':'median', })
demo_df = demo_df.melt(id_vars=['ID','SEX'], value_vars=['AGE', 'BMI','WT','HT'], var_name='ITEM', value_name="VALUE",)

# p_value
results = []
for item, subdf in demo_df.groupby('ITEM'):
    group0 = subdf.loc[subdf['SEX'] == 0, 'VALUE']
    group1 = subdf.loc[subdf['SEX'] == 1, 'VALUE']

    # Mann-Whitney U test
    stat, p = mannwhitneyu(group0, group1, alternative='two-sided')

    results.append({'ITEM': item, 'p_value': p})

demo_pval_df = pd.DataFrame(results)

# est. vals

demo_df['VALUE2'] = demo_df['VALUE'].copy()
demo_df = demo_df.groupby(['SEX','ITEM'],as_index=False).agg({'VALUE':'mean','VALUE2':'std'})
demo_df['VALUE_STR'] = demo_df.apply(lambda x:f"{round(x['VALUE'],1)} ({round(x['VALUE2'],1)})", axis=1)
demo_df['SEX'] = demo_df['SEX'].map({0.0:'Male',1.0:'Female'})
demo_df = demo_df[['SEX','ITEM','VALUE_STR']].pivot(index='ITEM', columns='SEX', values='VALUE_STR').copy()
demo_df.columns.name = None
demo_df.index.name = 'ITEM'
demo_df = demo_df.reset_index(drop=False)

## DOSE

dose_df = df[df['MDV']==1].copy()
dose_df['AMT'] = dose_df['AMT'].astype(float)
dose_df['WT'] = dose_df['WT'].astype(float)
dose_df = dose_df.groupby(['ID','SEX','ROUTE'],as_index=False).agg({'AMT':'mean','WT':'mean'})
dose_df['DOSE_MEAN'] = dose_df['AMT'].copy()
dose_df['DOSE_SD'] =dose_df['DOSE_MEAN'].copy()
dose_df['DOSEpWT_MEAN'] = dose_df['AMT']/dose_df['WT']
dose_df['DOSEpWT_SD'] =dose_df['DOSEpWT_MEAN'].copy()

dose_pval_df = dose_df[['ID','SEX','ROUTE','DOSE_MEAN','DOSEpWT_MEAN']].melt(id_vars=['ID','SEX','ROUTE'], value_vars=['DOSE_MEAN','DOSEpWT_MEAN'], var_name='ITEM', value_name="VALUE",)
dose_pval_df['ROUTE'] = dose_pval_df['ROUTE'].map({'1':'IV','2':'SC'})
dose_pval_df['SEX'] = dose_pval_df['SEX'].map({0:'Male',1:'Female'})
# dose_pval_df = dose_pval_df.groupby(['ID','SEX','ROUTE','ITEM'],as_index=False).agg({'VALUE':'mean'})

results = []
for route in ['IV','SC']:
    for item, subdf in dose_pval_df.groupby('ITEM'):
        group0 = subdf.loc[(subdf['SEX'] == 'Male')&(subdf['ROUTE'] == route), 'VALUE'].copy()
        group1 = subdf.loc[(subdf['SEX'] == 'Female')&(subdf['ROUTE'] == route), 'VALUE'].copy()

        # Mann-Whitney U test
        stat, p = mannwhitneyu(group0, group1, alternative='two-sided')

        results.append({'ITEM': item, 'ROUTE':route,'p_value': p})

dose_pval_df = pd.DataFrame(results)
dose_pval_df['ITEM'] = dose_pval_df.apply(lambda x:f"{x['ITEM'].replace('_MEAN','')} ({x['ROUTE']})", axis=1)
dose_pval_df = dose_pval_df.drop(['ROUTE'],axis=1)


dose_df = dose_df.groupby(['SEX','ROUTE'],as_index=False).agg({'DOSE_MEAN':'mean','DOSEpWT_MEAN':'mean','DOSE_SD':'std','DOSEpWT_SD':'std'})
dose_df['DOSE'] = dose_df.apply(lambda x:f"{round(x['DOSE_MEAN'],1)} ({round(x['DOSE_SD'],1)})", axis=1)
dose_df['DOSEpWT'] = dose_df.apply(lambda x:f"{round(x['DOSEpWT_MEAN'],1)} ({round(x['DOSEpWT_SD'],1)})", axis=1)
dose_df['ROUTE'] = dose_df['ROUTE'].map({'1':'IV','2':'SC'})
dose_df['SEX'] = dose_df['SEX'].map({0:'Male',1:'Female'})
dose_df = dose_df[['SEX','ROUTE','DOSE','DOSEpWT']].melt(id_vars=['SEX','ROUTE'], value_vars=['DOSE', 'DOSEpWT'], var_name='ITEM', value_name="VALUE_STR",)
dose_df = dose_df.pivot(index=['ROUTE','ITEM'], columns='SEX', values='VALUE_STR').reset_index(drop=False)
dose_df['ITEM'] = dose_df.apply(lambda x:f"{x['ITEM']} ({x['ROUTE']})", axis=1)
dose_df = dose_df[['ITEM','Female','Male']].copy()

## Lab
# df['ADA']
lab_df = df.groupby(['ID','SEX'],as_index=False).agg({'ALB':'first', 'AST':'first', 'ALT':'first', 'CRP':'first', 'FCAL':'first', 'CREATININE':'first', 'ADA':'max'}).copy()
lab_df['SEX'] = lab_df['SEX'].map({0:'Male',1:'Female'})
co_lab_df = lab_df.drop(['ADA'],axis=1)


# 분석할 변수 리스트
vars_to_test = ['ALB', 'AST', 'ALT', 'CRP', 'FCAL', 'CREATININE']

# 결과 저장 리스트
lab_results = []
lab_pval_df = list()
for var in vars_to_test:
    # 그룹별 데이터
    g0 = co_lab_df.loc[lab_df['SEX'] == 'Male', var].dropna()
    g1 = co_lab_df.loc[lab_df['SEX'] == 'Female', var].dropna()

    # 평균(SD)
    mean0, sd0 = np.mean(g0), np.std(g0, ddof=1)
    mean1, sd1 = np.mean(g1), np.std(g1, ddof=1)

    # 비모수 검정 (Mann–Whitney U test)
    stat, pval = mannwhitneyu(g0, g1, alternative='two-sided')

    lab_results.append({
        'ITEM': var,
        'Female': f"{mean1:.1f} ({sd1:.1f})",
        'Male': f"{mean0:.1f} ({sd0:.1f})",
    })
    lab_pval_df.append({
        'ITEM': var,
        'p_value': pval
    })

# 결과 DataFrame
lab_results = pd.DataFrame(lab_results)
lab_pval_df = pd.DataFrame(lab_pval_df)
# print(result_df)
## ADA

ca_lab_df = lab_df[['SEX','ADA']].copy()
# 교차표 생성
ct = pd.crosstab(ca_lab_df['SEX'], ca_lab_df['ADA'])

# Fisher exact test (2x2 테이블 전제)
# 만약 ADA 값이 0/1 외에 다른 경우가 있다면 사전에 처리 필요
oddsratio, pval = fisher_exact(ct)

# 각 그룹별 n과 % 계산
summary = (
    ca_lab_df.groupby('SEX')['ADA']
    .agg(['sum', 'count'])
    .rename(columns={'sum': 'ADA_Pos', 'count': 'Total'})
)

summary['Rate'] = summary['ADA_Pos'] / summary['Total'] * 100

# 결과 정리
ca_result = {
    'ITEM': 'ADA, n (%)',
    'Male': f"{int(summary.loc['Male', 'ADA_Pos'])} ({summary.loc['Male', 'Rate']:.1f}%)",
    'Female': f"{int(summary.loc['Female', 'ADA_Pos'])} ({summary.loc['Female', 'Rate']:.1f}%)",
}
ca_pval_df = {'ITEM': 'ADA, n (%)', 'p_value': pval}

ca_res_df = pd.DataFrame([ca_result])
ca_pval_df = pd.DataFrame([ca_pval_df])


total_df = pd.concat([count_df, demo_df, lab_results, dose_df, ca_res_df, ], ignore_index=True)
total_pval_df = pd.concat([count_pval_df,demo_pval_df,lab_pval_df,dose_pval_df, ca_pval_df, ], ignore_index=True)
total_df = total_df.merge(total_pval_df, on='ITEM', how='left')
total_df['sig'] = (total_df['p_value']<0.05).map({True:'*',False:''})
total_df['p_value'] = total_df['p_value'].map(lambda x:round(x,4)).replace(0.0000,'<0.0001')
total_df['ITEM'] = total_df['ITEM'].replace('ALL_IBD','IBD, n (%)').replace('CD','CD, n (%)').replace('UC','UC, n (%)')
total_df.to_csv(f"{output_dir}/sex_diff_res.csv", index=False, encoding='utf-8-sig')
