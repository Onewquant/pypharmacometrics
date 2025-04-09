from tools import *
from pynca.tools import *
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor


resource_dir_path = "/resource/KSCPTSPRWS25/resource"
results_dir_path = "/resource/KSCPTSPRWS25/results"

########################################
## MAD 분석
########################################

m_pg_df = pd.read_csv(f'{resource_dir_path}/healthy_data/KSCPTSPRWS2025(healthy)_PG(MAD).csv')
m_auc_df = pd.read_csv(f'{resource_dir_path}/healthy_data/KSCPTSPRWS25(healthy)_AUC24.csv')
m_uge_df = pd.read_excel(f'{resource_dir_path}/healthy_data/KSCPTSPRWS25(healthy)_PD_prep_MAD.xlsx',sheet_name='Urine PD')

m_pg_df['ID'] = m_pg_df['Randomization No.'].map(lambda x:x[-5:])
m_pg_df = m_pg_df[['ID'] + [c for c in m_pg_df.columns if c.split('D')[0] in ['1','2','14','15','16']]].copy()
m_pg_df = m_pg_df.melt(id_vars=['ID'],var_name='DT',value_name='PG')
m_pg_df['TIME'] = m_pg_df['DT'].map(lambda x:(float(x.split('D')[0])-1)*24+float(x.split('D')[-1][:-1]))
m_pg_df = m_pg_df[['ID','TIME','PG']].sort_values(['ID','TIME'],ignore_index=True)

nss_mpg_df = m_pg_df[m_pg_df['TIME']<=24].copy()
ss_mpg_df = m_pg_df[m_pg_df['TIME']>=333].copy()
ss_mpg_df['TIME']=ss_mpg_df['TIME']-336

# iAUC_df=pd.DataFrame({'Name':['T[0-6h]','T[6-12h]','T[12-24h]','T[0-24h]'], 'Start':[0,6,12,0], 'End':[6,12,24,24]})
# result = tblNCA(concData=mdprep_conc_df,key=['ID','DAY'],colTime="A_TIME",colConc='CONC',down = "Linear",iAUC=iAUC_df,dose=0.5,slopeMode='BEST',colStyle='pw')
ss_gluauc_result = tblNCA(concData=ss_mpg_df,key=['ID'],colTime="TIME",colConc='PG',down = "Log",dose=1,slopeMode='BEST',colStyle='pw')
nss_gluauc_result = tblNCA(concData=nss_mpg_df,key=['ID'],colTime="TIME",colConc='PG',down = "Log",dose=1,slopeMode='BEST',colStyle='pw')

ss_gluauc_result = ss_gluauc_result.rename(columns={'AUClast':'AUC_GLU'}).iloc[1:,:]
nss_gluauc_result = nss_gluauc_result.rename(columns={'AUClast':'AUC_GLU'}).iloc[1:,:]

ss_gluauc_result['PG_AVG'] = ss_gluauc_result['AUC_GLU'].astype(float)/24
nss_gluauc_result['PG_AVG'] = nss_gluauc_result['AUC_GLU'].astype(float)/24

ss_gluauc_result = ss_gluauc_result[['ID','AUC_GLU','PG_AVG']]
nss_gluauc_result = nss_gluauc_result[['ID','AUC_GLU','PG_AVG']]

# ss_mpg_df = ss_mpg_df.merge(ss_gluauc_result[['ID','AUC_GLU','PG_AVG']], on=['ID'], how='left')
# nss_mpg_df = nss_mpg_df.merge(nss_gluauc_result[['ID','AUC_GLU','PG_AVG']], on=['ID'], how='left')

m_uge_df = m_uge_df[['Subject','C_TAD','Cumul_C','steady']].iloc[1:].dropna()
m_uge_df = m_uge_df[m_uge_df['C_TAD']==24].reset_index(drop=True)
m_uge_df = m_uge_df.pivot(index=['Subject','C_TAD'],columns='steady',values='Cumul_C').reset_index(drop=False)
m_uge_df.columns.name = None
m_uge_df['ID'] = m_uge_df['Subject'].map(lambda x:str(int(x)))
m_uge_df['TIME'] = m_uge_df['C_TAD']
m_uge_df = m_uge_df[['ID','TIME','s','ss']].rename(columns={'s':'NSS','ss':'SS'})
m_uge_df['NSS'] = m_uge_df['NSS'].astype(float)
m_uge_df['SS'] = m_uge_df['SS'].astype(float)
nss_muge_df = m_uge_df[['ID','TIME','NSS']].copy()
nss_muge_df = nss_muge_df.rename(columns={'NSS':'NSS_UGE24'})
ss_muge_df = m_uge_df[['ID','TIME','SS']].copy()
ss_muge_df = ss_muge_df.rename(columns={'SS':'SS_UGE24'})

ss_pd_df = ss_muge_df.merge(ss_gluauc_result, on=['ID'], how='left')
ss_pd_df['SS_EFFECT1'] = ss_pd_df['SS_UGE24']/ss_pd_df['PG_AVG']
nss_pd_df = nss_muge_df.merge(nss_gluauc_result, on=['ID'], how='left')
nss_pd_df['NSS_EFFECT1'] = nss_pd_df['NSS_UGE24']/nss_pd_df['PG_AVG']

mpd_df = ss_pd_df[['ID','TIME','SS_UGE24','SS_EFFECT1']].merge(nss_pd_df[['ID','TIME','NSS_UGE24','NSS_EFFECT1']], on=['ID','TIME'], how='left')
# 상관계수 및 p-value 계산

gdf = mpd_df
x = 'NSS_EFFECT1'
y = 'SS_EFFECT1'

r, p_value = stats.pearsonr(mpd_df[x], mpd_df[y])
r_squared = r**2

# 산점도 + 회귀선 그래프 그리기
plt.figure(figsize=(20, 10))
sns.regplot(data=mpd_df, x=x, y=y, ci=95, line_kws={'color': 'royalblue'}, color='black')

# R²와 p-value 텍스트로 그래프에 추가
plt.title(f'Correlation Plot (N={len(mpd_df)})\nR = {r:.3f}, $R^2$ = {r_squared:.3f}, p-value = {p_value:.3g}', fontsize=14)
plt.xlabel(x)
plt.ylabel(y)
plt.grid(True)
plt.tight_layout()

plt.savefig(f"{results_dir_path}/[Healthy (MAD)] {x} vs {y}.png")  # PNG 파일로 저장

plt.cla()
plt.clf()
plt.close()

########################################
## SAD 분석
########################################

rand_df = pd.read_excel(f'{resource_dir_path}/healthy_data/KSCPTSPRWS25_RDZ_LIST_SAD.XLSX')
auc_df = pd.read_csv(f'{resource_dir_path}/healthy_data/KSCPTSPRWS25(healthy)_AUC24.csv')
pg_df = pd.read_csv(f'{resource_dir_path}/healthy_data/KSCPTSPRWS2025(healthy)_PG(SAD).csv')
uge_df = pd.read_excel(f'{resource_dir_path}/healthy_data/KSCPTSPRWS25(healthy)_PD_prep_SAD.xlsx',sheet_name='Urine')


rand_df = rand_df[rand_df['TREAT']==1].copy()
rand_df['ID'] = rand_df['RANDNO'].map(lambda x:x[-4:])
rand_df['DOSE'] = rand_df['TREATMENT'].map(lambda x:x.split(' ')[1])
rand_df = rand_df[['ID','DOSE']].copy()

auc_df['ID'] = auc_df['ID'].astype(str)
auc_df = auc_df[auc_df['Dose'].map(lambda x:x[-1])=='S'].copy()
auc_df['DOSE'] = auc_df['Dose'].map(lambda x:float(x[:-1]))
auc_df = auc_df[['ID','AUC_tau','F4','F5']].copy()
auc_df = auc_df.rename(columns={'AUC_tau':'AUClast'})

pg_df['ID'] = pg_df['Randomization No.'].map(lambda x:x[-4:])
pg_df = pg_df[['ID'] + [c for c in pg_df.columns if c[:2]in ['1D','2D']]].copy()
pg_df = pg_df.melt(id_vars=['ID'],var_name='DT',value_name='PG')
pg_df['TIME'] = pg_df['DT'].map(lambda x:(float(x[0])-1)*24+float(x[2:-1]))
pg_df = pg_df[['ID','TIME','PG']].sort_values(['ID','TIME'],ignore_index=True)

gluauc_result = tblNCA(concData=pg_df,key=['ID'],colTime="TIME",colConc='PG',down = "Log",dose=1,slopeMode='BEST',colStyle='pw')
gluauc_result = gluauc_result.rename(columns={'AUClast':'AUC_GLU'}).iloc[1:,:]
gluauc_result['PG_AVG'] = gluauc_result['AUC_GLU'].astype(float)/24

uge_df = uge_df.iloc[1:].copy()
uge_df['ID'] = uge_df['Subject'].map(lambda x:str(int(x)))
uge_df = uge_df[uge_df['TRT'].map(lambda x:x[:3]=='DWP')]
uge_df = uge_df[uge_df['N_Time']==24][['ID','Cumul_C', 'N_Time']]
uge_df = uge_df.rename(columns={'Cumul_C':'UGE24', 'N_Time':'TIME'}).reset_index(drop=True)


hpd_df = uge_df.merge(gluauc_result, on=['ID'], how='left')
hpd_df = hpd_df.merge(auc_df, on=['ID'], how='left')
hpd_df = hpd_df.merge(rand_df, on=['ID'], how='left')

hpd_df['EFFECT1'] = hpd_df['UGE24']/hpd_df['PG_AVG']
hpd_df['CL'] =(hpd_df['DOSE'].astype(float))/hpd_df['AUClast']
hpd_df['CL_DOSE'] =hpd_df['CL']/(hpd_df['DOSE'].astype(float))
hpd_df['AUClast_DOSE'] =hpd_df['AUClast']/(hpd_df['DOSE'].astype(float))
hpd_df = hpd_df.sort_values(['DOSE'])

hue_order = list(hpd_df['DOSE'].sort_values().unique())
# sns.histplot(hpd_df, x='EFFECT1', bins=30)
# 평균과 표준편차 계산
mean = hpd_df['EFFECT1'].mean()
std = hpd_df['EFFECT1'].std()
q25 = hpd_df['EFFECT1'].quantile(0.25)
q75 = hpd_df['EFFECT1'].quantile(0.75)

# ±2 standard deviation 범위

lower_bound = mean - 2 * std
upper_bound = mean + 2 * std

for x in ['AUClast']:
    for y in ['EFFECT1', 'UGE24']:

        plt.figure(figsize=(20, 10))
        # x = 'AUClast'
        # y = 'EFFECT1'
        # # y = 'UGE24'
        # # y = 'EFFECTdelta'
        hue = 'DOSE'

        # sns.scatterplot(data=ss_df, x='AUClast', y='EFFECT0', hue='GRP')
        # sns.scatterplot(data=ss_df, x='AUClast', y='EFFECT1', hue='GRP')
        # sns.scatterplot(data=nss_df, x='AUClast', y='EFFECT0', hue='GRP')
        # sns.scatterplot(data=nss_df, x='AUClast', y='EFFECT1', hue='GRP')
        sns.scatterplot(data=hpd_df, x=x, y=y, hue=hue,  hue_order=hue_order)
        plt.title(f'{x} vs {y} by {hue}')
        plt.xlabel(x)
        plt.ylabel(y)
        plt.grid(True)
        plt.tight_layout()
        # plt.show()

        plt.savefig(f"{results_dir_path}/[Healthy (SAD)] {x} vs {y}.png")  # PNG 파일로 저장

        plt.cla()
        plt.clf()
        plt.close()


categories = hpd_df["DOSE"].unique()

# 원하는 컬러맵에서 n개의 색상 추출
cmap = plt.get_cmap("Greys")
colors = [cmap(i / len(categories)) for i in range(len(categories))]

# 범주에 색상 매핑
palette = dict(zip(categories, colors))

for x in ['DOSE']:
    for y in ['AUClast', 'AUClast_DOSE', 'CL', 'CL_DOSE']:

        plt.figure(figsize=(20, 10))
        # x = 'AUClast'
        # y = 'EFFECT1'
        # # y = 'UGE24'
        # # y = 'EFFECTdelta'
        hue = 'DOSE'

        # sns.scatterplot(data=ss_df, x='AUClast', y='EFFECT0', hue='GRP')
        # sns.scatterplot(data=ss_df, x='AUClast', y='EFFECT1', hue='GRP')
        # sns.scatterplot(data=nss_df, x='AUClast', y='EFFECT0', hue='GRP')
        # sns.scatterplot(data=nss_df, x='AUClast', y='EFFECT1', hue='GRP')
        # sns.scatterplot(data=hpd_df, x=x, y=y, hue=hue,  hue_order=hue_order)
        sns.boxplot(data=hpd_df, x=x, y=y, linewidth=2, hue=hue,  hue_order=hue_order, palette=palette)
        plt.title(f'{x} vs {y}')
        plt.xlabel(x)
        plt.ylabel(y)
        plt.grid(True)
        plt.tight_layout()
        # plt.show()

        plt.savefig(f"{results_dir_path}/[Healthy (SAD)] {x} vs {y}.png")  # PNG 파일로 저장

        plt.cla()
        plt.clf()
        plt.close()