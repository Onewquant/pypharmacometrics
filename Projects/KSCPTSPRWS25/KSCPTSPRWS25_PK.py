import os.path

from tools import *
from pynca.tools import *

# df = pd.read_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/prep_data/KSCPT2025SPRWS_ConcPrep_SGLT2INH(R).csv")
# df[df['GRP'].isin([1,3])].reset_index(drop=True).to_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/prep_data/KSCPTSPRWS25MULTI_ConcPrep_SGLT2INH(R).csv", index=False, encoding='utf-8-sig')
# df[df['GRP'].isin([2,4])].reset_index(drop=True).to_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/prep_data/KSCPTSPRWS25SINGLE_ConcPrep_SGLT2INH(R).csv", index=False, encoding='utf-8-sig')

results_dir_path = "/resource/KSCPTSPRWS25/results"
ws_dataset_dir_path = "/resource/KSCPTSPRWS25/sglt2i_dataset"
if not os.path.exists(ws_dataset_dir_path):
    os.mkdir(ws_dataset_dir_path)

## 기본정보입력

project_dict = {'KSCPTSPRWS25MULTI':['SGLT2INH'], 'KSCPTSPRWS25SINGLE':['SGLT2INH']}
modeling_dir_path = "/resource/KSCPTSPRWS25"
prepconc_dir_path = "/resource/KSCPTSPRWS25/prep_data"

## Modeling and Dosing Policy 파일 불러오기

mdpolicy_file_path = f'{modeling_dir_path}/Modeling_and_Dosing_Policy - MDP.csv'
dspol_df = modeling_dosing_policy(mdpolicy_file_path, selected_models=[], model_colname='MODEL')

## Project, Drug 별로 Conc data 모으기

drugconc_dict=get_drug_conc_data_dict_of_multiple_projects(project_dict, prepconc_dir_path=prepconc_dir_path, conc_filename_format="[project]_ConcPrep_[drug](R).csv")
# drugconc_dict['SGLT2INH'].columns


# 개별 농도그래프

conc_df = drugconc_dict['SGLT2INH'].copy()
conc24_df = conc_df[conc_df['NTIME'] <= 24].copy()
# conc24_df.columns

# 대상자 ID별로 고유한 마커 생성
gdf = conc24_df.copy()
gdf[['ID', 'GRP', 'NTIME', 'ATIME', 'CONC']].to_csv(f'{ws_dataset_dir_path}/KSCPTSPRWS25_SGLT2i_PK.csv',index=False)
unique_ids = gdf['ID'].unique()
# markers = ['o', 's', 'D', '^', 'v', '<', '>', 'P', 'X', '*', '+', '1', '2', '3', '4', '8']
# marker_map = {uid: markers[i % len(markers)] for i, uid in enumerate(unique_ids)}

# GRP에 따른 색상 구분을 위해 palette 설정
palette = sns.color_palette("tab10", gdf['GRP'].nunique())

# 플롯 그리기
plt.figure(figsize=(15, 12))

# 대상자별로 라인을 그리고 마커는 다르게, 색은 GRP 기준
for uid in unique_ids:
    subset = gdf[gdf['ID'] == uid]
    grp = subset['GRP'].iloc[0]
    sns.lineplot(
        data=subset,
        x='ATIME',
        y='CONC',
        color=palette[int(grp) - 1],  # GRP가 1~4라고 가정
        # marker=marker_map[uid],
        marker = 'o',
        markers=True,
        markersize=10,
        dashes=False,
        legend=False,
    )

# 범례는 GRP 기준으로 수동 생성
from matplotlib.lines import Line2D

legend_elements = [
    Line2D([0], [0], color=palette[i], lw=2, label=f'GRP {i+1}')
    for i in range(gdf['GRP'].nunique())
]

plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=15)
fig_title = '[WSCT] Individual Time-Concentration Profiles by GRP'
plt.title(fig_title, fontsize=15)
plt.xlabel('Time (ATIME)', fontsize=15)
plt.ylabel('Concentration (CONC)', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(True)
plt.tight_layout()
# plt.show()

plt.savefig(f"{results_dir_path}/{fig_title}.png")  # PNG 파일로 저장

plt.cla()
plt.clf()
plt.close()

## Population

plt.figure(figsize=(15, 12))

# 2. 그룹별 Mean ± SD 계산
summary_df = gdf.groupby(['GRP', 'NTIME']).agg(
    mean_CONC=('CONC', 'mean'),
    sd_CONC=('CONC', 'std')
).reset_index()

# 3. 그룹별 평균 라인 & SD 범위 그리기
for grp in summary_df['GRP'].unique():
    grp_df = summary_df[summary_df['GRP'] == grp]

    # Mean line
    plt.plot(
        grp_df['NTIME'],
        grp_df['mean_CONC'],
        color=palette[int(grp) - 1],
        label=f'GRP {grp} Mean',
        linewidth=2.5,
        marker='o',            # ★ marker 추가
        markersize=10,          # ★ marker 크기
    )

    # SD 범위 (음영)
    plt.fill_between(
        grp_df['NTIME'],
        grp_df['mean_CONC'] - grp_df['sd_CONC'],
        grp_df['mean_CONC'] + grp_df['sd_CONC'],
        color=palette[int(grp) - 1],
        alpha=0.2
    )

# 범례
legend_elements = [
    Line2D([0], [0], color=palette[i], lw=2, label=f'GRP {i + 1} Mean ± SD')
    for i in range(gdf['GRP'].nunique())
]
plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=13)

# 기타 설정
fig_title = '[WSCT] Population Time-Concentration Profiles by GRP'
plt.title(fig_title, fontsize=15)
plt.xlabel('Time (NTIME)', fontsize=15)
plt.ylabel('Concentration (CONC)', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

# 저장
plt.savefig(f"{results_dir_path}/{fig_title}.png", dpi=300)

plt.cla()
plt.clf()
plt.close()


## NCA

nca_result = tblNCA(concData=conc24_df ,key=['ID','GRP'],colTime="ATIME",colConc='CONC',down = "Log", dose=0.5,slopeMode='BEST',colStyle='pw')
nca_result = nca_result.iloc[1:].copy()
nca_result[['ID','GRP','N_Samples','Tmax', 'Cmax', 'AUClast']].to_excel(f"{results_dir_path}/[WSCT] NCARes_PK.xlsx", index=False)

# iauc_res_df[iauc_res_df['DAY']==1][['ID', 'N_Samples', 'Tmax', 'Cmax', 'AUClast', 'Vz_F_obs', 'Cl_F_obs']]

## ANOVA 비교

import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import t


# 결과 저장용 리스트
results = []

# 분석할 파라미터 리스트
columns = ['AUClast', 'Cmax', 'Tmax']

# GRP가 범주형인지 확인 및 변환
nca_result['GRP'] = nca_result['GRP'].astype(int).astype('category')

for col in columns:
    # 로그 변환
    nca_result[f'log_{col}'] = np.log(nca_result[col].astype(float))

    # ANOVA
    model = ols(f'log_{col} ~ C(GRP)', data=nca_result).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    pval = anova_table['PR(>F)'][0]

    coef = model.params
    se = model.bse
    df_resid = model.df_resid
    t_val = t.ppf(0.975, df_resid)

    # GMR 및 CI 계산
    summary = []
    for grp in nca_result['GRP'].cat.categories:
        if grp == 1:
            summary.append(f'GRP 1: reference')
            continue
        coef_name = f'C(GRP)[T.{grp}]'
        log_gmr = coef.get(coef_name, 0)
        log_se = se.get(coef_name, 0)
        ci_lower = log_gmr - t_val * log_se
        ci_upper = log_gmr + t_val * log_se

        gmr = np.exp(log_gmr)
        gmr_ci_lower = np.exp(ci_lower)
        gmr_ci_upper = np.exp(ci_upper)

        summary.append(f'GRP {grp}: {gmr:.3f} ({gmr_ci_lower:.3f}–{gmr_ci_upper:.3f})')

    results.append({
        'Parameter': col,
        'ANOVA p-value': pval,
        'GMR (95% CI) by GRP (vs GRP 1)': '; '.join(summary)
    })

# 결과 테이블 생성
result_df = pd.DataFrame(results)
result_df.columns = ['Parameter', 'ANOVA p-value', 'GMR (95% CI) by GRP (vs GRP 1)']
result_df.to_excel(f"{results_dir_path}/[WSCT] ANOVA_PK.xlsx", index=False)

# 출력
print(result_df)
