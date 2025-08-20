## PK parameter와 PD marker 관계 그려보기 - CD/UC 및 IND/MAINT 서브그룹 나눠서


from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr



prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## DEMO Covariates Loading

# datacheck_df = pd.read_csv(f"{output_dir}/infliximab_integrated_pdeda_df_dayscale.csv")
simulation_df = pd.read_csv(f"{output_dir}/modeling_df_covar/infliximab_integrated_pdeda_df_dayscale.csv")
final_sim_df = pd.read_csv(f"{nonmem_dir}/run/sim60",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')

# sim_conc_df = final_sim_df[final_sim_df['MDV']==1].copy()
# realworld_df = simulation_df[simulation_df['MDV']==0].copy()

ibd_type_dict = {1:'CD',2:'UC'}

# final_sim_df.columns
# # final_sim_df
# final_sim_df['ROUTE_SHIFT1'] = final_sim_df['ROUTE'].shift(1)
# final_sim_df['ROUTE_SHIFT2'] = final_sim_df['ROUTE'].shift(2)
# final_sim_df['AMT_SHIFT1'] = final_sim_df['AMT'].shift(1)
# final_sim_df['AMT_SHIFT2'] = final_sim_df['AMT'].shift(2)
# final_sim_df['']

final_sim_df['AMT'] = final_sim_df['AMT'].replace(0.0,np.nan)





# DOSING INTERVAL 채워넣기

sim_dosing_df = final_sim_df[(final_sim_df['MDV']==1)&(final_sim_df['REALDATA']==1)].copy()

dosing_df_list = list()
for id, frag_sim_dosing_df in sim_dosing_df.groupby('ID'):
    frag_sim_dosing_df['DOSING_INTERVAL'] = frag_sim_dosing_df['TIME'].diff().fillna(method='bfill')
    dosing_df_list.append(frag_sim_dosing_df)

sim_dosing_df = pd.concat(dosing_df_list)
# sim_dosing_df['DOSING_INTERVAL'].describe()
sim_dosing_df = sim_dosing_df[sim_dosing_df['DOSING_INTERVAL'] >= 3].copy()

# sim_dosing_df[sim_dosing_df['ID']==74]

dosing_interval_df = sim_dosing_df[['ID','TIME','MDV','DOSING_INTERVAL']].copy()

final_sim_df = final_sim_df.merge(dosing_interval_df, on=['ID','TIME','MDV'], how='left')

# 채혈데이터에 AMOUNT 채워넣기

sim_df_list = list()
for id, frag_sim_df in final_sim_df.groupby('ID'):
    frag_sim_df['AMT'] = frag_sim_df['AMT'].fillna(method='ffill')
    frag_sim_df['DOSING_INTERVAL'] = frag_sim_df['DOSING_INTERVAL'].fillna(method='ffill').fillna(method='bfill')
    sim_df_list.append(frag_sim_df)

final_sim_df = pd.concat(sim_df_list, ignore_index=True)


# NCA로 Cumulative AUC 삽입
# final_sim_sampling_df = final_sim_df[final_sim_df['MDV']==0].copy()
# final_sim_dosing_df = final_sim_df[final_sim_df['MDV']==1].copy()
# nca_result = tblNCA(final_sim_sampling_df, key=["ID"], colTime="TIME", colConc="IPRED", dose=1, adm="Extravascular", dur=0, doseUnit="mg", timeUnit="day", concUnit="ug/mL", down="Log", R2ADJ=0, MW=0, SS=False, iAUC="", excludeDelta=1, slopeMode='SNUHCPT', colStyle='pw')
# nca_result = tblNCA(final_sim_sampling_df, key=["ID"], colTime="TIME", colConc="IPRED", dose=1, adm="Extravascular", dur=0, doseUnit="mg", timeUnit="day", concUnit="ug/mL", down="Log", R2ADJ=0, MW=0, SS=False, iAUC="", excludeDelta=1, slopeMode='BEST', colStyle='pw')


# sim_conc_df['DOSING_INTERVAL'].min()

sim_conc_df = final_sim_df[(final_sim_df['TIME']!=0)&(final_sim_df['MDV']==0)&(final_sim_df['REALDATA']==1)].copy()
sim_conc_df['DAILY_DOSE'] = sim_conc_df['AMT']/sim_conc_df['DOSING_INTERVAL']
sim_conc_df['AUC24'] = sim_conc_df['DAILY_DOSE']/sim_conc_df['CL']

ibd_type_dict = {1:'CD',2:'UC'}
sim_conc_df['IBD_TYPE'] = sim_conc_df['IBD_TYPE'].map(ibd_type_dict)
# sim_conc_df[sim_conc_df['DOSING_INTERVAL']!=0]['DOSING_INTERVAL'].min()
# final_sim_df[(final_sim_df['MDV']==1)&(final_sim_df['REALDATA']==1)].copy()
# for i in range(1,8):
#     # print(i)
#     sim_conc_df[f'DV{i}DS'] =

# final_sim_df[['ID','TIME','DV','MDV','PD_PRO2','REALDATA']]
# sim_conc_df = sim_conc_df[['ID','IBD_TYPE','TIME','DV','PD_PRO2']].copy()

# cov_cand = ['AGE', 'SEX', 'HT', 'WT', 'BMI', 'ALB', 'eGFR', 'AST', 'ALT', 'SODIUM', 'TBIL', 'AUClast']
# corr_matrix = pd_df[cov_cand].corr().abs()
#
# # 히트맵 그리기
# plt.figure(figsize=(18, 15))
# ax = sns.heatmap(
#     corr_matrix,
#     annot=True,
#     fmt=".2f",
#     cmap='coolwarm',
#     vmin=-1,
#     vmax=1,
#     annot_kws={"size": 12}
# )
# fig_title = "[WSCT] Correlation matrix heatmap of the covariates"
# plt.title(fig_title, fontsize=15)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
#
# # colorbar 글씨 크기 조정
# ax.collections[0].colorbar.ax.tick_params(labelsize=15)
#
# plt.savefig(f"{results_dir_path}/{fig_title}.png")
# plt.cla()
# plt.clf()
# plt.close()



# 변수 간 spearman 상관계수 및 p-value 계산
# corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['PD_PRO2'])
# print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")
#
# corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['CRP'])
# print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")
#
# corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['CALPRTSTL'])
# print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")
#
# corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['SEX'])
# print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")
#
# corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['AGE'])
# print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")
#
# corr, pval = spearmanr(sim_conc_df['DV'], sim_conc_df['HT'])
# print(f"Spearman correlation: {corr:.3f}, p-value: {pval:.4f}")

if not os.path.exists(f"{output_dir}/PKPD_EDA"):
    os.mkdir(f"{output_dir}/PKPD_EDA")
if not os.path.exists(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups"):
    os.mkdir(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups")

pk_col_list = ['DV','DAILY_DOSE','AUC24','IPRED']
pd_col_list = ['PD_CRP','PD_CRP_DELT','PD_FCAL','PD_FCAL_DELT','PD_PRO2','PD_PRO2_DELT','PD_CR','PD_CR_DELT']

# scatter plot 그리기
plt.figure(figsize=(15, 12))

# mgdf['IBD_TYPE']
for x in pk_col_list:
    for y in pd_col_list:
        mgdf = sim_conc_df.copy()
        for ibd_type, ibd_mgdf in mgdf.groupby('IBD_TYPE'):

            ibd_meddf = ibd_mgdf[ibd_mgdf['PD_INDEXISTS']==1].copy()
            ibd_ind_df = ibd_meddf[ibd_meddf['TIME'] < 128].copy()
            ibd_maint_df = pd.concat([ibd_meddf[ibd_meddf['TIME'] >= 128].copy(), ibd_mgdf[ibd_mgdf['PD_INDEXISTS']==0]])
            phase_df_dict = {'IND':ibd_ind_df, 'MAINT':ibd_maint_df}
            # ibd_mgdf[ibd_mgdf['PD_INDEXISTS']==1]
            # ibd_mgdf[ibd_mgdf['PD_INDEXISTS']==0]

            for ind_phase, gdf in phase_df_dict.items(): #break

                gdf = gdf.replace([np.inf, -np.inf], np.nan).dropna(subset=[x, y])
                X = gdf[[x]].copy()
                X_const = sm.add_constant(X).applymap(float)
                y_vals = gdf[y].map(float)

                model = sm.OLS(y_vals, X_const).fit()

                intercept, slope = model.params
                r_squared = model.rsquared
                p_value = model.pvalues[x]

                corr, corr_pval = spearmanr(gdf[x], gdf[y])

                # sns.scatterplot(data=gdf, x=x, y=y, marker='o', hue=hue)
                sns.scatterplot(data=gdf, x=x, y=y, marker='o',legend=False)
                fig_title = f'[PKPD_EDA ({ibd_type}_{ind_phase})] {x} vs {y}\ncorr_coef: {corr:.4f}, p-value: {corr_pval:.4f}'
                # fig_title = f'[PKPD_EDA] {x} vs {y}\nR-squared:{r_squared:.4f}, p-value: {p_value:.4f}\nbeta: {slope:.4f}, intercept: {intercept:.4f}'
                plt.title(fig_title, fontsize=14)
                plt.xlabel(x, fontsize=14)
                plt.ylabel(y, fontsize=14)
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                # plt.legend(fontsize=14)
                # plt.legend().remove()
                plt.grid(True)
                plt.tight_layout()
                # plt.show()

                if not os.path.exists(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups/{ind_phase}"):
                    os.mkdir(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups/{ind_phase}")
                if not os.path.exists(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups/{ind_phase}/{ibd_type}"):
                    os.mkdir(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups/{ind_phase}/{ibd_type}")
                if not os.path.exists(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups/{ind_phase}/{ibd_type}/{x}"):
                    os.mkdir(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups/{ind_phase}/{ibd_type}/{x}")

                # plt.savefig(f"{output_dir}/PKPD_EDA/{fig_title.split('R-squared')[0].strip()}.png")  # PNG 파일로 저장
                plt.savefig(f"{output_dir}/PKPD_EDA/PKvsPD_Corr_Subgroups/{ind_phase}/{ibd_type}/{x}/{fig_title.split('corr_coef')[0].strip()}.png")  # PNG 파일로 저장

                plt.cla()
                plt.clf()
                plt.close()
