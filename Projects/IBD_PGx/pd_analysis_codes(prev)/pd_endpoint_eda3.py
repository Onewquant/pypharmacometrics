## PK parameter와 PD marker 관계 그려보기


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

sim_conc_df = add_iAUC(final_sim_df[(final_sim_df['MDV']==0)].copy(), time_col="TIME", conc_col="DV", id_col="ID")  # ID가 없으면 자동으로 전체 처리
sim_dose_df = final_sim_df[(final_sim_df['MDV']==1)].copy()
sim_dose_df['iAUC']	= np.nan
sim_dose_df['cumAUC'] = np.nan
final_sim_df = pd.concat([sim_conc_df,sim_dose_df]).sort_values(['ID','TIME','MDV'])

# sim_conc_df = final_sim_df[final_sim_df['MDV']==1].copy()
# realworld_df = simulation_df[simulation_df['MDV']==0].copy()
# final_sim_df[(final_sim_df['ID']==0)&(final_sim_df['MDV']==0)].to_csv(f"{output_dir}/modeling_df_covar/conc_test.csv", index=False,encoding='utf-8-sig' )

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


if not os.path.exists(f"{output_dir}/PKPD_EDA"):
    os.mkdir(f"{output_dir}/PKPD_EDA")
if not os.path.exists(f"{output_dir}/PKPD_EDA/PKvsPD_Corr"):
    os.mkdir(f"{output_dir}/PKPD_EDA/PKvsPD_Corr")

pk_col_list = ['DV','DAILY_DOSE','AUC24','IPRED','cumAUC']
pd_col_list = ['PD_CRP','PD_CRP_DELT','PD_FCAL','PD_FCAL_DELT','PD_PRO2','PD_PRO2_DELT','PD_CR']

# scatter plot 그리기
# plt.figure(figsize=(15, 12))


for x in pk_col_list:
    # x='ALT'
    # x='cumAUC'
    for y in pd_col_list:
        # y='PD_CRP_DELT'
        hue = 'IBD_TYPE'
        # sim_conc_df[[x,y]].dropna()
        gdf = sim_conc_df.copy()
        # gdf = sim_conc_df.copy()
        y_base = '_'.join(y.split('_')[:2])
        if y_base=='PD_CRP':
            gdf = gdf[gdf[f"{y_base}_BL"] > 0.5].copy()
        elif y_base=='PD_FCAL':
            gdf = gdf[gdf[f"{y_base}_BL"] > 50].copy()
        elif y_base=='PD_PRO2':
            gdf = gdf[gdf[f"{y_base}_BL"] >= 1].copy()
        elif y_base=='PD_CR':
            gdf = gdf[gdf[f"{y_base}_BL"] == 0].copy()
        #     # raise ValueError
        # gdf.columns
        gdf = gdf.replace([np.inf, -np.inf], np.nan).dropna(subset=[x, y])
        X = gdf[[x]].copy()
        X_const = sm.add_constant(X).applymap(float)
        y_vals = gdf[y].map(float)

        model = sm.OLS(y_vals, X_const).fit()

        intercept, slope = model.params
        r_squared = model.rsquared
        p_value = model.pvalues[x]

        corr, corr_pval = spearmanr(gdf[x], gdf[y])

        if y_base == 'PD_CR':
            sns.scatterplot(data=gdf, x=x, y=y, marker='o', hue=hue)
        else:
            sns.lmplot(
                data=gdf, x=x, y=y,
                height=5, aspect=1.2,
                scatter_kws={"s": 40},  # 점 크기
                line_kws={"linewidth": 2},  # 회귀선 두께
                hue = hue
            )


        # fig_title = f'[PKPD_EDA] {x} vs {y}\ncorr_coef: {corr:.4f}, p-value: {corr_pval:.4f}'
        fig_title = f'[PKPD_EDA] {x} vs {y}\nR-squared:{r_squared:.4f}, p-value: {p_value:.4f}\nbeta: {slope:.4f}, intercept: {intercept:.4f}'
        plt.title(fig_title, fontsize=14)
        plt.xlabel(x, fontsize=14)
        plt.ylabel(y, fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=14)
        # plt.legend().remove()
        plt.grid(True)
        plt.tight_layout()
        # plt.show()

        # plt.savefig(f"{output_dir}/PKPD_EDA/{fig_title.split('R-squared')[0].strip()}.png")  # PNG 파일로 저장

        if not os.path.exists(f"{output_dir}/PKPD_EDA/PKvsPD_Corr/{x}"):
            os.mkdir(f"{output_dir}/PKPD_EDA/PKvsPD_Corr/{x}")
        # plt.savefig(f"{output_dir}/PKPD_EDA/PKvsPD_Corr/{x}/{fig_title.split('corr_coef')[0].strip()}.png")  # PNG 파일로 저장
        plt.savefig(f"{output_dir}/PKPD_EDA/PKvsPD_Corr/{x}/{fig_title.split('R-squared')[0].strip()}.png")  # PNG 파일로 저장

        plt.cla()
        plt.clf()
        plt.close()
