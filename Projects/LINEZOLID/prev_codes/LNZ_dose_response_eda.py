from tools import *
from pynca.tools import *
from datetime import datetime, timedelta
from scipy.stats import spearmanr



prj_name = 'LINEZOLID'
prj_dir = './Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
# nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

dose_df = pd.read_csv(f"{output_dir}/lnz_final_dose_df.csv")
lab_df = pd.read_csv(f"{output_dir}/lnz_final_lab_df.csv")

# lab_df.columns

dose_df['INTERVAL_HOUR'] = dose_df['INTERVAL'].map(lambda x:int(x.replace('q','').replace('h','')))
dose_df['DOSEpHOUR'] = dose_df['DOSE']/dose_df['INTERVAL_HOUR']

dose_res_df = dose_df.merge(lab_df, on=['UID','DATE'], how='left')
dose_res_df = dose_res_df.dropna(axis=1, how='all').drop('Blast', axis=1)

pd_list = list(dose_res_df.loc[:,'ANC':].columns)

pk_col_list = ['DOSE','DOSEpHOUR']
pd_col_list = pd_list

# scatter plot 그리기
plt.figure(figsize=(15, 12))

if not os.path.exists(f"{output_dir}/PKPD_EDA"):
    os.mkdir(f"{output_dir}/PKPD_EDA")
for x in pk_col_list:
    # x='ALT'
    for y in pd_col_list:
        hue = 'IBD_TYPE'
        # sim_conc_df[[x,y]].dropna()
        gdf = dose_res_df.copy()
        gdf = gdf.replace([np.inf, -np.inf], np.nan).dropna(subset=[x, y])
        X = gdf[[x]].copy()
        X_const = sm.add_constant(X, has_constant='add').applymap(float)
        y_vals = gdf[y].map(float)

        model = sm.OLS(y_vals, X_const).fit()

        intercept, slope = model.params
        r_squared = model.rsquared
        p_value = model.pvalues[x]

        corr, corr_pval = spearmanr(dose_res_df[x], dose_res_df[y])

        sns.scatterplot(data=gdf, x=x, y=y, marker='o', legend=False)
        fig_title = f'[PKPD_EDA] {x} vs {y}\ncorr_coef: {corr:.4f}, p-value: {corr_pval:.4f}'
        # fig_title = f'[PKPD_EDA] {x} vs {y}\nR-squared:{r_squared:.4f}, p-value: {p_value:.4f}\nbeta: {slope:.4f}, intercept: {intercept:.4f}'
        plt.title(fig_title, fontsize=14)
        plt.xlabel(x, fontsize=14)
        plt.ylabel(y, fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # plt.legend(fontsize=14)
        plt.legend().remove()
        plt.grid(True)
        plt.tight_layout()
        # plt.show()

        # plt.savefig(f"{output_dir}/PKPD_EDA/{fig_title.split('R-squared')[0].strip()}.png")  # PNG 파일로 저장
        plt.savefig(f"{output_dir}/PKPD_EDA/{fig_title.split('corr_coef')[0].strip()}.png")  # PNG 파일로 저장

        plt.cla()
        plt.clf()
        plt.close()

# dose_df.columns