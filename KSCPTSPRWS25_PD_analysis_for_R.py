from tools import *
from pynca.tools import *
from statsmodels.stats.outliers_influence import variance_inflation_factor

wd_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25"
input_dir_path = f"{wd_path}/sglt2i_dataset"
results_dir_path = f"{wd_path}/results_r"

pk_nca_df = pd.read_csv(results_dir_path+'/[WSCT] NCARes_PK.csv')
pd_df = pd.read_csv(input_dir_path+'/KSCPTSPRWS25_SGLT2i_PD.csv')
pd_df = pd_df.merge(pk_nca_df, on=['ID','GRP'], how='left')

# scatter plot 그리기

plt.figure(figsize=(15, 12))

for x in ['AUClast', 'eGFR', 'TBIL', 'ALT']:
    # x='ALT'
    y='EFFECT1'
    hue = 'GRP'

    X = pd_df[[x]].copy()
    X_const = sm.add_constant(X).applymap(float)
    y_vals = pd_df[y].map(float)

    model = sm.OLS(y_vals, X_const).fit()

    intercept, slope = model.params
    r_squared = model.rsquared
    p_value = model.pvalues[x]


    sns.scatterplot(data=pd_df, x=x, y=y, hue=hue, marker='o')
    fig_title = f'[WSCT] {x} vs {y} by {hue}\nR-squared:{r_squared:.4f}, p-value: {p_value:.4f}\nbeta: {slope:.4f}, intercept: {intercept:.4f} '
    plt.title(fig_title, fontsize=14)
    plt.xlabel(x, fontsize=14)
    plt.ylabel(y, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(f"{results_dir_path}/{fig_title.split('R-squared')[0].strip()}.png")  # PNG 파일로 저장

    plt.cla()
    plt.clf()
    plt.close()