## 시간에 따른 개별 농도 및 PD 경향성 그려보기

from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

## DEMO Covariates Loading

# simulation_df = pd.read_csv(f"{output_dir}/infliximab_integrated_modeling_df_dayscale.csv")
simulation_df = pd.read_csv(f"{output_dir}/modeling_df_covar/infliximab_integrated_pdeda_df_dayscale.csv")
# simulation_df['TIME_DIFF'] = simulation_df['TIME'].diff()
# simulation_df[simulation_df['TIME_DIFF']>0]['TIME_DIFF'].sort_values()
# interval = 3/24
interval = 1
for_sim_df = get_model_population_sim_df(df=simulation_df, interval=interval, add_on_period=0)
# for_sim_df = for_sim_df[for_sim_df['PD_CALPRTSTL'].isna()].replace(np.nan,'.')
for_sim_df = for_sim_df.replace(np.nan,'.')

# for_sim_df = for_sim_df[['ID', 'TIME', 'DV', 'MDV', 'AMT', 'DUR', 'CMT', 'ROUTE', 'IBD_TYPE', 'ALB', 'ADA', 'AGE', 'SEX', 'WT', 'HT', 'BMI', 'REALDATA']].copy()
for_sim_df = for_sim_df[['ID', 'TIME', 'DV', 'MDV', 'AMT', 'DUR', 'CMT', 'ROUTE', 'IBD_TYPE', 'ALB', 'ADA', 'AGE', 'SEX', 'WT', 'HT', 'BMI', 'PD_INDEXISTS', 'PD_PRO2', 'PD_PRO2_BL', 'PD_PRO2_DELT', 'PD_CR', 'PD_CR_BL', 'PD_CR_DELT', 'PD_CRP', 'PD_CRP_BL', 'PD_CRP_DELT', 'PD_FCAL', 'PD_FCAL_BL', 'PD_FCAL_DELT', 'REALDATA', 'RATE', 'TAD']].copy()
# for_sim_df = for_sim_df[['ID', 'TIME', 'DV', 'MDV', 'AMT', 'DUR', 'CMT', 'ROUTE', 'IBD_TYPE', 'ALB', 'ADA', 'AGE', 'SEX', 'WT', 'HT', 'BMI', 'PD_CRP', 'PD_CALPRTSTL','PD_PRO2', 'REALDATA', 'RATE', 'TAD']].copy()
for_sim_df.to_csv(f"{output_dir}/modeling_df_covar/infliximab_integrated_simulation_df.csv",index=False, encoding='utf-8-sig')


# for_sim_df.columns
# str(for_sim_df.columns).replace("', '"," ").replace("',\n       '"," ")
# final_sim_df = pd.read_csv(f"{nonmem_dir}/run/sim57",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
final_sim_df = pd.read_csv(f"{nonmem_dir}/run/sim60",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
realworld_df = simulation_df[simulation_df['MDV']==0].copy()
# final_sim_df['PD_CR']
# final_sim_df['IPRED']
# final_sim_df['DV']
# final_sim_df.columns
ibd_type_dict = {1:'CD',2:'UC'}
for uid, id_sim_df in final_sim_df.groupby('ID'): #break

    # raise ValueError
    #
    # gdf[(gdf['REALDATA']==1)&(gdf['AMT']==0)]

    gdf = id_sim_df[id_sim_df['MDV']==0][['ID', 'TIME', 'DV', 'IPRED', 'AMT', 'IBD_TYPE', 'PD_FCAL','PD_CRP', 'PD_PRO2', 'PD_CR','REALDATA']].copy()
    adf = id_sim_df[id_sim_df['MDV']==1][['ID', 'TIME', 'DV', 'IPRED', 'AMT', 'IBD_TYPE', 'ROUTE', 'REALDATA']].copy()
    adf['ROUTE'] = adf['ROUTE'].map({1.0:'IV',2.0:'SC'})
    adf['AMT'] = adf['AMT'].astype(int).astype(str)
    adf['DOSE_INFO'] = adf['ROUTE']+adf['AMT']
    adf['DOSE_REGIMEN_CHANGE'] = adf['DOSE_INFO']!=adf['DOSE_INFO'].shift(1)
    adf['DOSE_INFO'] = adf.apply(lambda x:x['DOSE_INFO'] if x['DOSE_REGIMEN_CHANGE'] else '*' ,axis=1)

    rdf = realworld_df[realworld_df['ID']==uid].copy()
    rdf['DV'] = rdf['DV'].astype(float)
    # gdf['DV'] = 1
    ibd_type = ibd_type_dict[gdf['IBD_TYPE'].iloc[0]]

    # TIME이 정렬되어 있어야 깔끔하게 보임
    gdf = gdf.sort_values(by='TIME')

    # subplot 준비 (2행 1열, sharex를 통해 TIME 축 공유)
    # fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 10), sharex=True)
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(15, 10), sharex=True)

    # Real data point 추가
    sns.scatterplot(data=rdf, x='TIME', y='DV', ax=ax1, color='red', s=50, marker='o')

    # 첫 번째 그래프: TIME vs DV
    y_col = 'IPRED'
    # y_col = 'DV'

    sns.lineplot(data=gdf, x='TIME', y=y_col, ax=ax1, color='steelblue')
    # sns.lineplot(data=gdf, x='TIME', y='DV', ax=ax1, color='steelblue')
    ax1.axhline(y=5, color='red', linestyle='--', linewidth=1.5)  # y=5 수평선
    ax1.set_title(f"[ID({int(uid)})_{ibd_type}] CONC over TIME")
    ax1.set_ylabel("CONC")
    ax1.grid(True,linestyle='--')

    # Dosing에 대한 Annotation
    for _, row in adf.iterrows():
        time = row['TIME']
        label = row['DOSE_INFO']

        # gdf에서 가장 가까운 TIME의 DV 값 찾기
        closest_idx = (gdf['TIME'] - time).abs().idxmin()
        # y = gdf.loc[closest_idx, 'DV']
        y = gdf.loc[closest_idx, y_col]

        # annotate (텍스트 + 화살표)
        ax1.annotate(
            text=label,
            xy=(time, y),  # 화살표 끝 (데이터 지점)
            xytext=(time, y + 2),  # 텍스트 위치 (조정 가능)
            arrowprops=dict(arrowstyle='->', color='black'),
            fontsize=9,
            ha='center'
        )


    pd_gdf = gdf[gdf['REALDATA']==1].copy()
    # pd_gdf = gdf.copy()

    # 두 번째 그래프: PD_PRO2 (왼쪽 y축)
    sns.scatterplot(data=pd_gdf, x='TIME', y='PD_PRO2', ax=ax2, color='darkorange', label='PD_PRO2', s=100)

    ax2.set_title(f"{ibd_type}_PRO2 and CRP over TIME")
    ax2.set_ylabel(f"{ibd_type}_PRO2")
    ax2.grid(True,linestyle='--')

    # ax2 범례 제거
    if ax2.get_legend():
        ax2.get_legend().remove()

    # 두 번째 그래프: CRP (오른쪽 y축)
    ax2_right = ax2.twinx()
    sns.lineplot(data=gdf, x='TIME', y='PD_CRP', ax=ax2_right, color='green', label='CRP')
    ax2_right.axhline(y=0.5, color='green', linestyle='--', linewidth=1.5)  # y=1 수평선 (CRP)
    ax2_right.set_ylabel("CRP")

    # 범례를 하나로 합치기 (ax2와 ax2_right)
    lines_1, labels_1 = ax2.get_legend_handles_labels()
    lines_2, labels_2 = ax2_right.get_legend_handles_labels()
    ax2_right.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

    # 세 번째 그래프: PD_CR (왼쪽 y축)
    sns.scatterplot(data=pd_gdf, x='TIME', y='PD_CR', ax=ax3, color='brown', label='CR', s=100)
    ax3.set_ylim(-0.3, 2.0)  # 원하는 최소값, 최대값으로 지정
    ax3.axhline(y=0, color='brown', linestyle='--', linewidth=1.5)  # y=1 수평선 (CRP)
    ax3.set_title(f"{ibd_type}_CR and FCAL over TIME")
    ax3.set_ylabel(f"{ibd_type}_CR")
    ax3.grid(True, linestyle='--')


    # ax2 범례 제거
    if ax3.get_legend():
        ax3.get_legend().remove()

    # 세 번째 그래프: FCAL (오른쪽 y축)
    ax3_right = ax3.twinx()
    sns.lineplot(data=gdf, x='TIME', y='PD_FCAL', ax=ax3_right, color='grey', label='FCAL')
    # ax3_right.set_ylim(-0.3)  # 원하는 최소값, 최대값으로 지정
    ax3_right.axhline(y=200, color='grey', linestyle='--', linewidth=1.5)  # y=1 수평선 (CRP)
    ax3_right.axhline(y=-0.3, color='grey', linestyle='--', linewidth=0.0)  # y=1 수평선 (CRP)
    ax3_right.set_ylabel("FCAL")

    # 범례를 하나로 합치기 (ax3와 ax3_right)
    lines_3, labels_3 = ax3.get_legend_handles_labels()
    lines_4, labels_4 = ax3_right.get_legend_handles_labels()
    ax3_right.legend(lines_3 + lines_4, labels_3 + labels_4, loc='upper right')


    plt.tight_layout()
    # plt.show()
    if not os.path.exists(f"{output_dir}/PKPD_EDA"):
        os.mkdir(f"{output_dir}/PKPD_EDA")
    if not os.path.exists(f"{output_dir}/PKPD_EDA/idv_trends_{y_col}"):
        os.mkdir(f"{output_dir}/PKPD_EDA/idv_trends_{y_col}")
    plt.savefig(f"{output_dir}/PKPD_EDA/idv_trends_{y_col}/IBDPGx_ID({int(uid)})IBD({ibd_type}).png")  # PNG 파일로 저장

    plt.cla()
    plt.clf()
    plt.close()


    # final_sim_df



# final_sim_df = pd.read_csv(f"{output_dir}/infliximab_integrated_simulation_df.csv")
# final_sim_df
# final_sim_df['PD_PRO2']

# uniq_ids = list(final_sim_df['ID'].unique())[:30]
# final_sim_df[final_sim_df['ID'].isin(uniq_ids)].to_csv(f"{output_dir}/amk_simulation_partial_df.csv",index=False, encoding='utf-8-sig')

# final_sim_df['ID'].drop_duplicates()
# final_sim_df['TAD'] = add_time_after_dosing_column(df=final_sim_df)



