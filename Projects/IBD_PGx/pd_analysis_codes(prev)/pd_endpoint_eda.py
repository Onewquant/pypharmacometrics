## 시간에 따른 모든 사람의 PD 경향성 그려보기

from tools import *
from pynca.tools import *
from datetime import datetime, timedelta

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

drug = 'infliximab'
mode_str = 'integrated'
# df = pd.read_csv(f'{output_dir}/{drug}_{mode_str}_pdeda_df_dayscale.csv')
df = pd.read_csv(f'{output_dir}/modeling_df_covar/{drug}_{mode_str}_datacheck(covar).csv')

# df[['ID','TIME','IBD_TYPE','PD_MARKER']]
df['DAY'] = df['TIME']/24
# df[['PD_MARKER']]
# pd_marker = 'PD_PRO2'

# df
# str(None)
# 그래프 그리기
# df[['PD_MARKER']]
# gdf['PRO2']
# pdmarker_dict = {'CD':'CDAI','UC':'MES'}
if not os.path.exists(f"{output_dir}/PKPD_EDA"):
    os.mkdir(f"{output_dir}/PKPD_EDA")
if not os.path.exists(f"{output_dir}/PKPD_EDA/PD_IDV_TRENDS"):
    os.mkdir(f"{output_dir}/PKPD_EDA/PD_IDV_TRENDS")

ibd_type_list = ['CD','UC']
pdmarker_list = ['PD_PRO2','PD_CR','PD_CRP','PD_FCAL']
for pdmarker in pdmarker_list:
    for ibd_type in ibd_type_list:
        ibdtype_gdf = df[df['IBD_TYPE']==ibd_type].copy()
        # for id_len in [10, None]:
        for id_len in [None,]:
            gdf = ibdtype_gdf[ibdtype_gdf['ID'].isin(list(ibdtype_gdf['ID'].unique())[:id_len])]

            plt.figure(figsize=(20, 10))
            sns.lineplot(data=gdf, x='DAY', y=pdmarker, hue='ID', marker='o')
            plt.title(f'[{ibd_type}] TIME vs {pdmarker} by ID')
            plt.xlabel('DAYS')
            # plt.ylabel('PD Marker')
            plt.ylabel(f"{pdmarker} ({ibd_type})")
            plt.legend().remove()
            plt.grid(True)



            plt.savefig(f'{output_dir}/PKPD_EDA/PD_IDV_TRENDS/[{ibd_type}] DAYS_VS_{pdmarker.replace("PD_","")}({str(id_len).replace("None","All")}).png', dpi=600, bbox_inches='tight')

            plt.cla()
            plt.clf()
            plt.close()

    # # 그래프 저장
    # plt.savefig('time_vs_pd_marker_by_id.png', dpi=300, bbox_inches='tight')



