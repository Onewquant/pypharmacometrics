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
df = pd.read_csv(f'{output_dir}/{drug}_{mode_str}_datacheck(covar).csv')

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
pdmarker_dict = {'CD':'PD_PRO2','UC':'PD_PRO2'}
for ibd_type, pdmarker in pdmarker_dict.items():
    ibdtype_gdf = df[df['IBD_TYPE']==ibd_type].copy()
    for id_len in [10, None]:
        gdf = ibdtype_gdf[ibdtype_gdf['ID'].isin(list(ibdtype_gdf['ID'].unique())[:id_len])]

        plt.figure(figsize=(20, 10))
        sns.lineplot(data=gdf, x='DAY', y=pdmarker, hue='ID', marker='o')
        plt.title(f'TIME vs {pdmarker} by ID')
        plt.xlabel('DAYS')
        # plt.ylabel('PD Marker')
        plt.ylabel(pdmarker)
        plt.legend().remove()
        plt.grid(True)

        if not os.path.exists(f"{output_dir}/PKPD_EDA"):
            os.mkdir(f"{output_dir}/PKPD_EDA")

        plt.savefig(f'{output_dir}/PKPD_EDA/[{ibd_type}] DAYS_VS_{pdmarker.replace("PD_","")}({str(id_len).replace("None","All")}).png', dpi=600, bbox_inches='tight')

        plt.cla()
        plt.clf()
        plt.close()

# # 그래프 저장
# plt.savefig('time_vs_pd_marker_by_id.png', dpi=300, bbox_inches='tight')



