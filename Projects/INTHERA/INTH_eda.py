import numpy as np
import pandas as pd
import seaborn as sns
from pynca.tools import *

prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/INTHERA'
df = pd.read_csv(f"{prj_dir}/in1011_id.csv")

# comp_dict = {1:'Blood',2:'cDRG',3:'lDRG',4:'tDRG',5:'Heart',6:"Kidney",7:'Brain',8:'Liver',9:'Lung',10:'Stomach',11:'Femur'}
# comp_dict = {1:'Blood',2:'cDRG',3:'lDRG',4:'tDRG',5:'Heart',6:'Intestine',7:"Kidney",8:'Brain',9:'Liver',10:'Lung',11:'Stomach',12:'Femur',13:'UB',14:'GB'}
comp_dict = {5:'Blood',2:'cDRG',3:'lDRG',4:'tDRG',1:'Heart',6:'Intestine',7:"Kidney",8:'Brain',9:'Liver',10:'Lung',11:'Stomach',12:'Femur',13:'UB',14:'GB'}

# rev_comp_dict =

gdf = df.copy()
gdf['CMT'] = gdf['CMT'].map(comp_dict)

# ## 그래프 그리기
#
# pdf = gdf[(gdf['MDV']!=1)&(gdf['DV']!='0')].copy()
# # pdf['CMT'] = pdf['CMT'].map(comp_dict)
# pdf['DV'] = pdf['DV'].map(float)
# pdf['TIME'] = pdf['TIME'].map(float)

# pdf.to_csv(f"{prj_dir}/in1011_graph.csv", encoding='utf-8-sig', index=False)

##

# fdf = gdf[~gdf['CMT'].isin(['Intestine', 'GB', 'UB', "Kidney"])]
fdf = gdf[~gdf['CMT'].isin(['Intestine', 'GB', 'UB', 'Blood'])]
# for id, id_df in gdf.groupby('ID'):
#     id_df.sort_values(['CMT'])

uniue_dict = {c:(inx+1) for inx, c in enumerate(list(fdf['CMT'].unique()))}
fdf['CMT'] = fdf['CMT'].map(uniue_dict).copy()
fdf.sort_values(['ID','CMT','TIME'])
fdf.to_csv(f"{prj_dir}/in1011_id4.csv", encoding='utf-8-sig', index=False)


wdf = fdf[(fdf['MDV']!=1)&(fdf['DV']!='0')].groupby('CMT', as_index=False).agg({'DV':np.median})
for inx, c in enumerate(wdf['DV']/20):
    print(f"(0, {round(c,3)})")
    print(f"(0, {round(c,3)})")


fnca = fdf[(fdf['MDV']==0)].copy()
fnca['DV'] = fnca['DV'].astype(float)
fnca = fnca.groupby(['CMT'],as_index=False).agg({'DV':'var'})
# fnca['DV'] = fnca['DV']/10**5

for inx, c in enumerate(fnca['DV']/10**5):
    print(f"(0, {round(c,3)})")
    print(f"(0, {round(c,3)})")

#
#
# result = tblNCA(fnca, key=["ID"], colTime="TIME", colConc="DV",
#                 dose=370, adm="BOLUS", dur=0, doseUnit="mg",
#                 timeUnit="h", concUnit="ug/mL", down="Log", R2ADJ=0,
#                 MW=0, SS=False, iAUC="", excludeDelta=1, slopeMode='SNUHCPT', colStyle='pw')
#
