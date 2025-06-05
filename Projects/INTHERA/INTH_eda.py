import numpy as np
import pandas as pd
import seaborn as sns

prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/INTHERA'
df = pd.read_csv(f"{prj_dir}/in1011_id.csv")

# comp_dict = {1:'Blood',2:'cDRG',3:'lDRG',4:'tDRG',5:'Heart',6:"Kidney",7:'Brain',8:'Liver',9:'Lung',10:'Stomach',11:'Femur'}
comp_dict = {1:'Blood',2:'cDRG',3:'lDRG',4:'tDRG',5:'Heart',6:'Intestine',7:"Kidney",8:'Brain',9:'Liver',10:'Lung',11:'Stomach',12:'Femur',13:'UB',14:'GB'}
# rev_comp_dict =

gdf = df.copy()
gdf['CMT'] = gdf['CMT'].map(comp_dict)

## 그래프 그리기

pdf = gdf[(gdf['MDV']!=1)&(gdf['DV']!='0')].copy()
# pdf['CMT'] = pdf['CMT'].map(comp_dict)
pdf['DV'] = pdf['DV'].map(float)
pdf['TIME'] = pdf['TIME'].map(float)

# pdf.to_csv(f"{prj_dir}/in1011_graph.csv", encoding='utf-8-sig', index=False)

##

fdf = gdf[~gdf['CMT'].isin(['Intestine', 'GB', 'UB', "Kidney"])]
# for id, id_df in gdf.groupby('ID'):
#     id_df.sort_values(['CMT'])

uniue_dict = {c:(inx+1) for inx, c in enumerate(list(fdf['CMT'].unique()))}
fdf['CMT'] = fdf['CMT'].map(uniue_dict).copy()
fdf.to_csv(f"{prj_dir}/in1011_id3.csv", encoding='utf-8-sig', index=False)


# wdf = gdf[(gdf['MDV']!=1)&(gdf['DV']!='0')].groupby('CMT', as_index=False).agg({'DV':np.median})
# for inx, c in enumerate(wdf['DV']/10):
#     print(f"(0, {round(c,2)})")
#     print(f"(0, {round(c,2)})")