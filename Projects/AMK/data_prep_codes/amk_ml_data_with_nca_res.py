from pynca.tools import *
import numpy as np
import pandas as pd
from tools import *


prj_name = 'AMK'
prj_dir = f'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/{prj_name}'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
nonmem_dir = f'C:/Users/ilma0/NONMEMProjects/{prj_name}'

########## NCA ###########
sim_df = pd.read_csv(f"{resource_dir}/sim053",encoding='utf-8-sig', skiprows=1, sep=r"\s+", engine='python')
conc_df = sim_df[sim_df['MDV']==0][['ID','TIME','DV','MDV','AMT']].copy()

# sim_df.drop_duplicates(['ID'], keep='last')
# conc_df.drop_duplicates(['ID'], keep='last')
total_auc_df = add_iAUC(conc_df, time_col="TIME", conc_col="DV", id_col="ID")
total_auc_df['ID'] = total_auc_df['ID'].astype(int)

# total_auc_df

auc_all_df = total_auc_df.drop_duplicates(['ID'], keep='last')
auc_all_df['TIME24BF'] = auc_all_df['TIME']-24
add_on_df = auc_all_df[['ID','TIME24BF','cumAUC']].rename(columns={'cumAUC':'cumAUC_ALL'})

auc_24bf_df = total_auc_df.merge(add_on_df, on='ID', how='left')
auc_24bf_df = auc_24bf_df[auc_24bf_df['TIME'] >= auc_24bf_df['TIME24BF']].copy()
auc_24bf_df['cumAUC_LAST24'] = auc_24bf_df['cumAUC_ALL'] - auc_24bf_df['cumAUC']
auc_24bf_df_1st_rows = auc_24bf_df.drop_duplicates(['ID'], keep='first')

AUC_F = auc_24bf_df_1st_rows[['ID','cumAUC_ALL','cumAUC_LAST24']].copy()

########## 24H 중에 maximum (Cmax,Ctrough) 구하기 ##########


def find_local_extrema(group: pd.DataFrame) -> pd.DataFrame:
    # TIME 기준 정렬
    g = group.sort_values("TIME").reset_index(drop=False)  # 기존 index 보존
    t = g["TIME"].to_numpy(dtype=float)
    dv = g["DV"].to_numpy(dtype=float)
    n = len(g)

    # 기본값은 'none'
    extreme_type = np.array(["none"] * n, dtype=object)

    if n >= 3:
        for i in range(1, n - 1):
            dv_prev = dv[i - 1]
            dv_curr = dv[i]
            dv_next = dv[i + 1]

            # ---- local maximum: 기존 정의 그대로 ----
            if (dv_curr > dv_prev) and (dv_curr > dv_next):
                extreme_type[i] = "max"
                continue

            # ---- local minimum: 확장 조건 적용 ----
            if (dv_curr < dv_prev) and (dv_curr < dv_next):
                # 기준 시간
                ti = t[i]

                # 앞쪽 2.167시간 이내의 시점들 (현재 지점 제외)
                prev_mask = (t < ti) & (t >= ti - 2.167)
                prev_idx = np.where(prev_mask)[0]

                # 뒤쪽 2.167시간 이내의 시점들 (현재 지점 제외)
                next_mask = (t > ti) & (t <= ti + 2.167)
                next_idx = np.where(next_mask)[0]

                # 각각 최소 2개 이상 존재하고, 해당 구간의 DV가 모두 현재 값보다 커야 함
                if (
                    prev_idx.size >= 2
                    and next_idx.size >= 2
                    and (dv[prev_idx] > dv_curr).all()
                    and (dv[next_idx] > dv_curr).all()
                ):
                    extreme_type[i] = "min"

    g["extreme_type"] = extreme_type

    # local minimum 또는 maximum만 남기기
    g = g[g["extreme_type"] != "none"]

    return g

# ID별로 local extrema 추출
local_total_auc_df_24H = (
    auc_24bf_df.groupby("ID", group_keys=False)[["ID", "TIME", "DV"]].apply(find_local_extrema)
)

# 확인
print(local_total_auc_df_24H.head())

def pick_final_extrema(group: pd.DataFrame) -> pd.DataFrame:
    rows = []

    # local minimum 들
    mins = group[group["extreme_type"] == "min"]
    # local maximum 들
    maxs = group[group["extreme_type"] == "max"]

    # local minimum 중 DV < 15 조건 추가 + 그 중 최대값 1개 선택
    mins_filtered = mins[mins["DV"] < 15]
    if not mins_filtered.empty:
        rows.append(mins_filtered.sort_values("DV", ascending=False).iloc[0])

    # local maximum 중 DV가 가장 큰 것 1개
    if not maxs.empty:
        rows.append(maxs.sort_values("DV", ascending=False).iloc[0])

    if rows:
        return pd.DataFrame(rows)
    else:
        # 해당 ID에 local extrema가 전혀 없으면 빈 DataFrame 반환
        return pd.DataFrame(columns=group.columns)

final_local_24H = (
    local_total_auc_df_24H
    .groupby("ID", group_keys=False)[["ID", "TIME", "DV", "extreme_type"]]
    .apply(pick_final_extrema)
    .reset_index(drop=True)
)
final_local_24H = final_local_24H[['ID','DV','extreme_type']].copy()
final_local_24H['extreme_type'] = final_local_24H['extreme_type'].map({'min':'TROUGH_LOCMAX','max':'PEAK_LOCMAX'})
final_local_24H = final_local_24H[['ID','DV','extreme_type']].pivot_table(values=['DV'], index=['ID'], columns=['extreme_type'])
# final_local_24H.columns.names[]
final_local_24H.columns = [c[1] for c in final_local_24H.columns]
# final_local_24H.index.name = 'ID'
final_local_24H = final_local_24H.reset_index(drop=False)
pk_df = AUC_F.merge(final_local_24H, on='ID', how='left')

ml_res_df = pd.read_csv(f"{output_dir}/final_mlres_data.csv")
ml_res_df = ml_res_df[~ml_res_df['UID'].isin({25524226, 24961411, 10617861, 15499525, 22006666, 25389067, 19551122, 18347926, 14009765, 24311845, 30514726, 10190892, 10451629, 11895215, 16574645, 34728899, 11584452, 11845188, 24913861, 28650698, 23809356, 10006221, 13115597, 21210190, 26599247, 13158741, 14783830, 18146774, 26367192, 28080857, 18991455, 10474848, 24925408, 26668770, 35441638, 21249383, 10963177, 25351536, 11675891, 24785011, 32592760, 26809209, 23084924, 26948351})]
# print(final_local_24H.head())
ml_res_df1 = ml_res_df.loc[:,:'Empirical'].copy()
ml_res_df2 = pd.concat([ml_res_df[['ID']],ml_res_df.loc[:,'AKI_OCCURRENCE':].copy()], axis=1)
ml_res_df = ml_res_df1.merge(pk_df, on='ID',how='left').merge(ml_res_df2, on='ID',how='left')
ml_res_df = ml_res_df.rename(columns={'γ-GT':'GGT'})
ml_res_df.to_csv(f"{output_dir}/final_mlres_data(with_pk).csv", index=False, encoding='utf-8')

# local_total_auc_df_24H.to_csv(r"C:/PYCHARM/AMK/results/local_total_auc_3_24H.csv", index=False)
# final_local_24H.to_csv(r"C:/PYCHARM/AMK/results/final_local_24H.csv", index=False)
