# EV ~ (나머지 covariates) 다변량 로지스틱 회귀 + 다중공선성 제거(correlation, VIF)
# - CSV 로딩
# - EV 이진화
# - 숫자: inf→NaN→median 대치
# - 범주: mode 대치 + 원-핫 인코딩(drop_first)
# - 0-분산 컬럼 제거, 상수항 추가
# - [추가] 상관계수 기반 제거(|r|>corr_threshold)
# - [추가] VIF 기반 반복 제거(VIF>vif_threshold)
# - Logit 적합(문제 시 GLM Binomial 폴백)
# - OR, 95% CI, p-value 테이블 저장/출력

import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.tools import add_constant
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.stats.multitest import multipletests
import os


def _to_binary_series(s: pd.Series) -> pd.Series:
    """EV를 이진(0/1)으로 변환: bool, {0,1}, 임의의 2값, 문자열 2수준 지원."""
    if s.dtype == bool:
        return s.astype(int)
    if np.issubdtype(s.dropna().dtype, np.number):
        uniq = pd.Series(sorted(s.dropna().unique()))
        if set(uniq) <= {0, 1}:
            return s.astype(int)
        if len(uniq) == 2:
            mapping = {uniq.iloc[0]: 0, uniq.iloc[1]: 1}
            return s.map(mapping).astype(int)
        raise ValueError("EV는 이진이어야 합니다(0/1 또는 2개 고유값).")
    uniq = s.dropna().unique()
    if len(uniq) == 2:
        mapping = {uniq[0]: 0, uniq[1]: 1}
        return s.map(mapping).astype(int)
    raise ValueError("EV는 이진이어야 합니다(문자열 2수준 가능).")

def drop_by_correlation(X: pd.DataFrame, threshold: float = 0.9):
    """
    상관계수 절대값이 threshold를 넘는 쌍들에 대해 한 변수를 제거.
    규칙: 두 변수의 '다른 변수들과의 평균|r|'이 더 큰 쪽을 제거(연결성이 큰 쪽을 드롭).
    """
    if X.shape[1] < 2:
        return X, []

    corr = X.corr().abs()
    # 상삼각만 사용하여 중복 제거
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
    to_drop = set()

    # 반복적으로 탐색하면서 가장 '연결성 높은' 변수 제거
    while True:
        # 현재 남아있는 쌍 중 threshold 초과가 존재하는지 확인
        mask = (upper > threshold)
        if not mask.any().any():
            break

        # 초과 쌍 리스트
        pairs = np.argwhere(mask.values)
        # 각 변수의 평균 연결 강도(본인 제외)
        mean_conn = (corr.sum() - 1) / (corr.shape[0] - 1)

        # 초과쌍마다 하나 선택해서 제거 후보에 추가
        local_drop = []
        for i, j in pairs:
            vi = upper.index[i]
            vj = upper.columns[j]
            # 이미 제거 예정이면 스킵
            if vi in to_drop or vj in to_drop:
                continue
            # 평균 연결 강도가 큰 것을 제거
            if mean_conn[vi] >= mean_conn[vj]:
                local_drop.append(vi)
            else:
                local_drop.append(vj)

        # 이번 라운드에 제거할 게 없다면 종료
        if not local_drop:
            break

        # 반영
        for v in local_drop:
            if v not in to_drop:
                to_drop.add(v)
                # 행/열 제거하여 이후 탐색 업데이트
                corr = corr.drop(index=v, columns=v)
                upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
                if corr.shape[0] < 2:
                    break

        if corr.shape[0] < 2:
            break

    kept_cols = [c for c in X.columns if c not in to_drop]
    return X[kept_cols].copy(), sorted(list(to_drop))

def compute_vif(X: pd.DataFrame) -> pd.DataFrame:
    """상수항(const) 제외하고 VIF 계산 후 DataFrame 반환."""
    cols = [c for c in X.columns if c != "const"]
    if len(cols) == 0:
        return pd.DataFrame(columns=["feature", "VIF"])
    X_ = X[cols].astype(float)
    vif_vals = []
    for i in range(X_.shape[1]):
        vif_vals.append(variance_inflation_factor(X_.values, i))
    return pd.DataFrame({"feature": cols, "VIF": vif_vals}).sort_values("VIF", ascending=False)

def drop_by_vif(X: pd.DataFrame, threshold: float = 10.0, max_iter: int = 50):
    """
    반복적으로 VIF가 threshold를 넘는 컬럼 중 최댓값을 가진 컬럼을 하나씩 제거.
    상수항(const)은 유지.
    """
    dropped = []
    X_work = X.copy()
    print(f"[Drop by VIF value]\n")
    for _ in range(max_iter):
        vif_df = compute_vif(X_work)
        if vif_df.empty:
            break
        max_vif = vif_df["VIF"].max()
        if not np.isfinite(max_vif) or max_vif <= threshold:
            break
        # 최댓값을 가진 feature 제거
        drop_feat = vif_df.iloc[0]["feature"]
        dropped.append(drop_feat)
        X_work = X_work.drop(columns=[drop_feat])
        # print(f"{drop_feat} / VIF: {max_vif} / Dropped")
        if X_work.shape[1] <= 1:  # const만 남는 상황
            break
    return X_work, dropped



# 사용 예시: 경로만 바꿔서 실행하세요.
output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'
multivar_totres_df = list()
multivar_res_df = list()
orthreshold_dict = {True: 0.5, False: 0.33}
for endpoint in ['PLT', 'Hb', 'WBC', 'ANC', 'Lactate']:

    csv_path = f"{output_dir}/b1da/mvlreg/b1da_lnz_mvlreg_{endpoint}_df.csv"     # <-- 본인 파일 경로

    if not os.path.exists(f'{output_dir}/b1da/mvlreg_output_glm'):
        os.mkdir(f'{output_dir}/b1da/mvlreg_output_glm')

    out_csv = f"{output_dir}/b1da/mvlreg_output_glm/b1da_lnz_mvlreg_{endpoint}_res.csv"

    # raise ValueError
    df_ori = pd.read_csv(csv_path)
    active_cols = [c for c in df_ori.columns if ('(TOTAL' not in c)]
    df_ori = df_ori[active_cols].copy()

    age_subgroup = 'Total_Adult'
    df = df_ori[df_ori['DOSE_PERIOD'] >= 1].copy()
    if age_subgroup=='Adult':
        df = df[(df['AGE'] >= 19)&(df['AGE'] < 65)].copy()
    elif age_subgroup=='Total_Adult':
        df = df[(df['AGE'] >= 19)].copy()
    elif age_subgroup=='Elderly':
        df = df[(df['AGE'] >= 65)].copy()
    elif age_subgroup=='Pediatric':
        df = df[df['AGE'] < 19].copy()
    else:
        raise ValueError




    # res, or_table, info = fit_ev_logistic(
    corr_threshold=0.8     # |r| 임계값 (보통 0.8~0.95)
    vif_threshold=10.0     # VIF 임계값 (보통 5 또는 10)
    enable_corr_filter=True
    enable_vif_filter=True
    verbose=True
    filename_add=age_subgroup

    # 1) 타깃 y 준비(이진화)
    y = _to_binary_series(df["EV"])

    # 2) 설명변수 X = EV 제외 전부
    X = df.drop(columns=["EV"]).copy()

    # 3) 숫자/범주 분리
    num_cols = X.select_dtypes(include=[np.number]).columns.tolist()

    # 4) 숫자 컬럼 정리: inf -> NaN -> median 대치
    if num_cols:
        X[num_cols] = X[num_cols].replace([np.inf, -np.inf], np.nan)
        X[num_cols] = X[num_cols].fillna(X[num_cols].median())

    # 5) 0-분산(상수) 컬럼 제거
    zero_var_cols = [c for c in X.columns if X[c].nunique(dropna=False) <= 1]
    if zero_var_cols:
        if verbose:
            pass
        X = X.drop(columns=zero_var_cols)

    # 6) [추가] 상관계수 기반 컬럼 제거
    corr_dropped = []
    if enable_corr_filter and X.shape[1] > 1:
        X, corr_dropped = drop_by_correlation(X, threshold=corr_threshold)
        if verbose and corr_dropped:
            pass

    # 7) 상수항 추가
    X = add_constant(X, has_constant='add')

    # 8) [추가] VIF 기반 반복 제거 (상수항 제외)
    vif_dropped = []
    if enable_vif_filter and X.shape[1] > 2:  # const + >=2 features
        X, vif_dropped = drop_by_vif(X, threshold=vif_threshold, max_iter=100)
        if verbose and vif_dropped:
            # print(f"[VIF filter > {vif_threshold}] dropped: {vif_dropped}")
            pass

    # 11) 형식 정리
    X = X.astype(float)
    y = y.astype(int)

    X = X.rename(columns={'const': 'PID'})
    X = X.reset_index(drop=True)
    X['PID'] = X.index
    X = add_constant(X, has_constant='add')
    X = X.rename(columns={'const':'y'})
    X['y'] = y
    X = add_constant(X, has_constant='add')
    X = X.rename(columns={'const': 'ENDPOINT'})
    X['ENDPOINT'] = endpoint
    if not os.path.exists(f'{output_dir}/b1da/mvlreg_output_glm/r'):
        os.mkdir(f'{output_dir}/b1da/mvlreg_output_glm/r')
    X.to_csv(f"{output_dir}/b1da/mvlreg_output_glm/r/b1da_lnz_mvlreg_datasubset({endpoint}).csv", index=False, encoding='utf-8-sig')

    raise ValueError


    if not os.path.exists(f'{output_dir}/b1da/mvlreg_output/datasubset'):
        os.mkdir(f'{output_dir}/b1da/mvlreg_output/datasubset')
    df.to_csv(f"{output_dir}/b1da/mvlreg_output/datasubset/b1da_lnz_mvlreg_datasubset({age_subgroup})({endpoint}).csv", index=False, encoding='utf-8-sig')

    print(sig_res_frag)
    # print(f"\nSaved OR table to: {out_csv}")

    print(f"\n[Correlation dropped ({endpoint})]")
    print(info["corr_dropped"])
    print(f"\n[VIF dropped ({endpoint})]")
    print(info["vif_dropped"])
    print(f"\n[Final features ({endpoint})]")
    print(info["final_features"])
    # print(f"\n[Final VIF ({endpoint})]")
    # print(info["final_vif"])

# df[df['DOSE_PERIOD']!=0]
# multivar_res_df.columns
# df[(df['EV']==1)]['DOSE_PERIOD'].mean()
# df[(df['EV']==0)]['DOSE_PERIOD'].mean()
multivar_totres_df = pd.concat(multivar_totres_df)
# multivar_totres_df.to_csv(f"{output_dir}/b1da/mvlreg_output/b1da_lnz_mvlreg_total_results({str(or_threshold)[-1]}).csv", index=False, encoding='utf-8-sig')
multivar_totres_df.to_csv(f"{output_dir}/b1da/mvlreg_output/b1da_lnz_mvlreg_total_results.csv", index=False, encoding='utf-8-sig')

sg_cond = (multivar_totres_df['subgroup'] == 'Total_Adult')
uor_cond = (np.abs(multivar_totres_df['uni_OR']-1) >= ((multivar_totres_df['uni_OR']>=1).map(orthreshold_dict)))
upv_cond = (multivar_totres_df['uni_pvalue'] < 0.05)
aor_cond = (np.abs(multivar_totres_df['aOR']-1) >= ((multivar_totres_df['aOR']>=1).map(orthreshold_dict)))
apv_cond = (multivar_totres_df['pvalue_adj'] < 0.05)

# multivar_totres_df.columns

sig_dig = 3

multivar_totres_df['EV_Count (%)'] = multivar_totres_df.apply(lambda x:f"{round(x['EVN'],sig_dig)} ({round(x['EVPct'],sig_dig)})",axis=1)

multivar_totres_df['OR (95% CI)'] = multivar_totres_df.apply(lambda x:f"{round(x['uni_OR'],sig_dig)} ({round(x['uni_CI2.5%'],sig_dig)}-{round(x['uni_CI97.5%'],sig_dig)})",axis=1)
multivar_totres_df['pval'] = multivar_totres_df['uni_pvalue'].map(lambda x:round(x,sig_dig))

multivar_totres_df['aOR (95% CI)'] = multivar_totres_df.apply(lambda x:f"{round(x['aOR'],sig_dig)} ({round(x['aOR CI2.5%'],sig_dig)}-{round(x['aOR CI97.5%'],sig_dig)})",axis=1)
multivar_totres_df['pval (adj)'] = multivar_totres_df['pvalue_adj'].map(lambda x:round(x,sig_dig))

# multivar_totres_df[(multivar_totres_df['subgroup']=='Total_Adult')&(multivar_totres_df['endpoint']=='ANC')&(multivar_totres_df['feature']=='SEX')][['OR (95% CI)','pval','aOR (95% CI)','pval (adj)']]
# multivar_totres_df[upv_cond&uor_cond]
# multivar_totres_df['pval']
uni_res_df = multivar_totres_df[upv_cond&uor_cond].sort_values(['subgroup','endpoint','aOR'],ascending=[False,True,False])[['subgroup','endpoint','feature','EV_Count (%)','OR (95% CI)','pval','aOR (95% CI)','pval (adj)']].reset_index(drop=True)
multi_res_df = multivar_totres_df[apv_cond&aor_cond].sort_values(['subgroup','endpoint','aOR'],ascending=[False,True,False])[['subgroup','endpoint','feature','EV_Count (%)','OR (95% CI)','pval','aOR (95% CI)','pval (adj)']].reset_index(drop=True)
uni_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output/b1da_lnz_mvlreg_univar_res_table.csv", index=False, encoding='utf-8-sig')
multi_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output/b1da_lnz_mvlreg_multivar_res_table.csv", index=False, encoding='utf-8-sig')

tot_res_df = multivar_totres_df.sort_values(['subgroup','endpoint','aOR'],ascending=[False,True,False])[['subgroup','endpoint','feature','EV_Count (%)','OR (95% CI)','pval','aOR (95% CI)','pval (adj)']].reset_index(drop=True)
tot_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output/b1da_lnz_mvlreg_total_res_table.csv", index=False, encoding='utf-8-sig')

# multivar_totres_df[pv_cond&sg_cond].sort_values(['subgroup','endpoint','OR'],ascending=[False,True,False])[['subgroup','endpoint','feature','n (%)','OR (95% CI)','pval','aOR (95% CI)','pval (adj)']]
# multivar_totres_df.sort_values(['subgroup','endpoint','OR'],ascending=[False,True,False])[['endpoint','feature','n (%)','OR (95% CI)','pval','aOR (95% CI)','pval (adj)']]

# or_table[or_table['feature']=='CUM_DOSE']
# multivar_res_df = pd.concat(multivar_res_df)
# multivar_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output/b1da_lnz_mvlreg_significant_results({str(or_threshold)[-1]}).csv", index=False, encoding='utf-8-sig')
# print(multivar_res_df[['subgroup','endpoint','feature','N','EVN','EVPct','OR','pvalue','pvalue_adj']].reset_index(drop=True))
######################################


            
# sns.relplot(data=df_ori, x='SCR', y='CUM_DOSE')
# sns.relplot(data=df_ori, x='SCR', y='DOSE_PERIOD(TOTAL)')
# sns.relplot(data=df_ori, x='DBIL', y='CUM_DOSE')
# sns.relplot(data=df_ori, x='DBIL', y='DOSE_PERIOD(TOTAL)')
# sns.relplot(data=df_ori, x='TBIL', y='DOSE_PERIOD(TOTAL)')
# sns.relplot(data=df_ori, x='AST', y='DOSE_PERIOD(TOTAL)')
# sns.relplot(data=df_ori, x='eGFR', y='DOSE_PERIOD(TOTAL)')
# # sns.relplot(data=df_ori, x='TBIL', y='DOSE_PERIOD(TOTAL)')
# sns.relplot(data=df_ori, x='AGE', y='DOSE_PERIOD(TOTAL)')
# sns.distplot(df_ori['DOSE_PERIOD(TOTAL)'])
# df_ori['DOSE_PERIOD(TOTAL)']