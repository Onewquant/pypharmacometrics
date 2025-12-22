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
    # upper = corr.where(np.tril(np.ones(corr.shape), k=-1).astype(bool))
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

def fit_ev_logistic(df,
                    out_csv: str = None,
                    corr_threshold: float = 0.9,
                    vif_threshold: float = 10.0,
                    enable_corr_filter: bool = True,
                    enable_vif_filter: bool = True,
                    verbose: bool = True,
                    drop_intercept_in_output: bool = True,
                    filename_add: str = '',):


    # csv_path=csv_path
    # out_csv=out_csv
    # corr_threshold=0.7     # |r| 임계값 (보통 0.8~0.95)
    # vif_threshold=10.0     # VIF 임계값 (보통 5 또는 10)
    # enable_corr_filter=True
    # enable_vif_filter=True
    # verbose=True

    # 1) 데이터 로드

    if "EV" not in df.columns:
        raise ValueError("데이터에 'EV' 컬럼이 필요합니다.")

    # 2) 타깃 y 준비(이진화)
    y = _to_binary_series(df["EV"])

    # 3) 설명변수 X = EV 제외 전부
    X = df.drop(columns=["EV"]).copy()

    # 4) 숫자/범주 분리
    num_cols = X.select_dtypes(include=[np.number]).columns.tolist()

    # 5) 숫자 컬럼 정리: inf -> NaN -> median 대치
    if num_cols:
        X[num_cols] = X[num_cols].replace([np.inf, -np.inf], np.nan)
        X[num_cols] = X[num_cols].fillna(X[num_cols].median())


    # cat_cols = [c for c in X.columns if c not in num_cols]
    # # 6) 범주 컬럼 정리: mode 대치 -> 원-핫 인코딩
    # if cat_cols:
    #     for c in cat_cols:
    #         mode_val = X[c].mode(dropna=True)
    #         fill_val = mode_val.iloc[0] if len(mode_val) else "MISSING"
    #         X[c] = X[c].fillna(fill_val).astype(str)
    #     X = pd.get_dummies(X, columns=cat_cols, drop_first=True)

    if X.shape[1] == 0:
        raise ValueError("EV 외 설명변수가 없습니다.")

    X_ori = X.copy()

    # 7) 0-분산(상수) 컬럼 제거
    zero_var_cols = [c for c in X.columns if X[c].nunique(dropna=False) <= 1]
    if zero_var_cols:
        if verbose:
            # print(f"[Zero-variance drop] {zero_var_cols}")
            pass
        X = X.drop(columns=zero_var_cols)

    # 8) [추가] 상관계수 기반 컬럼 제거
    corr_dropped = []
    if enable_corr_filter and X.shape[1] > 1:
        X, corr_dropped = drop_by_correlation(X, threshold=corr_threshold)
        if verbose and corr_dropped:
            # print(f"[Correlation filter |r|>{corr_threshold}] dropped: {corr_dropped}")
            pass

    # 9) 상수항 추가
    X = add_constant(X, has_constant='add')

    # 10) 반드시 포함/제거될 항 강제 지정
    must_inclusive_cols = ['DOSE24PERWT','SEX']
    must_drop_cols = []
    # must_drop_cols = ['TF_CRYO', 'TF_FFP', 'TF_PLT', 'TF_RBC']

    X_mincls_df = X_ori[must_inclusive_cols].copy()

    # 10) [추가] VIF 기반 반복 제거 (상수항 제외)
    vif_dropped = []
    if enable_vif_filter and X.shape[1] > 2:  # const + >=2 features
        X, vif_dropped = drop_by_vif(X, threshold=vif_threshold, max_iter=100)
        if verbose and vif_dropped:
            # print(f"[VIF filter > {vif_threshold}] dropped: {vif_dropped}")
            pass

    if len(must_inclusive_cols) > 0:
        for mincls_col in must_inclusive_cols:
            if mincls_col in X.columns:
                continue
            else:
                X[mincls_col]=X_mincls_df[mincls_col].copy()
    if len(must_drop_cols)>0:
        for mdrop_col in must_drop_cols:
            if mdrop_col in X.columns:
                X = X.drop([mdrop_col],axis=1)
            else:
                continue

    # 11) 형식 정리
    X = X.astype(float)
    y = y.astype(int)

    # 12) 적합 (Logit → 문제시 GLM Binomial 폴백)
    try:
        # print('test1')
        model = sm.Logit(y, X)
        res = model.fit(disp=False, maxiter=200)
        fitted_with = "Logit"
    except (PerfectSeparationError, np.linalg.LinAlgError):
        # print('test2')
        model = sm.GLM(y, X, family=sm.families.Binomial())
        res = model.fit()
        fitted_with = "GLM Binomial (fallback)"

    # 추가처리
    # if 'CUM_DOSE' in X.columns:
    #     X['CUM_DOSE'] = X['CUM_DOSE']/600

    # model = sm.Logit(y, X)
    # res = model.fit(disp=False, maxiter=200)
    # fitted_with = "Logit"

    # model = sm.GLM(y, X, family=sm.families.Binomial())
    # res = model.fit()
    # fitted_with = "GLM Binomial (fallback)"


    # except (PerfectSeparationError, np.linalg.LinAlgError):
    #     # print('test2')
    #     model = sm.Logit(y, X)
    #     res = model.fit(disp=False, maxiter=200)
    #     fitted_with = "Logit"

    # except Exception as e:
    #     raise

    # 13) OR 테이블 구성
    params = res.params
    conf = res.conf_int()
    conf.columns = ["2.5%", "97.5%"]
    or_table = pd.DataFrame({
        "feature": params.index,
        "beta": params.values,
        "aOR": np.exp(params.values),
        "aOR CI2.5%": np.exp(conf["2.5%"].values),
        "aOR CI97.5%": np.exp(conf["97.5%"].values),
        "pvalue": res.pvalues.values,
    })

    # 14) 다중비교 보정: Benjamini–Hochberg FDR
    _, p_adj, _, _ = multipletests(or_table["pvalue"].values, method="fdr_bh")
    or_table["pvalue_adj"] = p_adj

    # 14) 다중비교 보정: Bonferroni
    # _, p_adj, _, _ = multipletests(or_table["pvalue"].values, method="bonferroni")
    # or_table["pvalue_adj"] = p_adj


    # 15) Univariable Logistic Regression 추가

    uni_results = []

    for col in X.columns:
        if col == "const":
            continue
        try:
            Xi = sm.add_constant(X[[col]], has_constant="add")
            model_uni = sm.Logit(y, Xi)
            res_uni = model_uni.fit(disp=False, maxiter=200)
            params_uni = res_uni.params
            conf_uni = res_uni.conf_int(alpha=0.05)
            conf_uni.columns = ["2.5%", "97.5%"]

            uni_results.append({
                "feature": col,
                "uni_beta": params_uni[col],
                "uni_OR": np.exp(params_uni[col]),
                "uni_CI2.5%": np.exp(conf_uni.loc[col, "2.5%"]),
                "uni_CI97.5%": np.exp(conf_uni.loc[col, "97.5%"]),
                "uni_pvalue": res_uni.pvalues[col],
            })
        except (PerfectSeparationError, np.linalg.LinAlgError):
            uni_results.append({
                "feature": col,
                "uni_beta": np.nan,
                "uni_OR": np.nan,
                "uni_CI2.5%": np.nan,
                "uni_CI97.5%": np.nan,
                "uni_pvalue": np.nan,
            })

    uni_table = pd.DataFrame(uni_results)

    # -------------------------------
    # 15) 다변량 + 단변량 병합
    # -------------------------------
    or_table = or_table.merge(uni_table, on="feature", how="left")


    if drop_intercept_in_output and "const" in or_table["feature"].values:
        or_table = or_table[or_table["feature"] != "const"].copy()

    # or_table = or_table.sort_values("pvalue").reset_index(drop=True)
    or_table = or_table.sort_values("pvalue_adj").reset_index(drop=True)

    # 14) 결과 저장
    if out_csv:
        rp_str = '.csv' if filename_add=='' else f"({filename_add}).csv"
        or_table.to_csv(out_csv.replace('.csv',rp_str), index=False)

    # 15) 진단 정보 리턴
    final_vif = compute_vif(add_constant(X.drop(columns=["const"]), has_constant='add')) \
                if "const" in X.columns else compute_vif(X)

    info = {
        "fitted_with": fitted_with,
        "dropped_by_corr": corr_dropped,
        "dropped_by_vif": vif_dropped,
        "final_features": [c for c in X.columns if c != "const"],
        "final_vif": final_vif.sort_values("VIF", ascending=False).reset_index(drop=True).set_index('feature')['VIF'].to_dict()
    }
    return res, or_table, info


# 사용 예시: 경로만 바꿔서 실행하세요.
output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'
multivar_totres_df = list()
multivar_res_df = list()
covar_info_df = list()
orthreshold_dict = {True: 0.5, False: 0.33}
for endpoint in ['PLT', 'Hb', 'WBC', 'ANC', 'Lactate']:
# for endpoint in ['ANC', ]:
    # endpoint = 'PLT'
    # endpoint = 'Hb'
    # endpoint = 'WBC'
    # endpoint = 'ANC'
    # endpoint = 'Lactate'

    csv_path = f"{output_dir}/b1da/mvlreg/b1da_lnz_mvlreg_{endpoint}_df.csv"     # <-- 본인 파일 경로

    if not os.path.exists(f'{output_dir}/b1da/mvlreg_output(forced)'):
        os.mkdir(f'{output_dir}/b1da/mvlreg_output(forced)')

    out_csv = f"{output_dir}/b1da/mvlreg_output(forced)/b1da_lnz_mvlreg_{endpoint}_res.csv"

    # raise ValueError
    df_ori = pd.read_csv(csv_path)
    active_cols = [c for c in df_ori.columns if ('(TOTAL' not in c)]
    df_ori = df_ori[active_cols].copy()

    # df_ori.columns

    # raise ValueError

    # for age_subgroup in ['Adult','Elderly','Total_Adult']:
    for age_subgroup in ['Total_Adult']:
        # df.columns
        df = df_ori[df_ori['DOSE_PERIOD'] >= 1].copy()
        # df = df.drop(['DOSE_PERIOD(TOTAL)'], axis=1)
        # df = df.drop(['DOSE_PERIOD(TOTAL)','CUM_DOSE', 'DOSE_PERIOD', 'DOSE24', 'DOSE24PERWT'], axis=1)
        # df = df.drop(['DOSE_PERIOD(TOTAL)', 'DOSE_PERIOD', 'DOSE24',], axis=1)

        # if endpoint == 'Lactate':
        #     df = df.drop(['Lactate', 'pH'], axis=1)
        # df.columns
    # for age_subgroup in ['Pediatric','Adult','Elderly','Total_Adult']:
    # for age_subgroup in ['Elderly','Total_Adult']:

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

        # 데이터 Subset 저장

        # if endpoint=='lactate':
        #     raise ValueError
        df = df.drop(['clozapine', 'cyclophosphamide', 'methimazole', 'propofol', 'propylthiouracil', 'TF_WBLD'],axis=1)

        res, or_table, info = fit_ev_logistic(
            df=df,
            out_csv=out_csv,
            corr_threshold=0.8,     # |r| 임계값 (보통 0.8~0.95)
            vif_threshold=10.0,     # VIF 임계값 (보통 5 또는 10)
            enable_corr_filter=True,
            enable_vif_filter=True,
            verbose=True,
            filename_add=age_subgroup,
        )

        # print(f"[Fitted with] / {endpoint} / {info['fitted_with']}")
        # print(f"\n[Model Summary ({endpoint})]")
        # print(res.summary())
        print(f"\n[Odds Ratio Table ({endpoint})]")
        # or_threshold = 0.5
        # or_threshold = 0.0
        or_table['aOR threshold'] = (or_table['aOR'] >= 1).map(orthreshold_dict)
        or_table['endpoint'] = endpoint
        or_table['subgroup'] = age_subgroup
        or_table['N'] = len(df)
        or_table['EVN'] = int(sum(df['EV']))
        or_table['EVPct'] = round(100*sum(df['EV'])/len(df),2)
        pv_cond = (or_table['pvalue_adj'] < 0.05)
        or_cond = ((np.abs(or_table['aOR']-1) >= or_table['aOR threshold']))
        # sig_res_frag = or_table[(or_table['pvalue'] < 0.05) & (np.abs(or_table['OR'] - 1) >= or_threshold)].copy()
        sig_res_frag = or_table[pv_cond&or_cond].copy()

        multivar_totres_df.append(or_table.copy())

        info['endpoint'] = endpoint
        info['subgroup'] = age_subgroup
        covar_info_df.append(info.copy())

        # sig_res_frag = sig_res_frag[['subgroup','endpoint','N','EVN','EVPct']+list(sig_res_frag.columns)[:-4]]
        # multivar_res_df.append(sig_res_frag.copy())
        if not os.path.exists(f'{output_dir}/b1da/mvlreg_output(forced)/datasubset'):
            os.mkdir(f'{output_dir}/b1da/mvlreg_output(forced)/datasubset')
        df.to_csv(f"{output_dir}/b1da/mvlreg_output(forced)/datasubset/b1da_lnz_mvlreg_datasubset({age_subgroup})({endpoint}).csv", index=False, encoding='utf-8-sig')

        print(sig_res_frag)
        # print(f"\nSaved OR table to: {out_csv}")

        print(f"\n[Correlation dropped ({endpoint})]")
        print(info["dropped_by_corr"])
        print(f"\n[VIF dropped ({endpoint})]")
        print(info["dropped_by_vif"])
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

multivar_totres_df.to_csv(f"{output_dir}/b1da/mvlreg_output(forced)/b1da_lnz_mvlreg_total_results.csv", index=False, encoding='utf-8-sig')

covar_info_df = pd.DataFrame(covar_info_df)
covar_info_df = covar_info_df[['subgroup','endpoint']+list(covar_info_df.columns)[0:-2]].copy()
covar_info_df_saving = pd.melt(covar_info_df,
                   id_vars=['subgroup','endpoint','fitted_with'],          # 녹이지 않고 그대로 둘 열
                   value_vars=['dropped_by_corr', 'dropped_by_vif', 'final_features', 'final_vif'],       # 녹일 열 (None이면 나머지 전부)
                   var_name='info_type',         # 녹인 열 이름
                   value_name="covar_list")    # 녹인 값의 열 이름)

# covar_info_df_saving.columns


covar_info_df_saving = covar_info_df_saving.sort_values(by=['subgroup','endpoint','info_type'],ascending=[False,True,True])
covar_info_df_saving = covar_info_df_saving[covar_info_df_saving['info_type']!='final_vif'].copy()
covar_info_df_saving = covar_info_df_saving.drop(['fitted_with'],axis=1)
covar_info_df_saving['covar_list'] = covar_info_df_saving['covar_list'].map(lambda x: str({c.replace('(ACTIVE)','') for c in x}).replace("', '",", ").replace("{'","").replace("'}",""))
covar_info_df_saving['info_type'] = covar_info_df_saving['info_type'].map({'dropped_by_corr':'(Dropped by correlation coefficient)','dropped_by_vif':'(Dropped by variance inflation factor)','final_features':'Covariates in the final model'})
total_covar_list = list(set((", ".join(list(covar_info_df_saving['covar_list']))).split(", ")))
total_covar_list.sort()
covar_primer_df = pd.DataFrame([{'subgroup':'Total Covariates',	'endpoint':'-',	'info_type':'Total Covariates',	'covar_list':str(total_covar_list).replace("', '",", ").replace("['","").replace("']","")}])
covar_info_df_saving = pd.concat([covar_primer_df,covar_info_df_saving])
# covar_info_df_saving.to_csv(f"{output_dir}/b1da/mvlreg_output/b1da_lnz_mvlreg_covar_info.csv", index=False, encoding='utf-8-sig')
covar_info_df_saving = covar_info_df_saving[(covar_info_df_saving['info_type'].isin(['Covariates in the final model','Total Covariates']))&(covar_info_df_saving['subgroup']!='Elderly')].copy()
covar_name_dict = {'TF_CRYO':'Cryoprecipitate transfusion', 'TF_FFP':'FFP transfusion', 'TF_PLT':'Platelet transfusion',
                   'TF_RBC':'RBC transfusion', 'TPRO':'Total protein', 'WBC':'WBC count', 'ALB':'Albumin', 'ALT':'ALT', 'ANC':'ANC', 'AST':'AST', 'BMI':'BMI', 'CRP':'CRP',
                   'CUM_DOSE':'Cumulative dose', 'DOSE24PERWT':'Weight-normalized daily dose', 'DOSE24': 'Daily dose',
                   'DOSE_INTERVAL':'Dosing interval', 'DOSE_PERIOD':'Dosing period', 'ELD':'Elderly', 'GLU':'Glucose',
                   'Hb':'Hemoglobin', 'Lactate':'Lactate', 'PLT':'Platelet count', 'SCR':'Serum creatinine',
                   'TBIL':'Total bilirubin',  'aspirinclopidogrel':'Aspirin/clopidogrel', 'clozapine':'Clozapine',
                   'cyclophosphamide':'Cyclophosphamide','eGFR':'eGFR', 'eGFR-CKD-EPI':'eGFR-CKD-EPI', 'heparin':'Heparin',
                   'isoniazid':'Isoniazid', 'metformin':'Metformin', 'methimazole':'Methimazole', 'SEX':'Sex', 'pH':'pH',
                   'piperacillintazobactam':'Piperacillin/tazobactam', 'propofol':'Propofol',
                   'tmpsmx':'Trimethoprim/sulfamethoxazole', 'valproate':'Valproate','WT':'Weight',
                   'AGE':'Age', 'HT':'Height', }

for inx, row in covar_info_df_saving.iterrows():
    for k, v in covar_name_dict.items():
        covar_info_df_saving.at[inx,'covar_list'] = covar_info_df_saving.at[inx,'covar_list'].replace(k,v)
    covar_info_df_saving.at[inx, 'covar_list'] = covar_info_df_saving.at[inx, 'covar_list'].replace('Clozapine, ', '').replace('Methimazole, ', '')


# raise ValueError
covar_info_df_saving.to_excel(f"{output_dir}/b1da/mvlreg_output(forced)/b1da_lnz_mvlreg_covar_info.xlsx", index=False, encoding='utf-8-sig')


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
# uni_res_df = uni_res_df.sort_values([''])
uni_res_df['pval_sig'] = (uni_res_df['pval (adj)'] < 0.05).map({True:' **',False:''})
uni_res_df['pval (adj)'] = uni_res_df['pval (adj)'].replace(0,'<0.001').map(str) + uni_res_df['pval_sig']
uni_res_df['pval'] = uni_res_df['pval'].replace(0,'<0.001').map(str)
uni_res_df = uni_res_df.drop(['pval_sig'],axis=1)


multi_res_df['pval_sig'] = (multi_res_df['pval (adj)'] < 0.05).map({True:' **',False:''})
multi_res_df['pval (adj)'] = multi_res_df['pval (adj)'].replace(0,'<0.001').map(str) + multi_res_df['pval_sig']
multi_res_df['pval'] = multi_res_df['pval'].replace(0,'<0.001').map(str)
multi_res_df = multi_res_df.drop(['pval_sig'],axis=1)
# multi_res_df['pval'] = multi_res_df['pval'].replace(0,'<0.001')
# multi_res_df['pval (adj)'] = multi_res_df['pval (adj)'].replace(0,'<0.001')

uni_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output(forced)/b1da_lnz_mvlreg_univar_res_table.csv", index=False, encoding='utf-8-sig')
multi_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output(forced)/b1da_lnz_mvlreg_multivar_res_table.csv", index=False, encoding='utf-8-sig')

# uni_res_df[uni_res_df['subgroup']=='Total_Adult'][['endpoint','feature','OR (95% CI)','pval']]

tot_res_df = multivar_totres_df.sort_values(['subgroup','endpoint','aOR'],ascending=[False,True,False])[['subgroup','endpoint','feature','EV_Count (%)','OR (95% CI)','pval','aOR (95% CI)','pval (adj)']].reset_index(drop=True)
# raise ValueError
# tot_res_df.columns
covar_name_dict = {'TF_CRYO':'Cryoprecipitate transfusion', 'TF_FFP':'FFP transfusion', 'TF_PLT':'Platelet transfusion',
                   'TF_RBC':'RBC transfusion', 'TPRO':'Total protein', 'WBC':'WBC count (/μL)', 'ALB':'Albumin (g/dL)',
                   'ALT':'ALT (IU/L)', 'ANC':'ANC (/μL)', 'AST':'AST (IU/L)', 'BMI':'BMI (kg/m2)', 'CRP':'CRP (mg/dL)',
                   'CUM_DOSE':'Cumulative dose (mg)', 'DOSE24PERWT':'Weight-normalized daily dose (mg/day·kg)',
                   'DOSE24': 'Daily dose (mg/day)', 'DOSE_INTERVAL':'Dosing interval (hour)',
                   'DOSE_PERIOD':'Dosing period (days)', 'ELD':'Elderly', 'GLU':'Glucose (mg/dL)',
                   'Hb':'Hemoglobin (g/dL)', 'Lactate':'Lactate (mmol/L)', 'PLT':'Platelet count (/μL)',
                   'SCR':'Serum creatinine (mg/dL)', 'TBIL':'Total bilirubin (mg/dL)',
                   'aspirinclopidogrel':'Aspirin/clopidogrel', 'clozapine':'Clozapine',
                   'cyclophosphamide':'Cyclophosphamide', 'eGFR':'eGFR (mL/min/1.73 m²)',
                   'eGFR-CKD-EPI':'eGFR-CKD-EPI (mL/min/1.73 m²)', 'heparin':'Heparin',
                   'isoniazid':'Isoniazid', 'metformin':'Metformin', 'methimazole':'Methimazole', 'SEX':'Sex', 'pH':'pH',
                   'piperacillintazobactam':'Piperacillin/tazobactam', 'propofol':'Propofol',
                   'tmpsmx':'Trimethoprim/sulfamethoxazole', 'valproate':'Valproate','WT':'Weight (kg)',
                   'AGE':'Age (years)', 'HT':'Height (cm)', }

for k, v in covar_name_dict.items():
    tot_res_df['feature'] = tot_res_df['feature'].replace(k, v)
tot_res_df = tot_res_df[tot_res_df['subgroup']=='Total_Adult'].copy()
# raise ValueError
# tot_res_df.columns
tot_res_df = tot_res_df[~tot_res_df['feature'].isin(['Clozapine','Methimazole'])].copy()
tot_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output(forced)/b1da_lnz_mvlreg_total_res_table.csv", index=False, encoding='utf-8-sig')

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