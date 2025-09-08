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
        print(f"{drop_feat} / VIF: {max_vif} / Dropped")
        if X_work.shape[1] <= 1:  # const만 남는 상황
            break
    return X_work, dropped

def fit_ev_logistic(csv_path: str,
                    out_csv: str = None,
                    corr_threshold: float = 0.9,
                    vif_threshold: float = 10.0,
                    enable_corr_filter: bool = True,
                    enable_vif_filter: bool = True,
                    verbose: bool = True,
                    drop_intercept_in_output: bool = True):


    # csv_path=csv_path
    # out_csv=out_csv
    # corr_threshold=0.7     # |r| 임계값 (보통 0.8~0.95)
    # vif_threshold=10.0     # VIF 임계값 (보통 5 또는 10)
    # enable_corr_filter=True
    # enable_vif_filter=True
    # verbose=True

    # 1) 데이터 로드
    df = pd.read_csv(csv_path)
    if "EV" not in df.columns:
        raise ValueError("데이터에 'EV' 컬럼이 필요합니다.")

    # 2) 타깃 y 준비(이진화)
    y = _to_binary_series(df["EV"])

    # 3) 설명변수 X = EV 제외 전부
    X = df.drop(columns=["EV"]).copy()

    # 4) 숫자/범주 분리
    num_cols = X.select_dtypes(include=[np.number]).columns.tolist()
    cat_cols = [c for c in X.columns if c not in num_cols]

    # 5) 숫자 컬럼 정리: inf -> NaN -> median 대치
    if num_cols:
        X[num_cols] = X[num_cols].replace([np.inf, -np.inf], np.nan)
        X[num_cols] = X[num_cols].fillna(X[num_cols].median())

    # 6) 범주 컬럼 정리: mode 대치 -> 원-핫 인코딩
    if cat_cols:
        for c in cat_cols:
            mode_val = X[c].mode(dropna=True)
            fill_val = mode_val.iloc[0] if len(mode_val) else "MISSING"
            X[c] = X[c].fillna(fill_val).astype(str)
        X = pd.get_dummies(X, columns=cat_cols, drop_first=True)

    if X.shape[1] == 0:
        raise ValueError("EV 외 설명변수가 없습니다.")

    # 7) 0-분산(상수) 컬럼 제거
    zero_var_cols = [c for c in X.columns if X[c].nunique(dropna=False) <= 1]
    if zero_var_cols:
        if verbose:
            print(f"[Zero-variance drop] {zero_var_cols}")
        X = X.drop(columns=zero_var_cols)

    # 8) [추가] 상관계수 기반 컬럼 제거
    corr_dropped = []
    if enable_corr_filter and X.shape[1] > 1:
        X, corr_dropped = drop_by_correlation(X, threshold=corr_threshold)
        if verbose and corr_dropped:
            print(f"[Correlation filter |r|>{corr_threshold}] dropped: {corr_dropped}")

    # 9) 상수항 추가
    X = add_constant(X, has_constant='add')

    # 10) [추가] VIF 기반 반복 제거 (상수항 제외)
    vif_dropped = []
    if enable_vif_filter and X.shape[1] > 2:  # const + >=2 features
        X, vif_dropped = drop_by_vif(X, threshold=vif_threshold, max_iter=100)
        if verbose and vif_dropped:
            print(f"[VIF filter > {vif_threshold}] dropped: {vif_dropped}")

    # 11) 형식 정리
    X = X.astype(float)
    y = y.astype(int)

    # 12) 적합 (Logit → 문제시 GLM Binomial 폴백)
    try:
        model = sm.Logit(y, X)
        res = model.fit(disp=False, maxiter=200)
        fitted_with = "Logit"
    except (PerfectSeparationError, np.linalg.LinAlgError):
        model = sm.GLM(y, X, family=sm.families.Binomial())
        res = model.fit()
        fitted_with = "GLM Binomial (fallback)"
    except Exception as e:
        raise

    # 13) OR 테이블 구성
    params = res.params
    conf = res.conf_int()
    conf.columns = ["2.5%", "97.5%"]
    or_table = pd.DataFrame({
        "feature": params.index,
        "beta": params.values,
        "OR": np.exp(params.values),
        "CI2.5%": np.exp(conf["2.5%"].values),
        "CI97.5%": np.exp(conf["97.5%"].values),
        "pvalue": res.pvalues.values,
    })

    if drop_intercept_in_output and "const" in or_table["feature"].values:
        or_table = or_table[or_table["feature"] != "const"].copy()

    or_table = or_table.sort_values("pvalue").reset_index(drop=True)

    # 14) 결과 저장
    if out_csv:
        or_table.to_csv(out_csv, index=False)

    # 15) 진단 정보 리턴
    final_vif = compute_vif(add_constant(X.drop(columns=["const"]), has_constant='add')) \
                if "const" in X.columns else compute_vif(X)

    info = {
        "fitted_with": fitted_with,
        "corr_dropped": corr_dropped,
        "vif_dropped": vif_dropped,
        "final_features": [c for c in X.columns if c != "const"],
        "final_vif": final_vif.sort_values("VIF", ascending=False).reset_index(drop=True)
    }
    return res, or_table, info


# 사용 예시: 경로만 바꿔서 실행하세요.
output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'
endpoint = 'PLT'
# endpoint = 'Hb'
# endpoint = 'WBC'
# endpoint = 'ANC'

csv_path = f"{output_dir}/lnz_mvlreg_{endpoint}_df.csv"     # <-- 본인 파일 경로
out_csv  = f"{output_dir}/lnz_mvlreg_{endpoint}_res.csv"




# raise ValueError
res, or_table, info = fit_ev_logistic(
    csv_path=csv_path,
    out_csv=out_csv,
    corr_threshold=0.7,     # |r| 임계값 (보통 0.8~0.95)
    vif_threshold=13.0,     # VIF 임계값 (보통 5 또는 10)
    enable_corr_filter=True,
    enable_vif_filter=True,
    verbose=True
)

print(f"[Fitted with] / {endpoint} / {info['fitted_with']}")
print(f"\n[Model Summary ({endpoint})]")
print(res.summary())
print(f"\n[Odds Ratio Table ({endpoint})]")
print(or_table)
print(f"\nSaved OR table to: {out_csv}")

print(f"\n[Correlation dropped ({endpoint})]")
print(info["corr_dropped"])
print("\n[VIF dropped ({endpoint})]")
print(info["vif_dropped"])
print("\n[Final features ({endpoint})]")
print(info["final_features"])
print("\n[Final VIF ({endpoint})]")
print(info["final_vif"])
