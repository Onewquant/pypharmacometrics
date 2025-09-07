import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.tools import add_constant
from statsmodels.tools.sm_exceptions import PerfectSeparationError

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
    # 문자열/범주형
    uniq = s.dropna().unique()
    if len(uniq) == 2:
        mapping = {uniq[0]: 0, uniq[1]: 1}
        return s.map(mapping).astype(int)
    raise ValueError("EV는 이진이어야 합니다(문자열 2수준 가능).")

def fit_ev_logistic(csv_path: str,
                    out_csv: str = None,
                    drop_intercept_in_output: bool = True):
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
        X = X.drop(columns=zero_var_cols)

    # 8) 상수항 추가 및 형식 정리
    X = add_constant(X, has_constant='add')
    X = X.astype(float)
    y = y.astype(int)

    # 9) 적합 (Logit → 문제시 GLM Binomial 폴백)
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

    # 10) OR 테이블 구성
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

    # 11) 저장 옵션
    if out_csv:
        or_table.to_csv(out_csv, index=False)

    return res, or_table, fitted_with


# 사용 예시: 경로만 바꿔서 실행하세요.
output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'
out_csv  = f"{output_dir}/lnz_mvlreg_PLT_res.csv"

res, or_table, fitted_with = fit_ev_logistic(f"{output_dir}/lnz_mvlreg_PLT_df.csv", out_csv=out_csv)
print(f"[Fitted with] {fitted_with}")
print("\n[Model Summary]")
print(res.summary())
print("\n[Odds Ratio Table]")
print(or_table)
print(f"\nSaved OR table to: {out_csv}")
