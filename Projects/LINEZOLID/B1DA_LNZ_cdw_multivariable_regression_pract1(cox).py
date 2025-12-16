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
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test


def _to_binary_series(s: pd.Series) -> pd.Series:
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
    if X.shape[1] < 2:
        return X, []

    corr = X.corr().abs()
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
    to_drop = set()

    while True:
        mask = (upper > threshold)
        if not mask.any().any():
            break

        pairs = np.argwhere(mask.values)
        mean_conn = (corr.sum() - 1) / (corr.shape[0] - 1)

        local_drop = []
        for i, j in pairs:
            vi = upper.index[i]
            vj = upper.columns[j]
            if vi in to_drop or vj in to_drop:
                continue
            local_drop.append(vi if mean_conn[vi] >= mean_conn[vj] else vj)

        if not local_drop:
            break

        for v in local_drop:
            if v not in to_drop:
                to_drop.add(v)
                corr = corr.drop(index=v, columns=v)
                upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
                if corr.shape[0] < 2:
                    break
        if corr.shape[0] < 2:
            break

    kept_cols = [c for c in X.columns if c not in to_drop]
    return X[kept_cols].copy(), sorted(list(to_drop))


def compute_vif_df(X: pd.DataFrame) -> pd.DataFrame:
    # const(상수항) 제외하고 계산
    cols = [c for c in X.columns if c != "const"]
    if len(cols) == 0:
        return pd.DataFrame(columns=["feature", "VIF"])
    X_ = X[cols].astype(float)
    vifs = [variance_inflation_factor(X_.values, i) for i in range(X_.shape[1])]
    return pd.DataFrame({"feature": cols, "VIF": vifs}).sort_values("VIF", ascending=False)


def drop_by_vif(X: pd.DataFrame, threshold: float = 10.0, max_iter: int = 50):
    dropped = []
    X_work = X.copy()
    for _ in range(max_iter):
        vif_df = compute_vif_df(X_work)
        if vif_df.empty:
            break
        max_vif = vif_df["VIF"].max()
        if not np.isfinite(max_vif) or max_vif <= threshold:
            break
        drop_feat = vif_df.iloc[0]["feature"]
        dropped.append(drop_feat)
        X_work = X_work.drop(columns=[drop_feat])
        if X_work.shape[1] <= 1:  # const만 남는 경우
            break
    return X_work, dropped


def fit_ev_cox(
    df: pd.DataFrame,
    duration_col: str,          # <-- Cox에서 필수
    event_col: str = "EV",      # <-- EV(0/1)
    out_csv: str = None,
    corr_threshold: float = 0.8,
    vif_threshold: float = 10.0,
    enable_corr_filter: bool = True,
    enable_vif_filter: bool = True,
    penalizer: float = 0.0,     # 수렴 불안정/준분리 있으면 0.01~0.1 권장
    filename_add: str = "",
    drop_intercept_in_output: bool = True,
):

    """
    duration_col='DOSE_PERIOD'
    event_col='EV'
    out_csv=out_csv
    corr_threshold=0.8     # |r| 임계값 (보통 0.8~0.95)
    vif_threshold=10.0     # VIF 임계값 (보통 5 또는 10)
    enable_vif_filter=True
    filename_add=age_subgroup
    enable_corr_filter = True
    penalizer = 0.05

    """
    if event_col not in df.columns:
        raise ValueError(f"'{event_col}' 컬럼이 필요합니다.")
    if duration_col not in df.columns:
        raise ValueError(f"'{duration_col}' 컬럼이 필요합니다 (time-to-event).")

    # 1) y(event), T(duration)
    E = _to_binary_series(df[event_col])
    T = df[duration_col].copy()

    # Cox는 T>0 필요(0이면 epsilon 처리)
    T = pd.to_numeric(T, errors="coerce")
    if T.isna().any():
        raise ValueError(f"{duration_col}에 NaN이 있습니다. 먼저 처리하세요.")
    T = T.clip(lower=1e-6)

    # 2) X: event/duration 제외 전부
    X = df.drop(columns=[event_col, duration_col]).copy()

    # 3) 숫자 처리: inf -> NaN -> median
    num_cols = X.select_dtypes(include=[np.number]).columns.tolist()
    if num_cols:
        X[num_cols] = X[num_cols].replace([np.inf, -np.inf], np.nan)
        X[num_cols] = X[num_cols].fillna(X[num_cols].median())

    # (원하면) 범주 처리: mode + get_dummies  (당신 코드에서 주석 처리돼 있어서 그대로 둠)
    # cat_cols = [c for c in X.columns if c not in num_cols]
    # if cat_cols:
    #     for c in cat_cols:
    #         mode_val = X[c].mode(dropna=True)
    #         fill_val = mode_val.iloc[0] if len(mode_val) else "MISSING"
    #         X[c] = X[c].fillna(fill_val).astype(str)
    #     X = pd.get_dummies(X, columns=cat_cols, drop_first=True)

    if X.shape[1] == 0:
        raise ValueError("EV, duration 외 설명변수가 없습니다.")

    # 4) 0-분산 제거
    zero_var_cols = [c for c in X.columns if X[c].nunique(dropna=False) <= 1]
    if zero_var_cols:
        X = X.drop(columns=zero_var_cols)

    # 5) 상관 기반 제거 (Cox도 동일하게 가능)
    corr_dropped = []
    if enable_corr_filter and X.shape[1] > 1:
        X, corr_dropped = drop_by_correlation(X, threshold=corr_threshold)

    # 6) VIF 제거용으로 const 추가해서 VIF 계산 (Cox 모델에는 const 넣지 않음)
    X_for_vif = add_constant(X, has_constant="add")
    vif_dropped = []
    if enable_vif_filter and X_for_vif.shape[1] > 2:
        X_for_vif, vif_dropped = drop_by_vif(X_for_vif, threshold=vif_threshold, max_iter=100)
    # const 제거하고 최종 X 확정
    X_final = X_for_vif.drop(columns=["const"], errors="ignore").astype(float)

    # 7) CoxPH 적합
    cox_df = pd.concat([T.rename(duration_col), E.rename(event_col), X_final], axis=1)
    # cox_df.columns
    # cox_df['WBC'] = ((cox_df['WBC'] - cox_df['WBC'].mean())/cox_df['WBC'].std(ddof=0))
    # cox_df['WBC'] = cox_df['WBC'].clip(-5, 5)
    # cox_df['WBC'].isinf().sum()
    # cox_df.drop(['WBC'], axis=1)
    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(cox_df, duration_col=duration_col, event_col=event_col)

    # 2) Schoenfeld residuals 기반 PH test (time_transform은 보통 'rank' 또는 'km'를 사용)
    ph_test = proportional_hazard_test(cph, cox_df, time_transform="rank")

    # 3) 테이블로 정리
    # ph_test.summary에는 covariate별 test statistic / p 값이 들어있습니다.
    ph_table = ph_test.summary.reset_index().rename(columns={
        "index": "feature",
        "test_statistic": "chi2",
        "p": "pvalue"
    }).sort_values("pvalue", ignore_index=True)

    # print(ph_table)

    # (옵션) 유의한 변수만
    alpha = 0.05
    cox_viol_table = ph_table[ph_table["pvalue"] < alpha].copy()
    print(f"\n[PH violation candidates (p < 0.05)] / {len(cox_viol_table)} (총: {len(ph_table)}) / {list(cox_viol_table['feature'])}")
    # print()


    # cph.check_assumptions(
    #     cox_df.drop(['WBC'], axis=1),
    #     p_value_threshold=0.05,
    #     show_plots=True  # 각 rcovariate별 잔차-시간 플롯
    # )
    # cox_df['WBC'].sort_values()
    # 8) 결과 테이블(HR/CI/p)
    summ = cph.summary.reset_index().rename(columns={"covariate": "feature"})
    # lifelines summary: coef, exp(coef), p, exp(coef) lower 95%, upper 95%, ...
    out = pd.DataFrame({
        "feature": summ["feature"],
        "beta": summ["coef"],
        "aHR": summ["exp(coef)"],
        "aHR CI2.5%": summ["exp(coef) lower 95%"],
        "aHR CI97.5%": summ["exp(coef) upper 95%"],
        "pvalue": summ["p"],
    })

    # 9) 다중비교 보정 (Bonferroni; 당신 코드와 동일)
    _, p_adj, _, _ = multipletests(out["pvalue"].values, method="bonferroni")
    out["pvalue_adj"] = p_adj

    # 10) 단변량 Cox
    uni_results = []
    for col in X_final.columns:
        # pass
        try:
            tmp = pd.concat([T.rename(duration_col), E.rename(event_col), X_final[[col]]], axis=1)
            cph_u = CoxPHFitter(penalizer=penalizer)
            cph_u.fit(tmp, duration_col=duration_col, event_col=event_col)
            s = cph_u.summary.iloc[0]
            uni_results.append({
                "feature": col,
                "uni_beta": s["coef"],
                "uni_HR": s["exp(coef)"],
                "uni_CI2.5%": s["exp(coef) lower 95%"],
                "uni_CI97.5%": s["exp(coef) upper 95%"],
                "uni_pvalue": s["p"],
            })
        except Exception:
            uni_results.append({
                "feature": col,
                "uni_beta": np.nan,
                "uni_HR": np.nan,
                "uni_CI2.5%": np.nan,
                "uni_CI97.5%": np.nan,
                "uni_pvalue": np.nan,
            })

    uni_table = pd.DataFrame(uni_results)
    out = out.merge(uni_table, on="feature", how="left")

    out = out.sort_values("pvalue_adj").reset_index(drop=True)

    out["cox_viol"] = f"{'|'.join(list(cox_viol_table['feature']))}"
    out["cox_viol_len"] = f"{len(cox_viol_table)}"
    out["total_covars"] = f"{len(ph_table)}"

    if out_csv:
        rp_str = ".csv" if filename_add == "" else f"({filename_add}).csv"
        out.to_csv(out_csv.replace(".csv", rp_str), index=False)

    # 최종 VIF 리포트
    final_vif = compute_vif_df(add_constant(X_final, has_constant="add"))
    info = {
        "dropped_by_corr": corr_dropped,
        "dropped_by_vif": vif_dropped,
        "final_features": list(X_final.columns),
        "final_vif": final_vif.set_index("feature")["VIF"].to_dict(),
        "penalizer": penalizer,
    }
    return cph, out, info


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

    if not os.path.exists(f'{output_dir}/b1da/mvlreg_output(cox)'):
        os.mkdir(f'{output_dir}/b1da/mvlreg_output(cox)')

    out_csv = f"{output_dir}/b1da/mvlreg_output(cox)/b1da_lnz_mvlreg_{endpoint}_res.csv"

    # raise ValueError
    df_ori = pd.read_csv(csv_path)
    active_cols = [c for c in df_ori.columns if ('(TOTAL' not in c)]
    df_ori = df_ori[active_cols].copy()

    # df_ori.columns

    # raise ValueError

    # for age_subgroup in ['Adult','Elderly','Total_Adult']:
    for age_subgroup in ['Elderly','Total_Adult']:
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
    #     df['clozapine']
        df = df.drop(['clozapine','methimazole'], axis=1)

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
        # df.columns
        res, hr_table, info = fit_ev_cox(
            df=df,
            duration_col='DOSE_PERIOD',
            event_col='EV',
            out_csv=out_csv,
            corr_threshold=0.8,     # |r| 임계값 (보통 0.8~0.95)
            vif_threshold=10.0,     # VIF 임계값 (보통 5 또는 10)
            enable_corr_filter=True,
            enable_vif_filter=True,
            penalizer = 0.05,
            filename_add=age_subgroup,
        )

        # print(f"[Fitted with] / {endpoint} / {info['fitted_with']}")
        # print(f"\n[Model Summary ({endpoint})]")
        # print(res.summary())
        print(f"\n[Hazard Ratio Table ({endpoint})]")
        # or_threshold = 0.5
        # or_threshold = 0.0
        hr_table['aHR threshold'] = (hr_table['aHR'] >= 1).map(orthreshold_dict)
        hr_table['endpoint'] = endpoint
        hr_table['subgroup'] = age_subgroup
        hr_table['N'] = len(df)
        hr_table['EVN'] = int(sum(df['EV']))
        hr_table['EVPct'] = round(100*sum(df['EV'])/len(df),2)
        pv_cond = (hr_table['pvalue_adj'] < 0.05)
        hr_cond = ((np.abs(hr_table['aHR']-1) >= hr_table['aHR threshold']))
        # sig_res_frag = or_table[(or_table['pvalue'] < 0.05) & (np.abs(or_table['OR'] - 1) >= or_threshold)].copy()
        sig_res_frag = hr_table[pv_cond&hr_cond].copy()

        multivar_totres_df.append(hr_table.copy())

        info['endpoint'] = endpoint
        info['subgroup'] = age_subgroup
        covar_info_df.append(info.copy())

        # sig_res_frag = sig_res_frag[['subgroup','endpoint','N','EVN','EVPct']+list(sig_res_frag.columns)[:-4]]
        # multivar_res_df.append(sig_res_frag.copy())
        if not os.path.exists(f'{output_dir}/b1da/mvlreg_output(cox)/datasubset'):
            os.mkdir(f'{output_dir}/b1da/mvlreg_output(cox)/datasubset')
        df.to_csv(f"{output_dir}/b1da/mvlreg_output(cox)/datasubset/b1da_lnz_mvlreg_datasubset({age_subgroup})({endpoint}).csv", index=False, encoding='utf-8-sig')

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

multivar_totres_df.to_csv(f"{output_dir}/b1da/mvlreg_output(cox)/b1da_lnz_mvlreg_total_results.csv", index=False, encoding='utf-8-sig')

covar_info_df = pd.DataFrame(covar_info_df)
covar_info_df = covar_info_df[['subgroup','endpoint']+list(covar_info_df.columns)[0:-2]].copy()
covar_info_df_saving = pd.melt(covar_info_df,
                   id_vars=['subgroup','endpoint'],          # 녹이지 않고 그대로 둘 열
                   value_vars=['dropped_by_corr', 'dropped_by_vif', 'final_features', 'final_vif'],       # 녹일 열 (None이면 나머지 전부)
                   var_name='info_type',         # 녹인 열 이름
                   value_name="covar_list")    # 녹인 값의 열 이름)

# covar_info_df_saving.columns


covar_info_df_saving = covar_info_df_saving.sort_values(by=['subgroup','endpoint','info_type'],ascending=[False,True,True])
covar_info_df_saving = covar_info_df_saving[covar_info_df_saving['info_type']!='final_vif'].copy()
# covar_info_df_saving = covar_info_df_saving.drop(['fitted_with'],axis=1)
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
covar_info_df_saving.to_excel(f"{output_dir}/b1da/mvlreg_output(cox)/b1da_lnz_mvlreg_covar_info.xlsx", index=False, encoding='utf-8-sig')


sg_cond = (multivar_totres_df['subgroup'] == 'Total_Adult')
uor_cond = (np.abs(multivar_totres_df['uni_HR']-1) >= ((multivar_totres_df['uni_HR']>=1).map(orthreshold_dict)))
upv_cond = (multivar_totres_df['uni_pvalue'] < 0.05)
aor_cond = (np.abs(multivar_totres_df['aHR']-1) >= ((multivar_totres_df['aHR']>=1).map(orthreshold_dict)))
apv_cond = (multivar_totres_df['pvalue_adj'] < 0.05)

# multivar_totres_df.columns

sig_dig = 3

multivar_totres_df['EV_Count (%)'] = multivar_totres_df.apply(lambda x:f"{round(x['EVN'],sig_dig)} ({round(x['EVPct'],sig_dig)})",axis=1)

multivar_totres_df['HR (95% CI)'] = multivar_totres_df.apply(lambda x:f"{round(x['uni_HR'],sig_dig)} ({round(x['uni_CI2.5%'],sig_dig)}-{round(x['uni_CI97.5%'],sig_dig)})",axis=1)
multivar_totres_df['pval'] = multivar_totres_df['uni_pvalue'].map(lambda x:round(x,sig_dig))

multivar_totres_df['aHR (95% CI)'] = multivar_totres_df.apply(lambda x:f"{round(x['aHR'],sig_dig)} ({round(x['aHR CI2.5%'],sig_dig)}-{round(x['aHR CI97.5%'],sig_dig)})",axis=1)
multivar_totres_df['pval (adj)'] = multivar_totres_df['pvalue_adj'].map(lambda x:round(x,sig_dig))

# multivar_totres_df[(multivar_totres_df['subgroup']=='Total_Adult')&(multivar_totres_df['endpoint']=='ANC')&(multivar_totres_df['feature']=='SEX')][['OR (95% CI)','pval','aOR (95% CI)','pval (adj)']]
# multivar_totres_df[upv_cond&uor_cond]
# multivar_totres_df['pval']
uni_res_df = multivar_totres_df[upv_cond&uor_cond].sort_values(['subgroup','endpoint','aHR'],ascending=[False,True,False])[['subgroup','endpoint','feature','EV_Count (%)','HR (95% CI)','pval','aHR (95% CI)','pval (adj)']].reset_index(drop=True)
multi_res_df = multivar_totres_df[apv_cond&aor_cond].sort_values(['subgroup','endpoint','aHR'],ascending=[False,True,False])[['subgroup','endpoint','feature','EV_Count (%)','HR (95% CI)','pval','aHR (95% CI)','pval (adj)']].reset_index(drop=True)
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

uni_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output(cox)/b1da_lnz_mvlreg_univar_res_table.csv", index=False, encoding='utf-8-sig')
multi_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output(cox)/b1da_lnz_mvlreg_multivar_res_table.csv", index=False, encoding='utf-8-sig')

# uni_res_df[uni_res_df['subgroup']=='Total_Adult'][['endpoint','feature','OR (95% CI)','pval']]

tot_res_df = multivar_totres_df.sort_values(['subgroup','endpoint','aHR'],ascending=[False,True,False])[['subgroup','endpoint','feature','EV_Count (%)','HR (95% CI)','pval','aHR (95% CI)','pval (adj)']].reset_index(drop=True)
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
tot_res_df.to_csv(f"{output_dir}/b1da/mvlreg_output(cox)/b1da_lnz_mvlreg_total_res_table.csv", index=False, encoding='utf-8-sig')

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