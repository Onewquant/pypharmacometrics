import pandas as pd
import numpy as np
from typing import Dict, Any

output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'
comed_df = pd.read_csv(f"{output_dir}/lnz_final_comed_df.csv")
transfusion_df = pd.read_csv(f"{output_dir}/lnz_final_transfusion_df.csv")
surv_res_df = pd.read_csv(f"{output_dir}/b1da/b1da_lnz_surv_res_df(365).csv")

multivar_res_df = list()

demo_cols = ['SUBTOTAL_N','EV = 1','SEX = 1','ELD = 1','AGE','HT','WT','BMI']
total_dose_cols = ['CUM_DOSE(TOTAL)','DOSE_PERIOD(TOTAL)','DOSE24(TOTAL_ACTIVE)','DOSE24(TOTAL_PERIOD)','DOSE24PERWT(TOTAL_ACTIVE)','DOSE24PERWT(TOTAL_PERIOD)','DOSE_INTERVAL(TOTAL)']
ev_dose_cols = ['CUM_DOSE', 'DOSE_PERIOD','DOSE_PERIOD(ACTIVE)', 'DOSE24','DOSE24(ACTIVE)', 'DOSE24PERWT','DOSE24PERWT(ACTIVE)', 'DOSE_INTERVAL']

comed_cols = list(comed_df['DRUG'].drop_duplicates())
transfusion_cols = list(transfusion_df['TF_TYPE'].drop_duplicates())


total_df = list()
for inx, endpoint in enumerate(['PLT', 'Hb', 'WBC', 'ANC', 'Lactate']):
    csv_path = f"{output_dir}/B1DA/datacheck/b1da_lnz_datacheck_{endpoint}_df.csv"  # <-- 본인 파일 경로
    # out_csv  = f"{output_dir}/mvlreg_output/lnz_mvlreg_{endpoint}_res.csv"
    df = pd.read_csv(csv_path)
    df['SEX'] = df['SEX'].map({'남': 0, '여': 1})
    total_df.append(df.copy())

total_df = pd.concat(total_df)
# str 컬럼을 모두 'TOTAL'로 치환
str_cols = total_df.select_dtypes(include="object").columns
total_df[str_cols] = "TOTAL"

# UID별로 median 집계
# total_df = total_df.groupby("UID").median(numeric_only=True).reset_index()
total_df = total_df.groupby("UID").mean(numeric_only=True).reset_index()
total_df['EV'] = (total_df['EV']>0)*1

# str 컬럼은 다시 TOTAL 붙여줌
for col in str_cols:
    total_df[col] = "TOTAL"


# total_df.columns
# total_df


demo_table_df = list()
for inx, endpoint in enumerate(['TOTAL','PLT', 'Hb', 'WBC', 'ANC', 'Lactate']):
# for endpoint in ['ANC', ]:
    # endpoint = 'PLT'
    # endpoint = 'Hb'
    # endpoint = 'WBC'
    # endpoint = 'ANC'
    # endpoint = 'Lactate'
    if endpoint=='TOTAL':
        df = total_df.copy()
        df['SEX'] = df['SEX'].map({0.0:0,1.0:1})
    else:
        csv_path = f"{output_dir}/B1DA/datacheck/b1da_lnz_datacheck_{endpoint}_df.csv"     # <-- 본인 파일 경로
        df = pd.read_csv(csv_path)
        df['SEX'] = df['SEX'].map({'남':0,'여':1})
    # total_df.append(df.copy())
    # df_ori.columns

    # 제외할 컬럼과 변수 유형 지정
    exclude_cols = ["UID", "BL_DATE", "ENDPOINT"]
    categorical_cols = ["EV", "SEX", "ELD"] + comed_cols + transfusion_cols      # 범주형
    positive_levels: Dict[str, Any] = {    # N(%)로 표시할 수준
        "EV": 1,
        "SEX": 1,
        "ELD": 1
        # "SEX": 'F'
    }

    # 분석 대상 컬럼만 선택
    cols = [c for c in df.columns if c not in exclude_cols]
    work = df[cols].copy()

    # 실제 존재하는 범주형만 유지
    categorical_cols = [c for c in categorical_cols if c in work.columns]

    # 연속형 = 나머지
    continuous_cols = [c for c in work.columns if c not in categorical_cols]

    # 연속형 숫자 변환(문자/결측 안전 처리)
    for c in continuous_cols:
        work[c] = pd.to_numeric(work[c], errors="coerce")

    rows = []
    N_total = len(work)
    rows.append({
                    "Variable": 'SUBTOTAL_N', "Type": "-",
                    "Summary": f"{N_total} ({round(100*N_total/len(total_df),1)})", "N (non-missing)": N_total
                })

    # 연속형: mean (SD), non-missing N
    for c in continuous_cols:
        s = work[c].dropna()
        if s.empty:
            rows.append({
                "Variable": c, "Type": "Continuous",
                "Summary": "NA (NA)", "N (non-missing)": 0
            })
        else:
            rows.append({
                "Variable": c, "Type": "Continuous",
                "Summary": f"{s.mean():.1f} ({s.std(ddof=1):.1f})",
                "N (non-missing)": int(s.shape[0])
            })

    # 범주형(EV, SEX): 지정 수준의 N(%) 표시
    for c in categorical_cols:
        level = positive_levels.get(c, 1)
        s = work[c]
        n_pos = int((s == level).sum())
        pct_pos = (n_pos / N_total * 100.0) if N_total > 0 else np.nan
        rows.append({
            "Variable": f"{c} = {level}" if c not in (comed_cols+transfusion_cols) else c,
            "Type": "Categorical",
            "Summary": f"{n_pos} ({pct_pos:.1f})",
            "N (non-missing)": int(s.notna().sum())
        })

    demo_table = pd.DataFrame(rows)
    # # 보기 좋게 정렬: 범주형 → 연속형
    # order_map = {"Categorical": 0, "Continuous": 1}
    # demo_table["__o__"] = demo_table["Type"].map(order_map).fillna(2)
    # demo_table = demo_table.sort_values(["__o__", "Variable"]).drop(columns="__o__").reset_index(drop=True)

    demo_table_frag = demo_table[['Variable','Type','Summary']].copy()
    demo_table_frag = demo_table_frag.rename(columns={'Summary':f'Group({endpoint})'})
    if inx==0:
        demo_table_df = demo_table_frag
    else:
        demo_table_df = demo_table_df.merge(demo_table_frag, on=['Variable','Type'], how='left')



# pd.concat(demo_table_df, axis=1)
#
# print(demo_table)
# 파일로 저장 (원하면)
rest_cols = list(set(demo_table_df['Variable'].unique()) - set(demo_cols) - set(total_dose_cols) - set(ev_dose_cols) - set(comed_cols) - set(transfusion_cols))
rest_cols.sort()
demo_table_df = demo_table_df.set_index("Variable").loc[demo_cols+rest_cols+total_dose_cols+ev_dose_cols+comed_cols+transfusion_cols].reset_index(drop=False)
demo_table_df.to_csv(f"{output_dir}/B1DA/b1da_lnz_demographic_df.csv", index=False)

