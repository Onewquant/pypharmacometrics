import pandas as pd
import numpy as np
from typing import Dict, Any

output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'
multivar_res_df = list()

demo_cols = ['SUBTOTAL_N','EV = 1','SEX = 1','AGE']
dose_cols = ['DOSE_PERIOD(TOTAL)','DOSE_PERIOD','DOSE24','DOSE_INTERVAL','CUM_DOSE']


for endpoint in ['PLT', 'Hb', 'WBC', 'ANC', 'Lactate']:
# for endpoint in ['ANC', ]:
    # endpoint = 'PLT'
    # endpoint = 'Hb'
    # endpoint = 'WBC'
    # endpoint = 'ANC'
    # endpoint = 'Lactate'

    csv_path = f"{output_dir}/B1DA/datacheck/b1da_lnz_datacheck_{endpoint}_df.csv"     # <-- 본인 파일 경로
    # out_csv  = f"{output_dir}/mvlreg_output/lnz_mvlreg_{endpoint}_res.csv"
    df = pd.read_csv(csv_path)
    df['SEX'] = df['SEX'].map({'남':0,'여':1})
    # df_ori.columns

    # 제외할 컬럼과 변수 유형 지정
    exclude_cols = ["UID", "BL_DATE", "ENDPOINT"]
    categorical_cols = ["EV", "SEX"]       # 범주형
    positive_levels: Dict[str, Any] = {    # N(%)로 표시할 수준
        "EV": 1,
        "SEX": 1
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
                    "Summary": f"{N_total} ({round(100*N_total/1260,1)}%)", "N (non-missing)": N_total
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
                "Summary": f"{s.mean():.2f} ({s.std(ddof=1):.2f})",
                "N (non-missing)": int(s.shape[0])
            })

    # 범주형(EV, SEX): 지정 수준의 N(%) 표시
    for c in categorical_cols:
        level = positive_levels.get(c, 1)
        s = work[c]
        n_pos = int((s == level).sum())
        pct_pos = (n_pos / N_total * 100.0) if N_total > 0 else np.nan
        rows.append({
            "Variable": f"{c} = {level}", "Type": "Categorical",
            "Summary": f"{n_pos} ({pct_pos:.1f}%)",
            "N (non-missing)": int(s.notna().sum())
        })

    demo_table = pd.DataFrame(rows)
    # # 보기 좋게 정렬: 범주형 → 연속형
    # order_map = {"Categorical": 0, "Continuous": 1}
    # demo_table["__o__"] = demo_table["Type"].map(order_map).fillna(2)
    # demo_table = demo_table.sort_values(["__o__", "Variable"]).drop(columns="__o__").reset_index(drop=True)

    demo_table

    print(demo_table)
    # 파일로 저장 (원하면)
    demo_table.to_csv("/mnt/data/demographic_table_overall.csv", index=False)



    raise ValueError