import pandas as pd
import numpy as np
import glob
from scipy.stats import mannwhitneyu, fisher_exact


# 성별 차이 Subgroup analysis

# 데이터 로드

output_dir = '/Projects/LINEZOLID/results'

subgroup_files = glob.glob(f"{output_dir}/mvlreg_output/lnz_mvlreg_datasubset(*)(*).csv")
for file_path in subgroup_files:
    # raise ValueError
    subset_group = file_path.split(')(')[0].split('(')[-1]
    pd_endpoint = file_path.split(')(')[-1].split(')')[0]

    df = pd.read_csv(file_path)

    group_col = 'SEX'      # 그룹 변수
    binary_col = 'EV'      # 이분형 변수
    numeric_cols = [c for c in df.columns if c not in [group_col, binary_col]]

    rows = []
    rows.append({'Covar':'N',
                 'M':len(df[df[group_col]==0]),
                 'F':len(df[df[group_col]==1]),
                 'p_value': 0,
                 })

    # 1️⃣ 연속형 변수: mean (SD)와 Mann-Whitney U p-value
    for col in numeric_cols:
        g0 = df.loc[df[group_col] == 0, col].dropna()
        g1 = df.loc[df[group_col] == 1, col].dropna()
        pval = mannwhitneyu(g0, g1, alternative='two-sided')[1] if len(g0) and len(g1) else np.nan
        rows.append({
            'Covar': col,
            'M': f"{g0.mean():.2f} ({g0.std():.2f})",
            'F': f"{g1.mean():.2f} ({g1.std():.2f})",
            'p_value': round(pval,4)
        })

    # 2️⃣ EV 컬럼: N (%)와 Fisher의 정확 검정
    g0_total = (df[group_col] == 0).sum()
    g1_total = (df[group_col] == 1).sum()
    g0_ev1 = df.loc[df[group_col] == 0, binary_col].sum()
    g1_ev1 = df.loc[df[group_col] == 1, binary_col].sum()

    contingency = [[g0_ev1, g0_total - g0_ev1],
                   [g1_ev1, g1_total - g1_ev1]]
    pval_ev = fisher_exact(contingency)[1]

    rows.append({
        'Covar': binary_col,
        'M': f"{g0_ev1} ({g0_ev1 / g0_total * 100:.1f}%)",
        'F': f"{g1_ev1} ({g1_ev1 / g1_total * 100:.1f}%)",
        'p_value': round(pval_ev,4)
    })

    # 3️⃣ 결과 테이블 생성
    demographic_table = pd.DataFrame(rows).sort_values('p_value').reset_index(drop=True)
    demographic_table.to_csv(f"{output_dir}/mvlreg_output/lnz_subgrpanalysis({subset_group})({pd_endpoint}).csv", index=False, encoding='utf-8-sig')
    print(demographic_table)
