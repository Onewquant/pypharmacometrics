import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.tools import add_constant
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from statsmodels.stats.outliers_influence import variance_inflation_factor
import os
from scipy.stats import kruskal
import scikit_posthocs as sp


output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'

def cliffs_delta(x, y):
    """
    Cliff's delta: 두 비모수 표본 간 효과크기
    해석 기준 (Romano 2006):
      |delta| < 0.147 : negligible
      0.147 ≤ |delta| < 0.33 : small
      0.33 ≤ |delta| < 0.474 : medium
      |delta| ≥ 0.474 : large
    """
    m, n = len(x), len(y)
    # x[i] > y[j] : +1 / x[i] < y[j] : -1 / 같으면 0
    gt = sum(xx > yy for xx in x for yy in y)
    lt = sum(xx < yy for xx in x for yy in y)
    delta = (gt - lt) / (m * n)
    return delta

# df['SEX'].sum()
# pd.crosstab(df['SEX'], df['EV'])

group_sum_df = list()
posthoc_df = list()
for endpoint in ['ANC']:

    csv_path = f"{output_dir}/mvlreg/lnz_mvlreg_{endpoint}_df.csv"     # <-- 본인 파일 경로
    df = pd.read_csv(csv_path)
    # df['DOSE24']

    for col,lab_range in {"SCR":(0.8, 1.25), "TBIL":(0.8, 1.2)}.items():
        for dose_endpoint in ['DOSE_PERIOD(TOTAL)', 'CUM_DOSE', 'DOSE24']:
            group_series = list()
            for lab_val in df[col]:
                if lab_val < lab_range[0]:
                    group_series.append(1)
                elif (lab_val >= lab_range[0]) and (lab_val < lab_range[1]):
                    group_series.append(2)
                elif (lab_val >= lab_range[1]):
                    group_series.append(3)
                else:
                    pass

            df['GROUP'] = group_series
            # raise ValueError
            # ------------------------------------------
            # 2. Kruskal-Wallis 전체 비교
            # ------------------------------------------
            grouped = [group[dose_endpoint].values for _, group in df.groupby('GROUP')]
            h_stat, p_val = kruskal(*grouped)

            # ------------------------------------------
            # 3. Dunn 사후분석 (Bonferroni 보정)
            # ------------------------------------------
            dunn_p = sp.posthoc_dunn(df, val_col=dose_endpoint, group_col='GROUP', p_adjust='bonferroni'
            )

            # ------------------------------------------
            # 4. 그룹별 요약 통계
            # ------------------------------------------
            summary = (
                df.groupby('GROUP')[dose_endpoint]
                    .agg(['median', 'mean', 'std', 'count',
                          lambda x: x.quantile(0.25),
                          lambda x: x.quantile(0.75)])
                    .rename(columns={'<lambda_0>': 'Q1', '<lambda_1>': 'Q3'})
            )
            summary['Kruskal_H'] = h_stat
            summary['Kruskal_p'] = p_val
            summary['Covar_type'] = col
            summary['Dose_ep'] = dose_endpoint
            summary = summary.reset_index(drop=False)

            # ------------------------------------------
            # 5. 사후분석 p-value + 효과크기(Cliff's delta)
            # ------------------------------------------
            rows = []
            for g1 in sorted(df['GROUP'].unique()):
                for g2 in sorted(df['GROUP'].unique()):
                    if g1 < g2:
                        x = df.loc[df['GROUP'] == g1, 'DOSE_PERIOD(TOTAL)']
                        y = df.loc[df['GROUP'] == g2, 'DOSE_PERIOD(TOTAL)']
                        delta = cliffs_delta(x, y)
                        rows.append({
                            "Group1": g1,
                            "Group2": g2,
                            "Dunn_p": dunn_p.loc[g1, g2],
                            "Cliffs_delta": delta
                        })

            pairwise_dunn = pd.DataFrame(rows)


            # 효과크기 해석 추가
            def interpret_delta(d):
                ad = abs(d)
                if ad < 0.147:
                    return "negligible"
                elif ad < 0.33:
                    return "small"
                elif ad < 0.474:
                    return "medium"
                else:
                    return "large"


            pairwise_dunn['Effect_size'] = pairwise_dunn['Cliffs_delta'].apply(interpret_delta)
            pairwise_dunn['Covar_type'] = col
            pairwise_dunn['Dose_ep'] = dose_endpoint
            pairwise_dunn = pairwise_dunn

            # ------------------------------------------
            # 6. 최종 result
            # ------------------------------------------
            result = {
                "group_summary": summary,
                "pairwise_dunn": pairwise_dunn
            }

            group_sum_df.append(summary.copy())
            posthoc_df.append(pairwise_dunn.copy())
            print("\n[그룹별 요약 + Kruskal-Wallis 결과]")
            print(result["group_summary"])
            print("\n[그룹간 Dunn 사후분석 + Cliff's delta 효과크기]")
            print(result["pairwise_dunn"])

            # raise ValueError

group_sum_df = pd.concat(group_sum_df, ignore_index=True)
kw_res_df = group_sum_df[['Covar_type', 'Dose_ep','GROUP','median','Kruskal_p']].copy()
kw_res_df.to_csv()
# group_sum_df.columns
# posthoc_df.columns
posthoc_df = pd.concat(posthoc_df, ignore_index=True)
dunn_res_df = posthoc_df[['Covar_type', 'Dose_ep', 'Group1', 'Group2', 'Dunn_p','Effect_size']].copy()


        # sns.boxplot(x='GROUP', y='DOSE_PERIOD(TOTAL)', data=df_ori)
        #
        # sns.swarmplot(x='GROUP', y='DOSE_PERIOD(TOTAL)', data=df, color=".25")
        # plt.show()
