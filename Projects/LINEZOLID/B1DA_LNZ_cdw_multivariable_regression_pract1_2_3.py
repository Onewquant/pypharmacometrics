import pandas as pd
import numpy as np
import glob
import os
from scipy.stats import mannwhitneyu, fisher_exact

os.environ['R_HOME'] = r'C:\Program Files\R\R-4.5.1'  # ← R 설치 경로로 변경

from pymer4.models import Lmer

# 성별 차이 Subgroup analysis

# 데이터 로드

output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'

# subgroup_files = glob.glob(f"{output_dir}/b1da/mvlreg_output/datasubset/b1da_lnz_mvlreg_datasubset(*)(*).csv")
subgroup_files = glob.glob(f"{output_dir}/b1da/mvlreg_output/datasubset/b1da_lnz_mvlreg_datasubset(Total_Adult)(Lactate).csv")
for file_path in subgroup_files:
    # raise ValueError
    subset_group = file_path.split(')(')[0].split('(')[-1]
    pd_endpoint = file_path.split(')(')[-1].split(')')[0]

    df = pd.read_csv(file_path)
    raise ValueError

    # 로짓 확률로부터 ADR 생성


    # rand_eff = {'A': -0.3, 'B': 0.1, 'C': 0.5, 'D': -0.1, 'E': 0.3}
    # linpred = -5 + 0.003 * df.dose + 0.5 * df.sex - 0.6 * df.albumin + 0.02 * df.age + df.site.map(rand_eff)
    # p = 1 / (1 + np.exp(-linpred))
    # df["ADR"] = np.random.binomial(1, p)
    df['PID'] = df.index
    str(df.columns).replace("', '"," + ").replace("',\n       '"," + ")
    # GLMM 적합
    model = Lmer("EV ~ SEX + AGE + ELD + Hb + PLT + WBC + SCR + ANC + ALB + TBIL + TPRO + AST + ALT + eGFR + CRP + GLU + Lactate + pH + HT + WT + BMI + CUM_DOSE + DOSE24PERWT + DOSE_INTERVAL + (1|PID)", data=df, family='binomial')
    print(model.fit())




    # # 3️⃣ 결과 테이블 생성
    # demographic_table = pd.DataFrame(rows).sort_values('p_value').reset_index(drop=True)
    # if not os.path.exists(f'{output_dir}/b1da/mvlreg_output/subgroup_analysis'):
    #     os.mkdir(f'{output_dir}/b1da/mvlreg_output/subgroup_analysis')
    # demographic_table.to_csv(f"{output_dir}/b1da/mvlreg_output/subgroup_analysis/b1da_lnz_subgrpanalysis({subset_group})({pd_endpoint}).csv", index=False, encoding='utf-8-sig')
    # print(demographic_table)
