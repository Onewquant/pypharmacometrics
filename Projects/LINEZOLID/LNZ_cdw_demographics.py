import pandas as pd
import numpy as np

output_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results'
multivar_res_df = list()
for endpoint in ['PLT', 'Hb', 'WBC', 'ANC', 'Lactate']:
# for endpoint in ['ANC', ]:
    # endpoint = 'PLT'
    # endpoint = 'Hb'
    # endpoint = 'WBC'
    # endpoint = 'ANC'
    # endpoint = 'Lactate'

    csv_path = f"{output_dir}/datacheck/lnz_datacheck_{endpoint}_df.csv"     # <-- 본인 파일 경로
    # out_csv  = f"{output_dir}/mvlreg_output/lnz_mvlreg_{endpoint}_res.csv"
    df_ori = pd.read_csv(csv_path)
    # df_ori.columns

    raise ValueError