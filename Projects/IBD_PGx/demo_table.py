import numpy as np
import pandas as pd

prj_name = 'IBDPGX'
prj_dir = './Projects/IBD_PGx'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

lab_df = pd.read_csv(f"{output_dir}/conc_df(lab).csv")
cumlab_df = pd.read_csv(f"{output_dir}/conc_df(cum_lab).csv")