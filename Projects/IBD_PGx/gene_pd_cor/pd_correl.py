import pandas as pd
import numpy as np
import glob
import os
from scipy.stats import mannwhitneyu, fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns

prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx'
resource_dir = f'{prj_dir}/gene_pd_cor'
output_dir = f"{prj_dir}/gene_pd_cor"

pid_df = pd.read_csv(f"{resource_dir}/pid_df.csv")
rsid_df = pd.read_csv(f"{resource_dir}/rsid_dosage_matrix_with_alleles.csv")
ep_df = pd.read_csv(f"{resource_dir}/for_genomics_df1(all_drugs)_prev.csv")

rsid_df.rename(columns={'sid':'s'}).merge(pid_df, on=['s'], how='left')
rsid_df
ep_df


