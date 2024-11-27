import os
import re
import math
import numpy as np
import pandas as pd
import seaborn as sns

def load_data_dict(drug_list, filename_format, input_file_dir_path):
    drug_prep_df_dict = dict()
    for drug in drug_list:
        result_file_path = f"{input_file_dir_path}/" + filename_format.replace('[drug]',drug)
        if filename_format.split('.')[-1]=='csv':
            drug_prep_df_dict[drug] = pd.read_csv(result_file_path)
        if filename_format.split('.')[-1] == 'xls':
            drug_prep_df_dict[drug] = pd.read_excel(result_file_path)
    return drug_prep_df_dict

def convert_to_numeric_value(value):
    try:
        if int(value)==float(value):
            return int(value)
        else:
            return float(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return str(value)