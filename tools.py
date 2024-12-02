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


def dosing_cmt_for_advan_type(advan=0, route='', forced_dosing_cmt=np.nan):
    # Specific ADVAN
    if (advan == 1):
        return 1
    elif (advan == 2) and (route == 'IV'):
        return 1
    elif (advan == 2) and (route != 'IV'):
        return 2
    elif (advan == 3) and (route == 'IV'):
        return 1
    elif (advan == 3) and (route != 'IV'):
        return 1
    elif (advan == 4) and (route == 'IV'):
        return 2
    elif (advan == 4) and (route != 'IV'):
        return 1
    elif (advan == 11) and (route == 'IV'):
        return 1
    elif (advan == 11) and (route != 'IV'):
        return 1
    elif (advan == 12) and (route == 'IV'):
        return 2
    elif (advan == 12) and (route != 'IV'):
        return 1
    # General ADVAN
    elif (advan in (5, 6, 7, 8, 9, 13)):
        if (type(forced_dosing_cmt) in (float, int)):
            if not np.isnan(forced_dosing_cmt):
                return forced_dosing_cmt
            else:
                raise ValueError("Dosing Compartment 결정시 / forced_dosing_cmt가 NAN입니다.")
        else:
            raise ValueError("Dosing Compartment 결정시 / forced_dosing_cmt를 숫자로 정확히 입력하세요.")
    else:
        raise ValueError("Dosing Compartment 결정시 / ADVAN과 ROUTE를 정확히 입력하세요.")


def sampling_cmt_for_specific_advan_type(advan=0, forced_sampling_cmt=np.nan):
    # Specific ADVAN
    if (advan == 1):
        return 1
    elif (advan == 2):
        return 2
    elif (advan == 3):
        return 1
    elif (advan == 4):
        return 2
    elif (advan == 11):
        return 1
    elif (advan == 12):
        return 2
    # General ADVAN
    elif (advan in (5, 6, 7, 8, 9, 13)):
        if (type(forced_sampling_cmt) in (float, int)):
            if not np.isnan(forced_sampling_cmt):
                return forced_sampling_cmt
            else:
                raise ValueError("Sampling Compartment 결정시 / forced_sampling_cmt가 NAN입니다.")
        else:
            raise ValueError("Sampling Compartment 결정시 / forced_sampling_cmt를 숫자로 정확히 입력하세요.")
    else:
        raise ValueError("Sampling Compartment 결정시 / ADVAN을 정확히 입력하세요.")
