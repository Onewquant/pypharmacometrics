import os
import re
import math
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import pynca
from datetime import datetime, timedelta

## OLS result report

import statsmodels.api as sm
from fpdf import FPDF
from sklearn.preprocessing import StandardScaler

# Standardized Coefficient 계산
def get_standardized_coefficients(X, y):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    y_scaled = (y - y.mean()) / y.std()
    model_std = sm.OLS(y_scaled, X_scaled).fit()
    return model_std.params

# 결과 DataFrame 만들기
def ols_result_df(result, X, y):
    # standardized_coefs = get_standardized_coefficients(X, y)
    df = pd.DataFrame({
        'Coefficient': result.params,
        'Std.Err': result.bse,
        't-value': result.tvalues,
        'P-value': result.pvalues,
        'CI Lower': result.conf_int()[0],
        'CI Upper': result.conf_int()[1],
        # 'Standardized Coef': standardized_coefs
    }).round(4)
    return df

# PDF 리포트 생성 함수
def create_pdf_report(result, result_df, filepath="OLS_Report.pdf"):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", "B", 16)
    pdf.cell(0, 10, "OLS Multivariable Linear Regression Report", ln=True, align="C")
    pdf.ln(10)

    # 모델 요약
    pdf.set_font("Arial", "", 12)
    pdf.multi_cell(0, 8,
        f"Dependent Variable: {result.model.endog_names}\n"
        f"Number of Observations: {int(result.nobs)}\n"
        f"Degrees of Freedom: {int(result.df_model)} (Model), {int(result.df_resid)} (Residual)\n"
        f"R-squared: {result.rsquared:.4f}\n"
        f"Adjusted R-squared: {result.rsquared_adj:.4f}\n"
        f"F-statistic: {result.fvalue:.4f}\n"
        f"Prob (F-statistic): {result.f_pvalue:.4g}\n"
        f"AIC: {result.aic:.2f}, BIC: {result.bic:.2f}"
    )
    pdf.ln(10)

    # Coefficient Table
    pdf.set_font("Arial", "B", 12)
    pdf.cell(0, 10, "Coefficients Summary", ln=True)
    pdf.set_font("Arial", "", 10)
    col_width = pdf.w / (len(result_df.columns)+2)

    # Header
    row_height = 8
    pdf.cell(col_width, row_height, ' ', border=1, align="C")
    for col in result_df.columns:
        pdf.cell(col_width, row_height, col, border=1, align="C")
    pdf.ln()

    # Rows
    for index, row in result_df.iterrows():
        pdf.cell(col_width, row_height, str(index), border=1, align="C")
        for item in row:
            pdf.cell(col_width, row_height, str(item), border=1, align="C")
        pdf.ln()

    # Save
    pdf.output(filepath)
    print(f"✅ PDF Report saved as: {filepath}")



## Basic functions

def load_data_dict(drug_list, filename_format, input_file_dir_path):
    drug_prep_df_dict = dict()
    for drug in drug_list:
        result_file_path = f"{input_file_dir_path}/" + filename_format.replace('[drug]',drug)
        if filename_format.split('.')[-1]=='csv':
            drug_prep_df_dict[drug] = pd.read_csv(result_file_path)
        if filename_format.split('.')[-1] == 'xls':
            drug_prep_df_dict[drug] = pd.read_excel(result_file_path)
    return drug_prep_df_dict


def dates_every_step_days(start_date, day_step, addl, fmt="%Y-%m-%d", include_start=True):
    """
    start_date: 'YYYY-MM-DD' 문자열
    n         : 날짜간격으로 더하는 횟수 (기본 10회 → 10개 날짜)
    fmt       : 반환 문자열 포맷. None이면 date 객체로 반환
    include_start: True면 시작일도 리스트의 첫 항목으로 포함
    """
    base = datetime.strptime(start_date, "%Y-%m-%d").date()

    # 구간 생성
    if include_start:
        dates = [base + timedelta(days=day_step * i) for i in range(0, addl)]
    else:
        dates = [base + timedelta(days=day_step * i) for i in range(1, addl + 1)]

    return [d.strftime(fmt) for d in dates] if fmt else dates

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

def dosing_duration_for_abs_policy(input_dur, input_abs_pol):
    if input_abs_pol.upper()=='ZERO':
        return input_dur
    elif input_abs_pol.upper() == 'FIRST':
        return '.'
    elif input_abs_pol.upper() == 'ELANG':
        return '.'
    elif input_abs_pol.upper() == 'TRANSIT':
        return '.'
    else:
        return '.'

    return None

def dosing_cmt_for_advan_type(advan=0, route='', forced_dosing_cmt=np.nan, abs_policy=''):
    # Specific ADVAN
    if (advan== 1):                          # 1 comp / IV or zero-order abs
        return 1
    elif (advan == 2) and (route == 'IV'):   # 1 comp / ExtraVascular (IV로 주므로 central comp에)
        return 2
    elif (advan == 2) and (route != 'IV'):   # 1 comp / ExtraVascular (IV아니므로 depot comp에)
        return 1
    elif (advan == 3) and (route == 'IV'):   # 2 comp / IV (IV 이므로 central comp에)
        return 1
    elif (advan == 3) and (route != 'IV'):   # 2 comp / IV - zero-order abs 만 가능
        return 1
    elif (advan == 4) and (route == 'IV'):   # 2 comp / ExtraVascular (IV로 주므로 central comp에)
        return 2
    elif (advan == 4) and (route != 'IV'):   # 2 comp / ExtraVascular (IV아니므로 depot comp에)
        return 1
    elif (advan == 11) and (route == 'IV'):   # 3 comp / IV (IV 이므로 central comp에)
        return 1
    elif (advan == 11) and (route != 'IV'):   # 3 comp / IV - zero-order abs 만 가능
        return 1
    elif (advan == 12) and (route == 'IV'):   # 3 comp / ExtraVascular (IV로 주므로 central comp에)
        return 2
    elif (advan == 12) and (route != 'IV'):   # 3 comp / ExtraVascular (IV아니므로 depot comp에)
        return 1
    # General ADVAN
    elif (advan in (5, 6, 7, 8, 9, 13)):

        ## Erlang model인 경우
        erlang_patterns = re.findall(r'^ERLANG[\d]+$',abs_policy.upper())
        if len(erlang_patterns)>0:
            erlang_num = int(erlang_patterns[0].replace('ERLANG',''))
            return 1

        ## 기타 다른 모델인 경우
        if (type(forced_dosing_cmt) in (float, int)):
            if not np.isnan(forced_dosing_cmt):
                return forced_dosing_cmt
            else:
                raise ValueError("Dosing Compartment 결정시 / forced_dosing_cmt가 NAN입니다.")
        else:
            raise ValueError("Dosing Compartment 결정시 / forced_dosing_cmt를 숫자로 정확히 입력하세요.")
    else:
        raise ValueError("Dosing Compartment 결정시 / ADVAN과 ROUTE를 정확히 입력하세요.")


def sampling_cmt_for_specific_advan_type(advan=0, forced_sampling_cmt=np.nan, abs_policy=''):
    # Specific ADVAN
    if (advan == 1):     # 1 comp / IV
        return 1
    elif (advan == 2):   # 1 comp / ExtraVascular
        return 2
    elif (advan == 3):   # 2 comp / IV
        return 1
    elif (advan == 4):   # 2 comp / ExtraVascular
        return 2
    elif (advan == 11):   # 3 comp / IV
        return 1
    elif (advan == 12):   # 3 comp / ExtraVascular
        return 2
    # General ADVAN
    elif (advan in (5, 6, 7, 8, 9, 13)):
        ## Erlang model인 경우
        erlang_patterns = re.findall(r'^ERLANG[\d]+$', abs_policy.upper())
        if len(erlang_patterns) > 0:
            erlang_num = int(erlang_patterns[0].replace('ERLANG', ''))
            erlang_sampling_cmt_num = 2+erlang_num
            return erlang_sampling_cmt_num

        ## 기타 다른 모델인 경우
        if (type(forced_sampling_cmt) in (float, int)):
            if not np.isnan(forced_sampling_cmt):
                return forced_sampling_cmt
            else:
                raise ValueError("Sampling Compartment 결정시 / forced_sampling_cmt가 NAN입니다.")
        else:
            raise ValueError("Sampling Compartment 결정시 / forced_sampling_cmt를 숫자로 정확히 입력하세요.")
    else:
        raise ValueError("Sampling Compartment 결정시 / ADVAN을 정확히 입력하세요.")

## DDI

## 겹치는 / 안 겹치는 Person_ID 및 Date 생성함수

def get_intsec_or_diff_pid_date(df1, df2, mode='INTSEC', dc_name='DC_NAME', id_col='UID', start_date_col='START_DATE',end_date_col='END_DATE',addl_col='ADDL'):
    """
    ## df1 및 df2 필요컬럼
    ['PERSON_ID', 'DRUG_EXPOSURE_START_DATETIME', 'DRUG_EXPOSURE_END_DATETIME', 'DAYS_SUPPLY']

    ## 배열을 이용한 Vector 연산
    # 배열(pid x date)만들어 (최소, 최대 date 날짜범위로) drug1, drug2 모두 사용한 경우의 pid_date만 추려냄
    # row : pid -> exist_get_intsec_pids
    # col : date -> total_dates
    # val : count -> 처음엔 zero

    # mode :
    # 'intsec': 두 df에서 겹치는 pid_date 반환
    # 'diff' : d1-d2의 차이를 나타내는 pid_date 반환
    """

    exist_gen_intsec_pids = pd.concat([df1[[id_col]], df2[[id_col]]]).drop_duplicates(keep='first').sort_values([id_col],ignore_index=True)[id_col]

    total_min_sdate = min(df1[start_date_col].min(), df2[start_date_col].min())[:10]
    total_max_edate = max(df1[end_date_col].max(), df2[end_date_col].max())[:10]

    total_pids_modi = exist_gen_intsec_pids.astype(str) + '_'
    total_dates = list(pd.date_range(start=total_min_sdate, end=total_max_edate, closed='left').strftime('%Y-%m-%d')) + [total_max_edate]

    # 배열의 index로의 translation

    pid_to_rinx_dict = {pid: rinx for rinx, pid in enumerate(exist_gen_intsec_pids)}
    date_to_cinx_dict = {pid: cinx for cinx, pid in enumerate(total_dates)}

    # pid x date 배열지도

    total_date_array = np.tile(total_dates, (len(exist_gen_intsec_pids), 1))
    total_pid_array = np.tile(total_pids_modi, (len(total_dates), 1)).T

    total_map_array = total_pid_array + total_date_array

    # hit 여부 대시보드

    pid_date_dashboard_dict = dict()

    count = 1
    for longform_df in df1, df2:
        pid_date_dashboard_dict[count] = np.zeros((len(exist_gen_intsec_pids), len(total_dates)))
        longform_df = longform_df[[id_col, start_date_col, addl_col]].copy()
        # longform_df['DRUG_EXPOSURE_START_DATETIME'] = longform_df['DRUG_EXPOSURE_START_DATETIME'].map(lambda x:x[:10])

        for tlf_inx, tlf_row in longform_df.iterrows():
            s_cinx = date_to_cinx_dict[tlf_row[start_date_col][:10]]
            rinx = pid_to_rinx_dict[tlf_row[id_col]]
            adays = tlf_row[addl_col]

            pid_date_dashboard_dict[count][rinx, s_cinx:s_cinx + adays] = 1
        count += 1

    if mode == 'INTSEC':
        pid_date_dashboard = pid_date_dashboard_dict[1] + pid_date_dashboard_dict[2]
        intsec_result = total_map_array[pid_date_dashboard >= 2]
    elif mode == 'DIFF':
        pid_date_dashboard = pid_date_dashboard_dict[1] - pid_date_dashboard_dict[2]
        intsec_result = total_map_array[pid_date_dashboard >= 1]
    else:
        raise ValueError('Mode 재설정 필요')

    # intersection pid date 정보 저장

    intsec_df = pd.DataFrame(intsec_result, columns=[f'{id_col}_DATE'])
    intsec_df['DRUG_COMB'] = dc_name
    intsec_df[id_col] = intsec_df[f'{id_col}_DATE'].map(lambda x: int(x.split('_')[0]))
    intsec_df[f'{mode}_DATE'] = intsec_df[f'{id_col}_DATE'].map(lambda x: x.split('_')[1])
    intsec_df = intsec_df[['DRUG_COMB', id_col, f'{mode}_DATE']].sort_values([id_col, f'{mode}_DATE'])

    return intsec_df


# intsec_df = get_intsec_or_diff_pid_date(df1, df2, mode='intsec', dc_name='dc_str')
# diff_df = get_intsec_or_diff_pid_date(df2, df1, mode='diff', dc_name='dc_str')


## DataReading

def modeling_dosing_policy(mdpolicy_file_path, selected_models=[], model_colname='MODEL'):

    dspol_df = pd.read_csv(mdpolicy_file_path)
    if len(selected_models)!=0:
        dspol_df = dspol_df[dspol_df[model_colname].isin(selected_models)].reset_index(drop=True)
    dspol_df = dspol_df.applymap(lambda x:convert_to_numeric_value(x))
    return dspol_df


def get_drug_conc_data_dict_of_multiple_projects(project_dict, prepconc_dir_path, conc_filename_format="[project]_ConcPrep_[drug](R).csv"):
    project_drugconc_dict = dict()
    drugconc_dict = dict()
    for projectname, drug_list in project_dict.items():
        filename_format = conc_filename_format.replace('[project]', projectname)

        # ct_project_dir_path = modeling_dir_path + '/' + projectname

        project_drugconc_dict[projectname] = load_data_dict(drug_list, filename_format=filename_format, input_file_dir_path=prepconc_dir_path)

        ## Drug별 데이터 따로 모으기

        for drug, drugconc_df in project_drugconc_dict[projectname].items():
            drugconc_df['PROJECT'] = projectname
            if drug not in drugconc_dict:
                drugconc_dict[drug] = drugconc_df.copy()
            else:
                drugconc_dict[drug] = pd.concat([drugconc_dict[drug], drugconc_df], ignore_index=True)

            # drugconc_df = drugconc_df.drop(['PROJECT'],axis=1)
            # drugconc_df['SEQUENCE'] = drugconc_df['ID'].map(lambda x:x[0]).map({'A':1,'B':2})
            # drugconc_df.to_csv(f'{modeling_dir_path}/{projectname}/{projectname}_ConcPrep_{drug}(R).csv', index=False, encoding='utf-8-sig')

    return drugconc_dict

def generate_nonmem_subject_id(df, sid_col, uid_cols):
    ## ID Column

    if len(uid_cols)==0:
        dcdf = df.copy()
    else:
        dcdf = df[list(set([sid_col]+uid_cols))].copy()
    dcdf['NONMEM_ID'] = '1'  # 처음에 0이 안나오도록 설정: 사람이면 1, 동물이면 2로 하면 될듯

    # uid (Project 내의 고유 대상자 번호를 포함한 ID) 생성
    # uid에 구분할 수 있도록 포함시키려는 Column의 Value 값들에 숫자로 구분하여 배치표(dictionary)를 생성해두고, Nonmem에서 사용할 ID 생성
    uid_inf_dict = dict()
    for uicol in uid_cols:
        if len(dcdf[uicol].unique()) > 10:
            raise ValueError(f"{uicol} 컬럼의 Unique Value가 10개 이상입니다. 새로운 구분법이 필요합니다.")
        else:
            if str in dcdf[uicol].map(lambda x: type(x)).unique():
                uid_inf_dict[uicol] = {c: i for i, c in enumerate(dcdf[uicol].unique())}
            else:
                uid_inf_dict[uicol] = {c: c for c in dcdf[uicol].unique()}

        dcdf['NONMEM_ID'] += dcdf[uicol].map(lambda x: str(uid_inf_dict[uicol][x]))

    dcdf['NONMEM_ID'] += dcdf[sid_col]
    # dcdf['NONMEM_ID'] = dcdf['NONMEM_ID'].map(int)
    # dcdf['NONMEM_UID'] = dcdf['NONMEM_ID'].copy()
    return dcdf['NONMEM_ID'].copy()

def get_basic_nonmem_code(prep_df, dspol_df, modeling_dir_path, output_dir_name='modeling_basic_codes'):
    # prep_df = drugconc_dict['SGLT2INH']

    ess_dspol_cols = ['MODEL','DRUG','ROUTE','RATE','DOSING_CMT','SAMPLING_CMT','ABS_POLICY','ELIM_POLICY','ADVAN','TRANS']
    dspol_df_for_codes = dspol_df[ess_dspol_cols].drop_duplicates(ess_dspol_cols, ignore_index=True)
    basic_code_dir = f"{modeling_dir_path}/{output_dir_name}"
    for inx, row in dspol_df_for_codes.iterrows(): #break
        basic_code_path = f"{basic_code_dir}/{row['MODEL']}_{row['DRUG']}.txt"
        basic_code = ''

        PROBLEM = f"$PROBLEM {row['MODEL']}_{row['DRUG']}\n\n"
        INPUT = ("$INPUT "+ str(list(prep_df.columns)).replace("['","").replace("']","").replace("', '"," ") + "\n\n").replace('UID','UID=DROP')
        DATA = f"$DATA ..//MDP_{row['MODEL']}_{row['DRUG']}.csv IGNORE=@\n\n"

        SUBROUTINES = ""
        PK = ""

        MODEL = ""
        DES = ""
        ERROR = "$ERROR\nIPRED = F\nW = SQRT(THETA(?)**2*IPRED**2 + THETA(?)**2)\nY = IPRED + W*EPS(1)\nIRES = DV-IPRED\nIWRES = IRES/W"

        if row['ADVAN'] in [1, 2, 3, 4, 11, 12]:
            SUBROUTINES = f"$SUBROUTINES ADVAN{row['ADVAN']} TRANS{row['TRANS']}\n\n"
            PK = "$PK\n\n"
            MODEL = ""
            DES = ""
        else:
            if row['ADVAN'] in [5, 7]:
                SUBROUTINES = f"$SUBROUTINES ADVAN{row['ADVAN']}\n\n"
                MODEL = "$MODEL\n\n\n\n"
                DES = ""

                ## ERLANG 모델인 경우
                # row['ABS_POLICY'] = 'Erlang6'
                erlang_patterns = re.findall(r'^ERLANG[\d]+$', row['ABS_POLICY'].upper())
                if len(erlang_patterns) > 0:
                    erlang_num = int(erlang_patterns[0].replace('ERLANG',''))
                    MODEL = "$MODEL\nCOMP=(DEPOT,DEFDOSE)\n"
                    PK = "$PK\nK12=THETA(1)*EXP(ETA(1))\n"
                    for erl_comp_num in range(1,erlang_num+1):
                        erlang_comp_frag = f"COMP=(DELA{erl_comp_num})\n"

                        if erl_comp_num==erlang_num:
                            erlang_rateconstant_frag = ''
                        else:
                            erlang_rateconstant_frag = f"K{erl_comp_num+1}{erl_comp_num+2}=K12\n"

                        MODEL+=erlang_comp_frag
                        PK += erlang_rateconstant_frag

                    MODEL+='COMP=(CENTRAL, DEFOBS)\n\n'
                    PK+=f'CL=THETA(2)*EXP(ETA(2))\nV{erlang_num+2}=THETA(3)*EXP(ETA(3))\n\n'


            else:
                SUBROUTINES = f"$SUBROUTINES ADVAN{row['ADVAN']} TOL=4\n\n"
                MODEL = "$MODEL\n\n\n\n"
                DES = "$DES\n\n\n\n"

        THETA = "$THETA\n\n\n\n"
        OMEGA = "$OMEGA\n\n\n\n"
        SIGMA = "$SIGMA\n\n\n\n"
        EST = "$EST METHOD=1 INTER MAXEVAL=9999 NOABORT SIG=3 PRINT=1 POSTHOC\n"
        COV = "$COV\n\n"

        Xpose = "; Xpose\n"
        SDTAB = "$TABLE ID AMT TAD TIME RATE DV MDV IPRED IWRES CWRES ONEHEADER NOPRINT FILE=sdtab004"
        PATAB = "$TABLE CL KA V2 V3 Q ONEHEADER NOPRINT  FIRSTONLY FILE=patab004"
        COTAB = ""
        CATAB = ""
        TABLE = f"{SDTAB}\n{PATAB}\n{COTAB}\n{CATAB}\n"

        basic_code += (PROBLEM + INPUT + DATA + SUBROUTINES + MODEL + PK + DES + THETA + OMEGA + SIGMA + EST + COV + Xpose + TABLE)

        if not os.path.exists(f"{modeling_dir_path}/{output_dir_name}"):
            os.mkdir(f"{modeling_dir_path}/{output_dir_name}")
        with open(basic_code_path, 'w', encoding='utf-8-sig') as f:
            f.write(basic_code)




def formatting_data_nca_to_nonmem(drugconc_dict, dspol_df, uid_cols, modeling_dir_path, covar_cols=[],
                                  add_covar_df=pd.DataFrame(columns=['UID']), uid_on=True,
                                  term_dict={'TIME': 'ATIME', 'TAD': 'NTIME', 'DV': 'CONC', 'ID': 'ID'}, output_dir_name='modeling_prep_data'):
    """
    # add_covar_df 를 추가할때는 uid_cols를 반드시 포함하는 primary key로서 사용
    """

    uid_cols = list(set(uid_cols))
    cat_covar_dict = dict()

    unique_models = dspol_df['MODEL'].unique()
    for model_name in unique_models:  #break
        for drug in drugconc_dict.keys():  #break

            model_drug_dspol = dspol_df[(dspol_df['MODEL'] == model_name) & (dspol_df['DRUG'] == drug)].copy()
            # if drug=='Metformin': break
            dcdf = drugconc_dict[drug].copy()

            ## conc data에 이미 covar이 존재하는 경우 covar col 남김 (But, categorical의 경우 numeric으로 바꿔서)

            drug_covar_cols = list()
            for covar_c in covar_cols:
                if covar_c in dcdf.columns:
                    numeric_covar = True
                    frag_covar_val_dict = dict()
                    for ucovar_c in dcdf[covar_c].unique():
                        if type(ucovar_c) in [str, object]:
                            numeric_covar=False
                    if not numeric_covar:
                        frag_covar_val_dict = {covar_val:covar_inx+1 for covar_inx, covar_val in enumerate(dcdf[covar_c].unique())}
                        dcdf[f'NONMEM_{covar_c}'] = dcdf[covar_c].map(frag_covar_val_dict)
                    else:
                        dcdf[f'NONMEM_{covar_c}'] = dcdf[covar_c].copy()
                    cat_covar_dict[covar_c] = frag_covar_val_dict
                    drug_covar_cols.append(covar_c)
                else:
                    print(f"{model_name} / {drug} / Raw data에 {covar_c}라는 컬럼이 존재하지 않습니다.")
                    continue

            ## Result columns 설정
            if uid_on:
                result_cols = ['ID'] + ['UID'] + ['TAD', 'TIME', 'DV', 'MDV', 'AMT', 'RATE', 'DUR', 'CMT'] + drug_covar_cols
            else:
                result_cols = ['ID'] + ['TAD', 'TIME', 'DV', 'MDV', 'AMT', 'RATE', 'DUR', 'CMT'] + drug_covar_cols

            ## add_covar_df 추가하는 경우 (uid_cols 제외하고) (####### numeric 아닌경우 numeric으로 바꾸고 cat_covar_dict에 해당 내용 추가하는 코드 수정해야함)

            add_covar_df_copy = add_covar_df.copy()
            new_add_covar_df_cols = list()
            if len(uid_cols) != 0:
                for add_var_df_col in add_covar_df_copy.columns:
                    if add_var_df_col in uid_cols:
                        new_add_covar_df_cols.append(add_var_df_col)
                    else:
                        new_add_covar_df_cols.append(f"NONMEM_{add_var_df_col}")
                add_covar_df_copy.columns = new_add_covar_df_cols

            if len(add_covar_df)!=0:
                dcdf = dcdf.merge(add_covar_df_copy, on=uid_cols, how='left')

            ## nonmem_subject_id 생성

            dcdf['NONMEM_ID'] = generate_nonmem_subject_id(df=dcdf, sid_col=term_dict['ID'], uid_cols=uid_cols)
            if uid_on:
                dcdf['NONMEM_UID'] = dcdf['NONMEM_ID'].copy()

            # ## ID Column
            #
            # dcdf['NONMEM_ID']='1'  # 처음에 0이 안나오도록 설정: 사람이면 1, 동물이면 2로 하면 될듯
            #
            # # uid (Project 내의 고유 대상자 번호를 포함한 ID) 생성
            # # uid에 구분할 수 있도록 포함시키려는 Column의 Value 값들에 숫자로 구분하여 배치표(dictionary)를 생성해두고, Nonmem에서 사용할 ID 생서
            # for uicol in uid_inf_cols:
            #     if len(dcdf[uicol].unique())>10:
            #         raise ValueError(f"{uicol} 컬럼의 Unique Value가 10개 이상입니다. 새로운 구분법이 필요합니다.")
            #     else:
            #         if str in dcdf[uicol].map(lambda x:type(x)).unique():
            #             uid_inf_dict[uicol] = {c:i for i, c in enumerate(dcdf[uicol].unique())}
            #         else:
            #             uid_inf_dict[uicol] = {c:c for c in dcdf[uicol].unique()}
            #
            #     dcdf['NONMEM_ID']+=dcdf[uicol].map(lambda x:str(uid_inf_dict[uicol][x]))
            # dcdf['NONMEM_ID']+=dcdf['ID'].map(lambda x:re.findall(r'[\d]+', x)[0])
            # dcdf['NONMEM_ID']=dcdf['NONMEM_ID'].map(int)
            # dcdf['NONMEM_UID']=dcdf['NONMEM_ID'].copy()

            ## TIME Columns

            dcdf['NONMEM_TIME'] = dcdf[term_dict['TIME']].copy()
            dcdf['NONMEM_TAD'] = dcdf[term_dict['TAD']].copy()

            ## DV, MDV Column

            dcdf['NONMEM_DV'] = dcdf[term_dict['DV']].copy()
            dcdf['NONMEM_MDV'] = 0

            ## AMT, RATE, DUR

            dcdf['NONMEM_AMT'] = '.'

            ## CMT / DUR / RATE

            dspol_first_row = model_drug_dspol.iloc[0]

            dcdf['NONMEM_RATE'] = '.'
            dcdf['NONMEM_DUR'] = '.'
            dcdf['NONMEM_CMT'] = sampling_cmt_for_specific_advan_type(advan=dspol_first_row['ADVAN'], forced_sampling_cmt=dspol_first_row['SAMPLING_CMT'], abs_policy=dspol_first_row['ABS_POLICY'])

            """
            # 주로 농도측정이 Systemic compartment에서 진행되므로 2로 표시하는게 나을듯. 나중에 수정필요할수도
            # dosing policy 파일에 이 정보가 반영되도록 만들어야. 
            # (dosing policy / modeling policy 가 모두 영향을 주는 듯)
            """

            # Dosing row 추가 : Project / Drug 종류에 따라 Dosing Policy 다를 수 있다고 가정. (Drug은 위에서 dcdf로 이미 구분되어있음)

            drug_nmprep_df = list()
            for projectname, prj_dcdf in dcdf.groupby(['PROJECT']): #break

                # 해당 dosing policy
                prj_dspol = model_drug_dspol[model_drug_dspol['PROJECT']==projectname].reset_index(drop=True)

                for nmid, nmid_df in prj_dcdf.groupby(['NONMEM_ID']):   #break

                    nmid_df = nmid_df.reset_index(drop=True)

                    for dspol_inx, dspol_row in prj_dspol.iterrows():  # break

                        # 추가할 Dosing row

                        add_row = nmid_df.iloc[0:1, :].copy()
                        add_row['NONMEM_TIME'] = dspol_row['RELTIME']
                        add_row['NONMEM_TAD'] = dspol_row['RELTIME']
                        add_row['NONMEM_DV'] = '.'
                        add_row['NONMEM_MDV'] = 1
                        add_row['NONMEM_CMT'] = dosing_cmt_for_advan_type(advan=dspol_row['ADVAN'], route=dspol_row['ROUTE'], forced_dosing_cmt=dspol_row['DOSING_CMT'],abs_policy=dspol_row['ABS_POLICY'])
                        add_row['NONMEM_AMT'] = dspol_row['DOSE']

                        # add_row.iloc[0]
                        """
                        # RATE=-1 : Rate 추정 코드 추가해야
                        # RATE=-2 : Dur 추정 코드 추가해야

                        # [ADVAN 1 - 1구획 iv] : 1CMT (TRANS1 - K / TRANS2 - CL, V)
                        # [ADVAN 2 - 1구획 PO] : 1CMT-depot, 2CMT-central (TRANS1 - K, KA / TRANS2 - CL, V, KA)
                        # [ADVAN 3 - 2구획 iv] : 1CMT-central, 2CMT-peripheral (TRANS1 - K, K12, K21 / TRANS3 - CL, V, Q, VSS / TRANS4 - CL, V1, Q, V2 / TRANS5 - AOB, ALPHA, BETA / TRANS6 - ALPHA, BETA, K21)
                        # [ADVAN 4 - 2구획 PO] : 1CMT-depot, 2CMT-central, 3CMT-peripheral (TRANS1 - K, K23, K32, KA / TRANS3 - CL, V, Q, VSS, KA / TRANS4 - CL, V2, Q, V3, KA / TRANS5 - AOB, ALPHA, BETA, KA / TRANS6 - ALPHA, BETA, K32, KA)
                        # [ADVAN 11 - 3구획 iv] : 1CMT-central, 2CMT-peripheral1, 3CMT-peripheral2
                        # [ADVAN 12 - 3구획 PO] : 1CMT-depot, 2CMT-central, 3CMT-peripheral1, 4CMT-peripheral2
                        """

                        """
                        [ADVAN1 TRANS2 (1-comp IV세팅, dosing : CMT 1)] 
                        - AMT(O), RATE( _ or point ), DUR( _ or point ) => 안 돌아감
                        - AMT(O), RATE(-2), DUR(point) => 돌아감

                        [ADVAN2 TRANS2 (1-comp PO세팅, dosing : CMT 2)] 
                        - AMT(O), RATE( _ or point ), DUR( _ or point ) => 돌아감
                        """

                        add_row['NONMEM_RATE'] = dspol_row['RATE']
                        # add_row['NONMEM_DUR'] = dspol_row['DUR'] if dspol_row['ABS_POLICY'].upper() in ['ZERO']
                        add_row['NONMEM_DUR'] = dosing_duration_for_abs_policy(input_dur=dspol_row['DUR'], input_abs_pol=dspol_row['ABS_POLICY'])
                        # add_row['NONMEM_CMT'] = dosing_cmt_for_advan_type(advan=dspol_row['ADVAN'], route=dspol_row['ROUTE'])

                        # DUR값을 음수인 경우 -> 에러처리
                        if (dspol_row['DUR'] not in ('.', 0)):
                            if dspol_row['DUR'] < 0:
                                raise ValueError(f"Dosing Policy에서 {dspol_row['ROUTE'].upper()} 투여인데 DUR 값이 음수로 기록됨.")
                            else:
                                pass

                        if (dspol_row['DOSE'] < 0):  ### 여기
                            raise ValueError(f"DOSE가 음수로 기록되었습니다.")
                        else:
                            if (dspol_row['DOSE'] > 0) and (dspol_row['RATE'] in (0, '.')):
                                # print('Bolus와 같이 추정합니다. / DUR이 필요한 경우일지도. / 아직 어떻게 진행될지 모름 -> 해보자')
                                add_row['NONMEM_RATE'] = dspol_row['RATE']
                                add_row['NONMEM_DUR'] = dspol_row['DUR']
                            else:
                                if (dspol_row['RATE'] == -1):
                                    # print('Absorption의 Rate을 추정합니다. / NONMEM 코드에 Rn = theta(m) 구문을 추가하세요')
                                    add_row['NONMEM_RATE'] = -1
                                    add_row['NONMEM_DUR'] = '.'
                                elif (dspol_row['RATE'] == -2):
                                    # print('Absorption의 Duration을 추정합니다. / NONMEM 코드에 Dn = theta(m) 구문을 추가하세요')
                                    add_row['NONMEM_DUR'] = '.'
                                elif (dspol_row['RATE'] > 0):
                                    # print('양수인 RATE을 가지고 있습니다. / 코드에 어떻게 추가해야하는지 아직 모름 -> 해보자')
                                    add_row['NONMEM_DUR'] = '.'
                                else:
                                    raise ValueError(f"RATE, DURATION 값을 확인하세요.")

                        # 추가할 Dosing row를 위치할 곳 설정

                        if dspol_row['RELPOSITION'] not in (0, 1):
                            raise ValueError(f"Dosing Policy에서 RELPOSITION 은 0 또는 1이어야함. 해당 row 기준으로 (0: 직전 / 1: 직후)에 추가 row 삽입")

                        ## 원하는 RELTIME 기준이 데이터상에 존재하지 않는 경우 처리

                        if dspol_row['RELPOSITION'] == 1:
                            rt_inx_ds = nmid_df[f"NONMEM_{dspol_row['RELTIMECOL']}"] <= dspol_row['RELTIME']
                            if len(nmid_df[rt_inx_ds])==0:
                                continue
                            rt_inx = nmid_df[rt_inx_ds].iloc[-1].name
                            add_pos_inx = rt_inx + dspol_row['RELPOSITION'] ###############################################

                        elif dspol_row['RELPOSITION']==0:
                            rt_inx_ds = nmid_df[f"NONMEM_{dspol_row['RELTIMECOL']}"] >= dspol_row['RELTIME']
                            if len(nmid_df[rt_inx_ds])==0:
                                continue
                            rt_inx = nmid_df[rt_inx_ds].iloc[0].name
                            add_pos_inx = rt_inx + dspol_row['RELPOSITION'] ###############################################


                        if add_pos_inx < 0:
                            raise ValueError(f"Dosing row를 추가하려는 위치가 0번째 row보다 작은값입니다. / {projectname}, {drug}, {nmid}")
                        elif add_pos_inx == 0:
                            nmid_df = pd.concat([add_row, nmid_df], ignore_index=True)
                        elif (add_pos_inx > 0) and (add_pos_inx < len(nmid_df)):
                            nmid_df = pd.concat([nmid_df.loc[:add_pos_inx - 1, :].copy(), add_row, nmid_df.loc[add_pos_inx:, :].copy()], ignore_index=True)
                        elif (add_pos_inx == len(nmid_df)):
                            nmid_df = pd.concat([nmid_df, add_row], ignore_index=True)
                        else:
                            raise ValueError(f"Dosing row를 추가하려는 위치가 마지막 row보다 큰 값입니다. / {projectname}, {drug}, {nmid}")

                    while nmid_df['NONMEM_MDV'].iloc[-1]==1:
                        nmid_df = nmid_df.iloc[:-1]

                    drug_nmprep_df.append(nmid_df)

            # Prep Data 생성

            drug_nmprep_df = pd.concat(drug_nmprep_df, ignore_index=True)

            # Prep Data 컬럼편집

            drug_nonmem_cols = [c for c in drug_nmprep_df.columns if c.split('_')[0] == 'NONMEM']
            drug_nmprep_df = drug_nmprep_df[drug_nonmem_cols].copy()
            drug_nmprep_df.columns = [c[7:] for c in drug_nonmem_cols]

            # ID 값을 작은 값으로 재조정 및 dictionary 생성

            short_id_dict = {id: short_id + 1 for short_id, id in enumerate(drug_nmprep_df['ID'].unique())}
            drug_nmprep_df['ID'] = drug_nmprep_df['ID'].map(short_id_dict)

            # Prep Data 에 Covariate정보 추가 (추후 구현)

            # Prep Data 저장

            if not os.path.exists(f"{modeling_dir_path}/{output_dir_name}"):
                os.mkdir(f"{modeling_dir_path}/{output_dir_name}")
            drug_nmprep_df[result_cols].to_csv(f"{modeling_dir_path}/{output_dir_name}/MDP_{model_name}_{drug}.csv", index=False, encoding='utf-8-sig')
            # cat_covar_dict

            print(f"{model_name} / {drug} / data-formatting completed")

            # Basic Code 저장

            get_basic_nonmem_code(prep_df=drug_nmprep_df[result_cols], dspol_df=model_drug_dspol, modeling_dir_path=modeling_dir_path)



            # CMT, RATE, DUR, EVID, SS, ADDL

            ## Covariates


## TAD 컬럼 추가 함수

def add_time_after_dosing_column(df):
    ## 최근 투약 기준 시간으로 TAD 설정
    tad_df = df.copy()
    tad_df['TAD'] = np.nan

    # 농도값 바로 전에 투약기록 있는 경우
    tad_cond0 = (tad_df['MDV'] == 0)
    tad_cond1 = (tad_df['MDV'].shift(1).fillna(1.0) == 1)
    tad_cond2 = (tad_df['ID'] == tad_df['ID'].shift(1))
    tad_cond3 = (tad_df['TAD'].isna())
    for inx in tad_df[tad_cond0 & tad_cond1 & tad_cond2 & tad_cond3].index:
        tad_df.at[inx, 'TAD'] = tad_df.at[inx, 'TIME'] - tad_df.at[inx - 1, 'TIME']

    # 농도값 바로 전에 농도기록 있는 경우
    tad_cond0 = (tad_df['MDV'] == 0)
    tad_cond1 = (tad_df['MDV'].shift(1).fillna(1.0) == 1)
    tad_cond2 = (tad_df['ID'] == tad_df['ID'].shift(1))
    tad_cond3 = (tad_df['TAD'].isna())
    tad_cond4 = (~(tad_df['TAD'].shift(1).isna()))
    while len(tad_df[tad_cond0 & (~tad_cond1) & tad_cond2 & tad_cond3 & tad_cond4]) != 0:

        # covar_modeling_df[covar_modeling_df['TAD']<0]
        # covar_modeling_df.to_csv(f'{output_dir}/{drug}_checkcheck.csv',index=False, encoding='utf-8-sig')
        for inx in tad_df[tad_cond0 & (~tad_cond1) & tad_cond2 & tad_cond3 & tad_cond4].index:
            tad_df.at[inx, 'TAD'] = tad_df.at[inx - 1, 'TAD'] + tad_df.at[inx, 'TIME'] - tad_df.at[inx - 1, 'TIME']

        tad_cond0 = (tad_df['MDV'] == 0)
        tad_cond1 = (tad_df['MDV'].shift(1).fillna(1.0) == 1)
        tad_cond2 = (tad_df['ID'] == tad_df['ID'].shift(1))
        tad_cond3 = (tad_df['TAD'].isna())
        tad_cond4 = (~(tad_df['TAD'].shift(1).isna()))

    tad_df['TAD'] = tad_df['TAD'].replace(np.nan, 0)
    return tad_df['TAD']

## Modeling에 사용된 population의 simulation data set 생성 함수

def get_model_population_sim_df(df, interval=10/60, add_on_period=4*24, id_col='ID',time_col='TIME',dv_col='DV',mdv_col='MDV',cmt_col='CMT',amt_col='AMT',rate_col='RATE', tad_col='TAD'):

    final_sim_df = list()
    loop_cnt = 0
    for id, id_sim_df in df.groupby(id_col): #break
        # raise ValueError
        print(f"({loop_cnt}) {id}")
        id_conc_rows = id_sim_df[id_sim_df[mdv_col]==0].copy()
        id_conc_rows[dv_col] = np.nan
        id_conc_rows['REALDATA'] = 1
        id_dose_rows = id_sim_df[id_sim_df[mdv_col]==1].copy()
        id_dose_rows['REALDATA'] = 1

        internal_start_time = id_sim_df[time_col].min()
        internal_end_time = id_sim_df[time_col].max()
        external_end_time = internal_end_time + add_on_period

        time_arr = np.arange(internal_start_time, external_end_time + interval, interval)  # stop보다 큰 수까지 포함하려면 stop + step
        time_arr = time_arr[time_arr <= external_end_time]  # 정확히 stop보다 작거나 같은 값만 남기기
        uniq_time_arr = list(set(time_arr)-set(id_conc_rows[time_col]))
        uniq_time_arr.sort()

        added_sim_df = pd.DataFrame(columns=id_conc_rows.columns)
        added_sim_df[time_col] = uniq_time_arr
        added_sim_df[id_col] = id
        # added_sim_df['DV'] = np.nan
        # added_sim_df['DV'] = np.nan
        added_sim_df[mdv_col] = 0
        added_sim_df[cmt_col] = np.int64(id_conc_rows[cmt_col].median())
        added_sim_df['AMT'] = '.'
        added_sim_df['RATE'] = '.'
        added_sim_df['REALDATA'] = 0

        final_id_sim_df = pd.concat([id_conc_rows, id_dose_rows, added_sim_df]).sort_values([time_col,mdv_col])
        final_id_sim_df = final_id_sim_df.fillna(method='ffill').fillna(method='bfill')
        # id_conc_rows['RATE'].iloc[0]
        # id_conc_rows['CMT']

        # added_sim_df.columns
        final_sim_df.append(final_id_sim_df)
        loop_cnt+=1

    final_sim_df = pd.concat(final_sim_df, ignore_index=True)
    final_sim_df[tad_col] = final_sim_df[time_col]
    return final_sim_df


## 구간 AUC 구하기

def auc_linear_up_log_down(t0, c0, t1, c1):
    """Linear-up / Log-down trapezoid between two points (t0,c0) -> (t1,c1)."""
    dt = float(t1 - t0)
    if dt <= 0:
        return 0.0

    # log trapezoid (감소 구간, 두 농도 모두 >0)
    if (c0 > 0) and (c1 > 0) and (c1 < c0):
        ratio = c0 / c1
        if np.isclose(ratio, 1.0):
            return 0.5 * (c0 + c1) * dt
        return dt * (c0 - c1) / np.log(ratio)

    # linear trapezoid
    return 0.5 * (c0 + c1) * dt


def add_iAUC(df, time_col="TIME", conc_col="DV", id_col="ID"):
    """
    Linear-up/Log-down 방식으로 각 interval AUC(iAUC)와 누적 AUC(cumAUC)를 추가.
    첫 행(i=0)은 이전 구간이 없으므로 iAUC는 NaN.
    """
    # 정렬
    if id_col in df.columns:
        df_sorted = df.sort_values([id_col, time_col]).reset_index(drop=True)
    else:
        df_sorted = df.sort_values([time_col]).reset_index(drop=True)

    def _compute_group(g):
        t = g[time_col].to_numpy(dtype=float)
        c = g[conc_col].to_numpy(dtype=float)
        iAUC = np.full_like(c, np.nan, dtype=float)

        for i in range(1, len(g)):
            iAUC[i] = auc_linear_up_log_down(t[i-1], c[i-1], t[i], c[i])

        out = g.copy()
        out["iAUC"] = iAUC
        out["cumAUC"] = np.nancumsum(iAUC)  # 누적 AUC
        return out

    if id_col in df_sorted.columns:
        out = df_sorted.groupby(id_col, group_keys=False).apply(_compute_group)
    else:
        out = _compute_group(df_sorted)

    return out

