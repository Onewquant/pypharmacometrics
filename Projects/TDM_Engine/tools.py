import numpy as np
import pandas as pd
from scipy.optimize import minimize
from numpy.linalg import inv


class EBEEnvironment:
    def __init__(self):
        self.PRED = None
        self.DATAi = None
        self.TH = None
        self.OM = None
        self.SG = None
        self.nEta = None
        self.invOM = None
        self.Fi = None
        self.Ri = None
        self.Hi = None
        self.Vi = None


e = EBEEnvironment()


def mGrad(func, x, nRec):

    # func = PRED
    # x = EBEi
    # nRec = len(y2)

    # func = PredVanco
    # x = np.array([ 0.45030954, -0.06429413, -0.70336727,  0.49293616])
    # nRec = 505

    n = len(x)
    gr = np.zeros((nRec, n))
    for i in range(n):
        hi = 1e-4 if abs(x[i]) <= 1 else 1e-4 * abs(x[i])
        ga = np.zeros((nRec, 4))
        for k in range(4):
            x1 = np.copy(x)
            x2 = np.copy(x)
            x1[i] -= hi
            x2[i] += hi
            ga[:, k] = (func(x2) - func(x1)) / (2 * hi)
            hi /= 2
        ga[:, 0] = (ga[:, 1] * 4 - ga[:, 0]) / 3
        ga[:, 1] = (ga[:, 2] * 4 - ga[:, 1]) / 3
        ga[:, 2] = (ga[:, 3] * 4 - ga[:, 2]) / 3
        ga[:, 0] = (ga[:, 1] * 16 - ga[:, 0]) / 15
        ga[:, 1] = (ga[:, 2] * 16 - ga[:, 1]) / 15
        gr[:, i] = (ga[:, 1] * 64 - ga[:, 0]) / 63
    return gr


def objEta(ETAi):
    DV = e.DATAi.loc[~e.DATAi['DV'].isna(), 'DV'].values
    pred = e.PRED(e.TH, ETAi, e.DATAi)
    e.Fi = pred[~e.DATAi['DV'].isna()]
    e.Ri = DV - e.Fi
    e.Hi = np.column_stack((e.Fi, np.ones_like(e.Fi)))
    e.Vi = np.diag(e.Hi @ e.SG @ e.Hi.T)
    penalty = ETAi @ e.invOM @ ETAi
    return np.sum(np.log(e.Vi) + (e.Ri ** 2 / e.Vi)) + penalty


def EBE(PRED, DATAi, TH, OM, SG):

    # PRED = PredVanco

    e.PRED = PRED
    e.DATAi = DATAi.copy()
    e.TH = TH
    e.OM = OM
    e.SG = SG
    e.nEta = OM.shape[0]
    e.invOM = inv(OM)

    res = minimize(objEta, x0=np.zeros(e.nEta), method="BFGS", options={'disp': False})  ########## 요라인만 약간 다른값
    EBEi = res.x
    COV = 2 * res.hess_inv
    SE = np.sqrt(np.diag(COV))
    Fi = PRED(TH, EBEi, DATAi)
    Ri = DATAi.loc[~DATAi['DV'].isna(), 'DV'].values - Fi[~DATAi['DV'].isna()]

    def PREDij(ETA): return PRED(TH, ETA, DATAi)

    gr1 = mGrad(PREDij, EBEi, len(Fi))
    VF = np.diag(gr1 @ COV @ gr1.T)
    SEy = np.sqrt(VF)
    SDy = np.sqrt(VF + VF * SG[0, 0] + Fi ** 2 * SG[0, 0] + SG[1, 1])

    return {
        "EBEi": EBEi,
        "SE": SE,
        "COV": COV,
        "IPRED": Fi,
        "SE.IPRED": SEy,
        "SD.IPRED": SDy,
        "IRES": Ri
    }


def calcPI(PRED, DATAi, TH, SG, rEBE, npoints=500):  ########## 요부분에서 result df의 각 row가 2배로 나옴. 에러수정 필요

    # DATAi_augmented = addDATAi(DATAi.copy(), TIME, AMT, RATE, II, ADDL)
    # rTab = calcPI(PRED, DATAi_augmented, TH, SG, rEBE, npoints)
    # PRED = PredVanco

    EBEi = rEBE["EBEi"]
    COV = rEBE["COV"]
    max_time = DATAi["TIME"].max()
    pred_times = np.linspace(0, max_time, npoints)
    DATAi2 = pd.merge(pd.DataFrame({"TIME": pred_times}), DATAi, on="TIME", how="outer").sort_values("TIME")
    y2 = PRED(TH, EBEi, DATAi2)  ############################### 여기가 제대로 안나오고 있음 !!!

    def PREDij(ETA): return PRED(TH, ETA, DATAi2)

    gr2 = mGrad(PREDij, EBEi, len(y2))
    if np.any(np.isnan(gr2)):
        raise ValueError("❌ Gradient (gr2) contains NA values.")

    VF2 = np.diag(gr2 @ COV @ gr2.T)
    SEy2 = np.sqrt(VF2)
    SDy2 = np.sqrt(VF2 + VF2 * SG[0, 0] + y2 ** 2 * SG[0, 0] + SG[1, 1])

    return pd.DataFrame({
        "x": DATAi2["TIME"],
        "y": DATAi2["DV"],
        "y2": y2,
        "yciLL": y2 - 1.96 * SEy2,
        "yciUL": y2 + 1.96 * SEy2,
        "ypiLL": y2 - 1.96 * SDy2,
        "ypiUL": y2 + 1.96 * SDy2,
    })


def convDT(DATAi):
    # 'DATE' 컬럼이 있는 경우만 변환
    if "DATE" in DATAi.columns:
        DATAi["SDT"] = pd.to_datetime(DATAi["DATE"].astype(str) + " " + DATAi["TIME"].astype(str))
        FDT = DATAi["SDT"].iloc[0]
        DATAi["TIME"] = (DATAi["SDT"] - FDT).dt.total_seconds() / 3600
        DATAi = DATAi.drop(columns=["DATE", "SDT"])
    return DATAi


def addDATAi(DATAi, TIME, AMT, RATE, II, ADDL):  ########## 확인필요
    lRow = DATAi.iloc[[-1]].copy()
    nADD = ADDL + 1
    for i in range(nADD):
        aRow = lRow.copy()
        aRow["TIME"] = TIME + i * II
        aRow["AMT"] = AMT
        aRow["RATE"] = RATE
        aRow["DV"] = np.nan
        aRow["MDV"] = 1
        DATAi = pd.concat([DATAi, aRow], ignore_index=True)
    return DATAi

def expandDATA(DATAo):  ########## 조금 다르게 수정해봤는데, 수정한 것이 맞을지 확인 필요
    eDATAi = pd.DataFrame(columns=DATAo.columns)
    Added_flags = []

    for i in range(len(DATAo)):
        row = DATAo.iloc[i]
        eDATAi = pd.concat([eDATAi, row.to_frame().T], ignore_index=True)

        if pd.notna(row.get("ADDL", np.nan)) and row["ADDL"] > 0:
            cTIME = row["TIME"]
            cII = row["II"]
            nADD = int(row["ADDL"])
            for j in range(1, nADD + 1):
                new_row = row.copy()
                new_row["TIME"] = cTIME + j * cII
                new_row['II'] = cII
                new_row['ADDL'] = 0
                eDATAi = pd.concat([eDATAi, new_row.to_frame().T], ignore_index=True)
                Added_flags.append(True)
            eDATAi.at[i,'ADDL'] = 0
        Added_flags.append(False)

    eDATAi = eDATAi.sort_values("TIME").reset_index(drop=True)

    CovCols = [col for col in eDATAi.columns if col not in ["ID", "TIME", "AMT", "RATE", "II", "ADDL", "DV", "MDV"]]
    for i in range(1, len(eDATAi)):
        if i < len(Added_flags) and Added_flags[i]:
            eDATAi.loc[i, CovCols] = eDATAi.loc[i - 1, CovCols]
    return eDATAi


def PredVanco(TH, ETA, DATAi):

    # ETA = EBEi
    # TH, ETA = EBEi, DATAi2  <- 2배로 복제되어 나온다.

    """
    # DATAi2.to_csv("./Projects/TDM_Engine/DATAi2.csv",index=False)


    import numpy as np
    import pandas as pd

    TH = np.array([
        3.8135955291021233,
        39.889510090195238,
        44.981835351176571,
        2.0055189192561507
    ])
    DATAi = pd.read_csv("./Projects/TDM_Engine/DATAi2.csv")
    ETA = np.array([0.45030954, -0.06429413, -0.70336727,  0.49293616])

    """

    V1 = TH[1] * np.exp(ETA[1])
    V2 = TH[2] * np.exp(ETA[2])
    Q = TH[3] * np.exp(ETA[3])
    K12 = Q / V1
    K21 = Q / V2

    Conc = pd.DataFrame([[0.0, 0.0]], columns=["C1", "C2"])
    IPRE = [0.0]
    inf = False

    pCLCR = min(DATAi.iloc[0]["CLCR"], 150)

    for i in range(1, len(DATAi)):
        # if i==2: break

        row = DATAi.iloc[i]
        prev = DATAi.iloc[i - 1]
        CLCR = min(row["CLCR"], 150) if not pd.isna(row["CLCR"]) else pCLCR
        pCLCR = CLCR

        TVCL = TH[0] * CLCR / 100
        CL = TVCL * np.exp(ETA[0])
        K10 = CL / V1

        AlpBe = K10 + K12 + K21
        AlmBe = K10 * K21
        Det4 = np.sqrt(AlpBe**2 - 4 * AlmBe)
        Alpha = (AlpBe + Det4) / 2
        Beta = (AlpBe - Det4) / 2
        Divisor = V1 * (Alpha - Beta)

        cTime = row["TIME"]
        dTime = cTime - prev["TIME"]

        if not pd.isna(prev["AMT"]) and prev["AMT"] > 0:
            pTime = prev["TIME"]
            pAMT = prev["AMT"]
            pRate = prev["RATE"]
            pDur = pAMT / pRate

            if dTime <= pDur:
                infC1 = pRate * (Alpha - K21) / Divisor * (1 - np.exp(-Alpha * dTime)) / Alpha
                infC2 = pRate * (K21 - Beta) / Divisor * (1 - np.exp(-Beta * dTime)) / Beta
                C1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * dTime) + infC1
                C2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * dTime) + infC2
                inf = True
            else:
                eC1 = pRate * (Alpha - K21) / Divisor * (1 - np.exp(-Alpha * pDur)) / Alpha
                eC2 = pRate * (K21 - Beta) / Divisor * (1 - np.exp(-Beta * pDur)) / Beta
                pC1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * pDur)
                pC2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * pDur)
                C1 = (eC1 + pC1) * np.exp(-Alpha * (dTime - pDur))
                C2 = (eC2 + pC2) * np.exp(-Beta * (dTime - pDur))
                inf = False
        elif inf:
            # pTime = prev["TIME"]
            # pAMT = prev["AMT"]
            # pRate = prev["RATE"]
            # pDur = pAMT / pRate
            if cTime <= pTime + pDur:
                infC1 = pRate * (Alpha - K21) / Divisor * (1 - np.exp(-Alpha * dTime)) / Alpha
                infC2 = pRate * (K21 - Beta) / Divisor * (1 - np.exp(-Beta * dTime)) / Beta
                C1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * dTime) + infC1
                C2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * dTime) + infC2
                inf = True
            else:
                rDur = pTime + pDur - prev["TIME"]
                eC1 = pRate * (Alpha - K21) / Divisor * (1 - np.exp(-Alpha * rDur)) / Alpha
                eC2 = pRate * (K21 - Beta) / Divisor * (1 - np.exp(-Beta * rDur)) / Beta
                pC1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * rDur)
                pC2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * rDur)
                C1 = (eC1 + pC1) * np.exp(-Alpha * (cTime - pTime - pDur))
                C2 = (eC2 + pC2) * np.exp(-Beta * (cTime - pTime - pDur))
                inf = False
        else:
            C1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * dTime)
            C2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * dTime)
            inf = False

        Conc = pd.concat([Conc, pd.DataFrame([[C1, C2]], columns=["C1", "C2"])], ignore_index=True)
        IPRE.append(C1 + C2)

    return np.array(IPRE)

def auto_calcTDM(PRED, DATAi, TH, SG, rEBE, TIME, AMT, RATE, II=24, ADDL=10, npoints=500, pks_params={'peak_sampling_time':1.0,}):  ########## 확인필요

    """
    PRED = PredVanco
    TIME=50
    AMT=1000
    RATE=1000
    II=12
    ADDL=10
    npoints=500
    """

    """
    Add future dosing events to DATAi and calculate prediction interval (PI).

    Parameters:
        PRED: function for PK prediction (e.g., PredVanco)
        DATAi: input DataFrame with subject data
        TH: typical parameters
        SG: residual error covariance matrix
        rEBE: output from EBE(), including EBEi and COV
        TIME: time of next dose
        AMT: dose amount
        RATE: infusion rate
        II: interdose interval
        ADDL: number of additional doses
        npoints: number of time points for interpolation

    Returns:
        DataFrame with time, observed DV, model predictions and intervals
    """

    """
    # (1) T_half 구하기 위한 함수
    # (2) 현 용법시 TDM result dict 반환함수 (Input: T_half + 위의 모든 input / Output: TDM 결과 params)
    # (3) 변경용법시 TDM result dict 반환함수 (Input: T_half + 위의 모든 input / Output: TDM 결과 params)
    # (2), (3) 은 같은 함수로 구현
    """

    # Deside time_point for steady-state
    DATAi_forss = DATAi.sort_values(['ID', 'TIME'])
    ld_time = DATAi_forss.iloc[-1]['TIME']
    ld_CLcr = DATAi_forss['CLCR'].dropna().iloc[-1]

    ETA = rEBE['EBEi']

    # 입력값
    V1 = TH[1] * np.exp(ETA[1])
    V2 = TH[2] * np.exp(ETA[2])
    Q = TH[3] * np.exp(ETA[3])
    CL = (TH[0] * ld_CLcr / 100) * np.exp(ETA[0])

    # 구획 간 이동속도
    K10 = CL / V1
    K12 = Q / V1
    K21 = Q / V2

    # Eigenvalue 계산
    A = K10 + K12 + K21
    B = K10 * K21
    lambda1 = 0.5 * (A + np.sqrt(A ** 2 - 4 * B))
    lambda2 = 0.5 * (A - np.sqrt(A ** 2 - 4 * B))

    # Terminal half-life
    t_half_terminal = np.log(2) / min(lambda1, lambda2)

    ## Interpretation (현용법 / 변경용법의 therapeutic target 도달 적정성 판단 (적절한 Regimen 추천 알고리즘 추가필요) ###########

    if PRED==PredVanco:
        cur_AMT = DATAi['AMT'].dropna().iloc[-1]
        cur_II = DATAi[DATAi['MDV']==1]['TIME'].diff().iloc[-1]
        cur_AUC24 = (cur_AMT * 24 / cur_II) / CL

        if (cur_AUC24 < 410) or (cur_AUC24 > 610):
            Regimen_changed = {'STATE':True}

            target_AUC24 = 500
            regimen_cand_df = list()
            for II_cand in [6, 8, 12, 24, 48]:
                for AMT_cand in list(range(50,1200,10)):
                    regimen_cand_df.append({"DOSE": AMT_cand, "INTERVAL": II_cand, "AUC24": (AMT_cand * 24 / II_cand) / CL})
            regimen_cand_df = pd.DataFrame(regimen_cand_df)
            regimen_cand_df = regimen_cand_df[(regimen_cand_df['AUC24'] > 430)&(regimen_cand_df['AUC24'] < 570)].copy()
            regimen_cand_df = regimen_cand_df[(regimen_cand_df['DOSE'] < 1000)].copy()
            regimen_cand_df = regimen_cand_df[(regimen_cand_df['INTERVAL'].isin([8,12,24]))].copy()
            regimen_cand_df = regimen_cand_df.sort_values(['INTERVAL', 'DOSE'],ascending=[False,True], ignore_index=True)

            Regimen_changed.update(dict(regimen_cand_df[regimen_cand_df['AUC24']>=regimen_cand_df['AUC24'].median()].iloc[0]))
            AMT = Regimen_changed['DOSE']
            II = Regimen_changed['INTERVAL']
        else:
            Regimen_changed = {'STATE': False, 'DOSE': np.nan, 'INTERVAL': np.nan, 'AUC24': np.nan}

    ## [현용법 유지시] 추가 Dosing 반영하여 Data augmentation



    ## [변경용법 사용시] 추가 Dosing 반영하여 Data augmentation

    # 입력한 추가시간 (interval x 추가 dose)가 steady-state를 나타내기에 부족할 때 -> ADDL 더 추가해서 그리도록
    if (5 * t_half_terminal) / (ADDL * II) < 0.4:
        ADDL_aug = int(round((5*t_half_terminal)/(0.4 * II), -1))
    else:
        ADDL_aug = ADDL

    # Add future dosing to the dataset
    DATAi_augmented = addDATAi(DATAi.copy(), TIME, AMT, RATE, II, ADDL_aug)

    # Calculate prediction interval after adding the new doses
    rTab = calcPI(PRED, DATAi_augmented, TH, SG, rEBE, npoints)

    # Steady-state
    ss_time = ld_time + 5 * t_half_terminal
    ss_dosing = DATAi_augmented[DATAi_augmented['TIME'] >= ss_time].copy()
    ss_rTab = rTab[rTab['x'] >= ss_time].copy()

    ## TDM 주요 결과 반환

    # Cmax_ss, Cavg_ss, Cmin_ss 계산
    # peak_sampling_time=0
    peak_sampling_time = pks_params['peak_sampling_time']
    dose_times = ss_dosing["TIME"].sort_values().values

    # 결과 저장용 리스트
    results = []

    # 각 투여 간격 구간에서 Cmax, Cmin 계산
    for i in range(len(dose_times) - 1):
        start_time = dose_times[i]
        end_time = dose_times[i + 1]

        # 해당 구간의 ss_rTab rows 추출
        if (start_time + peak_sampling_time) > end_time:
            raise ValueError("SS 구간 dosing 후 peak_sampling_time이 다음 dosing 시간을 넘습니다.")
        mask = (ss_rTab["x"] >= start_time+peak_sampling_time) & (ss_rTab["x"] < end_time)
        subset = ss_rTab.loc[mask, "y2"]

        # NaN 제외하고 Cmax, Cmin 계산
        subset_clean = subset.dropna()
        if not subset_clean.empty:
            cmax = subset_clean.max()
            cmin = subset_clean.min()
        else:
            cmax = None
            cmin = None

        results.append({
            "interval_start": start_time,
            "peak_sampling": start_time+peak_sampling_time,
            "interval_end": end_time,
            "Cmax": cmax,
            "Cmin": cmin
        })

    cmax_cmin_df = pd.DataFrame(results)
    # from ace_tools import display_dataframe_to_user
    # display_dataframe_to_user("Steady-State Cmax/Cmin by Dosing Interval", cmax_cmin_df)

    peak_conc = cmax_cmin_df['Cmax'].median()
    trough_conc = cmax_cmin_df['Cmin'].median()

    pk_params_dict = {'Vd_ss':V1+V2, 'Total_CL':CL, 'T_half_terminal':t_half_terminal,}
    cur_regimen_dict = dict()
    fut_regimen_dict = dict()

    tdm_summary_dict = {'Dose':AMT,'Interval':II, 'Cmax_ss':peak_conc, 'Cmin_ss':trough_conc, 'AUC24':AUC24}

    return rTab


def calcTDM(PRED, DATAi, TH, SG, rEBE, TIME, AMT, RATE, II, ADDL, npoints=500, pks_params={'peak_sampling_time':1.0,}):  ########## 확인필요

    """
    PRED = PredVanco
    TIME=50
    AMT=1000
    RATE=1000
    II=12
    ADDL=10
    npoints=500
    """

    """
    Add future dosing events to DATAi and calculate prediction interval (PI).

    Parameters:
        PRED: function for PK prediction (e.g., PredVanco)
        DATAi: input DataFrame with subject data
        TH: typical parameters
        SG: residual error covariance matrix
        rEBE: output from EBE(), including EBEi and COV
        TIME: time of next dose
        AMT: dose amount
        RATE: infusion rate
        II: interdose interval
        ADDL: number of additional doses
        npoints: number of time points for interpolation

    Returns:
        DataFrame with time, observed DV, model predictions and intervals
    """

    # # Deside time_point for steady-state
    # DATAi_forss = DATAi.sort_values(['ID', 'TIME'])
    # ld_time = DATAi_forss.iloc[-1]['TIME']
    # ld_CLcr = DATAi_forss['CLCR'].dropna().iloc[-1]
    #
    # ETA = rEBE['EBEi']
    #
    # # 입력값
    # V1 = TH[1] * np.exp(ETA[1])
    # V2 = TH[2] * np.exp(ETA[2])
    # Q = TH[3] * np.exp(ETA[3])
    # CL = (TH[0] * ld_CLcr / 100) * np.exp(ETA[0])
    #
    # # 구획 간 이동속도
    # K10 = CL / V1
    # K12 = Q / V1
    # K21 = Q / V2
    #
    # # Eigenvalue 계산
    # A = K10 + K12 + K21
    # B = K10 * K21
    # lambda1 = 0.5 * (A + np.sqrt(A ** 2 - 4 * B))
    # lambda2 = 0.5 * (A - np.sqrt(A ** 2 - 4 * B))
    #
    # # Terminal half-life
    # t_half_terminal = np.log(2) / min(lambda1, lambda2)
    #
    # ## Interpretation (현용법 / 변경용법의 therapeutic target 도달 적정성 판단 (적절한 Regimen 추천 알고리즘 추가필요) ###########
    #
    # cur_AMT = DATAi['AMT'].dropna().iloc[-1]
    # cur_II = DATAi[DATAi['MDV']==1]['TIME'].diff().iloc[-1]
    # cur_AUC24 = (cur_AMT * 24 / cur_II) / CL
    #
    # ## [현용법 유지시] 추가 Dosing 반영하여 Data augmentation
    #
    #
    #
    # ## [변경용법 사용시] 추가 Dosing 반영하여 Data augmentation
    #
    # # 입력한 추가시간 (interval x 추가 dose)가 steady-state를 나타내기에 부족할 때 -> ADDL 더 추가해서 그리도록
    # if (5 * t_half_terminal) / (ADDL * II) < 0.4:
    #     ADDL_aug = int(round((5*t_half_terminal)/(0.4 * II), -1))
    # else:
    #     ADDL_aug = ADDL

    # Add future dosing to the dataset
    DATAi_augmented = addDATAi(DATAi.copy(), TIME, AMT, RATE, II, ADDL)
    # DATAi_augmented = addDATAi(DATAi.copy(), TIME, AMT, RATE, II, ADDL_aug)

    # Calculate prediction interval after adding the new doses
    rTab = calcPI(PRED, DATAi_augmented, TH, SG, rEBE, npoints)
    #
    # # Steady-state
    # ss_time = ld_time + 5 * t_half_terminal
    # ss_dosing = DATAi_augmented[DATAi_augmented['TIME'] >= ss_time].copy()
    # ss_rTab = rTab[rTab['x'] >= ss_time].copy()
    #
    # ## TDM 주요 결과 반환
    #
    # # Cmax_ss, Cavg_ss, Cmin_ss 계산
    # # peak_sampling_time=0
    # peak_sampling_time = pks_params['peak_sampling_time']
    # dose_times = ss_dosing["TIME"].sort_values().values
    #
    # # 결과 저장용 리스트
    # results = []
    #
    # # 각 투여 간격 구간에서 Cmax, Cmin 계산
    # for i in range(len(dose_times) - 1):
    #     start_time = dose_times[i]
    #     end_time = dose_times[i + 1]
    #
    #     # 해당 구간의 ss_rTab rows 추출
    #     if (start_time + peak_sampling_time) > end_time:
    #         raise ValueError("SS 구간 dosing 후 peak_sampling_time이 다음 dosing 시간을 넘습니다.")
    #     mask = (ss_rTab["x"] >= start_time+peak_sampling_time) & (ss_rTab["x"] < end_time)
    #     subset = ss_rTab.loc[mask, "y2"]
    #
    #     # NaN 제외하고 Cmax, Cmin 계산
    #     subset_clean = subset.dropna()
    #     if not subset_clean.empty:
    #         cmax = subset_clean.max()
    #         cmin = subset_clean.min()
    #     else:
    #         cmax = None
    #         cmin = None
    #
    #     results.append({
    #         "interval_start": start_time,
    #         "peak_sampling": start_time+peak_sampling_time,
    #         "interval_end": end_time,
    #         "Cmax": cmax,
    #         "Cmin": cmin
    #     })
    #
    # cmax_cmin_df = pd.DataFrame(results)
    # # from ace_tools import display_dataframe_to_user
    # # display_dataframe_to_user("Steady-State Cmax/Cmin by Dosing Interval", cmax_cmin_df)
    #
    # peak_conc = cmax_cmin_df['Cmax'].median()
    # trough_conc = cmax_cmin_df['Cmin'].median()
    #
    # pk_params_dict = {'Vd_ss':V1+V2, 'Total_CL':CL, 'T_half_terminal':t_half_terminal,}
    # cur_regimen_dict = dict()
    # fut_regimen_dict = dict()
    #
    # tdm_summary_dict = {'Dose':AMT,'Interval':II, 'Cmax_ss':peak_conc, 'Cmin_ss':trough_conc, 'AUC24':AUC24}

    return rTab
