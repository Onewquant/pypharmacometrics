import numpy as np
import pandas as pd
from scipy.optimize import minimize, least_squares
from numpy.linalg import inv

######## 기본 수치 계산 ######## 

def crcl_cockroft_gault(sex, age, wt, scr):
    sex_coef = 0.85 if sex in [1, 'F'] else 1
    crcl = (sex_coef * (140-age) * wt / (72*scr))
    return crcl

######## EBE 수행 관련 ######## 

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
