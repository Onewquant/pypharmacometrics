import numpy as np
import pandas as pd
from scipy.optimize import minimize
from numpy.linalg import inv
from datetime import datetime
import json

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


def objEta(ETAi, e):

    DV = e.DATAi.loc[~e.DATAi['DV'].isna(), 'DV'].values
    pred = e.PRED(e.TH, ETAi, e.DATAi)
    e.Fi = pred[~e.DATAi['DV'].isna()]
    e.Ri = DV - e.Fi
    e.Hi = np.column_stack((e.Fi, np.ones_like(e.Fi)))
    e.Vi = np.diag(e.Hi @ e.SG @ e.Hi.T)
    penalty = ETAi @ e.invOM @ ETAi
    return np.sum(np.log(e.Vi) + (e.Ri ** 2 / e.Vi)) + penalty


def EBE(PRED, DATAi, TH, OM, SG, e):

    # PRED = dmodel.Pred
    # TH = dmodel.TH
    # OM = dmodel.OM
    # SG = dmodel.SG

    e.PRED = PRED
    e.DATAi = DATAi.copy()
    e.TH = TH
    e.OM = OM
    e.SG = SG
    e.nEta = OM.shape[0]
    e.invOM = inv(OM)

    res = minimize(objEta, x0=np.zeros(e.nEta), args=(e,), method="BFGS", options={'disp': False})  ########## 요라인만 약간 다른값
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

def convDT(DATAi):
    # 'DATE' 컬럼이 있는 경우만 변환
    if "DATE" in DATAi.columns:
        DATAi["SDT"] = pd.to_datetime(DATAi["DATE"].astype(str) + " " + DATAi["TIME"].astype(str))
        FDT = DATAi["SDT"].iloc[0]
        DATAi["TIME"] = (DATAi["SDT"] - FDT).dt.total_seconds() / 3600
        DATAi = DATAi.drop(columns=["DATE", "SDT"])
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

    CovCols = [col for col in eDATAi.columns if col not in ["UID", "TIME", "AMT", "RATE", "II", "ADDL", "DV", "MDV"]]
    for i in range(1, len(eDATAi)):
        if i < len(Added_flags) and Added_flags[i]:
            eDATAi.loc[i, CovCols] = eDATAi.loc[i - 1, CovCols]
    return eDATAi