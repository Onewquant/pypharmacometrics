from Projects.TDM_Engine.engine_codes.engine_core_tools import *

def _p_to_eta(P, TH):
    pop = np.array([TH[0], TH[1], TH[2], TH[4], TH[5]])
    return np.log(P / pop)


def _get_pop_params(TH):
    return np.array([
        TH[0],  # CLcoef
        TH[1],  # CLslope
        TH[2],  # Vccoef
        TH[4],  # k12
        TH[5],  # k21
    ])


def _get_cv_from_omega(OM):
    return np.sqrt(np.exp(np.diag(OM)) - 1)


def residualP(P):
    # negative parameter 방지
    if np.any(P <= 0):
        return np.ones(100) * 1e6

    ETA = _p_to_eta(P, e.TH)

    obs_mask = ~e.DATAi["DV"].isna()
    DV = e.DATAi.loc[obs_mask, "DV"].values

    pred = e.PRED(e.TH, ETA, e.DATAi)
    Fi = pred[obs_mask]

    CV_assay = e.SG[0, 0]
    S_assay = e.SG[1, 1]

    sigma = CV_assay * Fi + S_assay
    sigma = np.maximum(sigma, 1e-12)

    obs_resid = (DV - Fi) / sigma

    P_pop = _get_pop_params(e.TH)
    P_cv = _get_cv_from_omega(e.OM)

    param_sigma = P_pop * P_cv

    # param_sigma = np.array([0.01, 0.248, 0.042, 0.28, 0.12])
    # param_sigma = np.maximum(param_sigma, 1e-12)
    # param_resid = (P - P_pop) / param_sigma

    param_sigma = np.maximum(param_sigma, 1e-12)
    param_resid = (P - P_pop) / param_sigma

    return np.concatenate([obs_resid, param_resid])


def PredVanco_PKS(TH, ETA, DATAi):
    """
    TH:
      0 CL_coef        mL/min/kg
      1 CLslope        dimensionless renal slope
      2 Vc_coef        L/kg
      3 Vc_slope       L per mL/min CrCL, usually 0
      4 k12            1/hr
      5 k21            1/hr

    ETA:
      0 ETA_CLcoef
      1 ETA_CLslope
      2 ETA_Vc
      3 ETA_k12
      4 ETA_k21
    """

    CL_coef_pop = TH[0]
    CLslope_pop = TH[1]
    Vc_coef_pop = TH[2]
    Vc_slope = TH[3]
    k12_pop = TH[4]
    k21_pop = TH[5]

    CL_coef_i = CL_coef_pop * np.exp(ETA[0])
    CLslope_i = CLslope_pop * np.exp(ETA[1])
    Vc_coef_i = Vc_coef_pop * np.exp(ETA[2])
    k12_i = k12_pop * np.exp(ETA[3])
    k21_i = k21_pop * np.exp(ETA[4])

    Conc = pd.DataFrame([[0.0, 0.0]], columns=["C1", "C2"])
    IPRE = [0.0]
    inf = False

    pCLCR = min(DATAi.iloc[0]["CLCR"], 150)

    for i in range(1, len(DATAi)):
        row = DATAi.iloc[i]
        prev = DATAi.iloc[i - 1]

        CLCR = min(row["CLCR"], 150) if not pd.isna(row["CLCR"]) else pCLCR
        pCLCR = CLCR

        WT = row["BWT"]

        CL = (CL_coef_i * WT + CLslope_i * CLCR) * 0.06
        V1 = Vc_coef_i * WT + Vc_slope * CLCR

        K10 = CL / V1
        K12 = k12_i
        K21 = k21_i

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

        Conc = pd.concat(
            [Conc, pd.DataFrame([[C1, C2]], columns=["C1", "C2"])],
            ignore_index=True
        )
        IPRE.append(C1 + C2)

    return np.array(IPRE)

def pks_lm_from_residual(
    x0,
    residual_func,
    max_iter=25,
    conv=1e-3,
    lam0=1.0,
    lam_up=10.0,
    lam_down=0.1,
    fd_eps=1e-4
):
    x = np.asarray(x0, dtype=float).copy()
    lam = lam0
    hist = []

    for it in range(1, max_iter + 1):
        r = residual_func(x)
        obj = np.sum(r ** 2)

        n = len(x)
        m = len(r)
        J = np.zeros((m, n))

        for j in range(n):
            h = fd_eps * max(abs(x[j]), 1.0)
            xp = x.copy()
            xm = x.copy()
            xp[j] += h
            xm[j] -= h
            J[:, j] = (residual_func(xp) - residual_func(xm)) / (2 * h)

        A = J.T @ J
        g = J.T @ r

        accepted = False

        for _ in range(20):
            A_mod = A.copy()
            for j in range(n):
                A_mod[j, j] = A[j, j] * (1 + lam)

            try:
                step = np.linalg.solve(A_mod, -g)
            except np.linalg.LinAlgError:
                lam *= lam_up
                continue

            x_next = x + step

            if np.any(x_next <= 0):
                lam *= lam_up
                continue

            r_next = residual_func(x_next)
            obj_next = np.sum(r_next ** 2)

            if obj_next < obj:
                p_delta = np.abs((x - x_next) / x_next)

                hist.append({
                    "iter": it,
                    "obj": obj_next,
                    "lambda": lam,
                    "P": x_next.copy(),
                    "p_delta": p_delta,
                    "max_delta": np.max(p_delta),
                    "accepted": True
                })

                x = x_next
                lam *= lam_down
                accepted = True

                if np.all(p_delta < conv):
                    return x, hist

                break

            lam *= lam_up

        if not accepted:
            hist.append({
                "iter": it,
                "obj": obj,
                "lambda": lam,
                "P": x.copy(),
                "p_delta": np.full(n, np.nan),
                "max_delta": np.nan,
                "accepted": False
            })

    return x, hist

def EBE(PRED, DATAi, TH, OM, SG, x0=None):

    e.PRED = PRED
    e.DATAi = DATAi.copy()
    e.TH = TH
    e.OM = OM
    e.SG = SG

    P_pop = _get_pop_params(TH)

    # if x0 is None:
    #     x0 = P_pop.copy()
    # else:
    #     # x0가 ETA로 들어오면 P로 변환
    #     x0 = P_pop * np.exp(x0)

    if x0 is None:
        x0 = P_pop.copy()
    else:
        x0 = np.asarray(x0, dtype=float)

    # res = least_squares(
    #     residualP,
    #     x0=x0,
    #     # x_scale=x0,
    #     method="lm",
    #     max_nfev=25,
    #     diff_step=1e-4,
    #     xtol=1e-3,
    #     ftol=1e-3,
    #     gtol=1e-3
    # )
    #
    # P_ind = res.x

    P_ind, hist = pks_lm_from_residual(
        x0=x0,
        residual_func=residualP,
        max_iter=25,
        conv=1e-3,
        lam0=1.0,
        lam_up=10.0,
        lam_down=0.1,
        fd_eps=1e-4
    )


    EBEi = _p_to_eta(P_ind, TH)

    Fi = PRED(TH, EBEi, DATAi)

    obs_mask = ~DATAi["DV"].isna()
    DV = DATAi.loc[obs_mask, "DV"].values
    IPRED_obs = Fi[obs_mask]

    Ri = DV - IPRED_obs

    CV_assay = SG[0, 0]
    S_assay = SG[1, 1]

    sigma = CV_assay * IPRED_obs + S_assay
    sigma = np.maximum(sigma, 1e-12)

    WSS = np.sum((Ri / sigma) ** 2)

    P_cv = _get_cv_from_omega(OM)
    param_sigma = P_pop * P_cv
    # param_sigma = np.array([0.01, 0.248, 0.042, 0.28, 0.12])
    PARAM_PENALTY = np.sum(((P_ind - P_pop) / param_sigma) ** 2)

    BAYES_SS = WSS + PARAM_PENALTY

    # try:
    #     J = res.jac
    #     COV_P = np.linalg.inv(J.T @ J)
    #     SE_P = np.sqrt(np.diag(COV_P))
    # except np.linalg.LinAlgError:
    #     COV_P = np.full((len(P_ind), len(P_ind)), np.nan)
    #     SE_P = np.full(len(P_ind), np.nan)

    return {
        "EBEi": EBEi,
        "P_ind": P_ind,
        # "SE_P": SE_P,
        # "COV_P": COV_P,
        "IPRED": Fi,
        "IRES": Ri,

        "ITERATION": len([h for h in hist if h["accepted"]]),
        "HISTORY": hist,
        "WSS": WSS,
        "PARAM_PENALTY": PARAM_PENALTY,
        "BAYES_SS": BAYES_SS,
        "SUM_OF_SQUARES": BAYES_SS,

        # "ITERATION": res.nfev,
        # "SUCCESS": res.success,
        # "MESSAGE": res.message
    }

###################################
## Parameter에 직접 penalty 가하는 방법이라고함 (이걸로 쓰려면 EBE 함수도 다시 구현해야)
# def obj_P(P):
#     if np.any(P <= 0):
#         return 1e12
#
#     P_pop = _get_pop_params(e.TH)
#     ETA = np.log(P / P_pop)
#
#     pred = PredVanco_PKS(e.TH, ETA, e.DATAi)
#
#     obs = e.DATAi["DV"].notna()
#     DV = e.DATAi.loc[obs, "DV"].values
#     Fi = pred[obs]
#
#     CV_assay = e.SG[0, 0]
#     S_assay = e.SG[1, 1]
#
#     sigma = CV_assay * Fi + S_assay
#
#     WSS = np.sum(((DV - Fi) / sigma) ** 2)
#
#     P_cv = _get_cv_from_omega(e.OM)
#     param_sigma = P_pop * P_cv
#
#     PARAM = np.sum(((P - P_pop) / param_sigma) ** 2)
#
#     return WSS + PARAM
#
# def numerical_grad_hess(func, x, eps=1e-4):
#     x = np.asarray(x, dtype=float)
#     n = len(x)
#
#     g = np.zeros(n)
#     H = np.zeros((n, n))
#
#     f0 = func(x)
#
#     for i in range(n):
#         hi = eps * max(abs(x[i]), 1.0)
#
#         xp = x.copy()
#         xm = x.copy()
#         xp[i] += hi
#         xm[i] -= hi
#
#         fp = func(xp)
#         fm = func(xm)
#
#         g[i] = (fp - fm) / (2 * hi)
#         H[i, i] = (fp - 2 * f0 + fm) / (hi ** 2)
#
#         for j in range(i + 1, n):
#             hj = eps * max(abs(x[j]), 1.0)
#
#             xpp = x.copy()
#             xpm = x.copy()
#             xmp = x.copy()
#             xmm = x.copy()
#
#             xpp[i] += hi
#             xpp[j] += hj
#
#             xpm[i] += hi
#             xpm[j] -= hj
#
#             xmp[i] -= hi
#             xmp[j] += hj
#
#             xmm[i] -= hi
#             xmm[j] -= hj
#
#             H_ij = (
#                 func(xpp)
#                 - func(xpm)
#                 - func(xmp)
#                 + func(xmm)
#             ) / (4 * hi * hj)
#
#             H[i, j] = H_ij
#             H[j, i] = H_ij
#
#     return g, H
#
# def pks_marquardt_levenberg(
#     func,
#     x0,
#     max_iter=25,
#     conv=1e-3,
#     lambda0=1.0,
#     lambda_up=10.0,
#     lambda_down=0.1,
#     eps=1e-4
# ):
#     x = np.asarray(x0, dtype=float).copy()
#     lam = lambda0
#     hist = []
#
#     for it in range(1, max_iter + 1):
#         f_cur = func(x)
#         g, H = numerical_grad_hess(func, x, eps=eps)
#
#         accepted = False
#
#         for _ in range(20):
#             H_mod = H.copy()
#
#             # PKS manual: H_jj = H_jj * (lambda + 1)
#             for j in range(len(x)):
#                 H_mod[j, j] = H[j, j] * (lam + 1)
#
#             try:
#                 step = np.linalg.solve(H_mod, -g)
#             except np.linalg.LinAlgError:
#                 lam *= lambda_up
#                 continue
#
#             x_next = x + step
#
#             if np.any(x_next <= 0):
#                 lam *= lambda_up
#                 continue
#
#             f_next = func(x_next)
#
#             if f_next < f_cur:
#                 p_delta = np.abs((x - x_next) / x_next)
#                 accepted = True
#
#                 hist.append({
#                     "iter": it,
#                     "obj": f_next,
#                     "lambda": lam,
#                     "P": x_next.copy(),
#                     "p_delta": p_delta,
#                     "max_delta": np.max(p_delta)
#                 })
#
#                 x = x_next
#                 lam *= lambda_down
#
#                 if np.all(p_delta < conv):
#                     return x, hist
#
#                 break
#
#             else:
#                 lam *= lambda_up
#
#         if not accepted:
#             hist.append({
#                 "iter": it,
#                 "obj": f_cur,
#                 "lambda": lam,
#                 "P": x.copy(),
#                 "p_delta": np.full(len(x), np.nan),
#                 "max_delta": np.nan
#             })
#
#     return x, hist


## INPUT DATA (AGE, SEX, WT, HT, SCR)

# AGE = 65
# SEX = 1            # 0: Male, 1: Female
# HT = 165
# WT = 50.1
# SCR = 0.56

# 1. CSV 파일 읽기
file_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/TDM_Engine/sample_vanco_loading_maintenance.csv"
# file_path = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/TDM_Engine/pks_vanco_test_data1.csv"
raw_data = pd.read_csv(file_path, na_values=["", ".", "NA"])

# 2. TIME 처리
data_prepped = convDT(raw_data)

# 3. ADDL 값이 유효하면 expand
if "ADDL" in data_prepped.columns and data_prepped["ADDL"].notna().any() and (data_prepped["ADDL"] > 0).any():
    data_prepped = expandDATA(data_prepped)

# 4. ID, TIME 기준 정렬
data_prepped = data_prepped.sort_values(by=["ID", "TIME"]).reset_index(drop=True)

# 5. 기본 수치 계산
# data_prepped['SEX'] = SEX   # 임시 (코드 수정 가능성)
# data_prepped['HT'] = HT   # 임시 (코드 수정 가능성)
# data_prepped['WT'] = WT   # 임시 (코드 수정 가능성)
# data_prepped['SCR'] = SCR # 임시 (코드 수정 가능성)
# CLCR = 144.82
# data_prepped['CLCR'] = CLCR
data_prepped['CLCR'] = data_prepped.apply(lambda row:crcl_cockroft_gault(sex=row['SEX'],age=row['AGE'],wt=row['BWT'],scr=row['SCR']),axis=1)
CLCR = data_prepped['CLCR'].iloc[-1] # 임시 (코드 수정 가능성)
WT = data_prepped['BWT'].iloc[-1] # 임시 (코드 수정 가능성)

# 5. 필요한 열만 추출
DATAi = data_prepped[["ID", "TIME", "AMT", "RATE", "DV", "MDV", "SEX", "AGE", "BWT", "SCR", "CLCR"]].copy()

# 6. Conc 측정 row의 AMT 및 Rate는 반드시 비어있어야함! (추정 잘못하게 됨)
obs_mask = DATAi["DV"].notna() | (DATAi["MDV"] == 0)
DATAi.loc[obs_mask, ["AMT", "RATE"]] = np.nan

## PRIOR VALUES

# BMI = WT/((HT/100)**2)

# PKS는 Total weight 사용 [참고: LBW = WT * 9270 / (8780+244 * BMI) if SEX==1 else WT * 9270 / (6680+216 * BMI)]

TVCLcoef = 0.05
TVCLslope = 0.75
TVVccoef = 0.21
TVVcslope = 0
TVK12 = 1.12
TVK21 = 0.48

TH = np.array([
    TVCLcoef,      # CL_coef, mL/min/kg
    TVCLslope,     # Renal CL slope
    TVVccoef,      # Vc_coef, L/kg
    TVVcslope,     # Renal Vc slope
    TVK12,         # k12, 1/hr
    TVK21          # k21, 1/hr
])

CV_CLcoef = 0.2
CV_CLslope = 0.33
CV_Vc = 0.20
CV_k12 = 0.25
CV_k21 = 0.25

OM = np.diag([
    np.log(1+CV_CLcoef**2),
    np.log(1+CV_CLslope**2),
    np.log(1+CV_Vc**2),
    np.log(1+CV_k12**2),
    np.log(1+CV_k21**2),
])

CV_assay = 0.15
S_assay = 0.25

SG = np.array([
    [CV_assay, 0.0],
    [0.0, S_assay]
])

# 7. EBE 추정
rEBE = EBE(PredVanco_PKS, DATAi, TH, OM, SG)

"""
CL_ind = TVCL*np.exp(rEBE['EBEi'][0])
Vc_ind = TVVc*np.exp(rEBE['EBEi'][1])
Vp_ind = TVVp*np.exp(rEBE['EBEi'][2])
Q_ind = TVQ*np.exp(rEBE['EBEi'][3])
k12_ind = Q_ind/Vc_ind    #k12
k21_ind = Q_ind/Vp_ind             #k21
k10_ind = CL_ind/Vc_ind
Vdss_ind = Vc_ind * (1 + k12_ind/k21_ind)

S = k10_ind + k12_ind + k21_ind
D = np.sqrt(S**2 - 4*k21_ind*k10_ind)
alpha = (S+D)/2
beta = (S-D)/2
T_half = 0.693/beta
"""

ETA = rEBE["EBEi"]

CL_coef_ind = TH[0] * np.exp(ETA[0])
CLslope_ind = TH[1] * np.exp(ETA[1])
Vc_coef_ind = TH[2] * np.exp(ETA[2])
k12_ind = TH[4] * np.exp(ETA[3])
k21_ind = TH[5] * np.exp(ETA[4])

CL_ind = (CL_coef_ind * WT + CLslope_ind * CLCR) * 0.06
Vc_ind = Vc_coef_ind * WT + TH[3] * CLCR

Q_ind = k12_ind * Vc_ind
Vp_ind = Q_ind / k21_ind

k10_ind = CL_ind / Vc_ind
Vdss_ind = Vc_ind + Vp_ind

S = k10_ind + k12_ind + k21_ind
D = np.sqrt(S**2 - 4 * k21_ind * k10_ind)

alpha = (S + D) / 2
beta = (S - D) / 2
T_half = 0.693 / beta

print("Iteration:", rEBE["ITERATION"])
print("Sum of squares:", round(rEBE["SUM_OF_SQUARES"],4))
print(f"Vc coef: {round(Vc_ind/WT,3)}")
print(f"Cl coef: {round(CL_coef_ind,4)}")
print(f"Renal Cl (slope): {round(CLslope_ind,3)}")
# print(f"Cl coef: {round(CL_ind*1000/60/WT,3)} ")
# print(f"Renal Cl (slope):  ")
print(f"K12: {round(k12_ind,3)}")
print(f"K21: {round(k21_ind,3)}")
print(f"Total Vc: {round(Vc_ind,3)}")
print(f"Total Cl: {round(CL_ind,3)}")
print(f"K10: {round(k10_ind,3)}")
print(f"Alpha: {round(alpha,3)}")
print(f"Beta: {round(beta,3)}")
print(f"Vdss: {round(Vdss_ind,3)}")
print(f"Half-life: {round(T_half,3)}")

# rEBE['WSS']
# rEBE['PARAM_PENALTY']


# cl_grid = np.linspace(0.90,1.05,40)
#
# objs = []
#
# for cls in cl_grid:
#
#     P = P_best.copy()
#     P[1] = cls
#
#     objs.append(obj_P(P))
#
# plt.plot(cl_grid,objs)
# plt.show()

# x0 = np.array([
#     0.05,  # CL coef
#     0.75,  # Renal CL slope
#     0.21,  # Vc coef
#     1.12,  # k12
#     0.48   # k21
# ])
#
# for n in [5, 8, 10, 12, 15, 19, 25]:
#     res = least_squares(
#         residualP,
#         x0=x0,
#         method="lm",
#         max_nfev=n,
#         xtol=1e-3,
#         ftol=1e-3,
#         gtol=1e-3
#     )
#
#     P = res.x
#     print(
#         n,
#         "CLcoef", round(P[0], 4),
#         "CLslope", round(P[1], 3),
#         "Vccoef", round(P[2], 3),
#         "k12", round(P[3], 3),
#         "k21", round(P[4], 3),
#     )
#
# ETA_pks = np.array([
#     np.log(0.0502 / 0.05),
#     np.log(0.978  / 0.75),
#     np.log(0.196  / 0.21),
#     np.log(1.04   / 1.12),
#     np.log(0.521  / 0.48),
# ])
#
# pred = PredVanco_PKS(TH, ETA_pks, DATAi)
# print(pd.DataFrame({
#     "TIME": DATAi.loc[DATAi["DV"].notna(), "TIME"],
#     "DV": DATAi.loc[DATAi["DV"].notna(), "DV"],
#     "IPRED": pred[DATAi["DV"].notna()]
# }))


# x0 = np.array([
#     0.05,  # CL coef
#     0.75,  # Renal CL slope
#     0.21,  # Vc coef
#     1.12,  # k12
#     0.48   # k21
# ])
#
# P_ind, hist = pks_marquardt_levenberg(
#     obj_P,
#     x0=x0,
#     max_iter=25,
#     conv=0.001,
#     lambda0=1.0,
#     eps=1e-4
# )
#
# for h in hist:
#     P = h["P"]
#     print(
#         h["iter"],
#         "OBJ", round(h["obj"], 4),
#         "lambda", round(h["lambda"], 6),
#         "max_delta", h["max_delta"],
#         "CLcoef", round(P[0], 4),
#         "CLslope", round(P[1], 3),
#         "Vccoef", round(P[2], 3),
#         "k12", round(P[3], 3),
#         "k21", round(P[4], 3),
#     )


# TVCL
#
# ## Estimation 쪽 문제인지 / Prediction 쪽 문제인지
#
# ETA_pks = np.array([
#     np.log(0.0502 / 0.05),   # CLcoef
#     np.log(0.978 / 0.75),   # CLslope
#     np.log(0.196 / 0.21),   # Vccoef
#     np.log(1.04 / 1.12),   # k12
#     np.log(0.521 / 0.48),   # k21
# ])
#
# pred = PredVanco_PKS_new(TH, ETA_pks, DATAi)
#
# obs = DATAi["DV"].notna()
#
# print(
#     pd.DataFrame({
#         "TIME": DATAi.loc[obs, "TIME"],
#         "DV": DATAi.loc[obs, "DV"],
#         "IPRED": pred[obs]
#     })
# )
#
# ## Peak 쪽 문제인지 / Trough쪽 문제인지 확인
#
# dose_rows = DATAi[DATAi["AMT"].notna()].copy()
#
# for t in [45.833333, 49.333333]:
#     prev_dose = dose_rows[dose_rows["TIME"] <= t]["TIME"].max()
#     print(
#         "obs", t,
#         "prev dose", prev_dose,
#         "delta", t - prev_dose
#     )
#

# ETA_pks = np.array([
#     np.log(0.0502 / 0.05),
#     np.log(0.978  / 0.75),
#     np.log(0.196  / 0.21),
#     np.log(1.04   / 1.12),
#     np.log(0.521  / 0.48),
# ])
#
# for dur in [0.5, 1, 1.5, 2, 3]:
#     DATAi_test = DATAi.copy()
#     dose_mask = DATAi_test["AMT"].notna() & (DATAi_test["AMT"] > 0)
#
#     # duration = AMT / RATE 이므로 RATE = AMT / duration
#     DATAi_test.loc[dose_mask, "RATE"] = DATAi_test.loc[dose_mask, "AMT"] / dur
#
#     pred = PredVanco_PKS_new(TH, ETA_pks, DATAi_test)
#     obs = DATAi_test["DV"].notna()
#
#     print("\nDUR =", dur)
#     print(pd.DataFrame({
#         "TIME": DATAi_test.loc[obs, "TIME"].values,
#         "DV": DATAi_test.loc[obs, "DV"].values,
#         "IPRED": pred[obs]
#     }))
#
#
# ETA_pks = np.array([
#     np.log(0.0502 / 0.05),
#     np.log(0.978  / 0.75),
#     np.log(0.196  / 0.21),
#     np.log(1.04   / 1.12),
#     np.log(0.521  / 0.48),
# ])
#
# obs_idx = DATAi[DATAi["DV"].notna()].index
#
# for shift in [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]:
#     DATAi_test = DATAi.copy()
#     DATAi_test.loc[obs_idx[1], "TIME"] += shift
#
#     pred = PredVanco_PKS_new(TH, ETA_pks, DATAi_test)
#
#     print(
#         shift,
#         round(pred[obs_idx[0]], 3),
#         round(pred[obs_idx[1]], 3)
#     )
#
# print(raw_data[['DATE','TIME','DV']])
#
#
#
# DATAi_test = DATAi.copy()
#
# dose_idx = DATAi_test[DATAi_test["AMT"].notna() & (DATAi_test["AMT"] > 0)].index
#
# # 첫 dose 기준 q12h로 강제
# first_dose_time = DATAi_test.loc[dose_idx[0], "TIME"]
# DATAi_test.loc[dose_idx, "TIME"] = first_dose_time + np.arange(len(dose_idx)) * 12
#
# ETA_pks = np.array([
#     np.log(0.0502 / 0.05),
#     np.log(0.978  / 0.75),
#     np.log(0.196  / 0.21),
#     np.log(1.04   / 1.12),
#     np.log(0.521  / 0.48),
# ])
#
# pred = PredVanco_PKS_new(TH, ETA_pks, DATAi_test)
#
# obs = DATAi_test["DV"].notna()
#
# print(DATAi_test[["TIME", "AMT", "RATE", "DV"]])
#
# print(pd.DataFrame({
#     "TIME": DATAi_test.loc[obs, "TIME"],
#     "DV": DATAi_test.loc[obs, "DV"],
#     "IPRED": pred[obs]
# }))


#
# # 7. 예측 인터벌 계산
# PI = calcTDM(PredVanco, DATAi, TH, SG, rEBE, TIME=50, AMT=1000, RATE=1000, II=12, ADDL=10)
#
# # 8. TDM 작성에 필수적인 결과정리
# # PI['y']
#
# # 9. 시각화
#
# # 스타일 지정 (선택)
# sns.set(style="whitegrid")
#
# # 그래프 그리기
# plt.figure(figsize=(10, 6))
#
# # 관측치
# mask_obs = ~PI["y"].isna()
# sns.scatterplot(x=PI["x"][mask_obs], y=PI["y"][mask_obs], color='black', s=40, label='Observed')
#
# # 예측치 및 신뢰구간
# sns.lineplot(x="x", y="y2", data=PI, color='dimgrey', label='Prediction (mean)', linewidth=2, alpha=0.5)
# sns.lineplot(x="x", y="yciLL", data=PI, linestyle='--', color='indianred', label='95% CI lower', alpha=0.5)
# sns.lineplot(x="x", y="yciUL", data=PI, linestyle='--', color='indianred', label='95% CI upper', alpha=0.5)
# sns.lineplot(x="x", y="ypiLL", data=PI, linestyle=':', color='royalblue', label='95% PI lower', alpha=0.5)
# sns.lineplot(x="x", y="ypiUL", data=PI, linestyle=':', color='royalblue', label='95% PI upper', alpha=0.5)
#
# # 기준선
# for y in [5, 15, 25, 35]:
#     plt.axhline(y=y, linestyle='--', color='gray', linewidth=0.8)
#
# # 축 및 레이블 설정
# plt.xlim(PI["x"].min(), PI["x"].max())
# plt.ylim([PI["ypiLL"].min(), PI["ypiUL"].max()])
# plt.xlabel("Time")
# plt.ylabel("Concentration ± 2SD")
# plt.title("Vancomycin Concentration Prediction ± 2SD")
# plt.legend()
# plt.tight_layout()
# plt.show()
# #
# # (CL_ind-CLnonrenal)/(CLCR * 0.06)
# #
# # Cl_coef = CL_nonrenal * 1000 / 60 / WT
# # Renal_CL_slope = (Total_CL - CL_nonrenal) / (CLCR * 0.06)
# #
