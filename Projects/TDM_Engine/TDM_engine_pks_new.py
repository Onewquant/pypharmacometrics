from Projects.TDM_Engine.tools import *
import matplotlib.pyplot as plt
import seaborn as sns


WT = 50.1
HT = 165
SCR = 0.56
AGE = 65
SEX = 1 # Female
SEX_Coef = 0.85 if SEX==1 else 1
BMI = WT/((HT/100)**2)
# LBW = WT * 9270 / (8780+244 * BMI) if SEX==1 else WT * 9270 / (6680+216 * BMI)
CLCR = (SEX_Coef * (140-AGE) * WT / (72*SCR))


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

# TH = np.array([
#     TVCL,
#     TVVc,
#     TVVp,
#     TVQ
# ])

# Inter-individual variability matrix (OMEGA)
# OM = np.array([
#     [ 0.0392,   0.0,  0.0,   0.0],
#     [0.0,    0.0392, 0.0, 0.0],
#     [0.0,      0.0,   0.1604,   0.0],
#     [0.0,     0.0,  0.0,   0.0998]
# ])

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


# Residual variability matrix (SIGMA)
# SG = np.array([
#     [0.01, 0.0],
#     [0.0,0.0001]
# ])

CV_assay = 0.15
S_assay = 0.25

SG = np.array([
    [CV_assay, 0.0],
    [0.0, S_assay]
])


# 1. CSV 파일 읽기
raw_data = pd.read_csv(
    "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/TDM_Engine/sample_vanco_loading_maintenance.csv",
    na_values=["", ".", "NA"]
)

raw_data['CLCR'] = CLCR

# 2. TIME 처리
data_prepped = convDT(raw_data)

# 3. ADDL 값이 유효하면 expand
if "ADDL" in data_prepped.columns and data_prepped["ADDL"].notna().any() and (data_prepped["ADDL"] > 0).any():
    data_prepped = expandDATA(data_prepped)

# 4. ID, TIME 기준 정렬
data_prepped = data_prepped.sort_values(by=["ID", "TIME"]).reset_index(drop=True)

# 5. 필요한 열만 추출
DATAi = data_prepped[["ID", "TIME", "AMT", "RATE", "DV", "MDV", "SEX", "AGE", "BWT", "SCR", "CLCR"]].copy()

# 6. Conc 측정 row의 AMT 및 Rate는 반드시 비어있어야함! (추정 잘못하게 됨)
obs_mask = DATAi["DV"].notna() | (DATAi["MDV"] == 0)
DATAi.loc[obs_mask, ["AMT", "RATE"]] = np.nan

# 7. EBE 추정
rEBE = EBEpks(PredVanco_PKS, DATAi, TH, OM, SG)

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
