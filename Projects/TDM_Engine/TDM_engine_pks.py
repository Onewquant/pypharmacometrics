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

## Typical parameters (THETA): CL, V1, V2, Q

# TVCL 구하기
NE=1
CLslope = 0.75

CLnonrenal = 0.05 * WT
CLrenal = CLslope * CLCR * 0.06
CL70kg = CLnonrenal + CLrenal
TVCL = CL70kg * ((WT/70)**NE)

# CL_ind * 1000 / 60 / WT

# TVCL / WT
# ((CL_ind / ((WT/70)**NE)) - 0.05 * WT)/ (CLCR * 0.06)
# CL_ind/WT

# TVVc 구하기
Vcnr = 0.21 * WT
Vc_slope = 0
TVVc = Vcnr + Vc_slope * CLCR

# TVVp 구하기
k12 = 1.12
k21 = 0.48
TVVp = TVVc * k12 / k21

# TVQ 구하기
TVQ = TVVc * k12

TH = np.array([
    TVCL,
    TVVc,
    TVVp,
    TVQ
])

# Inter-individual variability matrix (OMEGA)
OM = np.array([
    [0.0392,   0.0,  0.0,   0.0],
    [0.0,    0.0392, 0.0, 0.0],
    [0.0,      0.0,   0.1604,   0.0],
    [0.0,     0.0,  0.0,   0.0998]
])

# Residual variability matrix (SIGMA)
SG = np.array([
    [0.01, 0.0],
    [0.0,0.0001]
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
# DATAi.replace(np.nan,'.').to_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/TDM_Engine/vanco_tdm_sample_prep_data.csv", index=False)
# 6. EBE 추정
# rEBE = EBE(PredVanco, DATAi, TH, OM, SG)
rEBE = EBE(PredVanco_PKS, DATAi, TH, OM, SG)

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


print(f"Vc coef: {round(Vc_ind/WT,3)}")
print(f"Cl coef: 0.05 (fixed)")
print(f"Renal Cl (slope): 0.75 (fixed)")

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
# TVCL


# 7. 예측 인터벌 계산
PI = calcTDM(PredVanco, DATAi, TH, SG, rEBE, TIME=50, AMT=1000, RATE=1000, II=12, ADDL=10)

# 8. TDM 작성에 필수적인 결과정리
# PI['y']

# 9. 시각화

# 스타일 지정 (선택)
sns.set(style="whitegrid")

# 그래프 그리기
plt.figure(figsize=(10, 6))

# 관측치
mask_obs = ~PI["y"].isna()
sns.scatterplot(x=PI["x"][mask_obs], y=PI["y"][mask_obs], color='black', s=40, label='Observed')

# 예측치 및 신뢰구간
sns.lineplot(x="x", y="y2", data=PI, color='dimgrey', label='Prediction (mean)', linewidth=2, alpha=0.5)
sns.lineplot(x="x", y="yciLL", data=PI, linestyle='--', color='indianred', label='95% CI lower', alpha=0.5)
sns.lineplot(x="x", y="yciUL", data=PI, linestyle='--', color='indianred', label='95% CI upper', alpha=0.5)
sns.lineplot(x="x", y="ypiLL", data=PI, linestyle=':', color='royalblue', label='95% PI lower', alpha=0.5)
sns.lineplot(x="x", y="ypiUL", data=PI, linestyle=':', color='royalblue', label='95% PI upper', alpha=0.5)

# 기준선
for y in [5, 15, 25, 35]:
    plt.axhline(y=y, linestyle='--', color='gray', linewidth=0.8)

# 축 및 레이블 설정
plt.xlim(PI["x"].min(), PI["x"].max())
plt.ylim([PI["ypiLL"].min(), PI["ypiUL"].max()])
plt.xlabel("Time")
plt.ylabel("Concentration ± 2SD")
plt.title("Vancomycin Concentration Prediction ± 2SD")
plt.legend()
plt.tight_layout()
plt.show()
#
# (CL_ind-CLnonrenal)/(CLCR * 0.06)
#
# Cl_coef = CL_nonrenal * 1000 / 60 / WT
# Renal_CL_slope = (Total_CL - CL_nonrenal) / (CLCR * 0.06)
#
