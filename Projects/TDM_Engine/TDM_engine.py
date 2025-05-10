from Projects.TDM_Engine.tools import *
import matplotlib.pyplot as plt
import seaborn as sns

# Typical parameters (THETA): CL, V1, V2, Q
TH = np.array([
    3.8135955291021233,
    39.889510090195238,
    44.981835351176571,
    2.0055189192561507
])

# Inter-individual variability matrix (OMEGA)
OM = np.array([
    [ 0.10855133849022583,     -0.0152445093837639736, -0.19698189309256298,   0.11914555131547180],
    [-0.0152445093837639736,    0.00301016276351715687, 0.023128525633882277, -0.006615975862019478],
    [-0.19698189309256298,      0.023128525633882277,   1.0420667710781930,   -0.19223488058085114],
    [ 0.11914555131547180,     -0.006615975862019478,  -0.19223488058085114,   0.416]
])

# Residual variability matrix (SIGMA)
SG = np.array([
    [0.14019731106615912 ** 2, 0.0],
    [0.0,                     1.8662549226759475 ** 2]
])


# 1. CSV 파일 읽기
raw_data = pd.read_csv(
    "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/TDM_Engine/sample_vanco_loading_maintenance.csv",
    na_values=["", ".", "NA"]
)

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
rEBE = EBE(PredVanco, DATAi, TH, OM, SG)

# 7. 예측 인터벌 계산
PI = calcTDM(PredVanco, DATAi, TH, SG, rEBE, TIME=50, AMT=1000, RATE=1000, II=12, ADDL=10)

# 8. TDM 작성에 필수적인 결과정리
PI['y']

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
