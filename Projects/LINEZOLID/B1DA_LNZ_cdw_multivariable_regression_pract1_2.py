import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency, fisher_exact

prj_name = 'LNZ'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"

############################
# Thrombocytopenia 발생에 있어서의 Dose period의 영향 분석

endpoint = 'PLT'
csv_path = f"{output_dir}/b1da/datacheck/b1da_lnz_datacheck_{endpoint}_df.csv"     # <-- 본인 파일 경로
df = pd.read_csv(csv_path)

# ---- 사용자 설정 ----
DOSE_COL = "DOSE_PERIOD"
EV_COL = "EV"
CUTOFF = 6

# ---- 0) 필수 컬럼 확인 & 결측 처리 ----
required_cols = [DOSE_COL, EV_COL]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise KeyError(f"Missing required columns: {missing}")

work = df.loc[:, required_cols].copy()
work = work.dropna(subset=[DOSE_COL, EV_COL])

# EV가 True/False, 0/1, '0'/'1' 등 섞여 있어도 0/1로 정리
# (EV가 이미 0/1 int라면 그대로)
work[EV_COL] = pd.to_numeric(work[EV_COL], errors="coerce")
work = work.dropna(subset=[EV_COL])
work[EV_COL] = work[EV_COL].astype(int)

# 0/1 이외 값이 있으면 에러
uniq_ev = set(work[EV_COL].unique())
if not uniq_ev.issubset({0, 1}):
    raise ValueError(f"{EV_COL} must be binary (0/1). Found: {sorted(uniq_ev)}")

# ---- 1) DOSE_PERIOD 그룹 생성 (<=6 vs >6) ----
group_col = f"{DOSE_COL}_GROUP"
work[group_col] = np.where(work[DOSE_COL] <= CUTOFF, f"<={CUTOFF}", f">{CUTOFF}")

# ---- 2) Contingency table (Count) ----
ct_count = pd.crosstab(work[group_col], work[EV_COL])

# EV=0,1 컬럼이 모두 나오도록 보정 (한쪽에 event가 0만 있거나 1만 있는 경우 대비)
ct_count = ct_count.reindex(columns=[0, 1], fill_value=0)
# 그룹도 <=cutoff, >cutoff 순서로
ct_count = ct_count.reindex(index=[f"<={CUTOFF}", f">{CUTOFF}"], fill_value=0)

print("\n[Contingency table: counts]")
print(ct_count)

# ---- 3) 비율 테이블 ----
ct_prop = ct_count.div(ct_count.sum(axis=1), axis=0)

print("\n[Row-wise proportion table]")
print(ct_prop)

print("\n[EV=1 proportion by group]")
print(ct_prop[1])

# ---- 4) 통계검정: 기대빈도 기반 자동 선택 (Chi-square vs Fisher) ----
table = ct_count.values  # shape (2,2)

# 기대빈도 계산 (correction=False로 expected만 산출용)
chi2, p_chi, dof, expected = chi2_contingency(table, correction=False)

if (expected < 5).any():
    # Fisher exact test
    odds_ratio, p_value = fisher_exact(table)
    test_used = "Fisher exact test"
else:
    # Chi-square test (+ OR 계산)
    # OR = (a/b) / (c/d) where:
    # group1(<=6): [EV0=c, EV1=d], group2(>6): [EV0=b, EV1=a]  (아래 방식으로 직접 계산)
    # 안전하게 0 division 체크
    a = table[1, 1]  # >6 & EV=1
    b = table[1, 0]  # >6 & EV=0
    c = table[0, 1]  # <=6 & EV=1
    d = table[0, 0]  # <=6 & EV=0

    # 0이 있어 OR가 무한/0 되는 경우를 그대로 반영 (원하면 Haldane-Anscombe 보정 가능)
    if b == 0 or c == 0:
        odds_ratio = np.inf if (a > 0 and d > 0 and (b == 0 or c == 0)) else np.nan
    else:
        odds_ratio = (a / b) / (c / d)

    p_value = p_chi
    test_used = "Chi-square test"

print("\n[Statistical comparison]")
print(f"Test used: {test_used}")
print(f"Odds Ratio (> {CUTOFF} vs <= {CUTOFF}): {odds_ratio}")
print(f"P-value: {p_value:.6g}")

# ---- 5) 논문/리뷰어 대응용 요약 테이블 ----
summary_df = (
    ct_count
    .assign(
        Total=lambda x: x.sum(axis=1),
        EV1_rate=lambda x: (x[1] / x.sum(axis=1)).replace([np.inf, -np.inf], np.nan)
    )
)

print("\n[Summary table]")
print(summary_df)

# ---- 6) 결과 문장(옵션) ----
rate_le = summary_df.loc[f"<={CUTOFF}", "EV1_rate"]
rate_gt = summary_df.loc[f">{CUTOFF}", "EV1_rate"]
n_le = summary_df.loc[f"<={CUTOFF}", "Total"]
n_gt = summary_df.loc[f">{CUTOFF}", "Total"]

print("\n[Auto sentence]")
print(
    f"EV=1 proportion was {rate_gt:.3f} (n={n_gt}) in DOSE_PERIOD>{CUTOFF} and "
    f"{rate_le:.3f} (n={n_le}) in DOSE_PERIOD<={CUTOFF}; "
    f"{test_used} p={p_value:.3e}, OR={odds_ratio}."
)


##############################

bl_timedelta_df = pd.read_csv(f"{output_dir}/b1da/mvlreg_output/b1da_lnz_mvlreg_baseline_timedelta_df.csv")

bl_timedelta_dist_df = bl_timedelta_df.groupby('UID',as_index=False).agg('min').drop(['UID','ENDPOINT'], axis=1)
bl_timedelta_dist_stats_df = pd.DataFrame([bl_timedelta_dist_df.mean().map(lambda x:round(x,1)), bl_timedelta_dist_df.std().map(lambda x:round(x,1))]).T.rename(columns={0:'Mean',1:'SD'})
bl_timedelta_dist_stats_df['Day-Duration'] = bl_timedelta_dist_stats_df.apply(lambda x:f"{x['Mean']} ({x['SD']})", axis=1)
bl_timedelta_dist_stats_df[['Day-Duration']].T.to_csv(f"{output_dir}/b1da/mvlreg_output/b1da_lnz_mvlreg_baseline_time_interval_stats.csv")
# bl_timedelta_dist_df

df_long = bl_timedelta_dist_df.melt(var_name="Variable", value_name="Value")

# boxplot
plt.figure(figsize=(14, 6))
sns.boxplot(
    data=df_long,
    x="Variable",
    y="Value",
    showfliers=False  # outlier 숨김 (원하면 True)
)

font_size = 16
plt.xticks(rotation=45, ha="right", fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.ylabel("Time interval (days)", fontsize=font_size)
plt.xlabel(None)
plt.title("Distribution of the Time Interval from Baseline Laboratory Testing to Initiation of Linezolid Therapy", fontsize=font_size)


plt.tight_layout()
plt.savefig(f"{output_dir}/b1da/mvlreg_output/time_interval(baseline to 1st_dose).png")  # PNG 파일로 저장

plt.cla()
plt.clf()
plt.close()
