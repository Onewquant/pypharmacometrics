import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

prj_name = 'LNZ'
prj_dir = 'C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID'
resource_dir = f'{prj_dir}/resource'
output_dir = f"{prj_dir}/results"
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
