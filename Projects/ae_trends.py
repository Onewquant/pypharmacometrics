import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mcolors
import numpy as np

df = pd.DataFrame({
    "year": list(range(2000, 2026)),
    "reports": [
        264, 117, 158, 199, 620, 1400, 6239, 14453, 12796, 27010,
        64143, 74657, 92375, 183260, 183554, 198037, 228939, 252611,
        257438, 262983, 259089, 539441, 315867, 268148, 253486, 277279
    ]
})

sns.set_theme(style="whitegrid")

# blue gradient
cmap = mcolors.LinearSegmentedColormap.from_list(
    "blue_gradient",
    ["#DBEAFE", "#60A5FA", "#2563EB", "#1E3A8A"]
)
colors = [cmap(i) for i in np.linspace(0.15, 0.95, len(df))]

plt.figure(figsize=(13, 6))

ax = sns.barplot(
    data=df,
    x="year",
    y="reports",
    palette=colors,
    edgecolor="none"
)

ax.yaxis.set_major_formatter(
    FuncFormatter(lambda x, pos: f"{int(x):,}")
)

# title / labels 제거
ax.set_title("")
ax.set_xlabel("")
ax.set_ylabel("")

# grid 최소화
ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)
ax.grid(axis="x", visible=False)

sns.despine(left=True, bottom=True)

plt.xticks(rotation=45, color="#6B7280", fontsize=10)
plt.yticks(color="#9CA3AF", fontsize=10)

plt.tight_layout()

plt.savefig(
    "C:/Users/ilma0/Downloads/adverse_event_reports_gradient_bar.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()