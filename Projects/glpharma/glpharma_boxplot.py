import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# CSV 파일 불러오기
df = pd.read_csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/glpharma/resource/Ka Box plot (Tvs.R).csv")

# Boxplot 그리기
plt.figure(figsize=(12, 10))
sns.boxplot(
    x="TRT",
    y="K01",
    data=df,
    order=["T", "R"],
    palette={"T": "#ff7f0e", "R": "#1f77b4"}  # T: orange, R: blue
)
plt.title("Boxplot of Ka by Treatment Group", fontsize=20)
plt.xlabel("Treatment Group", fontsize=20, labelpad=8)
plt.ylabel("Ka", fontsize=20, labelpad=8)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True, axis='y')
plt.tight_layout(pad=3.5)
plt.show()
