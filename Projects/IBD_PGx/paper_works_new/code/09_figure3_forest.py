"""Figure 3 - Forest plot of candidate-variant effects on infliximab CL
(overall treatment period, log-scale ANCOVA).

Two panels (recessive model / dominant model). Each row is one variant:
adjusted geometric mean ratio (GMR) with 95% CI for the variant group vs
the reference group, sorted by P value within the panel; the numeric
GMR (95% CI), P and FDR q are printed as a text column. The dashed line
marks GMR = 1 (no association). The plot makes the main result visible
at a glance: no q value is below 0.05.

Run 03_pgx_ancova_fdr.py first.

Outputs -> paper_works_new/core_fig_tab/
  Figure3_forest_CL_overall.png (300 dpi) / .pdf
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter, NullLocator

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
output_dir = f"{prj_dir}/paper_works_new/output"
cft_dir = f"{prj_dir}/paper_works_new/core_fig_tab"

INK = "#1F2937"
MARK = "#0072B2"
REF_LINE = "#9CA3AF"
GRID = "#E5E7EB"

PANELS = [("HOM_vs_OTHERS", "Recessive model (homozygotes vs others)"),
          ("CARRIER_vs_NONCARRIER", "Dominant model (carriers vs non-carriers)")]


def genotype_labels(rsid_full):
    body = rsid_full[rsid_full.find("(") + 1: rsid_full.rfind(")")]
    ref = alt = "?"
    for part in body.split(","):
        if part.startswith("0="):
            ref = part[2:]
        elif part.startswith("1="):
            alt = part[2:]
    return ref + ref, ref + alt, alt + alt


def fmt_p(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}"


res = pd.read_csv(f"{output_dir}/Table_pgx_ancova_fdr_results.csv")
res = res[(res["END_POINT"] == "CL") & (res["MODEL_SCALE"] == "log(CL)")
          & (res["PHASE"] == "OVERALL")].copy()
res["RS"] = res["RSID"].str.split("(").str[0]
res["GENE_SHORT"] = res["GENE"].str.replace("TNF\u03b1 (TNF)", "TNF", regex=False)

n_rows = [len(res[res["COMPARISON"] == c]) for c, _ in PANELS]
fig, axes = plt.subplots(
    2, 1, figsize=(7.2, 0.27 * sum(n_rows) + 1.9),
    gridspec_kw={"height_ratios": n_rows, "hspace": 0.35},
)

X_TICKS = [0.7, 0.8, 1.0, 1.25, 1.5, 1.8]
X_LIM = (0.66, 1.9)

for ax, (comp, title) in zip(axes, PANELS):
    d = res[res["COMPARISON"] == comp].sort_values("P_VALUE", ascending=True)
    d = d.iloc[::-1].reset_index(drop=True)      # smallest P at the top
    y = np.arange(len(d))

    ax.axvline(1.0, color=REF_LINE, linewidth=0.9, linestyle="--", zorder=1)
    ax.hlines(y, d["CI_LOWER"], d["CI_UPPER"], color=MARK, linewidth=1.4, zorder=2)
    ax.scatter(d["EFFECT"], y, s=30, facecolor=MARK, edgecolor="white",
               linewidth=0.8, zorder=3)

    labels = []
    for _, r in d.iterrows():
        hom0, het, hom2 = genotype_labels(r["RSID"])
        grp = (f"{hom2} vs {hom0}+{het}" if comp == "HOM_vs_OTHERS"
               else f"{het}+{hom2} vs {hom0}")
        labels.append(f"{r['GENE_SHORT']} {r['RS']}   {grp} "
                      f"({int(r['VARIANT_N'])}/{int(r['REFERENCE_N'])})")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.4, color=INK)
    ax.set_ylim(-0.7, len(d) - 0.3)

    # text column: GMR (95% CI), P, q  - text tokens, not series color
    for yi, (_, r) in zip(y, d.iterrows()):
        ax.text(1.05, yi,
                f"{r['EFFECT']:.2f} ({r['CI_LOWER']:.2f}\u2013{r['CI_UPPER']:.2f})",
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=7.2, color=INK, clip_on=False)
        ax.text(1.57, yi, fmt_p(r["P_VALUE"]),
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=7.2, color=INK, clip_on=False)
        ax.text(1.80, yi, fmt_p(r["P_VALUE_FDR"]),
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=7.2, color=INK, clip_on=False)
    head_y = len(d) - 0.3 + 0.55
    for xpos, txt in [(1.05, "GMR (95% CI)"), (1.57, "P"), (1.80, "q")]:
        ax.text(xpos, head_y, txt, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=7.4, color=INK,
                fontweight="bold", clip_on=False)

    ax.set_xscale("log")
    ax.set_xlim(*X_LIM)
    ax.xaxis.set_major_locator(FixedLocator(X_TICKS))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in X_TICKS]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="x", labelsize=7.6, colors=INK)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, linewidth=0.6, zorder=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(INK)
    ax.set_title(title, fontsize=9, color=INK, loc="left", pad=16)

axes[-1].set_xlabel("Adjusted geometric mean ratio of infliximab clearance (95% CI)",
                    fontsize=8.5, color=INK)

fig.subplots_adjust(left=0.40, right=0.66, top=0.93, bottom=0.08)
for ext in ["png", "pdf"]:
    fig.savefig(f"{cft_dir}/Figure3_forest_CL_overall.{ext}", dpi=300,
                bbox_inches="tight")
print(f"saved: Figure3_forest_CL_overall.png / .pdf "
      f"({n_rows[0]} recessive + {n_rows[1]} dominant rows)")
