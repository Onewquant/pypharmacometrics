"""Figure - Infliximab CL by rs1061622 (TNFRSF1B) genotype.

Two panels (maintenance phase / whole phase): individual model-estimated
CL by genotype (TT / TG / GG), log-scaled y-axis, geometric mean with
95% CI per genotype, and the ANCOVA GMR + FDR-adjusted p for the
GG vs TT+TG comparison annotated from the main analysis output.

Run 03_pgx_ancova_fdr.py first.

Outputs:
  - paper_works_new/output/Figure_CL_by_rs1061622.png (300 dpi)
  - paper_works_new/output/Figure_CL_by_rs1061622.pdf
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter, NullLocator

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
data_dir = f"{prj_dir}/paper_works_new/data"
output_dir = f"{prj_dir}/paper_works_new/output"

RSID = "rs1061622(0=T,1=G)"
GENO_LABELS = {0: "TT", 1: "TG", 2: "GG"}

# categorical palette validated for CVD safety (dataviz six-checks):
# reference genotypes (TT, TG) vs homozygote variant (GG)
COLOR_REF = "#B25C00"
COLOR_HOM = "#0072B2"
COLOR_SUMMARY = "#1F2937"

PANELS = [("MAINT", "Maintenance phase"), ("ALL", "Whole phase")]


def build_phase_df(ep_df, phase):
    cond = (
        ep_df["PHASE"].isin(["IND", "MAINT"])
        if phase == "ALL" else ep_df["PHASE"] == phase
    )
    uid_df = (
        ep_df[cond]
        .groupby("UID", as_index=False)
        .agg(CL=("CL", "mean"))
    )
    return uid_df[uid_df["CL"] > 0]


def geomean_ci(x):
    log_x = np.log(x)
    m, se = log_x.mean(), log_x.std(ddof=1) / np.sqrt(len(x))
    return np.exp(m), np.exp(m - 1.96 * se), np.exp(m + 1.96 * se)


rsid_df = pd.read_csv(f"{data_dir}/rsid_dosage_matrix_with_alleles.csv")
rsid_df = rsid_df[~rsid_df["UID"].isna()].copy()
rsid_df["UID"] = rsid_df["UID"].map(lambda x: str(x).split(".")[0])
geno_df = rsid_df[["UID", RSID]].rename(columns={RSID: "DOS"})

ep_df = pd.read_csv(f"{data_dir}/for_genomics_df(all_drugs).csv")
ep_df["UID"] = ep_df["UID"].astype(str)
ep_df["PHASE"] = ep_df["PHASE"].map(lambda x: x.split("_")[0])
ep_df = ep_df[ep_df["DRUG"] == "infliximab"].copy()

main_df = pd.read_csv(f"{output_dir}/Table_pgx_ancova_fdr_results.csv")
annot_df = main_df[
    (main_df["RSID"] == RSID)
    & (main_df["COMPARISON"] == "HOM_vs_OTHERS")
    & (main_df["END_POINT"] == "CL")
    & (main_df["MODEL_SCALE"] == "log(CL)")
].set_index("PHASE")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
rng = np.random.default_rng(1)

for ax, (phase, title) in zip(axes, PANELS):
    d = build_phase_df(ep_df, phase).merge(geno_df, on="UID")
    d = d[d["DOS"].isin([0, 1, 2])]

    for dos in [0, 1, 2]:
        vals = d.loc[d["DOS"] == dos, "CL"].values
        if len(vals) == 0:
            continue

        color = COLOR_HOM if dos == 2 else COLOR_REF
        x_jit = dos + rng.uniform(-0.14, 0.14, len(vals))
        ax.scatter(
            x_jit, vals, s=22, facecolor=color, edgecolor="white",
            linewidth=0.6, alpha=0.85, zorder=3,
        )

        gm, lcl, ucl = geomean_ci(vals)
        ax.hlines(gm, dos - 0.26, dos + 0.26, color=COLOR_SUMMARY,
                  linewidth=1.8, zorder=4)
        ax.vlines(dos + 0.26, lcl, ucl, color=COLOR_SUMMARY,
                  linewidth=1.1, zorder=4)
        ax.hlines([lcl, ucl], dos + 0.22, dos + 0.30, color=COLOR_SUMMARY,
                  linewidth=1.1, zorder=4)

    counts = d["DOS"].value_counts()
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels([
        f"{GENO_LABELS[k]}\n(n={counts.get(k, 0)})" for k in [0, 1, 2]
    ], fontsize=8)

    ax.set_yscale("log")
    ax.set_title(title, fontsize=9.5, pad=30)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", which="major", color="0.9", linewidth=0.6, zorder=0)

    # plain decimal ticks instead of scientific notation on the log axis
    y_ticks = [0.15, 0.2, 0.3, 0.4, 0.5, 0.6]
    ax.yaxis.set_major_locator(FixedLocator(y_ticks))
    ax.yaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in y_ticks]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(0.12, 0.75)

    # bracket: GG vs TT+TG (ANCOVA on log CL, FDR-adjusted p)
    if phase in annot_df.index:
        row = annot_df.loc[phase]
        y_br = 0.55
        ax.plot([0, 0, 2, 2],
                [y_br, y_br * 1.05, y_br * 1.05, y_br],
                color=COLOR_SUMMARY, linewidth=0.9, clip_on=False)
        ax.text(
            1.0, y_br * 1.10,
            f"GMR {row['EFFECT']:.2f} "
            f"({row['CI_LOWER']:.2f}\u2013{row['CI_UPPER']:.2f}), "
            f"q = {row['P_VALUE_FDR']:.3f}",
            ha="center", va="bottom", fontsize=7.5, color=COLOR_SUMMARY,
        )

axes[0].set_ylabel("Infliximab clearance (L/day)", fontsize=9)
fig.suptitle("rs1061622 (TNFRSF1B) genotype and infliximab clearance",
             fontsize=10, y=1.02)
fig.tight_layout()

for ext in ["png", "pdf"]:
    fig.savefig(
        f"{output_dir}/Figure_CL_by_rs1061622.{ext}",
        dpi=300, bbox_inches="tight",
    )

print("saved: Figure_CL_by_rs1061622.png / .pdf")
