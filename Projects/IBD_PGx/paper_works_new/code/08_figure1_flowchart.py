"""Figure 1 - Eligibility flow chart (with PGx branch, current numbers).

Regenerates the cohort flow chart programmatically so the numbers stay
reproducible from Table_pgx_attrition.csv. Screening-stage counts that
come from the EMR extraction (total treated, exclusion 1/2 counts) are
placeholders ("n = X,XXX") to be filled in manually, as in the previous
PowerPoint version.

Outputs:
  - paper_works_new/core_fig_tab/Figure1_eligibility_flowchart.png (300 dpi)
  - paper_works_new/core_fig_tab/Figure1_eligibility_flowchart.pdf
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
out_dir = f"{prj_dir}/paper_works_new/output"
cft_dir = f"{prj_dir}/paper_works_new/core_fig_tab"

# reproducible counts
attr = pd.read_csv(f"{out_dir}/Table_pgx_attrition.csv")


def attr_n(step_substr):
    row = attr[attr["STEP"].str.contains(step_substr, regex=False)]
    return int(row["N"].iloc[0]) if len(row) else None


N_IFX = attr_n("Infliximab cohort")
N_EXCL_GENO = attr_n("genotype data removed")
N_PGX = attr_n("PGx analysis cohort")
N_IND = attr_n("[IND] covariate-complete")
N_MAINT = attr_n("[MAINT] covariate-complete")
N_OVERALL = attr_n("[OVERALL] covariate-complete")
N_ANALYTIC = 139
N_EXCL_NO_IFX = N_ANALYTIC - N_IFX

INK = "#1F2937"
BOX = dict(boxstyle="square,pad=0.55", facecolor="white",
           edgecolor=INK, linewidth=1.0)

fig, ax = plt.subplots(figsize=(8.6, 8.2))
ax.set_xlim(0, 10)
ax.set_ylim(0, 12)
ax.axis("off")

CX = 4.55   # main column center
EX = 8.15   # exclusion column center


def box(x, y, text, fontsize=9.5):
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            color=INK, bbox=BOX, linespacing=1.6)


def down_arrow(y_from, y_to):
    ax.annotate("", xy=(CX, y_to), xytext=(CX, y_from),
                arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.0))


def side_arrow(y):
    ax.annotate("", xy=(EX - 1.55, y), xytext=(CX, y),
                arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.0))


# --- main column -----------------------------------------------------------
box(CX, 11.2, "Treatment with anti-TNF inhibitors*\n($\\geq$ 1 day)\nn = X,XXX")
box(CX, 8.85, f"Analytic cohort\nn = {N_ANALYTIC}")
box(CX, 6.55, f"Infliximab PopPK modeling cohort\nn = {N_IFX}")
box(CX, 4.25, f"Pharmacogenomic analysis cohort\nn = {N_PGX}")

down_arrow(10.55, 9.35)
down_arrow(8.35, 7.05)
down_arrow(6.05, 4.75)

# --- exclusion boxes -------------------------------------------------------
box(EX, 9.95, "Exclusion 1:\nAnti-TNF inhibitor concentration\nnot available\nn = X,XXX\n\n"
              "Exclusion 2:\nWhole genome sequencing\nnot performed\nn = X,XXX", fontsize=8.5)
side_arrow(9.95)

box(EX, 7.7, f"Exclusion 3:\nNo infliximab records\nn = {N_EXCL_NO_IFX}", fontsize=8.5)
side_arrow(7.7)

box(EX, 5.4, "Exclusion 4:\nGenotype data removed during\n"
             f"quality control\nn = {N_EXCL_GENO}", fontsize=8.5)
side_arrow(5.4)

# --- phase-specific analysis sets ------------------------------------------
down_arrow(3.72, 2.9)
ax.plot([1.35, 7.75], [2.9, 2.9], color=INK, linewidth=1.0)
for x in [1.35, CX, 7.75]:
    ax.annotate("", xy=(x, 2.42), xytext=(x, 2.9),
                arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.0))

box(1.35, 1.7, f"Overall treatment\nn = {N_OVERALL}")
box(CX, 1.7, f"Induction phase\nn = {N_IND}")
box(7.75, 1.7, f"Maintenance phase\nn = {N_MAINT}")

# --- footnote ---------------------------------------------------------------
ax.text(0.05, 11.6, "* Anti-TNF inhibitors\n- Infliximab\n- Adalimumab\n- Ustekinumab",
        ha="left", va="top", fontsize=8.5, color=INK, linespacing=1.5)

fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(f"{cft_dir}/Figure1_eligibility_flowchart.{ext}",
                dpi=300, bbox_inches="tight")

print(f"saved: Figure1 (IFX {N_IFX} -> excl {N_EXCL_GENO} -> PGx {N_PGX}; "
      f"OVERALL {N_OVERALL} / IND {N_IND} / MAINT {N_MAINT})")
