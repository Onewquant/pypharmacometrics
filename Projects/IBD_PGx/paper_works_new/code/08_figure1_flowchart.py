"""Figure 1 - Eligibility flow chart (with PopPK and PGx branches).

Regenerates the cohort flow chart programmatically so the numbers stay
reproducible from Table_pgx_attrition.csv. Screening-stage counts that
come from the EMR extraction (total treated, exclusion 1/2 counts) are
placeholders ("n = X,XXX") to be filled in manually, as in the previous
PowerPoint version.

Flow (2026-09-06, EBE-based frame; former exclusion 4 merged into 3):
  analytic cohort 139
    -> infliximab PopPK modeling cohort 97     (excl. 3, n = 42: no infliximab
                                                records 41 + maintenance-
                                                starting patient with a single
                                                concentration, no EBE 1)
    -> pharmacogenomic analysis cohort 96      (excl. 4: sample removed in
                                                genotype QC)
    -> overall 96 / induction 83 / maintenance 96

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


N_ANALYTIC = 139
N_IFX = attr_n("Infliximab cohort")
N_EXCL_GENO = attr_n("genotype data removed")
N_GENOTYPED = attr_n("PGx analysis cohort")          # 97 = IFX - genotype QC
N_IND = attr_n("[IND] covariate-complete")
N_MAINT = attr_n("[MAINT] covariate-complete")
N_OVERALL = attr_n("[OVERALL] covariate-complete")
# patients in the infliximab cohort without an individual CL estimate
# (maintenance-starting, single concentration -> not in estimation dataset)
N_EXCL_NO_CL = N_GENOTYPED - N_OVERALL
N_POPPK = N_IFX - N_EXCL_NO_CL
# exclusion 3 = no infliximab records + no individual CL estimate (merged)
N_EXCL_3 = N_ANALYTIC - N_POPPK
N_PGX = N_POPPK - N_EXCL_GENO
assert N_PGX == N_OVERALL, (N_PGX, N_OVERALL)

INK = "#1F2937"
BOX = dict(boxstyle="square,pad=0.55", facecolor="white",
           edgecolor=INK, linewidth=1.0)

fig, ax = plt.subplots(figsize=(8.6, 8.4))
ax.set_xlim(0, 10)
ax.set_ylim(2.2, 14.2)
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
Y_TOP, Y_AN, Y_PK, Y_PGX = 13.3, 10.95, 8.55, 6.15
box(CX, Y_TOP, "Treatment with anti-TNF inhibitors*\n($\\geq$ 1 day)\nn = X,XXX")
box(CX, Y_AN, f"Analytic cohort\nn = {N_ANALYTIC}")
box(CX, Y_PK, f"Infliximab PopPK modeling cohort\nn = {N_POPPK}")
box(CX, Y_PGX, f"Pharmacogenomic analysis cohort\nn = {N_PGX}")

down_arrow(Y_TOP - 0.65, Y_AN + 0.5)
down_arrow(Y_AN - 0.5, Y_PK + 0.5)
down_arrow(Y_PK - 0.5, Y_PGX + 0.5)

# --- exclusion boxes -------------------------------------------------------
y = (Y_TOP + Y_AN) / 2 - 0.1
box(EX, y, "Exclusion 1:\nAnti-TNF inhibitor concentration\nnot available\nn = X,XXX\n\n"
           "Exclusion 2:\nWhole genome sequencing\nnot performed\nn = X,XXX", fontsize=8.5)
side_arrow(y)

y = (Y_AN + Y_PK) / 2
box(EX, y, "Exclusion 3:\nNo infliximab records or no\n"
           "evaluable infliximab concentration\n"
           f"n = {N_EXCL_3}", fontsize=8.5)
side_arrow(y)

y = (Y_PK + Y_PGX) / 2
box(EX, y, "Exclusion 4:\nGenotype data removed during\n"
           f"quality control\nn = {N_EXCL_GENO}", fontsize=8.5)
side_arrow(y)

# --- phase-specific analysis sets ------------------------------------------
Y_SPLIT, Y_PHASE = 4.7, 3.5
down_arrow(Y_PGX - 0.5, Y_SPLIT)
ax.plot([1.35, 7.75], [Y_SPLIT, Y_SPLIT], color=INK, linewidth=1.0)
for x in [1.35, CX, 7.75]:
    ax.annotate("", xy=(x, Y_PHASE + 0.52), xytext=(x, Y_SPLIT),
                arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.0))

box(1.35, Y_PHASE, f"Overall treatment\nn = {N_OVERALL}")
box(CX, Y_PHASE, f"Induction phase\nn = {N_IND}")
box(7.75, Y_PHASE, f"Maintenance phase\nn = {N_MAINT}")

# --- footnote ---------------------------------------------------------------
ax.text(0.05, 13.7, "* Anti-TNF inhibitors\n- Infliximab\n- Adalimumab\n- Ustekinumab",
        ha="left", va="top", fontsize=8.5, color=INK, linespacing=1.5)

fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(f"{cft_dir}/Figure1_eligibility_flowchart.{ext}",
                dpi=300, bbox_inches="tight")

print(f"saved: Figure1 (analytic {N_ANALYTIC} -excl3 {N_EXCL_3}-> PopPK {N_POPPK} -> PGx {N_PGX}; "
      f"OVERALL {N_OVERALL} / IND {N_IND} / MAINT {N_MAINT})")
