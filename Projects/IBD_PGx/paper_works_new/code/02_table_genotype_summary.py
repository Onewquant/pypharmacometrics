"""Table - Candidate variant genotype summary (counts, MAF, HWE).

For each candidate rsID: genotype counts (0/1/2 of the coded allele),
coded-allele frequency, MAF, and Hardy-Weinberg equilibrium exact test
p-value (Wigginton et al., 2005), in the genotyped cohort and in the
infliximab PGx cohort.

Inputs:
  - paper_works_new/data/rsid_dosage_matrix_with_alleles.csv
  - paper_works_new/data/for_genomics_df(all_drugs).csv

Output:
  - paper_works_new/output/Table_genotype_summary.csv
"""

import numpy as np
import pandas as pd

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
data_dir = f"{prj_dir}/paper_works_new/data"
output_dir = f"{prj_dir}/paper_works_new/output"

rsid_gene_dict = {
    "rs9828223": "CD96",
    "rs2097432": "HLA-DQA1",
    "rs396991": "FCGR3A",
    "rs1800629": "TNFα (TNF)",
    "rs4149570": "TNFRSF1A",
    "rs3397": "TNFRSF1B",
    "rs1061624": "TNFRSF1B",
    "rs5030728": "TLR4",
    "rs3804099": "TLR2",
    "rs10499563": "IL6",
    "rs2275913": "IL17A",
    "rs1800872": "IL10",
    "rs3024505": "IL10",
    "rs361525": "TNFα (TNF)",
    "rs767455": "TNFRSF1A",
    "rs1061622": "TNFRSF1B",
    "rs765249238": "SLCO2A1",
    "rs776813259": "SLCO2A1",
}


def hwe_exact_p(n_het, n_hom_rare, n_hom_common):
    """SNP-HWE exact test (Wigginton, Cutler & Abecasis, 2005)."""
    n_het = int(n_het)
    rare = 2 * int(n_hom_rare) + n_het
    n = int(n_het + n_hom_rare + n_hom_common)

    if n == 0 or rare == 0:
        return np.nan

    probs = np.zeros(rare + 1)
    mid = rare * (2 * n - rare) // (2 * n)
    if mid % 2 != rare % 2:
        mid += 1

    probs[mid] = 1.0

    het, hom_r, hom_c = mid, (rare - mid) // 2, n - mid - (rare - mid) // 2
    while het > 1:
        probs[het - 2] = (
            probs[het] * het * (het - 1.0)
            / (4.0 * (hom_r + 1.0) * (hom_c + 1.0))
        )
        het -= 2
        hom_r += 1
        hom_c += 1

    het, hom_r, hom_c = mid, (rare - mid) // 2, n - mid - (rare - mid) // 2
    while het <= rare - 2:
        probs[het + 2] = (
            probs[het] * 4.0 * hom_r * hom_c
            / ((het + 2.0) * (het + 1.0))
        )
        het += 2
        hom_r -= 1
        hom_c -= 1

    probs /= probs.sum()
    return float(min(1.0, probs[probs <= probs[n_het] * (1 + 1e-9)].sum()))


def summarize_cohort(df, rsid, label):
    dos = df[rsid].dropna()
    dos = dos[dos.isin([0, 1, 2])]
    n = len(dos)

    n0 = int((dos == 0).sum())
    n1 = int((dos == 1).sum())
    n2 = int((dos == 2).sum())

    coded_af = (n1 + 2 * n2) / (2 * n) if n > 0 else np.nan
    maf = min(coded_af, 1 - coded_af) if n > 0 else np.nan

    if coded_af <= 0.5:
        hwe_p = hwe_exact_p(n1, n2, n0)
    else:
        hwe_p = hwe_exact_p(n1, n0, n2)

    return {
        f"{label}_N": n,
        f"{label}_GENO_0/1/2": f"{n0}/{n1}/{n2}",
        f"{label}_CODED_ALLELE_FREQ": round(coded_af, 4),
        f"{label}_MAF": round(maf, 4),
        f"{label}_HWE_P": round(hwe_p, 4) if pd.notna(hwe_p) else np.nan,
    }


rsid_df = pd.read_csv(f"{data_dir}/rsid_dosage_matrix_with_alleles.csv")
rsid_df = rsid_df[~rsid_df["UID"].isna()].copy()
rsid_df["UID"] = rsid_df["UID"].map(lambda x: str(x).split(".")[0])
rsid_list = list(rsid_df.loc[:, "genomics_group":].columns)[1:]

ep_df = pd.read_csv(f"{data_dir}/for_genomics_df(all_drugs).csv")
ep_df["UID"] = ep_df["UID"].astype(str)
ifx_uids = set(ep_df.loc[ep_df["DRUG"] == "infliximab", "UID"])

ifx_rsid_df = rsid_df[rsid_df["UID"].isin(ifx_uids)].copy()

rows = []
for rsid in rsid_list:
    rs_base = rsid.split("(")[0]
    allele_info = rsid[rsid.find("("):] if "(" in rsid else ""

    row = {
        "RSID": rs_base,
        "GENE": rsid_gene_dict.get(rs_base, ""),
        "ALLELE_CODING": allele_info,
    }
    row.update(summarize_cohort(rsid_df, rsid, "GENOTYPED"))
    row.update(summarize_cohort(ifx_rsid_df, rsid, "IFX_COHORT"))
    rows.append(row)

result_df = pd.DataFrame(rows)
result_df.to_csv(
    f"{output_dir}/Table_genotype_summary.csv", index=False, encoding="utf-8-sig"
)

print(result_df.to_string(index=False))
