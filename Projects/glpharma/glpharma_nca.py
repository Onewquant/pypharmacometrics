from pynca.tools import *

## 기본정보입력

project_dict = {'GLPHARMA':['W2406']}
resource_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/glpharma/resource"

df = pd.read_csv(f'{resource_dir}/prep_data/GLPHARMA_ConcPrep_W2406_R.csv')

result = tblNCA(df, key=["DRUG", "ID"], colTime="ATIME", colConc="CONC",
                dose='DOSE', adm="Extravascular", dur=0, doseUnit="mg",
                timeUnit="h", concUnit="ug/mL", down="Log", R2ADJ=0,
                MW=0, SS=False, iAUC="", excludeDelta=1, slopeMode='SNUHCPT', colStyle='pw')


result.columns
result['Cl_F_obs']
result['Vz_F_obs']
result.to_excel(f'{resource_dir}/Final Parameters Pivoted.xlsx', index=False, encoding='utf-8')


# print(result)