from Projects.TDM_Engine.tools_core import *
from Projects.TDM_Engine.drug_models import *

# 0. JSON request data 생성
dose_conc_data = pd.read_csv(
    "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/TDM_Engine/sample_vanco_pks.csv",
    na_values=["", ".", "NA"]
)
comp_doseconc_data = dose_conc_data.drop(['ID','SEX','AGE','BWT',],axis=1)

comp_dose_conc_json = comp_doseconc_data.to_json(orient="records") # 행 단위로 json 변환
req_raw_dict = {
    "meta": {"TDM_DT": datetime.now().isoformat(),
             "HOSPITAL": 'SNUBH',
             "DEP": 'CPT',
             "USERNAME": 'BSH',
             "DESC": "샘플 사용자 데이터"},
    "patient": {"UID": 'test123',
                "NAME": '홍길동',
                "SEX": 'M',
                "AGE": 66,
                "HT":165.6,
                "WT":56
                },
    "data": comp_dose_conc_json
}

req_data = json.dumps(req_raw_dict, indent=2, ensure_ascii=False)

# 1. JSON request data를 dataframe으로 변환

pinfo_dict = json.loads(req_data)['patient']
raw_data = pd.DataFrame(json.loads(json.loads(req_data)['data']))
for pcol, pval in pinfo_dict.items():
    raw_data[pcol] = pval

# 2. TIME 처리
data_prepped = convDT(raw_data)

# 3. ADDL 값이 유효하면 expand
if "ADDL" in data_prepped.columns and data_prepped["ADDL"].notna().any() and (data_prepped["ADDL"] > 0).any():
    data_prepped = expandDATA(data_prepped)

# 4. ID, TIME 기준 정렬
DATAi = data_prepped.sort_values(by=["UID", "TIME"]).reset_index(drop=True)

# 6. EBE 추정
dmodel = vanco_adult()
e = EBEEnvironment()
rEBE = EBE(dmodel.Pred, DATAi, dmodel.TH, dmodel.OM, dmodel.SG, e)

json.dumps(rEBE, default=lambda x: x.tolist())


