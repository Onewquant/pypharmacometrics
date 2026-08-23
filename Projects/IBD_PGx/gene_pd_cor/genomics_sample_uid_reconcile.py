"""샘플 ID <-> 환자 UID 매칭 정규화 및 최종 dosage matrix 생성.

유전체 분석 스크립트가 만든 dosage matrix(샘플 ID 기준)에 pid_df의 UID를
붙이는 단계를 코드로 명시화한 것. 기존에는 이 병합 과정이 코드에 남아있지
않아 미매칭 사유를 추적할 수 없었음.

매칭 규칙:
  샘플 ID 표기가 두 가지 관례로 섞여 있음
    - VCF/matrix : 'sample_18-00544' (sample_ 접두사, EB 접미사 없음)
                   '17-00006_EB-01'  (접두사 없음, EB 접미사 있음)
    - pid_df     : '23-B02281_EB-01' 등
  -> sample_ 접두사와 _EB-nn 접미사를 제거한 정규화 key로 매칭하고,
     그래도 실패하면 숫자만 남긴 key로 2차 매칭.

출력:
  - rsid_dosage_matrix_with_uid.csv   : UID/genomics_group이 부여된 최종 matrix
  - sample_uid_audit.csv              : 샘플별 매칭 방식 및 미매칭 사유
  - (콘솔) 미매칭 요약 및 타깃 UID 판정
"""

import re

import numpy as np
import pandas as pd

prj_dir = "C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/IBD_PGx"
resource_dir = f"{prj_dir}/gene_pd_cor"
output_dir = f"{prj_dir}/gene_pd_cor"

MATRIX_IN = f"{resource_dir}/rsid_dosage_matrix_with_alleles.csv"
PID_IN = f"{resource_dir}/pid_df.csv"

MATRIX_OUT = f"{output_dir}/rsid_dosage_matrix_with_uid.csv"
AUDIT_OUT = f"{output_dir}/sample_uid_audit.csv"

# 추적 대상 UID (ID 매칭 실패로 PGx 코호트에서 빠졌던 환자)
TARGET_UIDS = ["17439372"]


def norm_key(x):
    x = str(x).strip().upper()
    x = re.sub(r"^SAMPLE[_-]", "", x)
    x = re.sub(r"[_-]EB[_-]?\d+$", "", x)
    x = re.sub(r"\s+", "", x)
    return x


def digit_key(x):
    # norm_key를 거친 뒤 숫자 추출 (EB 접미사 숫자 오염 방지)
    return re.sub(r"\D", "", norm_key(x))


mat_df = pd.read_csv(MATRIX_IN)
pid_df = pd.read_csv(PID_IN)
pid_df["UID"] = pid_df["UID"].map(lambda x: str(x).split(".")[0])

rsid_cols = [c for c in mat_df.columns if c.startswith("rs")]

mat_df["NORM_KEY"] = mat_df["s"].map(norm_key)
mat_df["DIGIT_KEY"] = mat_df["s"].map(digit_key)
pid_df["NORM_KEY"] = pid_df["s"].map(norm_key)
pid_df["DIGIT_KEY"] = pid_df["s"].map(digit_key)

# 정규화 key 충돌 확인 (매칭 신뢰도 전제조건)
for name, df in [("matrix", mat_df), ("pid_df", pid_df)]:
    dup = df[df["NORM_KEY"].duplicated(keep=False)]
    if len(dup):
        print(f"[경고] {name}에 정규화 key 중복 {len(dup)}건:")
        print(dup[["s", "NORM_KEY"]].to_string(index=False))

norm_lookup = dict(zip(pid_df["NORM_KEY"], pid_df["UID"]))
digit_lookup = (
    pid_df.drop_duplicates("DIGIT_KEY").set_index("DIGIT_KEY")["UID"].to_dict()
)

matched_uid = []
match_type = []
for _, row in mat_df.iterrows():
    if row["NORM_KEY"] in norm_lookup:
        matched_uid.append(norm_lookup[row["NORM_KEY"]])
        match_type.append("normalized")
    elif row["DIGIT_KEY"] in digit_lookup:
        matched_uid.append(digit_lookup[row["DIGIT_KEY"]])
        match_type.append("digit_only")
    else:
        matched_uid.append(np.nan)
        match_type.append("unmatched")

mat_df["UID_NEW"] = matched_uid
mat_df["MATCH_TYPE"] = match_type

# genomics_group: 연구 대상 환자(UID 부여됨)=1, 그 외 샘플=0
mat_df["GENOMICS_GROUP_NEW"] = np.where(mat_df["UID_NEW"].notna(), 1, 0)

# 기존 UID 컬럼과의 차이 확인 (새로 살아난 샘플 / 바뀐 샘플)
old_uid = mat_df["UID"].map(
    lambda x: np.nan if pd.isna(x) else str(x).split(".")[0]
)
newly_matched = mat_df[old_uid.isna() & mat_df["UID_NEW"].notna()]
changed = mat_df[
    old_uid.notna()
    & mat_df["UID_NEW"].notna()
    & (old_uid != mat_df["UID_NEW"])
]

print(f"matrix 샘플 {len(mat_df)}건 / pid_df UID {pid_df['UID'].nunique()}건")
print(f"  매칭 성공 : {int(mat_df['UID_NEW'].notna().sum())}건 "
      f"(정규화 {(mat_df['MATCH_TYPE']=='normalized').sum()}, "
      f"숫자매칭 {(mat_df['MATCH_TYPE']=='digit_only').sum()})")
print(f"  미매칭    : {int(mat_df['UID_NEW'].isna().sum())}건")
print(f"  기존 대비 새로 매칭된 샘플 : {len(newly_matched)}건")
print(f"  기존과 UID가 달라진 샘플   : {len(changed)}건")
if len(changed):
    print(changed[["s", "MATCH_TYPE"]].assign(OLD=old_uid[changed.index],
                                              NEW=changed["UID_NEW"]).to_string(index=False))

# pid_df에는 있으나 matrix(=VCF)에 없는 UID -> 유전형 확보 불가
unresolved = pid_df[~pid_df["NORM_KEY"].isin(set(mat_df["NORM_KEY"]))
                    & ~pid_df["DIGIT_KEY"].isin(set(mat_df["DIGIT_KEY"]))]
print()
print(f"pid_df에 있으나 VCF/matrix에 없는 샘플: {len(unresolved)}건")
if len(unresolved):
    print(unresolved[["UID", "s"]].to_string(index=False))

print()
for uid in TARGET_UIDS:
    hit = mat_df[mat_df["UID_NEW"] == uid]
    src = pid_df[pid_df["UID"] == uid]
    label = f"[타깃 UID {uid}]"
    if len(hit):
        print(f"{label} 매칭 성공 -> 샘플 {hit['s'].tolist()} "
              f"({hit['MATCH_TYPE'].iloc[0]})")
    elif len(src):
        print(f"{label} pid_df 상 샘플 = {src['s'].iloc[0]} 이지만 "
              f"VCF/matrix에 해당 샘플이 없어 유전형 확보 불가")
    else:
        print(f"{label} pid_df에 UID 자체가 없음")

audit_df = mat_df[["s", "NORM_KEY", "DIGIT_KEY", "MATCH_TYPE", "UID_NEW",
                   "GENOMICS_GROUP_NEW"]].rename(
    columns={"UID_NEW": "UID", "GENOMICS_GROUP_NEW": "genomics_group"})
audit_df.to_csv(AUDIT_OUT, index=False, encoding="utf-8-sig")

final_df = mat_df[["s", "UID_NEW", "GENOMICS_GROUP_NEW"] + rsid_cols].rename(
    columns={"UID_NEW": "UID", "GENOMICS_GROUP_NEW": "genomics_group"})
final_df.to_csv(MATRIX_OUT, index=False, encoding="utf-8-sig")

print()
print(f"[저장] {MATRIX_OUT}")
print(f"[저장] {AUDIT_OUT}")
