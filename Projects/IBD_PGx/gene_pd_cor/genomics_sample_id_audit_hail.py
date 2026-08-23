"""[Hail 환경용] VCF 샘플 목록 추출 + 미매칭 샘플 추적.

유전체 분석 스크립트(hl.import_vcf ~ dosage matrix 생성)의 **뒤에 이어서**
실행하거나, 해당 스크립트에 이 블록을 붙여 사용.

목적:
  1) 현재 VCF에 실제로 들어있는 전체 샘플 ID를 파일로 확정 저장
     -> pid_df와 매칭되지 않는 샘플이 "이름 표기 차이" 때문인지,
        "샘플 자체가 VCF에 없음" 때문인지 구분하기 위함
  2) 찾고자 하는 특정 샘플(예: UID 17439372 = 23-B02281_EB-01)이
     어떤 표기 변형으로든 VCF에 존재하는지 확정 판정

출력:
  - vcf_sample_list.csv          : VCF의 전체 샘플 ID (s)
  - vcf_sample_search_report.txt : 타깃 샘플 탐색 결과
"""

import re

import pandas as pd

# 찾고자 하는 샘플: pid_df.csv의 s 값 기준 (UID -> 시퀀싱 샘플 ID)
TARGET_SAMPLES = {
    "17439372": "23-B02281_EB-01",
}

OUT_SAMPLE_LIST = "vcf_sample_list.csv"
OUT_REPORT = "vcf_sample_search_report.txt"


def norm_key(x):
    """샘플 ID 표기 정규화: sample_ 접두사, _EB-01 접미사, 대소문자/공백 제거."""
    x = str(x).strip().upper()
    x = re.sub(r"^SAMPLE[_-]", "", x)
    x = re.sub(r"[_-]EB[_-]?\d+$", "", x)
    x = re.sub(r"\s+", "", x)
    return x


def digit_key(x):
    """숫자만 남긴 키 (구분자/접두사 차이를 무시한 2차 매칭용).

    반드시 norm_key를 거친 뒤 숫자를 추출한다. 원본에서 바로 추출하면
    '_EB-01' 접미사의 숫자까지 섞여 표기 관례가 다른 두 목록 간
    비교가 아예 성립하지 않는다.
    """
    return re.sub(r"\D", "", norm_key(x))


# mt은 유전체 분석 스크립트에서 이미 생성된 MatrixTable
vcf_samples = mt.s.collect()

pd.DataFrame({"s": vcf_samples}).to_csv(OUT_SAMPLE_LIST, index=False)
print(f"[저장] {OUT_SAMPLE_LIST} (샘플 {len(vcf_samples)}건)")

norm_map = {}
digit_map = {}
for s in vcf_samples:
    norm_map.setdefault(norm_key(s), []).append(s)
    digit_map.setdefault(digit_key(s), []).append(s)

lines = [f"VCF 샘플 수: {len(vcf_samples)}", ""]

for uid, target in TARGET_SAMPLES.items():
    nk, dk = norm_key(target), digit_key(target)
    exact = [s for s in vcf_samples if s == target]
    by_norm = norm_map.get(nk, [])
    by_digit = digit_map.get(dk, [])
    # 부분 일치(숫자 뒷부분)로도 훑어봄
    tail = digit_key(target)[-5:]
    by_part = [s for s in vcf_samples if tail and tail in digit_key(s)]

    lines += [
        f"[UID {uid}] 타깃 샘플: {target}",
        f"  정규화 key : {nk} / 숫자 key: {dk}",
        f"  완전일치   : {exact}",
        f"  정규화일치 : {by_norm}",
        f"  숫자일치   : {by_digit}",
        f"  부분일치({tail}) : {by_part}",
    ]
    if exact or by_norm or by_digit:
        lines.append("  => VCF에 존재함. 아래 reconcile 스크립트로 UID 부여 가능.")
    else:
        lines.append(
            "  => VCF에 존재하지 않음. 표기 문제가 아니라 샘플 자체가 "
            "이 VCF에 포함되지 않은 것 (시퀀싱 미실시 / 상위 QC 탈락 / "
            "다른 배치 여부 확인 필요)."
        )
    lines.append("")

report = "\n".join(lines)
with open(OUT_REPORT, "w", encoding="utf-8") as fh:
    fh.write(report)

print(report)
