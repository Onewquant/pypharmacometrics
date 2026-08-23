# =====================================================================
# [Hail 세션에 바로 붙여 실행] mt의 샘플 목록 vs VCF 파일 헤더 비교
#   목적: 타깃 샘플이 (a) 분석 코드 필터로 빠진 것인지
#                    (b) VCF 파일 자체에 없는 것인지 확정 판정
#   원리: filter_rows(변이 필터)는 샘플을 제거하지 않으므로,
#         mt의 샘플 집합 == VCF 헤더의 샘플 집합 이면 필터 영향 없음.
# =====================================================================

import gzip
import re

VCF_PATH = f"{resource_dir}/IBD_SNP.hailQC.annotated.tagged.vcf.gz"
TARGETS = {"17439372": "23-B02281_EB-01"}


def norm_key(x):
    x = str(x).strip().upper()
    x = re.sub(r"^SAMPLE[_-]", "", x)
    x = re.sub(r"[_-]EB[_-]?\d+$", "", x)
    return re.sub(r"\s+", "", x)


def digit_key(x):
    return re.sub(r"\D", "", norm_key(x))


def vcf_header_samples(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    return []


mt_samples = mt.cols().s.collect()
file_samples = vcf_header_samples(VCF_PATH)

print(f"mt 샘플 수        : {len(mt_samples)}  (count_cols={mt.count_cols()})")
print(f"VCF 헤더 샘플 수  : {len(file_samples)}")

only_file = sorted(set(file_samples) - set(mt_samples))
only_mt = sorted(set(mt_samples) - set(file_samples))
print(f"파일에만 있음(=필터로 빠짐) : {len(only_file)}건 {only_file[:10]}")
print(f"mt에만 있음(비정상)         : {len(only_mt)}건 {only_mt[:10]}")

if not only_file:
    print("=> mt와 VCF의 샘플 집합이 동일. 분석 코드에서 제거된 샘플 없음.")
else:
    print("=> 분석 코드 경로에서 샘플이 제거되었음. 위 목록 확인 필요.")

for uid, target in TARGETS.items():
    nk, dk = norm_key(target), digit_key(target)
    in_mt = [s for s in mt_samples if norm_key(s) == nk or digit_key(s) == dk]
    in_file = [s for s in file_samples if norm_key(s) == nk or digit_key(s) == dk]
    print()
    print(f"[UID {uid}] 타깃 {target} (정규화 {nk} / 숫자 {dk})")
    print(f"  mt   내 존재: {in_mt}")
    print(f"  파일 내 존재: {in_file}")
    if in_file and not in_mt:
        print("  => 파일엔 있으나 mt에 없음: 분석 코드 필터에서 제외된 것")
    elif in_file and in_mt:
        print("  => 존재함. UID 부여 가능 (reconcile 스크립트 실행)")
    else:
        print("  => VCF 파일 자체에 없음. 시퀀싱 미실시 또는 센터의 "
              "상위 QC/VCF 생성 단계에서 제외된 것 (본 분석 코드와 무관)")
