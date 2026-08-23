"""VCF 헤더에서 샘플 목록을 직접 읽어 특정 샘플의 존재를 확정 판정.

Hail/Spark를 거치지 않고 VCF의 #CHROM 라인(샘플명 나열)만 파싱하므로,
분석 코드의 필터와 무관한 '파일 원본 기준' 판정이 된다.
bgzip은 gzip 호환 스트림이라 gzip.open으로 헤더를 읽을 수 있다.

사용:
  python check_vcf_header_samples.py <VCF경로> [찾을샘플ID ...]
  (샘플ID 생략 시 기본값 23-B02281_EB-01 = UID 17439372)
"""

import gzip
import re
import sys

DEFAULT_TARGETS = ["23-B02281_EB-01"]


def norm_key(x):
    x = str(x).strip().upper()
    x = re.sub(r"^SAMPLE[_-]", "", x)
    x = re.sub(r"[_-]EB[_-]?\d+$", "", x)
    x = re.sub(r"\s+", "", x)
    return x


def digit_key(x):
    return re.sub(r"\D", "", norm_key(x))


def read_vcf_samples(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    return []


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return

    vcf_path = sys.argv[1]
    targets = sys.argv[2:] or DEFAULT_TARGETS

    samples = read_vcf_samples(vcf_path)
    print(f"VCF: {vcf_path}")
    print(f"헤더 기준 샘플 수: {len(samples)}")

    if not samples:
        print("[오류] #CHROM 라인을 찾지 못했습니다.")
        return

    print(f"샘플 예시: {samples[:5]}")

    norm_set = {norm_key(s): s for s in samples}
    digit_set = {digit_key(s): s for s in samples}

    for t in targets:
        nk, dk = norm_key(t), digit_key(t)
        print()
        print(f"[찾는 샘플] {t}  (정규화 {nk} / 숫자 {dk})")
        print(f"  완전일치   : {[s for s in samples if s == t]}")
        print(f"  정규화일치 : {[norm_set[nk]] if nk in norm_set else []}")
        print(f"  숫자일치   : {[digit_set[dk]] if dk in digit_set else []}")
        part = [s for s in samples if dk[-5:] and dk[-5:] in digit_key(s)]
        print(f"  부분일치   : {part}")
        found = (nk in norm_set) or (dk in digit_set)
        print("  => " + ("VCF에 존재함" if found
                         else "VCF 파일 자체에 없음 (분석 코드 필터와 무관)"))


if __name__ == "__main__":
    main()
