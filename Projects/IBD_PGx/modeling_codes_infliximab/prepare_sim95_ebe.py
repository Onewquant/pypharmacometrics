"""재추정 완료 후 (run 89 = 최종모델 구조 재실행, 2026-08-29), EBE 기반 시뮬레이션(run 95) 자동 구성.

배경 (2026-08):
  기존 시뮬레이션(90.mod)은 `$SIM ... ONLYSIM` + `$OMEGA 0.0711`로
  개인 ETA를 난수로 새로 뽑았음 → sim90의 개인 CL이 환자 EBE와 무관
  (sim ETA vs 93.phi EBE 상관 -0.08). PGx 분석 endpoint로 부적합.

  run 95는 같은 구조(시변 공변량 반영)를 유지하되, 개인 ETA를
  재추정(run 94)의 EBE에서 데이터 컬럼(EBE1)으로 받아 CL을 계산한다:
    CL_record = TVCL(공변량_t) x exp(EBE1_i)
  CL의 OMEGA는 0 FIX → CL은 결정론적. (DV의 잔차 시뮬은 그대로 유지되나
  PGx CL 추출에는 사용되지 않음)

절차:
  1. NONMEM으로 89.mod(또는 94.mod) 실행 (재추정, POSTHOC) - 완료됨
  2. 이 스크립트 실행 → simulation_df_ebe.csv + 95.mod 생성
  3. NONMEM으로 95.mod 실행 → sim95 생성
  4. pd_endpoint_data_for_genomics2.py 실행 (sim95를 읽도록 수정됨)

실행:
  venv python으로:  python prepare_sim95_ebe.py
"""

import io
import re

import pandas as pd

nonmem_dir = "C:/Users/ilma0/NONMEMProjects/IBDPGX"
run_dir = f"{nonmem_dir}/run"

THETA_NAMES = [
    "CL", "V2", "V3", "Q", "KA", "F1", "Prop.RE (sd)", "Add.RE (sd)",
    "CL_ADA", "CL_SEX", "CL_WT", "V2_ALB", "V2_WT",
]
OLD_SIM_THETAS = [0.294, 4.49, 0.405, 0.0646, 0.0553, 0.667,
                  0.38, 0.0, 0.711, -0.107, 0.515, 1.42, 0.789]


def read_final_thetas(ext_path, n=13):
    """ext 파일의 final estimate 행(-1000000000)에서 THETA 추출."""
    for line in open(ext_path):
        parts = line.split()
        if parts and parts[0] == "-1000000000":
            return [float(x) for x in parts[1:n + 1]]
    raise RuntimeError(f"final estimates not found in {ext_path}")


def read_phi_eta1(phi_path):
    phi = pd.read_csv(phi_path, skiprows=1, sep=r"\s+", engine="python")
    phi.columns = [c.strip() for c in phi.columns]
    phi["ID"] = phi["ID"].astype(float).astype(int)
    return phi.set_index("ID")["ETA(1)"]


# ---------------------------------------------------------------- 1) inputs
EST_RUN = "89"  # 최종모델 재추정 run 번호 (89: 2026-08-29 재실행)
thetas = read_final_thetas(f"{run_dir}/{EST_RUN}.ext")
eta1 = read_phi_eta1(f"{run_dir}/{EST_RUN}.phi")

print(f"=== run {EST_RUN} 최종 THETA (기존 sim(90.mod) 값과 비교) ===")
for name, new, old in zip(THETA_NAMES, thetas, OLD_SIM_THETAS):
    d = "" if old == 0 else f"  ({(new - old) / abs(old) * 100:+.1f}%)"
    print(f"  {name:14s} {new:10.4f}   (구 sim: {old:g}){d}")
print(f"\nEBE ETA(1): {len(eta1)}명, 범위 {eta1.min():.3f} ~ {eta1.max():.3f}")

# ------------------------------------------- 2) EBE 컬럼 붙인 시뮬 데이터셋
sim_df = pd.read_csv(f"{nonmem_dir}/infliximab_integrated_simulation_df.csv")
sim_ids = set(sim_df["ID"].astype(int))
missing = sorted(sim_ids - set(eta1.index))
if missing:
    # 추정 데이터셋에 없는 환자 (maintenance 시작 + 관찰농도 1건뿐이라
    # 추정에서 제외된 환자; 예: ID 97 = UID 38339532). EBE가 존재하지
    # 않으므로 개인 CL을 부여할 수 없어 시뮬레이션에서도 제외한다.
    # -> PGx 코호트는 96명 (98 - 유전형QC 1 - CL추정불가 1)
    n_rows = len(sim_df[sim_df["ID"].astype(int).isin(missing)])
    print(f"[제외] 추정모델에 없는 시뮬 ID {missing} ({n_rows}행) - EBE 없음")
    sim_df = sim_df[~sim_df["ID"].astype(int).isin(missing)].copy()

sim_df["EBE1"] = sim_df["ID"].astype(int).map(eta1).round(6)
out_csv = f"{nonmem_dir}/infliximab_integrated_simulation_df_ebe.csv"
sim_df.to_csv(out_csv, index=False, encoding="utf-8")
print(f"\n[저장] {out_csv}  ({len(sim_df)}행, EBE1 컬럼 추가)")

# ---------------------------------------------------------------- 3) 95.mod
src = io.open(f"{run_dir}/90.mod", encoding="utf-8", errors="replace").read()
mod = src.replace("\r\n", "\n")

mod = mod.replace(";; 1. Based on: 89", ";; 1. Based on: " + EST_RUN)
mod = mod.replace(
    "Covar CL_V2_Q, SIM",
    f"Covar CL_V2_Q, SIM with EBE etas from run {EST_RUN} (weight-corrected)",
)
mod = mod.replace(
    "$INPUT ID TIME DV MDV AMT DUR CMT ROUTE IBD_TYPE ALB ADA AGE SEX "
    "WT HT BMI REALDATA RATE TAD DT_YEAR DT_MONTH DT_DAY",
    "$INPUT ID TIME DV MDV AMT DUR CMT ROUTE IBD_TYPE ALB ADA AGE SEX "
    "WT HT BMI REALDATA RATE TAD DT_YEAR DT_MONTH DT_DAY EBE1",
)
mod = mod.replace(
    "$DATA ../infliximab_integrated_simulation_df.csv  IGNORE=@",
    "$DATA ../infliximab_integrated_simulation_df_ebe.csv  IGNORE=@",
)
# 개인 CL: 난수 ETA(1) 대신 데이터로 받은 EBE 사용 (OMEGA(1,1)=0 FIX)
mod = mod.replace(
    "CL = TVCL * EXP(ETA(1))",
    f"CL = TVCL * EXP(EBE1 + ETA(1)) ; EBE1 = run {EST_RUN} posthoc eta, OMEGA(1,1)=0",
)
mod = mod.replace("$SIM (12345) (54321) ONLYSIM",
                  "$SIM (12345) (54321) ONLYSIM")

# THETA 블록 교체 (run 94 추정치로 FIX)
t = thetas
theta_block = (
    "$THETA\n"
    f"({t[0]:.6g}) FIX ;CL\n"
    f"({t[1]:.6g}) FIX ;V2\n"
    f"({t[2]:.6g}) FIX ;V3\n"
    f"({t[3]:.6g}) FIX ;Q\n"
    f"({t[4]:.6g}) FIX ;KA\n"
    f"({t[5]:.6g}) FIX ; F1\n"
    f"(0, {t[6]:.6g}) FIX ; Prop.RE (sd)\n"
    f"(0) FIX ; Add.RE (sd)\n"
    f"({t[8]:.6g}) FIX ; CL_ADA\n"
    f"({t[9]:.6g}) FIX ; CL_SEX\n"
    f"({t[10]:.6g}) FIX ; CL_WT\n"
    f"({t[11]:.6g}) FIX ; V2_ALB\n"
    f"({t[12]:.6g}) FIX ; V2_WT\n"
)
mod = re.sub(r"\$THETA\n(?:.*\n)+?\n\$OMEGA", theta_block + "\n$OMEGA", mod)

# OMEGA CL -> 0 FIX (개인차는 EBE1 데이터 컬럼으로 반영)
mod = mod.replace(" 0.0711 FIX  ;CL", " 0 FIX  ;CL (EBE1 as data column)")

# 출력 테이블
mod = mod.replace("FILE=sim90", "FILE=sim95")
mod = mod.replace("FILE=sdtab90", "FILE=sdtab95")
mod = mod.replace("FILE=cotab90", "FILE=cotab95")
mod = mod.replace("FILE=catab90", "FILE=catab95")
mod = mod.replace("FILE=patab90", "FILE=patab95")
# EBE1도 테이블에 포함 (검증용)
mod = mod.replace(
    "PRED IPRED CL V2 V3 Q KA F1 F2 ONEHEADER",
    "PRED IPRED CL V2 V3 Q KA F1 F2 EBE1 ONEHEADER",
)

io.open(f"{run_dir}/95.mod", "w", encoding="utf-8",
        newline="").write(mod.replace("\n", "\r\n"))
print(f"[저장] {run_dir}/95.mod")
print("\n다음: NONMEM으로 95.mod 실행 -> sim95 생성 후 "
      "pd_endpoint_data_for_genomics2.py 재실행")
