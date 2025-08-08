
"""
Fit a latent K-PD + indirect response model (marrow suppression) for linezolid using dose & lab data.

Model:
  Latent exposure x(t):
    dx/dt = -k * x,  with instantaneous jumps at dose times: x(t_dose+) = x(t_dose-) + alpha * Dose

  Platelet turnover R(t):
    dR/dt = kin * (1 - Imax * x/(IC50 + x)) - kout * R
  Baseline link: kin = kout * R0

Assumptions:
- Doses are instantaneous boluses (reasonable for short infusions vs sampling granularity).
- Use per-subject non-linear least squares to get quick, interpretable estimates (no pooling).
- You can later move to hierarchical Bayes / nlmixr2 / NONMEM.

Required input CSVs:
- /mnt/data/lnz_final_dose_df.csv (columns: UID, DATE, TIME, DOSE, ...)
- /mnt/data/lnz_final_lab_df.csv  (columns: UID, DATE, PLT, ...)

Outputs:
- /mnt/data/lnz_kpd_indirect_platelet_params.csv : parameter estimates per UID
- /mnt/data/lnz_kpd_indirect_platelet_diagnostics.csv : fit stats
"""

import pandas as pd
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp

# -------------------- Config --------------------
prj_dir = './Projects/LINEZOLID'
output_dir = f"{prj_dir}/results"
DOSE_CSV = f"{output_dir}/lnz_final_dose_df.csv"
LAB_CSV  = f"{output_dir}/lnz_final_lab_df.csv"
LAB_NAME = "PLT"   # change to ANC/WBC if preferred

# Time unit: days
# Convert datetime to days from subject's first dose or first observation (whichever earlier)
TIME_UNIT = "D"

# Bounds and initial guesses
INIT = dict(
    k=0.3,          # 1/day  (effective elimination rate; half-life ~ ln2/k)
    alpha=0.01,     # exposure scale per mg
    Imax=0.6,       # max inhibition fraction
    IC50=0.5,       # exposure units
    kout=0.1,       # 1/day  (platelet turnover rate; half-life ~ 6.9 days if kout=0.1)
    R0=None         # baseline platelet (estimated from first labs if None)
)
BOUNDS = dict(
    k=(1e-4, 3.0),
    alpha=(1e-6, 10.0),
    Imax=(0.0, 0.99),
    IC50=(1e-4, 100.0),
    kout=(1e-3, 2.0),
    R0=(10.0, 1000.0),
)

# Integration settings
ATOL = 1e-6
RTOL = 1e-6
MAX_STEP = 0.25  # days

# ------------------------------------------------

def _parse_datetime(df, date_col="DATE", time_col=None):
    if time_col and time_col in df.columns:
        # DATE assumed like YYYY-MM-DD, TIME like HH:MM (or HH:MM:SS)
        dt = pd.to_datetime(df[date_col].astype(str) + " " + df[time_col].astype(str), errors="coerce")
    else:
        dt = pd.to_datetime(df[date_col], errors="coerce")
    return dt

def _make_subject_times(dose_sub, lab_sub):
    # t0 is earliest of first dose or first lab
    t0 = min(dose_sub["DT"].min(), lab_sub["DT"].min())
    dose_sub = dose_sub.assign(t=(dose_sub["DT"] - t0) / np.timedelta64(1, "D"))
    lab_sub  = lab_sub.assign(t=(lab_sub["DT"] - t0) / np.timedelta64(1, "D"))
    return dose_sub.sort_values("t"), lab_sub.sort_values("t")

def build_exposure_fun(dose_times: np.ndarray, dose_amts: np.ndarray, k: float, alpha: float):
    # Return a function that integrates x(t) between events; use piecewise solution with jumps.
    # Solution between events: x(t) = x(t_i) * exp(-k*(t - t_i)), and at dose: x += alpha * Dose.
    events = list(zip(dose_times, dose_amts))
    events.sort(key=lambda z: z[0])

    def x_of_t(query_times: np.ndarray):
        # Compute x(t) at arbitrary query times by stepping through events.
        # For efficiency, vectorized over sorted unique times.
        qs = np.asarray(query_times).copy()
        order = np.argsort(qs)
        qs_sorted = qs[order]

        x = 0.0
        xi = 0.0
        last_t = 0.0
        eidx = 0
        out = np.empty_like(qs_sorted, dtype=float)

        # Merge event times and query times timeline
        # process sorted times; at each step, decay from last_t, then if event at that time, add jump
        evt_times = [t for t, a in events]
        evt_dict = {}
        for t, a in events:
            evt_dict.setdefault(t, 0.0)
            evt_dict[t] += alpha * a

        for i, t in enumerate(qs_sorted):
            # progress event pointer up to current time
            # But we need to process all events between last_t and t
            # We'll iterate over all event times in (last_t, t], applying decay to the event time then jump
            # For performance, we could pre-index; here clarity first.
            mid_events = [et for et in evt_times if (et > last_t and et <= t)]
            # decay to each mid-event and apply jump
            for et in mid_events:
                xi *= np.exp(-k * (et - last_t))
                xi += evt_dict.get(et, 0.0)
                last_t = et
            # finally decay to t
            xi *= np.exp(-k * (t - last_t))
            last_t = t
            out[i] = xi

        # unsort to original order
        out_unsorted = np.empty_like(out)
        out_unsorted[order] = out
        return out_unsorted

    return x_of_t

def simulate_platelet(times: np.ndarray, x_t: np.ndarray, params: Dict[str, float], R0_guess: float):
    # Integrate dR/dt = kin*(1 - Imax*x/(IC50 + x)) - kout*R with kin=kout*R0
    Imax = params["Imax"]; IC50 = params["IC50"]; kout = params["kout"]
    R0 = params.get("R0", None) or R0_guess
    kin = kout * R0

    def rhs(t, y):
        # y[0] = R
        # interpolate x(t) by linear interpolation on provided time grid
        # assume times are sorted; np.interp does 1D linear
        x = np.interp(t, times, x_t)
        inhib = Imax * x / (IC50 + x)
        dR = kin * (1.0 - inhib) - kout * y[0]
        return [dR]

    sol = solve_ivp(rhs, t_span=(times.min(), times.max()), y0=[R0], t_eval=times,
                    method="RK45", atol=ATOL, rtol=RTOL, max_step=MAX_STEP)
    return sol.y[0]

def fit_subject(dose_sub: pd.DataFrame, lab_sub: pd.DataFrame, uid: str) -> Tuple[Dict[str, float], Dict[str, float]]:
    # Prepare arrays
    dose_times = dose_sub["t"].to_numpy(dtype=float)
    dose_amts  = dose_sub["DOSE"].to_numpy(dtype=float)
    lab_times  = lab_sub["t"].to_numpy(dtype=float)
    lab_vals   = lab_sub[LAB_NAME].to_numpy(dtype=float)

    # Initial R0: median of first 3 observations (or first)
    R0_init = float(np.nanmedian(lab_vals[:min(3, len(lab_vals))])) if len(lab_vals)>0 else 150.0

    # Parameter vector: [k, alpha, Imax, IC50, kout, R0]
    p0 = np.array([
        INIT["k"],
        INIT["alpha"],
        INIT["Imax"],
        INIT["IC50"],
        INIT["kout"],
        R0_init if INIT["R0"] is None else INIT["R0"]
    ])

    lb = np.array([BOUNDS["k"][0], BOUNDS["alpha"][0], BOUNDS["Imax"][0], BOUNDS["IC50"][0], BOUNDS["kout"][0], BOUNDS["R0"][0]])
    ub = np.array([BOUNDS["k"][1], BOUNDS["alpha"][1], BOUNDS["Imax"][1], BOUNDS["IC50"][1], BOUNDS["kout"][1], BOUNDS["R0"][1]])

    # Residual function
    def resid(p):
        k, alpha, Imax, IC50, kout, R0 = p
        x_fun = build_exposure_fun(dose_times, dose_amts, k=k, alpha=alpha)
        x_grid = x_fun(lab_times)
        y_hat  = simulate_platelet(lab_times, x_grid, {"Imax":Imax, "IC50":IC50, "kout":kout, "R0":R0}, R0_guess=R0)
        # proportional + additive error weight (rough)
        w = 1.0 / np.maximum(1.0, np.abs(lab_vals))
        return (y_hat - lab_vals) * w

    res = least_squares(resid, p0, bounds=(lb, ub), xtol=1e-6, ftol=1e-6, gtol=1e-6, max_nfev=1000, verbose=0)

    k, alpha, Imax, IC50, kout, R0 = res.x
    # Derived half-lives
    t12_eff = np.log(2.0)/k
    t12_pl  = np.log(2.0)/kout

    params = dict(UID=uid, k=k, alpha=alpha, Imax=Imax, IC50=IC50, kout=kout, R0=R0, t12_eff=t12_eff, t12_pl=t12_pl)
    diag = dict(UID=uid, n_obs=int(len(lab_vals)), n_dose=int(len(dose_times)), cost=float(0.5*np.sum(res.fun**2)), success=bool(res.success))

    return params, diag

def main():
    dose_df = pd.read_csv(DOSE_CSV)
    lab_df  = pd.read_csv(LAB_CSV)

    # Parse datetime
    dose_df = dose_df.copy()
    lab_df  = lab_df.copy()

    dose_df["DT"] = pd.to_datetime(dose_df["DATE"].astype(str) + " " + dose_df.get("TIME", pd.Series(["00:00"]*len(dose_df))).astype(str), errors="coerce")
    lab_df["DT"]  = pd.to_datetime(lab_df["DATE"], errors="coerce")

    # Filter/clean
    dose_df = dose_df.dropna(subset=["UID","DT","DOSE"])
    lab_df  = lab_df.dropna(subset=["UID","DT", LAB_NAME])

    # Work per UID present in both
    uids = sorted(set(dose_df["UID"]).intersection(set(lab_df["UID"])))

    all_params = []
    all_diag = []

    for uid in uids:
        dsub = dose_df.loc[dose_df["UID"]==uid].copy()
        lsub = lab_df.loc[lab_df["UID"]==uid].copy()

        # require at least 1 dose and 2 labs
        if len(dsub)==0 or len(lsub)<2:
            continue

        dsub, lsub = _make_subject_times(dsub, lsub)

        try:
            p, dg = fit_subject(dsub, lsub, uid=str(uid))
            all_params.append(p)
            all_diag.append(dg)
        except Exception as e:
            all_diag.append(dict(UID=str(uid), n_obs=int(len(lsub)), n_dose=int(len(dsub)), cost=np.nan, success=False, error=str(e)))

    params_df = pd.DataFrame(all_params)
    diag_df   = pd.DataFrame(all_diag)

    params_out = f"{output_dir}/lnz_kpd_indirect_platelet_params.csv"
    diag_out   = f"{output_dir}/lnz_kpd_indirect_platelet_diagnostics.csv"
    params_df.to_csv(params_out, index=False)
    diag_df.to_csv(diag_out, index=False)

    print(f"Saved parameter estimates to: {params_out}")
    print(f"Saved diagnostics to: {diag_out}")
    if len(params_df):
        print(params_df.describe(include='all'))


main()
