import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# -----------------------
# 1. Population parameters
# -----------------------
pop_params = {
    'Vc': 0.25,      # L/kg
    'Cl': 3.5,       # L/hr
    'Q': 6.0,        # L/hr (intercompartmental clearance)
    'Vp': 18.0       # L (peripheral volume)
}
param_sd = {
    'Vc': 0.05,
    'Cl': 0.7,
    'Q': 1.0,
    'Vp': 4.0
}

# -----------------------
# 2. Sample observed data (time in hours, conc in mg/L)
# -----------------------
obs_df = pd.DataFrame({
    'time': [1.5, 4.0, 8.0],
    'conc': [25.2, 12.3, 6.1]
})
dose = 1000  # mg

# -----------------------
# 3. PK model: 2-compartment IV bolus
# -----------------------
def pk_conc(t, Vc, Cl, Q, Vp, dose):
    k10 = Cl / Vc
    k12 = Q / Vc
    k21 = Q / Vp

    A = 0.5 * ((k10 + k12 + k21) + np.sqrt((k10 + k12 + k21)**2 - 4 * k10 * k21))
    B = 0.5 * ((k10 + k12 + k21) - np.sqrt((k10 + k12 + k21)**2 - 4 * k10 * k21))
    alpha = A
    beta = B

    A1 = dose / Vc * (alpha - k21) / (alpha - beta)
    A2 = dose / Vc * (k21 - beta) / (alpha - beta)
    return A1 * np.exp(-alpha * t) + A2 * np.exp(-beta * t)

# -----------------------
# 4. Bayesian Objective Function
# -----------------------
def bayesian_objective(theta):
    Vc, Cl, Q, Vp = theta
    pred = pk_conc(obs_df['time'].values, Vc, Cl, Q, Vp, dose)
    error = np.sum((obs_df['conc'].values - pred) ** 2)

    # Bayesian penalty (prior)
    penalty = 0
    for i, name in enumerate(['Vc', 'Cl', 'Q', 'Vp']):
        penalty += ((theta[i] - pop_params[name]) / param_sd[name]) ** 2

    return error + penalty

# -----------------------
# 5. Optimization
# -----------------------
initial_guess = [pop_params[k] for k in ['Vc', 'Cl', 'Q', 'Vp']]
bounds = [(0.1, 5), (0.1, 20), (0.1, 20), (1, 50)]
result = minimize(bayesian_objective, initial_guess, bounds=bounds)

# -----------------------
# 6. Result summary
# -----------------------
estimated_params = result.x
for name, val in zip(['Vc', 'Cl', 'Q', 'Vp'], estimated_params):
    print(f"{name}: {val:.3f}")

# -----------------------
# 7. Plotting
# -----------------------
time_plot = np.linspace(0, 12, 100)
pred_plot = pk_conc(time_plot, *estimated_params, dose)

plt.plot(time_plot, pred_plot, label='Predicted')
plt.scatter(obs_df['time'], obs_df['conc'], color='red', label='Observed')
plt.xlabel("Time (hr)")
plt.ylabel("Concentration (mg/L)")
plt.legend()
plt.title("Bayesian Estimated PK Curve (Vancomycin)")
plt.grid(True)

plt.savefig("sine_wave.png")  # 파일로 저장
print("그래프가 'sine_wave.png'로 저장되었습니다.")
