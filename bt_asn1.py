import numpy as np
import pandas as pd
from scipy.integrate import simpson

data = pd.read_csv('cleaned_data.csv')  
tau_exp = data['tau'].values
variables_exp = data.drop(columns='tau')  
variable_names = list(variables_exp.columns)

tau_start = 1.4 #Change
tau_end = 6.8 #Change


def mmkinetics_sim(eta, eps, kappa, dtau, tau_points):
    Ptau = 0.0 #Change
    Xtau = 0.0 #Change
    Vtau = 0.0 #Change
    Stau = 1.0 #Change
    Etau = 1.0 #Change
    tau_sim = 0.0 #Change

    tau_vals = [0.0]
    S_vals = [Stau]
    P_vals = [Ptau]
    E_vals = [Etau]
    V_vals = [Vtau]
    X_vals = [Xtau]  

    while tau_sim < max(tau_points):

        dX = ((1 - Xtau) * (1 - eps * Xtau - Ptau) - (eta + kappa) * Xtau) / eta
        Xtau += dtau * dX
        Ptau += dtau * eps * Xtau
        Etau = 1 - Xtau
        Stau = 1 - eps * Xtau - Ptau
        if 'V' in variable_names:
            Vtau = eps * Xtau

        Xtau = max(Xtau, 0)
        Stau = max(Stau, 0)
        Ptau = min(max(Ptau, 0), 1)
        Etau = max(min(Etau, 1), 0)
        Vtau = max(Vtau, 0)

        tau_sim += dtau
        tau_vals.append(tau_sim)
        S_vals.append(Stau)
        P_vals.append(Ptau)
        E_vals.append(Etau)
        V_vals.append(Vtau)
        X_vals.append(Xtau)

    interp_dict = {}
    if 'S' in variable_names:
        interp_dict['S'] = np.interp(tau_points, tau_vals, S_vals)
    if 'P' in variable_names:
        interp_dict['P'] = np.interp(tau_points, tau_vals, P_vals)
    if 'E' in variable_names:
        interp_dict['E'] = np.interp(tau_points, tau_vals, E_vals)
    if 'V' in variable_names:
        interp_dict['V'] = np.interp(tau_points, tau_vals, V_vals)
    if 'X' in variable_names:
        interp_dict['X'] = np.interp(tau_points, tau_vals, X_vals)

    sim_array = np.stack([tau_points] + [interp_dict[name] for name in variable_names], axis=1)
    return sim_array


def loss(params):
    eta, eps, kappa = params
    sim = mmkinetics_sim(eta, eps, kappa, dtau=0.1, tau_points=tau_exp)
    err = 0.0
    for i, col in enumerate(variable_names, start=1):
        err += np.mean((sim[:, i] - variables_exp[col].values)**2)
    return err


scan_range = np.arange(0.1, 1.5 + 0.05, 0.05)
results = []
best_loss = np.inf
best_params = None

for eta in scan_range:
    for eps in scan_range:
        for kappa in scan_range:
            l = loss([eta, eps, kappa])
            results.append({'eta': eta, 'eps': eps, 'kappa': kappa, 'loss': l})
            if l < best_loss:
                best_loss = l
                best_params = (eta, eps, kappa)

print("Most probable parameters (min loss):")
print("eta =", best_params[0])
print("eps =", best_params[1])
print("kappa =", best_params[2])


scan_df = pd.DataFrame(results)
scan_df.to_csv('parameter_scan_results.csv', index=False)


def arc_length_auto(eta, eps, kappa, dtau, tau_start, tau_end):
    tau_dense = np.arange(tau_start, tau_end + dtau, dtau)
    sim = mmkinetics_sim(eta, eps, kappa, dtau, tau_dense)
    col_indices = range(1, 1 + len(variable_names)) 
    derivatives = [np.gradient(sim[:, i], tau_dense) for i in col_indices]
    integrand = np.sqrt(np.sum([d**2 for d in derivatives], axis=0))
    return simpson(integrand, x=tau_dense)

arc_len = arc_length_auto(best_params[0], best_params[1], best_params[2],
                          dtau=0.1, tau_start=tau_start, tau_end=tau_end)

print(f"Arc length from tau={tau_start} to tau={tau_end} (over CSV variables): {arc_len}")
