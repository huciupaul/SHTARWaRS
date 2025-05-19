import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve



t1 = 0.0028 # m
t_mli = 50 * 1e-3 # m
L_in_max = 0.5 # m
k_mli = 0.035 # W/mK
k_vac = 10.44 * 1e-3 # W/mK
emis_mli = 0.03 # emissivity of MLI
Q_str = 0.4 # W
P0 = 101000 # Pa
T_tank = CP.PropsSI('T', 'P', P0, 'Q', 0.1, 'ParaHydrogen')
T_amb = 300 # K
P_amb = 101325 # Pa
r_in = 0.8
t2 = t1
Q_in_max = 294.5
k_1 = 134
k_2 = 134
eps2 = 0.048
SF=0.9
# Conduction...
def Q_cond(dv):
    # Conduction resistance
    R_cond = 0.0
    R_cond = np.log((r_in + t1) / r_in) / (2 * np.pi * L_in_max * k_1)
    R_cond += np.log((r_in + t1 + t_mli) / (r_in + t1)) / (2 * np.pi * L_in_max * k_mli)
    R_cond += np.log((r_in + t1 + t_mli + dv) / (r_in + t1 + t_mli)) / (2 * np.pi * L_in_max * k_vac)
    R_cond += np.log((r_in + t1 + dv + t_mli + t2) / (r_in + t1 + dv + t_mli)) / (2 * np.pi * L_in_max * k_2)
    return (T_amb - T_tank) / R_cond

# Radiation...
def Q_rad(dv):
    r1 = r_in + t1 + t_mli        # Outer surface of inner tank
    r2 = r1 + dv          # Inner surface of outer tank

    # Surface areas (cylinder + two hemispherical caps)
    A1 = 2 * np.pi * r1 * (L_in_max + 4 * r1 / 3)
    A2 = 2 * np.pi * r2 * (L_in_max + 4 * r2 / 3)

    # Radiation heat transfer
    denom = (1 / emis_mli) + (A1 / A2) * (1 / eps2 - 1)
    return 5.670374419e-8 * A1 * (T_amb**4 - T_tank**4) / denom

# Total heat influx...
def total_heat_influx(dv):
    Q_cond_value = Q_cond(dv)
    Q_rad_value = Q_rad(dv)
    return Q_cond_value + Q_rad_value + Q_str

# Optimization eq...
def equation(dv):
    return total_heat_influx(dv) - Q_in_max

# Initial guess for dv
dv_initial_guess = 0.01  # Initial guess for the vacuum gap thickness (m)

# Solve for dv...
dv_solution = fsolve(equation, dv_initial_guess)

dv = dv_solution[0]  # Extract the solution from the array


t2_min = P_amb * (r_in+t1+dv+t_mli) / 495*1e6*SF

print(t1, dv, t2)
print(Q_cond(dv), Q_rad(dv), total_heat_influx(dv))
