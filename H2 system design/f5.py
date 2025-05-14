import numpy as np
from scipy.optimize import fsolve

def heat_influx(Q_in_max, L_in_max, t1, T_amb, T_tank, r_in, k_1, k_vac, k_2, eps1, eps2, Q_structure):
    t2 = t1
    # Conduction...
    def Q_cond(dv):
        # Conduction resistance
        R_cond = 0.0
        R_cond = np.log((r_in + t1) / r_in) / (2 * np.pi * L_in_max * k_1)
        R_cond += np.log((r_in + t1 + dv) / (r_in + t1)) / (2 * np.pi * L_in_max * k_vac)
        R_cond += np.log((r_in + t1 + dv + t2) / (r_in + t1 + dv)) / (2 * np.pi * L_in_max * k_2)
        return (T_amb - T_tank) / R_cond

    # Radiation...
    def Q_rad(dv):
        r1 = r_in + t1        # Outer surface of inner tank
        r2 = r1 + dv          # Inner surface of outer tank

        # Surface areas (cylinder + two hemispherical caps)
        A1 = 2 * np.pi * r1 * (L_in_max + 4 * r1 / 3)
        A2 = 2 * np.pi * r2 * (L_in_max + 4 * r2 / 3)

        # Radiation heat transfer
        denom = (1 / eps1) + (A1 / A2) * (1 / eps2 - 1)
        return 5.670374419e-8 * A1 * (T_amb**4 - T_tank**4) / denom

    # Total heat influx...
    def total_heat_influx(dv):
        Q_cond_value = Q_cond(dv)
        Q_rad_value = Q_rad(dv)
        return Q_cond_value + Q_rad_value + Q_structure

    # Optimization eq...
    def equation(dv):
        return total_heat_influx(dv) - Q_in_max

    # Initial guess for dv
    dv_initial_guess = 0.01  # Initial guess for the vacuum gap thickness (m)

    # Solve for dv...
    dv_solution = fsolve(equation, dv_initial_guess)

    dv_sol = dv_solution[0]  # Extract the solution from the array

    return dv_sol