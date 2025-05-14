def heat_influx(self, L_in_max, Q_str,t1,emis_mli,k_vac,t_mli, k_mli):
        T_tank = self.T0
        T_amb = 300 # K
        P_amb = 101325 # Pa
        r_in = self.R_in
        t2 = t1
        Q_in_max = self.Q_leak_max
        k_1 = self.mat_property[2]
        k_2 = self.mat_property[2]
        eps1 = emis_mli
        eps2 = self.mat_property[3]

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

        t2 = max(t1, P_amb * (r_in+t1+dv) / self.mat_property[1])

        return dv, t2