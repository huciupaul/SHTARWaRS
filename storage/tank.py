import CoolProp.CoolProp as CP # Used to get properties of Hydrogen
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy.optimize as opt
from scipy.optimize import fsolve
import csv
import ast
import matplotlib.colors as mcolors
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import global_constants
from global_constants import * # Import global constants

# Constants
# dormancy = 24       #[hours]
# stratification_factor = 2
# k_str = 1.9         #[W/mK]
# MAWP_global = 600000       #[Pa]
# p_vent_global = 600000     #[Pa]
# P_sl = 101325             #[Pa]
# r_in = 0.75                #[m]
# Q_leak_min = 10 #W (determined from Nicolas's thesis)
# Q_leak_max = 1000 #W (determined from Nicolas's thesis)
# T_hot = 40      #[Celsius]
# gf = 'S-Glass fiber'
# fos = 2/3
# gf_density, gf_tensile_strength, gf_thermal_cond, gf_thermal_emis, gf_co2, gf_ee, gf_fvf = 1905,1730*1e6*fos,0.745,0.95,7.22,116.5,0.675
# Thesis values
# Q_original_str = 0.4 #W
# kevlar_thermal_cond = 1.9 #W/mK
# mass_original_str = 2.1 #kg
# mass_originalg_lh2 = 6.2 #kg
# gravimetric_index = 0.35 #from NASA report
# kevlar_co2 = 13.1 #kg/kg (Kevlar 149)
# kevlar_emb_energy = 257 #MJ/kg (Embodied Energy for Kevlar 149)
# t_min = 0.001 #m (minimum thickness of the tank wall)
# mli_density = 7900 #kg/m^3 https://www.sciencedirect.com/science/article/pii/S135943112200391X
# mli_emis = 0.21  #https://www.thermalengineer.com/library/effective_emittance.htm
# vacuum_thermal_cond = 0.015*1e-1#3 # W/mK https://www.researchgate.net/publication/321219004_Cylindrical_Cryogenic_Calorimeter_Testing_of_Six_Types_of_Multilayer_Insulation_Systems
# mli_thermal_cond = 17.4 # W/mK  https://www.sciencedirect.com/science/article/pii/S135943112200391X
# mli_ss_co2 = 3 #kg/kg (for SS)
# mli_ss_ee = 42.74 #MJ/kg (Embodied Energy for SS)
# mli_layers = 40
# mli_thickness = 0.03 *1e-3 * mli_layers #https://www.sciencedirect.com/science/article/pii/S135943112200391X



def main_storage(m_h2):
    
    class Tank:
        def __init__(self, MAWP, material, material2, mat_property,mass_h2,mat2_property,fill,V_in,p_vent):
            # ---------------------------------------------Database ------------------------------------------------------------------
            self.material = material
            self.material2 = material2
            self.mat_property = mat_property
            self.mat_density = mat_property[0] 
            self.mat_yield_strength = mat_property[1]
            self.mat_thermal_conductivity = mat_property[2]
            self.mat_emissivity = mat_property[3]
            self.mat_co2 = mat_property[4]
            self.mat_ee = mat_property[5]
            self.mat_fibre_ratio = mat_property[6]

            self.mat2_property = mat2_property
            self.mat2_density = mat2_property[0]
            self.mat2_yield_strength = mat2_property[1]
            self.mat2_thermal_conductivity = mat2_property[2]
            self.mat2_emissivity = mat2_property[3]
            self.mat2_co2 = mat2_property[4]
            self.mat2_ee = mat2_property[5]
            self.mat2_fibre_ratio = mat2_property[6]

            # ---------------------------------------------Constants ----------------------------------------------------------------
            self.dormancy = dormancy #hours
            self.fill_ratio = fill #percentage of LH2 at the start
            self.mass_h2 = mass_h2 #kg
            self.stratification_factor = stratification_factor
            self.k_str = k_str #[W/mK]


            # --------------------------------------------- Inputs ------------------------------------------------------------------

            self.MAWP = MAWP # Convert from Bar to Pa for the venting pressure
            self.Pvent = p_vent #Assumption

            # Hydrogen properties
            self.P0 = P_sl # Initial pressure in Pa
            self.T0 = CP.PropsSI('T', 'P', self.P0, 'Q', 0.1, 'ParaHydrogen') # Initial temperature in K
            self.rhol0 = CP.PropsSI('D', 'T', self.T0, 'Q', 0, 'ParaHydrogen') # density of liquid hydrogen
            self.rhog0 = CP.PropsSI('D', 'T', self.T0, 'Q', 1, 'ParaHydrogen') # density of gas hydrogen
            self.V_in = V_in
            self.mg0 = CP.PropsSI('D', 'T', self.T0, 'Q', 1, 'ParaHydrogen')*self.V_in*(1-self.fill_ratio) # kg
            self.ml0 = CP.PropsSI('D', 'T', self.T0, 'Q', 0, 'ParaHydrogen')*self.V_in*(self.fill_ratio) # kg
            self.hg0 = CP.PropsSI('H', 'T', self.T0, 'Q', 1, 'ParaHydrogen')#kJ/kg
            self.hl0 = CP.PropsSI('H', 'T', self.T0, 'Q', 0, 'ParaHydrogen')#kJ/kg
            self.x0 = self.mg0/(self.mg0+self.ml0) #unitless (vapor quality)

            # ---------------------------------------------- Constraint ---------------------------------------------------
            self.R_in = r_in #m
            self.Q_leak_min = Q_leak_min #W (determined from Nicolas's thesis)
            self.Q_leak_max = Q_leak_max #W (determined from Nicolas's thesis)
            

        # ----------------------------------------------- Maximum Heat Load ---------------------------------------------------
        def maximum_Qin(self,Qleak):
            P = [self.P0]
            T = [self.T0]
            mg = [self.mg0]
            ml = [self.ml0]
            rhog = [self.rhog0]
            rhol = [self.rhol0]
            hg = [self.hg0]
            hl = [self.hl0]
            x = [self.x0]
            Ps = [self.P0]
            time = [0]
            dt = 250 # Time step in seconds
            # Run the simulation until the pressure reaches the venting pressure [code section extracted from Nicolas's thesis]
            while P[-1]<self.Pvent:
                # Compute derivatives
                delrhog_delP = CP.PropsSI('d(D)/d(P)|T', 'T', T[-1], 'Q', 1, 'ParaHydrogen')
                delrhol_delP = CP.PropsSI('d(D)/d(P)|T', 'T', T[-1], 'Q', 0, 'ParaHydrogen')
                delrhog_delT = CP.PropsSI('d(D)/d(T)|P', 'P', P[-1], 'Q', 1, 'ParaHydrogen')
                delrhol_delT = CP.PropsSI('d(D)/d(T)|P', 'P', P[-1], 'Q', 0, 'ParaHydrogen')
                delhg_delP = CP.PropsSI('d(H)/d(P)|T', 'T', T[-1], 'Q', 1, 'ParaHydrogen')
                delhl_delP = CP.PropsSI('d(H)/d(P)|T', 'T', T[-1], 'Q', 0, 'ParaHydrogen')
                delhg_delT = CP.PropsSI('d(H)/d(T)|P', 'P', P[-1], 'Q', 1, 'ParaHydrogen')
                delhl_delT = CP.PropsSI('d(H)/d(T)|P', 'P', P[-1], 'Q', 0, 'ParaHydrogen')
                # Small temperature increment to calculate the rate of pressure change with temperature
                T_increment = 6.5e-06
                T_new = T[-1]+T_increment
                P_new = CP.PropsSI('P', 'T', T_new, 'Q', 0.1, 'ParaHydrogen')
                dPs_dT = (P_new-P[-1])/(T_new-T[-1])
                # Construct the coefficient matrix A
                A21 = -(((mg[-1]/(rhog[-1]**2))*delrhog_delP)+
                ((ml[-1]/(rhol[-1]**2))*delrhol_delP))
                A22 = -(((mg[-1]/(rhog[-1]**2))*delrhog_delT)+
                ((ml[-1]/(rhol[-1]**2))*delrhol_delT))
                A42 = (ml[-1]*delhl_delT)+(mg[-1]*delhg_delT)\
                +(((ml[-1]*delhl_delP)+(mg[-1]*delhg_delP)-self.V_in)*dPs_dT)
                A = np.array([[1, -dPs_dT, 0, 0],
                [A21, A22, (1/rhog[-1])-(1/rhol[-1]), 0],
                [0, 0, 1, 1],
                [0, A42, -(hl[-1]-hg[-1]), 0]])
                # Construct the B matrix
                B = np.array([[0],
                [0],
                [0],
                [Qleak]])
                # Solve the system of linear equations
                differential_matrix = np.linalg.solve(A, B)
                dP_dt = differential_matrix[0, 0]*self.stratification_factor
                dT_dt = differential_matrix[1, 0]
                dmg_dt = differential_matrix[2, 0]
                dml_dt = differential_matrix[3, 0]
                # Update the state variables
                P.append(P[-1]+(dt*dP_dt))
                T.append(CP.PropsSI('T', 'P', P[-1], 'Q', 0.1, 'ParaHydrogen'))
                if mg[-1]+(dt*dmg_dt)<=0:
                    mg.append(0)
                else:
                    mg.append(mg[-1]+(dt*dmg_dt))
                    ml.append(ml[-1]+(dt*dml_dt))
                    x.append(mg[-1]/(mg[-1]+ml[-1]))
                    rhog.append(CP.PropsSI('D', 'P', P[-1], 'Q', 1, 'ParaHydrogen'))
                    rhol.append(CP.PropsSI('D', 'P', P[-1], 'Q', 0, 'ParaHydrogen'))
                    hg.append(CP.PropsSI('H', 'P', P[-1], 'Q', 1, 'ParaHydrogen'))
                    hl.append(CP.PropsSI('H', 'P', P[-1], 'Q', 0, 'ParaHydrogen'))
                    Ps.append(CP.PropsSI('P', 'T', T[-1], 'Q', 0.5, 'ParaHydrogen'))
                    time.append(time[-1]+dt)
            # Convert time to hours and pressure to bars for plotting
            time_hours = np.array(time)/3600
            P_bars = np.array(P)/100000 # Convert pressure from Pa to bars
            return time_hours[-1] - self.dormancy

        # ----------------------------------------------- Inner tank Parameters ---------------------------------------------------
        def volume_equation(self,l):
            L_cyl = l - 2*self.R_in  # cylindrical length
            if L_cyl < 0:
                return 1e6  # Penalize invalid cases
            V = np.pi * self.R_in**2 * L_cyl + (4/3) * np.pi * self.R_in**3  # cylinder + 2 hemispheres = 1 sphere
            return V - self.V_in  # we want this to be zero

        
        def inner_tank_thickness(self):
            P_test = (self.MAWP) * 1.5 #Safety factor
            alpha = 0
            if self.material == 'S-Glass Fiber' or self.material == 'Carbon Fiber':
                #Calc membrane loads
                N_theta_p = P_test * self.R_in
                N_phi_p = P_test * self.R_in / 2
                self.mat_property[1] = self.mat_property[1] #may add safety factors later
                angle = np.arange(1, 54, 0.1) #helical winding angle in degrees
                t_min = 1000
                #thicknesses = []
                #angles_deg = []
                for ang in angle:
                    # Calc netting thickness in hoop and helical direction
                    t_heli = N_phi_p / (self.mat_property[1] * np.cos(ang*np.pi/180)**2)
                    t_hoop = (N_theta_p - N_phi_p * np.tan(ang*np.pi/180)**2) / (self.mat_property[1])
                    # Calc minimum thickness based on FVF
                    t = (t_heli + t_hoop) / self.mat_property[6]
                    #thicknesses.append(t)
                    #angles_deg.append(np.rad2deg(ang))
                    if t < t_min:
                        t_min = t
                        alpha = ang
                '''
                # Plot angle vs thickness
                plt.figure()
                plt.plot(angles_deg, thicknesses, label='Thickness vs Angle')
                plt.xlabel('Helical Winding Angle (deg)')
                plt.ylabel('Minimum Thickness (m)')
                plt.title('Angle vs Minimum Thickness for Composite Tank')
                plt.grid(True)
                plt.legend()
                plt.show()
                '''
                return t_min, alpha
            else:
                t = (P_test * self.R_in / self.mat_property[1]) * np.sqrt(3)/2
                return t, alpha

        # ------------------------------------------------ Vacuum and Outer tank  ---------------------------------------------------

        def heat_influx(self, L_in_max, Q_str,t1,emis_mli,k_vac,t_mli, k_mli,Qmax,n_mli,t2):
            T_tank = self.T0
            T_amb = 273.15 + T_hot # K
            P_amb = P_sl # Pa
            r_in = self.R_in
            Q_in_max = Qmax -Q_str
            k_1 = self.mat_property[2]
            k_2 = self.mat2_property[2]
            eps1 = self.mat_property[3]
            eps2 = self.mat2_property[3]
            if n_mli == 0:
                t_mli = 0
            # Conduction...
            def Q_cond(dv):
                # Conduction resistance
                R_cond = 0.0
                R_cond = np.log((r_in + t1) / r_in) / (2 * np.pi * L_in_max * k_1)
                R_cond += np.log((r_in + t1 + t_mli + dv) / (r_in + t1 + t_mli)) / (2 * np.pi * L_in_max * k_vac)
                R_cond += np.log((r_in + t1 + dv + t_mli + t2) / (r_in + t1 + dv + t_mli)) / (2 * np.pi * L_in_max * k_2)
                if n_mli != 0:
                    R_cond += np.log((r_in + t1 + t_mli) / (r_in + t1)) / (2 * np.pi * L_in_max * k_mli) 
                
                

                return (T_amb - T_tank) / R_cond

            # Radiation...
            def Q_rad(dv):
                r1 = r_in + t1 + t_mli        # Outer surface of inner tank
                r2 = r1 + dv          # Inner surface of outer tank

                # Surface areas (cylinder + two hemispherical caps)
                A1 = 2 * np.pi * r1 * (L_in_max + 2 * r1)
                A2 = 2 * np.pi * r2 * (L_in_max + 2 * r2)

                # Radiation heat transfer
                
                #denom = (1 / (eps2)) + (A2 / A1) * (1 / emis_mli - 1)
                #return 5.670374419e-8 * A2 * (T_amb**4 - T_tank**4) / denom
                #return num * A2 / denom
                #
                A_rat = A1/A2
                if n_mli == 0:
                    B21 = A_rat * eps1 / (1-(1-eps1)*(1-eps2)*A_rat-(1-eps2)*(1-A_rat))
                    return (5.670374419e-8 * eps2 * B21 * (T_amb**4 - T_tank**4))
                else:
                    #return (5.670374419e-8 * (T_amb**4 - T_tank**4))/((n_mli+1)/emis_mli) * A1
                    e_mli = (1 / (2 / emis_mli - 1)) * (1 / (N_MLI+1))
                    B21 = A_rat * e_mli / (1-(1-e_mli)*(1-eps2)*A_rat-(1-eps2)*(1-A_rat))
                    return (5.670374419e-8 * eps2 * B21 * (T_amb**4 - T_tank**4))

            def total_heat_influx(dv):
                Q_cond_value = Q_cond(dv)
                Q_rad_value = Q_rad(dv)
                return Q_cond_value + Q_rad_value + Q_str

            def equation(dv):
                return total_heat_influx(dv) - Q_in_max

            dv_initial_guess = 0.1  # Initial guess for the vacuum gap thickness (m)

            dv_solution = fsolve(equation, dv_initial_guess)

            dv = dv_solution[0]  # Extract the solution from the array
            t2_min = 1000
            P_test = (P_amb)*1.5 #Safety factor
            alpha = 0
            if self.material2 == 'S-Glass Fiber' or self.material2 == 'Carbon Fiber':
                #Calc membrane loads
                N_theta_p = P_test * (r_in+t1+dv+t_mli)
                N_phi_p = P_test * (r_in+t1+dv+t_mli) / 2
                #self.mat2_property[1] = self.mat2_property[1] #may add safety factors
                angle = np.arange(0, 54, 0.1) # helical winding angle in radians
                
                for ang in angle:
                    #Calc netting thickness in hoop and helical direction
                    t_heli = N_phi_p / (self.mat2_property[1] * np.cos(ang*np.pi/180)**2)
                    t_hoop = (N_theta_p - N_phi_p * np.tan(ang*np.pi/180)**2) / (self.mat2_property[1])
                    #Calc minimum thickness based on FVF
                    t = (t_heli + t_hoop) / self.mat2_property[6]
                    if t < t2_min:
                        t2_min = t
                        alpha = ang
                t2_final = t2_min
            else:
                t2_final = P_amb * (r_in+t1+dv+t_mli) / self.mat2_property[1] * np.sqrt(3)/2
            
            t2 = t2_final

            return dv, t2, alpha, Q_in_max, Q_cond(dv),  Q_rad(dv)

        # ------------------------------------------------ Tank Dimensions ---------------------------------------------------
        def total_volume(self, l,dv,t1,t2,t_mli):
            R_out = self.R_in + dv +t1 +t2 +t_mli
            L_out = 2*R_out + l
            V = np.pi * R_out**2 * (l) + (4/3) * np.pi * R_out**3
            return V, L_out,R_out
        
        def total_mass(self, l, dv, t1, t2,t_mli,dens_mli):
            R_out = self.R_in + dv +t1 +t2 +t_mli
            L_cyl = l - 2 * self.R_in  
            surface_inner = 2 * np.pi * self.R_in * L_cyl + 4 * np.pi * self.R_in**2
            surface_outer = 2 * np.pi * R_out * L_cyl + 4 * np.pi * R_out**2
            mass_inner = surface_inner * t1 * self.mat_density
            mass_outer = surface_outer * t2 * self.mat2_density
            mass_mli = surface_inner * t_mli * dens_mli
            return mass_inner,mass_outer,mass_mli

        def kg_co2(self,mass_inner,mass_outer,mass_mli,mass_str,mli_co2,co2_kevlar):
            lca_inner = mass_inner * self.mat_property[4]
            lca_outer = mass_outer * self.mat2_property[4]
            lca_mli = mass_mli * mli_co2
            lca_str = mass_str * co2_kevlar
            lca_total = lca_inner + lca_outer + lca_mli + lca_str
            return lca_total
        
        def embodied_energy(self,mass_inner,mass_outer,mass_mli,mass_str,kevlar_ee,mli_ee):
            ee_inner = mass_inner * self.mat_ee
            ee_outer = mass_outer * self.mat2_ee
            ee_mli = mass_mli * mli_ee
            ee_str = mass_str * kevlar_ee
            ee_total = ee_inner + ee_outer + ee_mli + ee_str
            return ee_total

    # -------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------- Tank Database ---------------------------------------------------
    materials = ['S-Glass Fiber'] #From Granta and Engineering Toolbox
    SF = fos #NASA s safety factor for the materials
    #density in kg/m^3, yield strength in Pa, thermal conductivity in W/mK, emissivity in [-], CO2 [kg/kg], Embodied Energy in MJ/kg, Fibre Ratio
    mat_properties = [ [gf_density, gf_tensile_strength, gf_thermal_cond, gf_thermal_emis, gf_co2, gf_ee, gf_fvf]]
    #MAWPS = [600000,650000,800000,1000000,1200000,1280000] #bar to Pa
    MAWP = MAWP_global
    P_vent = p_vent_global

    # ------------------------------------------------- Structure ------------------------------------------------------
    # Thesis values
    Q_og_str = Q_original_str #W
    k_kevlar = kevlar_thermal_cond #W/mK
    og_str_mass = mass_original_str #kg
    og_lh2 = mass_originalg_lh2 #kg
    grav_idx = gravimetric_index #from NASA report
    co2_kevlar = kevlar_co2 #kg/kg (Kevlar 149)
    kevlar_ee = kevlar_emb_energy #MJ/kg (Embodied Energy for Kevlar 149)

    #Our values
    mass_h2 = m_h2#278.9577  #kg
    estimated_mass = mass_h2/grav_idx - mass_h2
    t_limit = t_min #m (minimum thickness of the tank wall)

    #Insulation

    dens_mli = mli_density #kg/m^3 https://www.sciencedirect.com/science/article/pii/S135943112200391X
    emis_mli = mli_emis  #https://www.thermalengineer.com/library/effective_emittance.htm
    k_vac = vacuum_thermal_cond#3 # W/mK https://www.researchgate.net/publication/321219004_Cylindrical_Cryogenic_Calorimeter_Testing_of_Six_Types_of_Multilayer_Insulation_Systems
    k_mli = mli_thermal_cond # W/mK  https://www.sciencedirect.com/science/article/pii/S135943112200391X
    mli_co2 = mli_ss_co2 #kg/kg (for SS)
    mli_ee = mli_ss_ee #MJ/kg (Embodied Energy for SS)
    N_MLI = mli_layers
    t_mli = mli_thickness #https://www.sciencedirect.com/science/article/pii/S135943112200391X


    def fA(mh2, P_vent, fl_final = 0.98):

        # From the heat influx get the maximum initial fill fraction
        rho_l_f = CP.PropsSI('D', 'P', P_vent, 'Q', 0, 'ParaHydrogen')  # Final liquid density
        # print(f"Final liquid density: {rho_l_f} kg/m^3")
        V_l_fin = mh2 / rho_l_f  # Final liquid volume
        V_tot   = V_l_fin / fl_final

        # Calculate the initial fill fraction
        rho_l_0 = CP.PropsSI('D', 'P', P_sl, 'Q', 0, 'ParaHydrogen')  # Initial liquid density
        # print(f"Initial liquid density: {rho_l_0} kg/m^3")
        V_l_0   = mh2 / rho_l_0  # Initial liquid volume
        fl_init = V_l_0 / V_tot  # Initial fill fraction

        return V_tot, fl_init

    def compute_Qleak(material, material2, mat_property, MAWP,mass_h2, Q_str,mat2_property,str_mass,fill_ratio,V_in,P_vent):
        tankh2 = Tank(MAWP, material, material2, mat_property,mass_h2,mat2_property,fill_ratio,V_in,P_vent)
        # -----------------------------Maximum allowable heat load -----------------------------------
        Q_solution = opt.root_scalar(tankh2.maximum_Qin, bracket=[tankh2.Q_leak_min, tankh2.Q_leak_max], method='brentq')

        if Q_solution.converged:
            Qmax = Q_solution.root
            # print(f"Calculated Qin: {Qmax:.4f} W")
        else:
            print("Failed to find an Qin_max. Adjust the bounds.")
        return Qmax

    def compute_tank(material, material2, mat_property, MAWP,mass_h2, Q_str,mat2_property,str_mass,fill_ratio,V_in,P_vent,Qmax):
        tankh2 = Tank(MAWP, material, material2, mat_property,mass_h2,mat2_property,fill_ratio,V_in,P_vent)
        dv_list = []
        # -----------------------------Inner tank sizing -----------------------------------
        L_solution = opt.root_scalar(tankh2.volume_equation, bracket=[2*tankh2.R_in, 10], method='brentq')

        if L_solution.converged:
            L_in = L_solution.root
            # print(f"Calculated length: {L_in:.4f} meters")
        else:
            print("Failed to find an inner length. Adjust the bounds.")

        t1, ang1_w = tankh2.inner_tank_thickness()
        if t1 < t_limit:
            t1 = t_limit
        L_cyl = L_in - 2 * tankh2.R_in
        dv, t2, ang2_w, Qleak, Qcond, Qrad = tankh2.heat_influx(L_cyl, Q_str,t1,emis_mli,k_vac,t_mli,k_mli,Qmax, N_MLI,t1)
        dv_list.append(dv)
        t2 = max(t2,t_limit)
        t20 = t2 + 1
        while abs(t2-t20) > 1e-20:
            t20 = t2
            dv,t2,ang2_w,Qleak,Qcond,Qrad = tankh2.heat_influx(L_cyl, Q_str,t1,emis_mli,k_vac,t_mli,k_mli,Qmax, N_MLI,t20)
            dv_list.append(dv)
            t2 = max(t2,t_limit)

        Vt,L_out,R_out = tankh2.total_volume(L_cyl, dv,t1,t2,t_mli)
        mass_inner,mass_outer,mass_mli = tankh2.total_mass(L_in, dv, t1, t2,t_mli,dens_mli) 
        Mt = mass_inner + mass_outer + mass_mli + mass_h2
        Mt += str_mass
        mass_error = Mt - estimated_mass - mass_h2
        Vh2 = tankh2.V_in

        co2_kg = tankh2.kg_co2(mass_inner,mass_outer,mass_mli,str_mass,mli_co2,co2_kevlar)
        embodied_energy = tankh2.embodied_energy(mass_inner,mass_outer,mass_mli,str_mass,kevlar_ee,mli_ee)

        return Vt, Mt, mass_error,Vh2,t1,t2,dv,L_in, Vh2,tankh2.R_in, co2_kg, embodied_energy, ang1_w, ang2_w, P_vent,Qleak, Qcond, Qrad, dv_list, L_out, R_out

    # ------------------------------------------------- Main ------------------------------------------------------

    V_in, fill_ratio = fA(mass_h2, P_vent)
    V_in = V_in
    Qmax = compute_Qleak(materials[0], materials[0], mat_properties[0], MAWP,mass_h2,0,mat_properties[0],0,fill_ratio,V_in, P_vent)
    for material, mat_property in zip(materials, mat_properties):
        for material2, mat2_property in zip(materials, mat_properties):
            if material == 'S-Glass Fiber' or  material == 'Carbon Fiber': #Composites
                og_tank_mass = 4.8+3.15
            else: #Metallic materials (compared to aluminium)
                og_tank_mass = 8.4+3.6
            ratio = og_str_mass / (og_tank_mass+og_lh2)
            str_mass = estimated_mass * ratio
            Q_str = Q_og_str * np.sqrt(ratio/ratio) #Assuming Volume increase with the same ratio as the mass
            Vt, Mt, mass_error,Vh2,t1,t2,dv,L_in,Vh2,R_in,co2_kg,emb_energy, ang1_w, ang2_w, p_vent, Qleak, Qcond, Qrad, dv_list, L_out, R_out = compute_tank(material, material2, mat_property, MAWP,mass_h2,Q_str,mat2_property,str_mass,fill_ratio,V_in, P_vent,Qmax)
            # print(f"Material In: {material}, Material Out: {material2}, MAWP: {MAWP} Pa, P_vent:{P_vent}")
            # print(f"Tank Volume: {Vt:.4f} m^3")
            # print(f"Tank Mass: {Mt:.4f} kg")
            # print(f"CO2 emissions: {co2_kg:.4f} kg")
            # print(f"Embodied Energy: {emb_energy:.4f} MJ")
    return Mt, Vt

print(main_storage(278.9577))  # Example call with a mass of 278.9577 kg of hydrogen