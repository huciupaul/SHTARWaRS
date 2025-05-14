import CoolProp.CoolProp as CP # Used to get properties of Hydrogen
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.optimize import fsolve

class Tank:
    def __init__(self, MAWP, material, mat_property,mass_h2):
        # ---------------------------------------------Database ------------------------------------------------------------------
        self.material = material
        self.mat_property = mat_property
        self.mat_density = mat_property[0] 
        self.mat_yield_strength = mat_property[1]
        self.mat_thermal_conductivity = mat_property[2]
        self.mat_emissivity = mat_property[3]

        # ---------------------------------------------Constants ----------------------------------------------------------------
        self.dormancy = 24 #hours
        self.fill_ratio = 0.8 #percentage of LH2 at the start
        self.mass_h2 = mass_h2 #kg
        self.stratification_factor = 2
        self.k_str = 1.9 #[W/mK]


        # --------------------------------------------- Inputs ------------------------------------------------------------------

        self.MAWP = MAWP # Convert from Bar to Pa for the venting pressure
        self.Pvent = MAWP #Assumption

        # Hydrogen properties
        self.P0 = 101000 # Initial pressure in Pa
        self.T0 = CP.PropsSI('T', 'P', self.P0, 'Q', 0.1, 'ParaHydrogen') # Initial temperature in K
        self.rhol0 = CP.PropsSI('D', 'T', self.T0, 'Q', 0, 'ParaHydrogen') # density of liquid hydrogen
        self.rhog0 = CP.PropsSI('D', 'T', self.T0, 'Q', 1, 'ParaHydrogen') # density of gas hydrogen
        self.V_in = self.mass_h2 / (self.fill_ratio * self.rhol0) # Volume in m^3
        self.mg0 = CP.PropsSI('D', 'T', self.T0, 'Q', 1, 'ParaHydrogen')*self.V_in*(1-self.fill_ratio) # kg
        self.ml0 = CP.PropsSI('D', 'T', self.T0, 'Q', 0, 'ParaHydrogen')*self.V_in*(self.fill_ratio) # kg
        self.hg0 = CP.PropsSI('H', 'T', self.T0, 'Q', 1, 'ParaHydrogen')#kJ/kg
        self.hl0 = CP.PropsSI('H', 'T', self.T0, 'Q', 0, 'ParaHydrogen')#kJ/kg
        self.x0 = self.mg0/(self.mg0+self.ml0) #unitless (vapor quality)

        # ---------------------------------------------- Constraint ---------------------------------------------------
        self.R_in = 0.8 #m
        self.Q_leak_min = 50 #W (determined from Nicolas's thesis)
        self.Q_leak_max = 600 #W (determined from Nicolas's thesis)

    
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
        # Run the simulation until the pressure reaches the venting pressure
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
        print(0)
        return time_hours[-1] - self.dormancy

    # ----------------------------------------------- Inner tank Parameters ---------------------------------------------------
    def volume_equation(self,l):
        L_cyl = l - 2*self.R_in  # cylindrical length
        if L_cyl < 0:
            return 1e6  # Penalize invalid cases
        V = np.pi * self.R_in**2 * L_cyl + (4/3) * np.pi * self.R_in**3  # cylinder + 2 hemispheres = 1 sphere
        return V - self.V_in  # we want this to be zero

    
    def inner_tank_thickness(self):
        return MAWP * self.R_in / self.mat_property[1]



    # ------------------------------------------------ Vacuum and Outer tank  ---------------------------------------------------


    def heat_influx(self, L_in_max, Q_str,t1,emis_mli,k_vac,t_mli):
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
            R_cond += np.log((r_in + t1 + dv) / (r_in + t1)) / (2 * np.pi * L_in_max * k_vac)
            R_cond += np.log((r_in + t1 + dv + t2) / (r_in + t1 + dv)) / (2 * np.pi * L_in_max * k_2)
            return (T_amb - T_tank) / R_cond

        # Radiation...
        def Q_rad(dv):
            r1 = r_in + t1 + t_mli        # Outer surface of inner tank
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



    # ------------------------------------------------ Tank Dimensions ---------------------------------------------------
    def total_volume(self, l,dv,t1,t2,t_mli):
        R_out = self.R_in + dv +t1 +t2 +t_mli
        V = np.pi * R_out**2 * (l-2*self.R_in) + (4/3) * np.pi * R_out**3
        return V
    
    def total_mass(self, l, dv, t1, t2,t_mli,dens_mli):
        R_out = self.R_in + dv +t1 +t2 +t_mli
        L_cyl = l - 2 * self.R_in  
        surface_inner = 2 * np.pi * self.R_in * L_cyl + 4 * np.pi * self.R_in**2
        surface_outer = 2 * np.pi * R_out * L_cyl + 4 * np.pi * R_out**2
        mass_inner = surface_inner * t1 * self.mat_density
        mass_outer = surface_outer * t2 * self.mat_density
        mass_mli = surface_inner * t_mli * dens_mli
        return mass_inner + mass_outer + self.mass_h2 + mass_mli



# -------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------- Tank Database ---------------------------------------------------
materials = ['Al-7075-T6','G10']
SF = 0.9
mat_properties = [[2800,495*1e6*SF,134,0.05],[1900,298*1e6*SF,0.5,0.95]]  #density in kg/m^3, yield strength in Pa, thermal conductivity in W/mK, emissivity in [-]
#MAWPS = [400000,650000,800000,1000000,1200000,1280000] #bar to Pa
MAWPS = [600000,850000,1280000] #test


# ------------------------------------------------- Structure ------------------------------------------------------
# Thesis values
Q_og_str = 0.4 #W
k_kevlar = 1.9 #W/mK
og_str_mass = 2.1 #kg
og_lh2 = 6.2
grav_idx = 0.35 #from NASA report

#Our values
mass_h2 = 250 #kg
estimated_mass = mass_h2/grav_idx - mass_h2

#Insulation
t_mli = 5 *1e-3
dens_mli = 20 #kg/m^3
emis_mli = 0.03
k_vac = 0.01 # W/mK



def compute_tank_volume(material, mat_property, MAWP,mass_h2, Q_str):
    tankh2 = Tank(MAWP, material, mat_property,mass_h2)

    # -----------------------------Maximum allowable heat load -----------------------------------
    Q_solution = opt.root_scalar(tankh2.maximum_Qin, bracket=[tankh2.Q_leak_min, tankh2.Q_leak_max], method='brentq')

    if Q_solution.converged:
        Qmax = Q_solution.root
        print(f"Calculated Qin: {Qmax:.4f} W")
    else:
        print("Failed to find an Qin_max. Adjust the bounds.")

    # -----------------------------Inner tank sizing -----------------------------------
    L_solution = opt.root_scalar(tankh2.volume_equation, bracket=[2*tankh2.R_in, 10], method='brentq')

    if L_solution.converged:
        L_in = L_solution.root
        print(f"Calculated length: {L_in:.4f} meters")
    else:
        print("Failed to find an inner length. Adjust the bounds.")

    t1 = tankh2.inner_tank_thickness()
    print(f"Inner tank thickness: {t1:.4f} m")
    print(f"MAWP:",MAWP)
    print(f"R_in:",tankh2.R_in)
    print(f"Strength:",tankh2.mat_property[1])

    dv, t2 = tankh2.heat_influx(L_in, Q_str,t1,emis_mli,k_vac,t_mli)
    
    Vt = tankh2.total_volume(L_in, dv,t1,t2,t_mli)
    Mt = tankh2.total_mass(L_in, dv, t1, t2,t_mli,dens_mli) 
    Mt += str_mass
    mass_error = Mt - estimated_mass - mass_h2

    return Vt, Mt, mass_error


# ------------------------------------------------- Main ------------------------------------------------------

plot_mats = []
plot_MAWPS = []
plot_volumes = []
plot_masses = []
plot_mass_errors = []


for MAWP in MAWPS:
    for material, mat_property in zip(materials, mat_properties):
        if material == 'G10': #Composite materials
            og_tank_mass = 4.8+3.15
        else: #Metallic materials
            og_tank_mass = 8.4+3.6
        ratio = og_str_mass / (og_tank_mass+og_lh2)
        str_mass = estimated_mass * ratio 
        Q_str = Q_og_str * np.sqrt(ratio/ratio) #Assuming Volume increase with the same ratio as the mass
        Vt, Mt, mass_error = compute_tank_volume(material, mat_property, MAWP,mass_h2,Q_str)
        print(f"Material: {material}, MAWP: {MAWP} Pa")
        print(f"Tank Volume: {Vt:.4f} m^3")
        print(f"Tank Mass: {Mt:.4f} kg")
        print(f"Mass Error: {mass_error:.4f} kg\n")

        plot_mats.append(material)
        plot_MAWPS.append(MAWP)
        plot_volumes.append(Vt)
        plot_masses.append(Mt)
        plot_mass_errors.append(mass_error)


# ------------------------------------------------- Plotting ------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: Volume vs MAWP
for material in materials:
    material_volumes = [plot_volumes[i] for i in range(len(plot_mats)) if plot_mats[i] == material]
    material_MAWPS = [plot_MAWPS[i] for i in range(len(plot_mats)) if plot_mats[i] == material]
    axs[0].plot(material_MAWPS, material_volumes, label=material, marker='o')

axs[0].set_xlabel('MAWP (Pa)')
axs[0].set_ylabel('Tank Volume (m^3)')
axs[0].set_title('Tank Volume vs MAWP')
axs[0].legend()
axs[0].grid(True)

# Plot 2: Mass vs MAWP
for material in materials:
    material_masses = [plot_masses[i] for i in range(len(plot_mats)) if plot_mats[i] == material]
    material_MAWPS = [plot_MAWPS[i] for i in range(len(plot_mats)) if plot_mats[i] == material]
    axs[1].plot(material_MAWPS, material_masses, label=material, marker='o')

axs[1].set_xlabel('MAWP (Pa)')
axs[1].set_ylabel('Tank Mass (kg)')
axs[1].set_title('Tank Mass vs MAWP')
axs[1].legend()
axs[1].grid(True)

# Plot 3: Mass Error vs MAWP
for material in materials:
    material_mass_errors = [plot_mass_errors[i] for i in range(len(plot_mats)) if plot_mats[i] == material]
    material_MAWPS = [plot_MAWPS[i] for i in range(len(plot_mats)) if plot_mats[i] == material]
    axs[2].plot(material_MAWPS, material_mass_errors, label=material, marker='o')

axs[2].set_xlabel('MAWP (Pa)')
axs[2].set_ylabel('Mass Error (kg)')
axs[2].set_title('Mass Error vs MAWP')
axs[2].legend()
axs[2].grid(True)

plt.tight_layout()
plt.show()




