import CoolProp.CoolProp as CP # Used to get properties of Hydrogen
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

# ---------------------------------------------Constants ----------------------------------------------------------------
dormancy = 24 #hours
fill_ratio = 0.8 #percentage of LH2 at the start
mass_h2 = 100 #kg
stratification_factor = 2


# --------------------------------------------- Inputs ------------------------------------------------------------------

MAWP = 6.5*100000 # Convert from Bar to Pa for the venting pressure
Pvent = MAWP #Assumption

# Hydrogen properties
P0 = 101000 # Initial pressure in Pa
T0 = CP.PropsSI('T', 'P', P0, 'Q', 0.1, 'ParaHydrogen') # Initial temperature in K
rhog0 = CP.PropsSI('D', 'T', T0, 'Q', 0, 'ParaHydrogen') # density of liquid hydrogen
rhol0 = CP.PropsSI('D', 'T', T0, 'Q', 1, 'ParaHydrogen') # density of gas hydrogen
V_in = mass_h2 / (fill_ratio * rhol0) # Volume in m^3
print("Volume in m^3: ", V_in)
mg0 = CP.PropsSI('D', 'T', T0, 'Q', 1, 'ParaHydrogen')*V_in*(1-fill_ratio) # kg
ml0 = CP.PropsSI('D', 'T', T0, 'Q', 0, 'ParaHydrogen')*V_in*(fill_ratio) # kg
hg0 = CP.PropsSI('H', 'T', T0, 'Q', 1, 'ParaHydrogen')#kJ/kg
hl0 = CP.PropsSI('H', 'T', T0, 'Q', 0, 'ParaHydrogen')#kJ/kg
x0 = mg0/(mg0+ml0) #unitless (vapor quality)

# ---------------------------------------------- Constraint ---------------------------------------------------
L_max = 2 #m
Q_leak_min = 2 #W (determined from Nicolas's thesis)
Q_leak_max = 4 #W (determined from Nicolas's thesis)


# ----------------------------------------------- Maximum Heat Load ---------------------------------------------------
def maximum_Qin(Qleak):
    P = [P0]
    T = [T0]
    mg = [mg0]
    ml = [ml0]
    rhog = [rhog0]
    rhol = [rhol0]
    hg = [hg0]
    hl = [hl0]
    x = [x0]
    Ps = [P0]
    time = [0]
    dt = 250 # Time step in seconds
    i=0
    # Run the simulation until the pressure reaches the venting pressure
    while P[-1]<Pvent:
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
        +(((ml[-1]*delhl_delP)+(mg[-1]*delhg_delP)-V_in)*dPs_dT)
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
        dP_dt = differential_matrix[0, 0]*stratification_factor
        dT_dt = differential_matrix[1, 0]
        dmg_dt = differential_matrix[2, 0]
        dml_dt = differential_matrix[3, 0]
        print(i)
        i+=1
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
    return time_hours[-1] - dormancy

Q_solution = opt.root_scalar(maximum_Qin, bracket=[Q_leak_min, Q_leak_max], method='brentq')

if Q_solution.converged:
    Qmax = Q_solution.root
    print(f"Calculated radius: {Qmax:.4f} meters")
else:
    print("Failed to find an inner radius.")

# ----------------------------------------------- Inner tank Parameters ---------------------------------------------------
def volume_equation(r):
    L_cyl = L_max - 2*r  # cylindrical length
    if L_cyl < 0:
        return 1e6  # Penalize invalid cases
    V = np.pi * r**2 * L_cyl + (4/3) * np.pi * r**3  # cylinder + 2 hemispheres = 1 sphere
    return V - V_in  # we want this to be zero

r_solution = opt.root_scalar(volume_equation, bracket=[0.001, L_max/2], method='brentq')

if r_solution.converged:
    R_in = r_solution.root
    print(f"Calculated radius: {R_in:.4f} meters")
else:
    print("Failed to find an inner radius.")

