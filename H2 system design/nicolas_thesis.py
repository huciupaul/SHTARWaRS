import CoolProp.CoolProp as CP # Used to get properties of Hydrogen
import numpy as np
import matplotlib.pyplot as plt
# Defining the Scenarios
# Q_leak is the amount of heat leakage entering the vessel.
# Fill ratio is the percentage of LH2 at the start.
scenarios = [
{'fill_ratio': 0.40, 'Q_leak': 36.5},
{'fill_ratio': 0.80, 'Q_leak': 36.5},
{'fill_ratio': 0.40, 'Q_leak': 1.5},
{'fill_ratio': 0.80, 'Q_leak': 1.5}
]
Pvent = 6.5*100000 # Convert from Bar to Pa for the venting pressure
stratification_factor = 2 #Doubles the pressure rate for each increment
fig, ax = plt.subplots(figsize=(10, 6))
# List of vertical offsets for annotations to avoid overlapping
offsets = [-0.1, -0.2, 0.1, 0.2]
# Iterate over each scenario
for idx, scenario in enumerate(scenarios):
    fill_ratio_initial = scenario['fill_ratio']
    Q_leak = scenario['Q_leak']
    # Initial conditions
    P0 = 101000 # Initial pressure in Pa
    T0 = CP.PropsSI('T', 'P', P0, 'Q', 0.1, 'ParaHydrogen') # Initial temperature in K
    V = 0.091 # Volume in m^3
    # Initial masses of gas and liquid hydrogen
    mg0 = CP.PropsSI('D', 'T', T0, 'Q', 1, 'ParaHydrogen')*V*(1-fill_ratio_initial) # kg
    ml0 = CP.PropsSI('D', 'T', T0, 'Q', 0, 'ParaHydrogen')*V*(fill_ratio_initial) # kg
    # Initial densities and enthalpies of gas and liquid hydrogen
    rhog0 = CP.PropsSI('D', 'T', T0, 'Q', 1, 'ParaHydrogen')# density of gas
    rhol0 = CP.PropsSI('D', 'T', T0, 'Q', 0, 'ParaHydrogen')#density
    hg0 = CP.PropsSI('H', 'T', T0, 'Q', 1, 'ParaHydrogen')#kJ/kg
    hl0 = CP.PropsSI('H', 'T', T0, 'Q', 0, 'ParaHydrogen')#kJ/kg
    # Initial vapor quality
    x0 = mg0/(mg0+ml0) #unitless
    # Store initial conditions
    initial_conditions = [P0, T0, mg0, ml0]
    # Initialize lists to store time series data
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
        +(((ml[-1]*delhl_delP)+(mg[-1]*delhg_delP)-V)*dPs_dT)
        A = np.array([[1, -dPs_dT, 0, 0],
        [A21, A22, (1/rhog[-1])-(1/rhol[-1]), 0],
        [0, 0, 1, 1],
        [0, A42, -(hl[-1]-hg[-1]), 0]])
        # Construct the B matrix
        B = np.array([[0],
        [0],
        [0],
        [Q_leak]])
        # Solve the system of linear equations
        differential_matrix = np.linalg.solve(A, B)
        dP_dt = differential_matrix[0, 0]*stratification_factor
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
    # Plot the pressure vs time for each scenario
    label = f'Fill Rate: {fill_ratio_initial*100}%, Heat Load: {Q_leak} W'
    ax.plot(time_hours, P_bars, label=label)
# Set labels with increased font size
ax.set_xlabel('Time (hours)', fontsize=14)
ax.set_ylabel('Pressure (Bar)', fontsize=14)
ax.set_title('Pressure Rise vs Time for Different Initial Fill Levels and Heat Leaks',
fontsize=14)
ax.legend()
ax.set_xlim(0, 50)
ax.set_ylim(0.5, 7.3)
plt.show()