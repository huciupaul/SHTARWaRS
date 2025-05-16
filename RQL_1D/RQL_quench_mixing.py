import cantera as ct
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm
import matplotlib.pyplot as plt
from mixture_properties import (
    mixture_properties
)

'''
This script is used to estimate the maximum mixing time in the quench zone of a RQL combustor 
to keep the NOx emissions under the specified limit. 
'''

#----------Inputs----------
phi = 1.0
MECH = 'gri30.yaml'  # Mechanism for H2 combustion with NOx

# Rich zone exhaust mix
mdot_rich_ex = 0.1  # kg/s
X_rich_ex = "N2:0.7, O2:0.1"
T_rich_ex = 703.0     # K
P_rich_ex = 12.1 * ct.one_atm     # Pa

# Air
mdot_air = 0.2      # kg/s
X_air = "N2:0.78084, O2:0.20946"
T_air = 603.0       #  K

# Water
mdot_h2o = np.arange(0, 0.26, 0.01)
X_h2o = "H2O:1"
T_h2o = 750.0  # K


#----------Calculations----------

X_mix1, mdot_tot1, _ = mixture_properties(
    mdot1=mdot_rich_ex, X1=X_rich_ex,
    mdot2=mdot_air, X2=X_air
)

# Compute enthalpy-weighted temperature for the mixture
gas = ct.Solution(MECH)

h_tot = 0.0
for mdot, X, T in [
    (mdot_rich_ex, X_rich_ex, T_rich_ex),
    (mdot_air, X_air, T_air)
]:
    if mdot > 0:
        gas.TPX = T, P_rich_ex, X
        h_tot += mdot * gas.enthalpy_mass

h_mix = h_tot / mdot_tot1
gas.HP = h_mix, P_rich_ex
gas.equilibrate('HP')

# Set up 0D reactor
reactor = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([reactor])

times = []
NOx = []
temps = []

time = 0.0
t_end = 0.05  # seconds
dt = 1e-3
num_steps = int(t_end / dt)

for step in tqdm(range(num_steps), desc="Simulation time"):
    time = (step + 1) * dt
    times.append(time)
    NOx_row = []
    temps_row = []
    for mdot_h2o_i in mdot_h2o:
        # Mix rich+air with current h2o
        X_mix, mdot_mix, _ = mixture_properties(
            mdot1=mdot_tot1, X1=X_mix1,
            mdot2=mdot_h2o_i, X2=X_h2o
        )

        # Compute enthalpy-weighted temperature for the new mixture
        h_tot = 0.0
        for mdot, X, T in [
            (mdot_tot1, X_mix1, gas.T),  # Use last gas.T as approx
            (mdot_h2o_i, X_h2o, T_h2o)
        ]:
            if mdot > 0:
                gas.TPX = T, P_rich_ex, X
                h_tot += mdot * gas.enthalpy_mass

        h_mix = h_tot / (mdot_tot1 + mdot_h2o_i)
        gas.HP = h_mix, P_rich_ex
        gas.equilibrate('HP', solver='vcs')
        
        # Set up 0D reactor
        reactor = ct.IdealGasConstPressureReactor(gas)
        sim = ct.ReactorNet([reactor])
        sim.advance(time)

        temps_row.append(reactor.T)
        NOx_row.append(reactor.thermo['NO'].X[0] + reactor.thermo['NO2'].X[0])

    temps.append(temps_row)
    NOx.append(NOx_row)

# Convert lists to numpy arrays for plotting
times = np.array(times)
NOx = np.array(NOx)
temps = np.array(temps)
mdot_h2o = np.array(mdot_h2o)

# Create meshgrid for time and mdot_h2o
T, M = np.meshgrid(times, mdot_h2o, indexing='ij')

# Plot NOx surface
fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(T, M, NOx, cmap='viridis')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('mdot_h2o [kg/s]')
ax1.set_zlabel('NOx mole fraction')
ax1.set_title('NOx vs Time and mdot_h2o')

# Plot temperature surface
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(T, M, temps, cmap='plasma')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('mdot_h2o [kg/s]')
ax2.set_zlabel('Temperature [K]')
ax2.set_title('Temperature vs Time and mdot_h2o')

plt.tight_layout()
plt.show()
