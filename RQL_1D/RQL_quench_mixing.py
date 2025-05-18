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

X_mix1, mdot_mix1, _ = mixture_properties(
    mdot1=mdot_rich_ex, X1=X_rich_ex,
    mdot2=mdot_air, X2=X_air
)
T_mix1 = (T_rich_ex*mdot_rich_ex  + T_air*mdot_air) / mdot_mix1

# Compute enthalpy-weighted temperature for the mixture
gas = ct.Solution(MECH)
gas.TPX = T_mix1, P_rich_ex, X_mix1

# Set up 0D reactor
reactor = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([reactor])

times = []
NOx = []
temps = []

time = 0.0
t_end = 0.05  # seconds
dt = 1e-6
num_steps = int(t_end / dt)

for mdot_h2o_i in mdot_h2o:
    NOx_row = []
    temps_row = []
    time = 0.0
    # Recompute mixture for this mdot_h2o_i
    X_mix, mdot_mix, _ = mixture_properties(
        mdot1=mdot_mix1, X1=X_mix1,
        mdot2=mdot_h2o_i, X2=X_h2o
    )
    T_mix = (T_mix1*mdot_mix1  + T_h2o*mdot_h2o_i) / mdot_mix
    gas.TPX = T_mix, P_rich_ex, X_mix

    # Set up 0D reactor
    reactor = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([reactor])

    for step in range(num_steps):
        time += dt
        if len(times) < num_steps:
            times.append(time)
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

# Plot NOx heatmap
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
  
im1 = ax1.imshow(
    NOx, 
    aspect='auto', 
    origin='lower', 
    extent=[times[0], times[-1], mdot_h2o[0], mdot_h2o[-1]],
    cmap='viridis'
)
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('mdot_h2o [kg/s]')
ax1.set_title('NOx mole fraction')
fig.colorbar(im1, ax=ax1, label='NOx mole fraction')

# Add contour for 1 ppm NOx (1e-6 mole fraction)
T, M = np.meshgrid(times, mdot_h2o)
contour = ax1.contour(
    times, mdot_h2o, NOx, levels=[1e-38], colors='red', linewidths=2, origin='lower'
)
ax1.clabel(contour, fmt={1e-38: '1e-38'}, colors='red')


# Plot temperature heatmap
im2 = ax2.imshow(
    temps, 
    aspect='auto', 
    origin='lower', 
    extent=[times[0], times[-1], mdot_h2o[0], mdot_h2o[-1]],
    cmap='plasma'
)
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('mdot_h2o [kg/s]')
ax2.set_title('Temperature [K]')
fig.colorbar(im2, ax=ax2, label='Temperature [K]')

plt.tight_layout()
plt.show()