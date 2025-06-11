import cantera as ct
ct.suppress_thermo_warnings()
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm
import matplotlib.pyplot as plt
from mixture_properties import (
    mix_streams_const_HP
)

'''
This script is used to estimate the maximum mixing time in the quench zone of a RQL combustor 
to keep the NOx emissions under the specified limit.
'''

#----------Inputs----------
MECH = 'gri30.yaml'  # Mechanism for H2 combustion with NOx

# Rich zone exhaust mix
mdot_rich_ex = 0.23071218802678672  # kg/s
X_rich_ex = "N2:0.27396853378378067, H2:0.5784925158354742, H:0.0005133907896591396, O2:7.393652373124469e-10, O:4.591482831748451e-09, H2O:0.14701811528905126, OH:7.430992453322465e-06, H2O2:1.1742925935728515e-10, HO2:3.6285427485688574e-12, NO:5.4364866743951985e-09, N2O:9.901157119223157e-12, N:4.953608429224339e-12, NH:1.4251848100599631e-11, NNH:1.977256219409285e-09, NH2:4.148258601540983e-10"
T_rich_ex = 734.22      # K
P_rich_ex = 12.1 * 1e5  # Pa

# Air
mdot_air_tot = 4.635545784891179     # kg/s
mdot_air1 = 0.20131953579998843 # kg/s
mdot_bled_min = (mdot_air_tot - mdot_air1)*0.0525 # kg/s
mdot_bled_max = (mdot_air_tot - mdot_air1)*0.73 # kg/s
mdot_air2_min = mdot_air_tot - mdot_air1 - mdot_bled_max # kg/s
mdot_air2_max = mdot_air_tot - mdot_air1 - mdot_bled_min # kg/s
mdot_air2 = np.arange(mdot_air2_min, mdot_air2_max, (mdot_air2_max - mdot_air2_min) * 0.1) # kg/s
X_air = "N2:0.78084, O2:0.20946"
T_air = 603.0  #  K

# Water
mdot_h2o = np.arange(0, 0.26, 1e-2) # kg/s
X_h2o = "H2O:1"
T_h2o = 333.0  # K

# Total time for mixing
t_mix_tot = np.arange(0.01, 1.01, 1e-2)  # seconds

#----------Outputs----------
times_to_mix = []
NOx = []
temps = []

NOx_limit = 25 * 1e-6 # mole fraction

#----------Iteration----------
for t_end in tqdm(t_mix_tot):
    for mdot_air2_i in mdot_air2:
        for mdot_h2o_i in mdot_h2o:
            t = 0.0
            dt = 1e-2
            NOx_tot = 0.0
            m_air2_step = mdot_air2_i * (dt / t_end)
            m_h2o_step = mdot_h2o_i * (dt / t_end)

            # Create three separate gas objects for each stream to avoid overwriting state
            gas_rich_ex = ct.Solution(MECH)
            gas_rich_ex.TPX = T_rich_ex, P_rich_ex, X_rich_ex
            stream1 = ct.Quantity(gas_rich_ex, mass=mdot_rich_ex)

            gas_air = ct.Solution(MECH)
            gas_air.TPX = T_air, P_rich_ex, X_air
            stream2 = ct.Quantity(gas_air, mass=m_air2_step)

            gas_h2o = ct.Solution(MECH)
            gas_h2o.TPX = T_h2o, P_rich_ex, X_h2o
            gas_h2o.HP = gas_h2o.enthalpy_mass + 1.9871*1e6, P_rich_ex
            stream3 = ct.Quantity(gas_h2o, mass=m_h2o_step)

            while t <= t_end and NOx_tot <= NOx_limit:
                # Set up gas object
                gas = ct.Solution(MECH)                

                # Mix the streams at constant pressure and enthalpy
                X_mix, T_mix, mdot_mix = mix_streams_const_HP(stream1, stream2, stream3, P_rich_ex)
                gas.TPX = T_mix, P_rich_ex, X_mix
                stream1 = ct.Quantity(gas, mass=mdot_mix)

                # Set up 0D reactor
                reactor = ct.IdealGasConstPressureReactor(gas)
                sim = ct.ReactorNet([reactor])
                sim.advance(t)
                NOx_tot += reactor.thermo['NO'].X[0] + reactor.thermo['NO2'].X[0]
                t += dt

                if t >= t_end or NOx_tot >= NOx_limit:
                    times_to_mix.append(t)
                    NOx.append(NOx_tot)
                    temps.append(reactor.T)
                    break
                try:
                    # Continue simulation as normal
                    pass
                except Exception as e:
                    times_to_mix.append(t)
                    NOx.append(NOx_tot)
                    temps.append(reactor.T if 'reactor' in locals() else np.nan)
                    break




# Convert lists to numpy arrays for plotting
times = np.array(times_to_mix)
NOx = np.array(NOx)
temps = np.array(temps)
mdot_h2o = np.array(mdot_h2o)

# Reshape results for heatmap plotting
# Find the maximum lengths for mdot_air2 and mdot_h2o
n_air2 = len(mdot_air2)
n_h2o = len(mdot_h2o)

# Take the last t_end (or any fixed t_end) for plotting
NOx_arr = np.array(NOx).reshape(-1, n_air2, n_h2o)[-1]
temps_arr = np.array(temps).reshape(-1, n_air2, n_h2o)[-1]

# If arrays are smaller than the max shape, pad with np.nan (will show as white in imshow)
max_n_air2 = n_air2
max_n_h2o = n_h2o

NOx_padded = np.full((max_n_air2, max_n_h2o), np.nan)
temps_padded = np.full((max_n_air2, max_n_h2o), np.nan)

shape = NOx_arr.shape
NOx_padded[:shape[0], :shape[1]] = NOx_arr
temps_padded[:shape[0], :shape[1]] = temps_arr

# Plot NOx and temperature heatmaps
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# NOx heatmap
im1 = ax1.imshow(
    NOx_padded.T,
    aspect='auto',
    origin='lower',
    extent=[mdot_air2[0], mdot_air2[-1], mdot_h2o[0], mdot_h2o[-1]],
    cmap='viridis'
)
ax1.set_xlabel('mdot_air2 [kg/s]')
ax1.set_ylabel('mdot_h2o [kg/s]')
ax1.set_title('NOx (mole fraction)')
fig.colorbar(im1, ax=ax1, label='NOx mole fraction')

# Temperature heatmap
im2 = ax2.imshow(
    temps_padded.T,
    aspect='auto',
    origin='lower',
    extent=[mdot_air2[0], mdot_air2[-1], mdot_h2o[0], mdot_h2o[-1]],
    cmap='plasma'
)
ax2.set_xlabel('mdot_air2 [kg/s]')
ax2.set_ylabel('mdot_h2o [kg/s]')
ax2.set_title('Temperature [K]')
fig.colorbar(im2, ax=ax2, label='Temperature [K]')

plt.tight_layout()
plt.show()












"""
# Compute enthalpy-weighted temperature for the mixture
gas = ct.Solution(MECH)
stream1 = ct.Quantity(gas, mass_flow_rate=mdot_rich_ex, T=T_rich_ex, P=P_rich_ex, X=X_rich_ex)
stream2 = ct.Quantity(gas, mass_flow_rate=mdot_air, T=T_air, P=P_rich_ex, X=X_air)
# Mix the streams at constant pressure and enthalpy
X_mix1, T_mix1, mdot_mix1 = mix_streams_const_HP(stream1, stream2, P_rich_ex)
gas.TPX = T_mix1, P_rich_ex, X_mix1
stream3 = ct.Quantity(gas, mass_flow_rate=mdot_mix1, T=T_mix1, P=P_rich_ex, X=X_mix1)

# Set up 0D reactor
reactor = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([reactor])

for mdot_h2o_i in mdot_h2o:
    NOx_row = []
    temps_row = []
    time = 0.0
    # Recompute mixture for this mdot_h2o_i
    stream3 = ct.Quantity(gas, mass_flow_rate=mdot_h2o_i, T=T_h2o, P=P_rich_ex, X=X_h2o)
    # Mix the streams at constant pressure and enthalpy
    X_mix, T_mix, mdot_mix = mix_streams_const_HP(stream1, stream3, P_rich_ex)
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
"""