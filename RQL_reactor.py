import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

"""
This script models spark ignition of H2 and air.
Tracked params: temperature, H2, O2, H2O, NO, NO2
"""

#------Mechanism Selection------
mechanism_1 = 'gri30.yaml'  # Mechanism for H2 combustion with NOx
mechanism_2 = 'CapursoMechanism.yaml'  # Mechanism for H2 combustion with NOx
mechanism_3 = 'H2_pNOX_15_94_TC.yaml'  # Mechanism for H2 combustion with NOx

#------Parameters------
P = 12.1 * ct.one_atm         # Pressure [Pa]
T_initial = 603.25            # Initial temperature [K]
gas = ct.Solution(mechanism_1)
gas.TP = T_initial, P
gas.set_equivalence_ratio(2.5, 'H2', 'O2:2.1, N2:7.9')
gas_copy = gas
species_to_track = ['H2', 'O2', 'H2O', 'NO', 'NO2']
species_data = {sp: [] for sp in species_to_track if sp in gas.species_names}

#------Reactor Setup------
reactor = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([reactor])

#------Iteration Parameters------
time = 0.0
t_end = 0.01
dt = 1e-6

times = []
temps = []

IDT = 0.002  # ignition delay time [s]
volume = reactor.volume

#------Ignition Temperature------
gas_copy.equilibrate('HP')
T_igition = gas_copy.T  # Set ignition temperature slightly above equilibrium temperature

#------Simulation Loop------
while time < t_end:
    time += dt

    # Simulate spark by increasing temperature (energy deposition)
    if abs(time - IDT) < dt:
        gas.TP = T_igition, reactor.thermo.P # Set the gas to the ignition temperature
        reactor.syncState()
        
    sim.advance(time)
    times.append(time)
    temps.append(reactor.T)
    for sp in species_data:
        species_data[sp].append(reactor.thermo[sp].X[0])


#------Results------
print(f'Ignition temperature: {T_igition:.2f} K')
print(f'Initial temp: ', temps[0], '; Final temp: ', temps[-1])
print("\nFinal Species Content (mole fractions):")
for sp, values in species_data.items():
    print(f"{sp}: {values[-1]:.5e}")

# Combined plot for temperature rise and tracked species content
plt.figure(figsize=(10, 10))

# Plot temperature rise (ignition)
plt.subplot(2, 1, 1)
plt.plot(times, temps, color='red', label='Temperature [K]')
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.title('Temperature Rise Over Time')
plt.grid(True)
plt.legend()

# Plot tracked species content
plt.subplot(2, 1, 2)
for sp, values in species_data.items():
    plt.plot(times, values, label=sp)
plt.xlabel('Time [s]')
plt.ylabel('Mole Fraction')
plt.title('Tracked Species Content Over Time')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
