import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


mechanism_1 = 'gri30.yaml'  # Mechanism for H2 combustion with NOx
mechanism_2 = 'CapursoMechanism.yaml'  # Mechanism for H2 combustion with NOx
gas = ct.Solution(mechanism_2)
gas.TP = 630.25, 12.5 * ct.one_atm
gas.set_equivalence_ratio(2.5, 'H2', 'O2:21.0, N2:79.0')
gas_copy = gas
species_to_track = ['H2', 'O2', 'H2O', 'NO', 'NO2']
species_data = {sp: [] for sp in species_to_track if sp in gas.species_names}

reactor = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([reactor])

# Set up for data collection
time = 0.0
t_end = 0.01
dt = 1e-6

times = []
temps = []

# Parameters for spark
IDT = 0.002  # Time when spark occurs [s]
volume = reactor.volume
cv = gas.cp_mass  # Approximate specific heat at constant pressure

gas_copy.equilibrate('HP')
T_igition = gas_copy.T  # Set ignition temperature slightly above equilibrium temperature
print(f'Ignition temperature: {T_igition:.2f} K')

# Loop over time
while time < t_end:
    time += dt

    # Simulate spark by increasing temperature (energy deposition)
    if abs(time - IDT) < 1e-6:
        T_current = reactor.T
        mass = reactor.mass
        gas.TP = T_igition, reactor.thermo.P # Set the gas to the ignition temperature
        reactor.syncState()
        

    print(f'Temp: ', gas.T)
    sim.advance(time)
    times.append(time)
    temps.append(reactor.T)
    for sp in species_data:
        species_data[sp].append(reactor.thermo[sp].X[0])


print(f'Initial temp: ', temps[0], '; Final temp: ', temps[-1])
print("\nFinal Species Content (mole fractions):")
for sp, values in species_data.items():
    print(f"{sp}: {values[-1]:.5e}")

# Plot temperature rise (ignition)
plt.plot(times, temps)
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.title('Spark Ignition via Energy Deposition')
plt.grid(True)
plt.show()

# Plot tracked species content
plt.figure(figsize=(10, 6))
for sp, values in species_data.items():
    plt.plot(times, values, label=sp)
plt.xlabel('Time [s]')
plt.ylabel('Mole Fraction')
plt.title('Tracked Species Content Over Time')
plt.grid(True)
plt.legend()
plt.show()







""""
# --------------------
# Parameters
# --------------------
P = 12.1 * ct.one_atm         # Pressure [Pa]
T_initial = 603.25            # Initial temperature [K]
T_spark = 4000.0              # Spark reservoir temperature [K]
phi = 2.5                     # Equivalence ratio
# mechanism = 'H2_pNOX_15_94_TC.yaml' 
mechanism = 'gri30.yaml'
# mechanism = 'H2_pNOX_15_94_TC.yaml'  # Mechanism for H2 combustion with NOx
# mechanism = 'CapursoMechanism.yaml'  # Mechanism for H2 combustion with NOx

# --------------------
# Setup gas mixture
# --------------------
gas = ct.Solution(mechanism)
oxidizer = 'N2:0.78084, O2:0.20946'  # Air composition
fuel = 'H2:1.0'  # Hydrogen composition
gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=oxidizer)
gas.TP = T_initial, P

# --------------------
# Reactor setup
# --------------------
# Main reacting gas reactor
reactor = ct.IdealGasConstPressureReactor(gas)

# Hot reservoir simulating the spark
spark_gas = ct.Solution(mechanism)
spark_gas.TP = T_spark, P
spark_reservoir = ct.Reservoir(spark_gas)

# Wall between spark reservoir and reactor — simulates heat input
wall = ct.Wall(left=reactor, right=spark_reservoir, A=1.0)

# Reactor network
sim = ct.ReactorNet([reactor])
# --------------------
# Time integration
# --------------------
time = 0.0
dt = 1e-6
max_time = 1  # seconds
IDT = 0.01  # Ignition delay time [s]

times = []
temperatures = []
species_to_track = ['H2', 'O2', 'H2O', 'NO', 'NO2']
species_data = {sp: [] for sp in species_to_track if sp in gas.species_names}

while time < max_time:
    time += dt
    wall.heat_transfer_coeff = 2e10 if (time >= IDT and time < IDT + dt) else 0.0  # High heat flux for a short time
    sim.advance(time)

    times.append(time)
    temperatures.append(reactor.T)
    for sp in species_data:
        species_data[sp].append(reactor.thermo[sp].X[0])

    # # Stop early if temperature stabilizes
    # if time > 0.01 and abs(temperatures[-1] - temperatures[-2]) < 1e-4:
    #     break

# --------------------
# Results
# --------------------
print(f"\nFinal Temperature: {reactor.T:.2f} K")
print("Final Composition (mole fractions):")
for sp in species_data:
    print(f"  {sp}: {reactor.thermo[sp].X[0]:.5f}")

# --------------------
# Plotting
# --------------------
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(times, temperatures, color='red', label='Temperature [K]')
plt.ylabel('Temperature [K]')
plt.grid()
plt.legend()

plt.subplot(2, 1, 2)
for sp, values in species_data.items():
    plt.plot(times, values, label=sp)
plt.xlabel('Time [s]')
plt.ylabel('Mole Fraction')
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()
"""