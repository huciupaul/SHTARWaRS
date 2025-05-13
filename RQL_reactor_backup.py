import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


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
