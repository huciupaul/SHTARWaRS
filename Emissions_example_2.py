import cantera as ct
import numpy as np

import matplotlib.pyplot as plt

# Operating conditions
Pi = 30 * 101325  # Pascal
Ti = 900.0  # Kelvin
AFR_st = 9.52  # Stoichiometric air-fuel ratio (CH4-air mixture)

# Initialize reservoirs and gases
# Fuel inlet to rich-burn combustor
fuel = ct.Solution('gri30.yaml')
fuel.TPX = Ti, Pi, 'CH4:1.0'
fuel_in = ct.Reservoir(fuel)

# Oxidizer inlet
oxi = ct.Solution('gri30.yaml')
oxi.TPX = Ti, Pi, 'O2:0.21,N2:0.79'
air_in = ct.Reservoir(oxi)

# Reacting gas that acts as igniter
gas = ct.Solution('gri30.yaml')
gas.TPX = 3500.0, Pi, 'CH4:1.0,O2:2'

# Rich-burn combustor
combustor1 = ct.IdealGasReactor(gas)
combustor1.volume = 0.02

# Lean-burn combustor
combustor2 = ct.IdealGasReactor(gas)
combustor2.volume = 0.02

# Final exhaust reservoir
exhaust = ct.Reservoir(gas)

# Initialize flow rates and flow directions
# Fuel and air mass flow rates, kg/s
mdot_fuel = 3.0  # Fuel rate to rich-burn combustor
mdot_air_total = 73.5  # Total air flow rate (will be split into two)

# Fuel flow direction to rich-burn PSR
m1 = ct.MassFlowController(fuel_in, combustor1, mdot=mdot_fuel)

# Air flow direction to rich-burn PSR
m2 = ct.MassFlowController(air_in, combustor1)

# Air flow direction to lean-burn PSR
m3 = ct.MassFlowController(air_in, combustor2)

# Initialize valves between reservoirs
between_combustors_valve = ct.Valve(combustor1, combustor2)
between_combustors_valve.valve_coeff = 40.0

to_exhaust_valve = ct.Valve(combustor2, exhaust)
to_exhaust_valve.valve_coeff = 1.0

# Initialize reactor network
network = ct.ReactorNet([combustor1, combustor2])

# Compute final temperature and composition for different equivalence ratios in rich zone
n_points = 41
vector_phi_rb = np.linspace(1.0, 3.0, n_points)

x_NO = []
x_CO = []
x_CH4 = []
vector_temperature = []

for phi_rb in vector_phi_rb:
    # Reset reactors' network to initial state
    network.set_initial_time(0.0)

    # Calculating variable air flow rates
    mdot_air_rb = AFR_st / phi_rb * mdot_fuel
    mdot_air_lb = mdot_air_total - mdot_air_rb

    # Set computed air-flow rates to combustors
    m2.mdot = mdot_air_rb
    m3.mdot = mdot_air_lb

    # Solve ODEs
    t = 0.0
    dt = 5.0e-5
    stop = 0
    temp = []

    while stop < 100:
        t += dt
        network.advance(t)
        temp.append(combustor2.T)

        # Iterates until the temperature changes by less than 1e-6 for 100 consecutive time steps
        if len(temp) > 1 and abs(temp[-1] - temp[-2]) < 1e-6:
            stop += 1

    # Output
    x_NO.append(combustor2.thermo['NO'].X[0])
    x_CO.append(combustor2.thermo['CO'].X[0])
    x_CH4.append(combustor2.thermo['CH4'].X[0])
    vector_temperature.append(combustor2.T)

# Plots
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(vector_phi_rb, vector_temperature, '--.', linewidth=1.2)
plt.xlabel('Equivalence Ratio')
plt.ylabel('Temperature (K)')
plt.grid()

plt.subplot(2, 2, 2)
plt.semilogy(vector_phi_rb, x_NO, '--.', linewidth=1.2)
plt.xlabel('Equivalence Ratio')
plt.ylabel('Mole Fraction NO')
plt.grid()

plt.subplot(2, 2, 3)
plt.semilogy(vector_phi_rb, x_CO, '--.', linewidth=1.2)
plt.xlabel('Equivalence Ratio')
plt.ylabel('Mole Fraction CO')
plt.grid()

plt.subplot(2, 2, 4)
plt.semilogy(vector_phi_rb, x_CH4, '--.', linewidth=1.2)
plt.xlabel('Equivalence Ratio')
plt.ylabel('Mole Fraction CH4')
plt.grid()

plt.tight_layout()
plt.show()
