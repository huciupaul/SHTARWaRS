import cantera as ct
import numpy as np
import time

import matplotlib.pyplot as plt

def reactor1(fuel):
    """
    REACTOR1 Zero-dimensional kinetics: adiabatic, constant pressure.

    This example illustrates how to use Cantera's Reactor class for
    zero-dimensional kinetics simulations. The reactor is adiabatic
    and very close to constant pressure.
    """

    # Initial parameters
    P = 101325  # Pressure in Pascal
    T = 1125  # Temperature in Kelvin
    dt = 0.5e-5  # Time step in seconds
    TotalTime = 0.05  # Total simulation time in seconds
    nSteps = int(np.ceil(TotalTime / dt))  # Number of steps

    # Select the gas model based on the fuel type
    if fuel == 'kerosene':
        gas = ct.Solution('kerosene.cti')
    elif fuel == 'ethylene':
        gas = ct.Solution('gri30.cti')
    else:
        gas = ct.Solution('gri30.cti')

    # Set the initial conditions for ethylene (phi = 0.42)
    gas.TPX = T, P, 'C2H4:1, O2:7.143, N2:26.86'

    # Create a reactor and insert the gas
    reactor = ct.IdealGasReactor(gas)

    # Create a reservoir to represent the environment
    air = ct.Solution('air.yaml')
    air.TP = T, P
    env = ct.Reservoir(air)

    # Define a wall between the reactor and the environment
    wall = ct.Wall(left=reactor, right=env)
    wall.expansion_rate_coeff = 1.0e6  # dV/dt = KA(P_1 - P_2)
    wall.area = 1.0

    # Create a reactor network and insert the reactor
    sim = ct.ReactorNet([reactor])

    # Initialize arrays to store results
    tim = np.zeros(nSteps)
    temp = np.zeros(nSteps)
    x = np.zeros((nSteps, 7))
    kero = np.zeros((nSteps, 3))

    # Time integration
    t = 0.0
    t0 = time.process_time()
    for n in range(nSteps):
        t += dt
        sim.advance(t)
        tim[n] = sim.time
        temp[n] = reactor.T
        x[n, :] = gas['CH4', 'CO', 'CO2', 'H2O', 'NO', 'NO2', 'H2'].Y
        kero[n, :] = gas['NC10H22', 'PHC3H7', 'CYC9H18'].Y

    print(f'CPU time = {time.process_time() - t0:.2f} seconds')

    # Plot results
    plt.figure(figsize=(10, 8))

    plt.subplot(2, 2, 1)
    plt.plot(tim, temp, 'r', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')

    plt.subplot(2, 2, 2)
    plt.plot(tim, np.sum(kero, axis=1), 'k', label='Kerosene')
    plt.plot(tim, x[:, 2], 'r', label='CO2')
    plt.plot(tim, x[:, 3], 'b', label='H2O')
    plt.plot(tim, x[:, -1], 'g', label='H2')
    plt.xlabel('Time (s)')
    plt.ylabel('Mass fraction')
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(tim, x[:, 1])
    plt.xlabel('Time (s)')
    plt.ylabel('CO Mass Fraction')

    plt.subplot(2, 2, 4)
    plt.plot(tim, (x[:, 4] + x[:, 5]) * 1e6)
    plt.xlabel('Time (s)')
    plt.ylabel('NOX Mass Fraction (ppm)')

    plt.tight_layout()
    plt.show()

# Example usage
reactor1('ethylene')
