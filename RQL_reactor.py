import cantera as ct
import numpy as np
import scipy.integrate


class ReactorOde:
    def __init__(self, gas):
        # Parameters of the ODE system and auxiliary data are stored in the
        # ReactorOde object.
        self.gas = gas
        self.P = gas.P

    def __call__(self, t, y):
        """the ODE function, y' = f(t,y) """
        self.gas.set_unnormalized_mass_fractions(y[1:])
        self.gas.TP = y[0], self.P
        rho = self.gas.density

        wdot = self.gas.net_production_rates
        heat_source = 1e6 if t < 1e-5 else 0  # Energy input during the first 10 microseconds
        dTdt = - (np.dot(self.gas.partial_molar_enthalpies, wdot) /
                  (rho * self.gas.cp)) + heat_source / (rho * self.gas.cp)
        dYdt = wdot * self.gas.molecular_weights / rho

        return np.hstack((dTdt, dYdt))


gas = ct.Solution('gri30.yaml')

# Initial condition
P = 12.5 * ct.one_atm
T_spark = 630  # K
# Hydrogen as fuel and air with O2:21, N2:79
gas.TPX = T_spark, P, 'H2:2,O2:1,N2:3.76'
y0 = np.hstack((gas.T, gas.Y))

# Set up objects representing the ODE and the solver
ode = ReactorOde(gas)
solver = scipy.integrate.ode(ode)
solver.set_integrator('vode', method='bdf', with_jacobian=True)
solver.set_initial_value(y0, 0.0)

# Integrate the equations, keeping T(t) and Y(k,t)
t_end = 1e-3
states = ct.SolutionArray(gas, 1, extra={'t': [0.0]})
dt = 1e-5
temperature_tolerance = 1e-3  # Steady-state temperature change threshold (K)
previous_temperature = gas.T

while solver.successful() and solver.t < t_end:
    solver.integrate(solver.t + dt)
    gas.TPY = solver.y[0], P, solver.y[1:]
    states.append(gas.state, t=solver.t)

    # Check for steady-state temperature
    #temperature_change = abs(gas.T - previous_temperature)
    #if temperature_change < temperature_tolerance:
    #    print(f"Steady-state temperature reached at t = {solver.t:.6f} s")
    #    break
    #previous_temperature = gas.T

# Plot the results
try:
    import matplotlib.pyplot as plt

    # Plot temperature vs. time
    plt.figure()
    plt.plot(states.t, states.T, color='r', label='Temperature', lw=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature vs. Time')
    plt.legend()
    plt.grid()
    plt.show()

    # Plot concentration of H2O, NOx, and other intermediate species vs. time
    plt.figure()
    plt.plot(states.t, states('H2O').Y, label='H2O', lw=2)
    plt.plot(states.t, states('NO').Y, label='NO', lw=2)
    plt.plot(states.t, states('NO2').Y, label='NO2', lw=2)
    plt.plot(states.t, states('OH').Y, label='OH', lw=2)
    plt.plot(states.t, states('H').Y, label='H', lw=2)
    plt.plot(states.t, states('O').Y, label='O', lw=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Mass Fraction')
    plt.title('Concentration of H2O, NOx, and Intermediate Species vs. Time')
    plt.legend()
    plt.grid()
    plt.show()

except ImportError:
    print('Matplotlib not found. Unable to plot results.')