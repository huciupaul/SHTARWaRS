import CoolProp.CoolProp as CP # Used to get properties of Hydrogen
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np

#constant inputs:

dormancy = 24 #hours
mass_h2 = 200 #kg
fill_ratio = 0.8 #percentage of LH2 at the start

MAWP = [100000] # Convert from Bar to Pa for the venting pressure

T0 = CP.PropsSI('T', 'P', MAWP, 'Q', 0.1, 'ParaHydrogen') # Initial temperature in K
rho_lh2 = CP.PropsSI('D', 'T', T0, 'Q', 0, 'ParaHydrogen') # density of liquid hydrogen
rho_gh2 = CP.PropsSI('D', 'T', T0, 'Q', 1, 'ParaHydrogen') # density of gas hydrogen

#inner volume calculation
V_in = mass_h2 / (fill_ratio * rho_lh2) # Volume in m^3

L_max = 2 #m

def volume_equation(r):
    L_cyl = L_max - 2*r  # cylindrical length
    if L_cyl < 0:
        return 1e6  # Penalize invalid cases
    V = np.pi * r**2 * L_cyl + (4/3) * np.pi * r**3  # cylinder + 2 hemispheres = 1 sphere
    return V - V_in  # we want this to be zero

# Solve for radius using a root-finding method
r_solution = opt.root_scalar(volume_equation, bracket=[0.001, L_max/2], method='brentq')

if r_solution.converged:
    R_in = r_solution.root
    print(f"Calculated radius: {R_in:.4f} meters")
else:
    print("Failed to find an inner radius.")

