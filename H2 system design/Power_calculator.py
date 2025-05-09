import numpy as np

# Constants:
V_TO = 54.0167
rho_c = 0.550798
V_c = 144.044
S = 28.79
C_D_0 = 0.021
AR = 10.80825287
e = 0.7358093584
k = 0.0400248141
b = 17.64
propulsive_eff = 0.85
TOGA_0 = 1.908 * 1e6

# Inputs:
W_0 = 7765 * 9.81  # N, initial weight of the aircraft
W_new = 7756 * 9.81  # N, new weight of the aircraft

a = 0.5 * rho_c * V_c**2 * S * C_D_0  # N/m^2, thrust required for takeoff
b = 1 / (0.5 * rho_c * V_c**2 * S * np.pi * e * AR)  # N/m^2, thrust required for takeoff

P_cr_0 = (a + b * W_0**2 ) * V_c / propulsive_eff
P_cr_new = (a + b * W_new**2) * V_c / propulsive_eff

print(f"Original cruise power: {P_cr_0 * 1e-6} MW")
print(f"New cruise power: {P_cr_new * 1e-6} MW")

