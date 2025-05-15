import numpy as np

R_AIR = 287.05                                      # [J/(kg·K)]
G_0 = 9.80665                                       # [m/s²]
LAPSE = 0.0065                                      # [K/m]
k_air = 1.4                                         # [-] Specific heat ratio of air
A_inlet = 2*0.015446456799602548*np.pi*2.78**2/4    # [m^2] Inlet area
TOGA = 1.98e6                                       # [W] Take-off Go-around power