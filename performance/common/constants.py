import numpy as np
from CoolProp.CoolProp import PropsSI

R_AIR = 287.05                                      # [J/(kg·K)]
G_0 = 9.80665                                       # [m/s²]
LAPSE = -0.0065                                      # [K/m]
k_air = 1.4                                         # [-] Specific heat ratio of air
A_inlet = 2*0.038410516                             # [m^2] Inlet area
MAXC = 2*904e3                                      # [W] Max Continuous Power