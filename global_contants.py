import numpy as np

"""
This file comprises global constants used throughout the project
"""

# General constants
R_AIR = 287.05                                      # [J/(kg·K)]
G_0 = 9.80665                                       # [m/s²]
LAPSE = -0.0065                                     # [K/m]
k_air = 1.4                                         # [-] Specific heat ratio of air

# Sea level atmosphere constans
T_sl = 288.15                                      # [K] Sea level temperature
P_sl = 101_325.0                                   # [Pa] Sea level pressure
rho_sl = 1.225                                     # [kg/m³] Sea level density

# Original Aircraft geometry constants
C_d_0_orig = 0.023
AR, e = 10.808, 0.735809
k = 1 / (np.pi * AR * e)
S = 28.79                                           # [m²]
MTOW_orig = 7_765   * G_0                           # [N] Original Maximum Take-Off Weight
#Modified Aircraft MTOW(requirement)
MTOW_mod  = 8_037.6 * G_0                           # [N] Modified Maximum Take-Off Weight  

# Orifinal aircraft performance constants
V_TO, V_L  = 54.0167, 54.0167
V_stall_L  = 46.3
V_cruise   = 144.044
rho_cruise = 0.550798
ROC, V_climb = 10.16, 73.43


#######################################################################################

###########
# From FC #
###########



####################
# From performance #
####################
A_inlet = 2*0.038410516                             # [m^2] Inlet area
MAXC = 2*906e3                                      # [W] Max Continuous Power

eff_prop   = 0.80
mu_TO      = 0.04
S_g, S_land = 1163, 851
LW_max     = 7_604 * G_0
TOGA = 1908*1000

# Dynamic pressures
q_climb = 0.5 * rho_sl     * V_climb**2
q_c     = 0.5 * rho_cruise * V_cruise**2


################
# From storage #
################
