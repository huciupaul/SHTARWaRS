import numpy as np
from typing import Tuple

"""
This file comprises global constants used throughout the project
"""

# General constants
R_AIR = 287.05                                          # [J/(kg·K)]
G_0 = 9.80665                                           # [m/s²]
LAPSE = -0.0065                                         # [K/m]
k_air = 1.4                                             # [-] Specific heat ratio of air

# Sea level atmosphere constants
T_sl = 288.15                                           # [K] Sea level temperature
P_sl = 101_325.0                                        # [Pa] Sea level pressure
rho_sl = 1.225                                          # [kg/m³] Sea level density

# ISA model
def isa_atmosphere(h: np.ndarray, T_sl: float = T_sl, P_sl: float = P_sl) -> Tuple[float, float, float]:
    """Thin-layer ISA (no tropopause).
    Returns T [K], P [Pa], rho [kg/m^3], and speed of sound."""
    T = T_sl + LAPSE * h
    P = P_sl * (T / T_sl) ** (-G_0 / (LAPSE * R_AIR))
    rho = P / (R_AIR * T)
    a = np.sqrt(k_air * R_AIR * T)  # speed of sound
    return T, P, rho, a

# Original Aircraft geometry constants
C_d_0_orig = 0.023
AR, e = 10.808, 0.735809
k = 1 / (np.pi * AR * e)
S = 28.79                                               # [m²]
MTOW_orig = 7_765   * G_0                               # [N] Original Maximum Take-Off Weight

#Modified Aircraft MTOW(requirement)
MTOW_mod  = 8_037.6 * G_0                               # [N] Modified Maximum Take-Off Weight  

# Original aircraft performance constants
V_TO, V_L  = 54.0167, 54.0167
V_stall_L  = 46.3
V_cruise   = 144.044
rho_cruise = 0.550798
ROC, V_climb = 10.16, 73.43


#######################################################################################
#################### Code specific constants from each part ###########################

###########
# From FC #
###########



####################
# From performance #
####################
A_inlet = 2*0.038410516                                 # [m^2] Inlet area
MAXC = 2*906e3                                          # [W] Max Continuous Power

eff_prop   = 0.80
mu_TO      = 0.04
S_g, S_land = 1163, 851
LW_max     = 7_604 * G_0
TOGA = 1908*1000 + 8.4 * 1000                           # [W] Take-Off/Go-Around power including auxiliary power(8.4 kW)

# Dynamic pressures
q_climb = 0.5 * rho_sl     * V_climb**2
q_c     = 0.5 * rho_cruise * V_cruise**2


################
# From storage #
################

dormancy = 24                                           #[hours]
stratification_factor = 2
k_str = 1.9                                             #[W/mK]
MAWP_global = 600000                                    #[Pa]
p_vent_global = 600000                                  #[Pa]
p_sl = 101325                                           #[Pa]

r_in = 0.75                                             #[m]

Q_leak_min = 10                                         #W (determined from Nicolas's thesis)
Q_leak_max = 1000                                       #W (determined from Nicolas's thesis)

T_hot = 40                                              #[Celsius] conservative for 

gf = 'S-Glass fiber'
fos = 2/3                                               #factor of safety for the tensile strength of the glass fiber
gf_density, gf_tensile_strength, gf_thermal_cond, gf_thermal_emis, gf_co2, gf_ee, gf_fvf = 1905, 1730*1e6*fos, 0.745, 0.95, 7.22, 116.5, 0.675

# Thesis values
Q_original_str = 0.4                                    #W
kevlar_thermal_cond = 1.9                               #W/mK
mass_original_str = 2.1                                 #kg
mass_originalg_lh2 = 6.2                                #kg
gravimetric_index = 0.35                                #from NASA report
kevlar_co2 = 13.1                                       #kg/kg (Kevlar 149)
kevlar_emb_energy = 257                                 #MJ/kg (Embodied Energy for Kevlar 149)

t_min = 0.001                                           #m (minimum thickness of the tank wall)

mli_density = 7900                                      #kg/m^3 https://www.sciencedirect.com/science/article/pii/S135943112200391X
mli_emis = 0.21                                         #https://www.thermalengineer.com/library/effective_emittance.htm
vacuum_thermal_cond = 0.015*1e-1                        #3 # W/mK https://www.researchgate.net/publication/321219004_Cylindrical_Cryogenic_Calorimeter_Testing_of_Six_Types_of_Multilayer_Insulation_Systems
mli_thermal_cond = 17.4                                 # W/mK  https://www.sciencedirect.com/science/article/pii/S135943112200391X
mli_ss_co2 = 3                                          #kg/kg (for SS)
mli_ss_ee = 42.74                                       #MJ/kg (Embodied Energy for SS)
mli_layers = 40
mli_thickness = 0.03 *1e-3 * mli_layers                 #https://www.sciencedirect.com/science/article/pii/S135943112200391X