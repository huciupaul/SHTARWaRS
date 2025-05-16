from ADRpy import constraintanalysis as ca
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

# Example: Constraint Analysis for a Small Turboprop Aircraft (based on G-OUAV ADRpy notebook)

# Constants and Parameters
g = 9.8056                              # gravity [m/s^2]
rho_SL = 1.225                        # sea level air density [kg/m^3]

S = 28.79                             # wing surface area [m^2]
MTOW = 7765*g                         # MTOW [N]
W_fuel = 2022*g
LW_max = 7604*g                           # max landing weight [N]  
LW = MTOW-W_fuel
V_TO = 54.0167                        # take-off speed [m/s]
V_L =  60.2                           # landing speed  [m/s]

rho_cruise = 0.550798                 # assumed cruise altitude density [kg/m^3]
V_stall = 46.3                    # stall speed [m/s]  -->  clean, max T-O weight
V_cruise = 144.044                    # cruise speed [m/s]
LD_cruise = 15.07                     # Lift-to-Drag ratio at cruise
ROC = 6.35                            # rate of climb [m/s]

mu_TO = 0.03                          # ground friction coeff  --> source??????
S_g = 1163                            # ground run length [m]  -->  T-O flap setting
S_land = 851                          # required landing distance [m]  -->  at max landing W

CL_max_TO = 1.8
CL_max_L = 2.4

# why CL max TO > CL max L ?
#CL_max_TO = MTOW/(0.5*rho_SL*V_TO**2*S)                 # max lift coefficient in takeoff 
#CL_max_L = LW_max/(0.5*rho_SL*V_L**2*S)                 # max lift coefficient in landing    
#print("CL max for take off is:", CL_max_TO)
#print("CL max for landing is:", CL_max_L)

# Wing loading range [N/m^2]
# --> S = cte, W decreases
# MAX = MTOW*g/S = 2645.87183
# MIN = (MTOW-MAXF)/S = (7765-2022)*g/S = 1956.8888
#WS = np.linspace(1950, 2646, 5) 
WS = np.linspace(500, 5000, 50)  

# Constraint 1: Takeoff field length
TW_takeoff = (WS / (rho_SL * CL_max_TO)) / (S_g * g) 

# Constraint 2: Stall speed (converted to max WS)
WS_stall = 0.5 * rho_SL * V_stall**2 * CL_max_L

# Constraint 3: Cruise condition
TW_cruise = WS / (rho_cruise * V_cruise**2 * LD_cruise)

# Constraint 4: Climb performance (using simplified excess power method)
TW_climb = ROC / V_cruise + 1 / LD_cruise

# From Roskam: WS < K * S_land * rho * CL_max_L
K_land = 0.3 * g  
WS_land_limit = rho_SL * CL_max_L * K_land * S_land  # Max allowed wing loading for landing


# Feasible region
TW_envelope = np.maximum(TW_takeoff, TW_cruise)
TW_envelope = np.maximum(TW_envelope, TW_climb)
WS_feasible = WS[WS <= WS_stall]
TW_feasible = TW_envelope[WS <= WS_stall]

# Plot
plt.figure(figsize=(10, 6))
plt.fill_between(WS_feasible, TW_feasible, y2=1.0, color='lightgreen', alpha=0.4, label="Feasible region")

plt.plot(WS, TW_takeoff, label="Takeoff")
plt.axvline(WS_stall, color='r', linestyle='--', label="Stall limit")
plt.plot(WS, TW_cruise, label="Cruise")
plt.axvline(WS_land_limit, color='purple', linestyle='--', label="Landing limit")
plt.axhline(TW_climb, color='g', linestyle='--', label="Climb requirement")

plt.xlabel("Wing Loading W/S [N/m²]")
plt.ylabel("Thrust-to-Weight T/W")
plt.title("Constraint Diagram for Small Turboprop Aircraft")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
