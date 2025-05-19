from ADRpy import constraintanalysis as ca       # you kept this for later use
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# 1. Constants and aircraft data
# ------------------------------------------------------------------
g        = 9.8056          # [m/s²]
rho_SL   = 1.225           # sea-level density [kg/m³]

# ---- Beechcraft 1900D parameters ---------------------------------
C_d_0    = 0.024
AR       = 10.808
e        = 0.735809
k        = 1 / (np.pi * AR * e)

TOGA     = 1_908_000       # total take-off power available [W]
S        = 28.79           # wing area [m²]

MTOW     = 7_765 * g       # [N]
W_fuel   = 2_022 * g
LW_max   = 7_604 * g       # max landing weight [N]
LW       = MTOW - W_fuel

V_TO     = 54.0167         # take-off speed [m/s]
V_L      = 54.0167         # approach / landing speed [m/s]

rho_cruise = 0.550798
V_stall_L  = 46.3          # stall speed (TO flaps, LW_max) [m/s]
V_stall_c  = 51.9589       # clean stall speed at MTOW [m/s]
V_cruise   = 144.044       # cruise TAS [m/s]
LD_cruise  = 15.07         # L/D at cruise (unused here)
ROC        = 6.35          # required climb rate [m/s]

V_climb   = 73.43
rho_climb = rho_SL
q_climb   = 0.5 * rho_climb * V_climb**2

sigma      = rho_cruise / rho_SL
q_c        = 0.5 * rho_cruise * V_cruise**2
eff_prop   = 0.80

mu_TO      = 0.04
S_g        = 1163          # take-off ground run [m]
S_land     = 851           # landing distance [m]

# ------------------------------------------------------------------
# 2. Aerodynamic coefficients
# ------------------------------------------------------------------
CL_max_TO = MTOW / (0.5*rho_SL*V_stall_L**2*S)
CL_TO     = MTOW / (0.5*rho_SL*V_TO**2       *S)

CL_max_L  = LW_max / (0.5*rho_SL*V_stall_L**2*S)
CL_max    = max(CL_max_TO, CL_max_L)

print(f"CL_max (take-off): {CL_max_TO:.3f}")
print(f"CL_max (landing): {CL_max_L :.3f}")

# ------------------------------------------------------------------
# 3. Wing-loading grid
# ------------------------------------------------------------------
WS = np.linspace(500, 6000, 200)        # [N/m²]

# ------------------------------------------------------------------
# 4. Constraint equations
# ------------------------------------------------------------------
# 4.1 Take-off field length
C_D_TO    = C_d_0 + k * CL_TO**2
TW_TO     = (1.21 / (g * rho_SL * CL_max_TO * S_g)) * WS \
          + (0.605 / CL_TO) * (C_D_TO - mu_TO * CL_TO) + mu_TO
PW_TO     = TW_TO * V_TO / eff_prop

# 4.2 Stall limit  (gives a *vertical* WS limit)
WS_stall  = 0.5 * rho_SL * V_stall_L**2 * CL_max_TO

# 4.3 Cruise
TW_cruise = (q_c * C_d_0) / WS + (k / q_c) * WS
PW_cruise = TW_cruise * V_cruise / eff_prop

# 4.4 Climb
TW_climb  = ROC / V_climb + (q_climb * C_d_0) / WS + k * WS / q_climb
PW_climb  = TW_climb * V_climb / eff_prop

# 4.5 Landing distance (vertical line)
WS_land_limit = 0.5 * rho_SL * V_L**2 * CL_max_L

# ------------------------------------------------------------------
# 5. Actual aircraft point (optional reference)
# ------------------------------------------------------------------
PW_beechcraft = TOGA / MTOW           # “available” P/W at MTOW
WS_beechcraft = MTOW / S
print(f"Beech 1900D  P/W: {PW_beechcraft:.3f}  W/S: {WS_beechcraft:.1f}")

# ------------------------------------------------------------------
# 6. DESIGN point (intersection of stall-WS with climb curve)
# ------------------------------------------------------------------
PW_design = np.interp(WS_stall, WS, PW_climb)

# ------------------------------------------------------------------
# 7. PLOT
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

# 7.1 Constraint curves
ax.plot(WS, PW_TO,      label='Take-off')
ax.plot(WS, PW_cruise,  label='Cruise')
ax.plot(WS, PW_climb,   label='Climb requirement', ls='--', c='green')

# 7.2 Vertical limits
ax.axvline(WS_stall,      ls='--', c='red',    label='Stall limit')
ax.axvline(WS_land_limit, ls='--', c='purple', label='Landing limit')

# 7.3 Feasible region shading  (only up to WS_stall)
mask        = WS <= WS_stall
P_required  = np.maximum.reduce([PW_TO[mask],
                                 PW_cruise[mask],
                                 PW_climb[mask]])
_, y_top = ax.get_ylim()        # current y-axis upper bound
ax.fill_between(WS[mask], P_required, y_top,
                color='lightgreen', alpha=0.35,
                zorder=0, label='Feasible region')

# 7.4 Design point
ax.scatter(WS_stall, PW_design, marker='x', s=100, c='k', lw=2,
           label='Design point')

# 7.5 Cosmetics
ax.set_xlim(WS.min(), WS.max())
ax.set_ylim(0)                      # y-axis starts at zero
ax.set_xlabel('Wing loading  W/S  [N/m²]')
ax.set_ylabel('Power-to-weight  P/W  [W/N]')
ax.set_title('Constraint diagram – Beechcraft 1900D')
ax.grid(True, which='both', ls=':')
ax.legend(loc='upper right')
plt.tight_layout()
plt.show()
