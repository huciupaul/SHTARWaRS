"""
Constraint diagram – Beechcraft 1900D
Baseline vs. modified (higher C_D0 and MTOW)

Colours follow your original palette.
New: green fill over the *extra* feasible W/S gained by the higher stall limit.
"""

# ------------------------------------------------------------------
# 1. Imports
# ------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# 2. Constants (aircraft-independent)
# ------------------------------------------------------------------
g        = 9.8056
rho_SL   = 1.225

# Geometry & aero
C_d_0_orig = 0.024
AR, e      = 10.808, 0.735809
k          = 1 / (np.pi * AR * e)
S          = 28.79                       # [m²]

# Speeds & mission data
V_TO, V_L  = 54.0167, 54.0167
V_stall_L  = 46.3
V_cruise   = 144.044
rho_cruise = 0.550798
ROC, V_climb = 10.16, 73.43

# Misc
eff_prop   = 0.80
mu_TO      = 0.04
S_g, S_land = 1163, 851
LW_max     = 7_604 * g

# Dynamic pressures
q_climb = 0.5 * rho_SL     * V_climb**2
q_c     = 0.5 * rho_cruise * V_cruise**2

# ------------------------------------------------------------------
# 3. Wing-loading grid
# ------------------------------------------------------------------
WS = np.linspace(500, 6000, 200)         # [N m⁻²]

# ------------------------------------------------------------------
# 4. Helper: curves & intersections
# ------------------------------------------------------------------
def constraint_curves(C_d_0, MTOW):
    CL_max_TO = MTOW / (0.5 * rho_SL * V_stall_L**2 * S)
    CL_TO     = MTOW / (0.5 * rho_SL * V_TO**2 * S)
    CL_max_L  = LW_max / (0.5 * rho_SL * V_stall_L**2 * S)

    # Take-off
    C_D_TO = C_d_0 + k * CL_TO**2
    TW_TO  = (1.21 / (g * rho_SL * CL_max_TO * S_g)) * WS \
           + (0.605 / CL_TO) * (C_D_TO - mu_TO * CL_TO) + mu_TO
    PW_TO  = TW_TO * V_TO / eff_prop

    # Cruise
    TW_cr  = (q_c * C_d_0) / WS + (k / q_c) * WS
    PW_cr  = TW_cr * V_cruise / eff_prop

    # Climb
    TW_cl  = ROC / V_climb + (q_climb * C_d_0) / WS + k * WS / q_climb
    PW_cl  = TW_cl * V_climb / eff_prop

    # Limits
    WS_stall = 0.5 * rho_SL * V_stall_L**2 * CL_max_TO
    WS_land  = 0.5 * rho_SL * V_L**2       * CL_max_L

    # Design pt - stall line on climb curve
    PW_stall_int = np.interp(WS_stall, WS, PW_cl)

    # Design pt - cruise ⋂ climb
    diff = PW_cl - PW_cr
    idx  = np.where(np.diff(np.sign(diff)))[0]
    if idx.size:
        i       = idx[0]
        ws1, ws2   = WS[i], WS[i+1]
        pwc1, pwc2 = PW_cl[i], PW_cl[i+1]
        pwr1, pwr2 = PW_cr[i], PW_cr[i+1]
        slope      = (pwc2 - pwc1) - (pwr2 - pwr1)
        ws_x       = ws1 + (pwr1 - pwc1) * (ws2 - ws1) / slope
        pw_x       = np.interp(ws_x, WS, PW_cl)
    else:
        ws_x, pw_x = np.nan, np.nan

    PW_req = np.maximum.reduce([PW_TO, PW_cr, PW_cl])

    return dict(PW_TO=PW_TO, PW_cr=PW_cr, PW_cl=PW_cl,
                WS_stall=WS_stall, WS_land=WS_land,
                PW_stall_int=PW_stall_int,
                WS_x=ws_x, PW_x=pw_x,
                PW_req=PW_req)

# ------------------------------------------------------------------
# 5. Baseline vs. modified
# ------------------------------------------------------------------
MTOW_orig = 7_765   * g
MTOW_mod  = 8_037.6 * g

cur_orig = constraint_curves(C_d_0_orig,            MTOW_orig)
cur_mod  = constraint_curves(C_d_0_orig * 1.17,     MTOW_mod)

# ------------------------------------------------------------------
# 6. Plot
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

# Baseline curves
ax.plot(WS, cur_orig["PW_TO"],  label="Take-off",          lw=2)          # blue
ax.plot(WS, cur_mod["PW_TO"],   label="Take-off (mod)",    lw=2,
        ls=":",   c="blue")
ax.plot(WS, cur_orig["PW_cr"],  label="Cruise",            lw=2)          # orange
ax.plot(WS, cur_mod["PW_cr"],   label="Cruise (mod)",      lw=2,
        ls=":",   c="orange")
ax.plot(WS, cur_orig["PW_cl"],  label="Climb requirement", lw=2,  c="green")
ax.plot(WS, cur_mod["PW_cl"],   label="Climb req. (mod)",  lw=2,
        ls=":",  c="green")

# Vertical limits
ax.axvline(cur_orig["WS_stall"], c="red",    label="Stall limit")
ax.axvline(cur_mod ["WS_stall"], ls=":", c="red",
           label="Stall limit (mod)")
ax.axvline(cur_orig["WS_land"], c="purple", label="Landing limit ")
ax.axvline(cur_mod ["WS_land"], ls=":",  c="purple")

# ---------------- Feasible regions ---------------------------------
y_max = ax.get_ylim()[1]

# Baseline feasible (WS ≤ original stall)
mask_base = WS <= cur_orig["WS_stall"]
ax.fill_between(WS[mask_base], cur_orig["PW_req"][mask_base], y_max,
                color="lightgreen", alpha=0.35, zorder=0,
                label="Feasible region ")

# Extra feasible W/S gained (between the two stall limits)
mask_gain = (WS >= cur_orig["WS_stall"]) & (WS <= cur_mod["WS_stall"])
ax.fill_between(WS[mask_gain], cur_mod["PW_req"][mask_gain], y_max,
                color="lightgreen", alpha=0.35, zorder=0)

# Additional P/W gap (amber) – only where both designs share W/S
mask_common = WS <= cur_orig["WS_stall"]
ax.fill_between(WS[mask_common],
                cur_orig["PW_req"][mask_common],
                cur_mod ["PW_req"][mask_common],
                color="#ffdf8d", alpha=0.45, zorder=0,
                label="Reduced P/W (mod)")

# ---------------- Design points ------------------------------------
# Baseline
ax.scatter(cur_orig["WS_stall"], cur_orig["PW_stall_int"],
           marker="x", s=100, c="k", lw=2, label="Design pt")
ax.scatter(cur_orig["WS_x"],     cur_orig["PW_x"],
           marker="x", s=100, c="k", lw=2)

# Modified
ax.scatter(cur_mod["WS_stall"],  cur_mod["PW_stall_int"],
           marker="o", s=90, facecolors="none", edgecolors="k", lw=1.8,
           label="Design pt (mod)")
ax.scatter(cur_mod["WS_x"],      cur_mod["PW_x"],
           marker="o", s=90, facecolors="none", edgecolors="k", lw=1.8)

print("Left point", cur_mod["WS_x"], cur_mod["PW_x"])
print("Right point", cur_mod["WS_stall"], cur_mod["PW_stall_int"])

# ---------------- Cosmetics ----------------------------------------
ax.set_xlim(WS.min(), WS.max())
ax.set_ylim(0, y_max)
ax.set_xlabel("Wing loading  W/S  [N m⁻²]")
ax.set_ylabel("Power-to-weight  P/W  [W N⁻¹]")
ax.grid(True, which="both", ls=":")
ax.legend(loc="upper right")
plt.tight_layout()
plt.show()

