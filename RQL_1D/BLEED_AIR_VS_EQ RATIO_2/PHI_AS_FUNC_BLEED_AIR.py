import json, numpy as np, cantera as ct
import matplotlib.pyplot as plt
from mixture_properties import mixture_properties   # your helper

# ---------------- USER SETTINGS -----------------
F_MIN   = 0.00          # 0   % of TOTAL_AIR
F_MAX   = 0.80          # 30  %
N_STEPS = 300            # resolution of the map
MAPFILE = "phi2_air_map.npz"
# -----------------------------------------------

# ---------- 1. Load data from rich run ---------
with open("rich_to_lean.json") as fp:
    J = json.load(fp)

MECH       = J["MECH"]
P_mix      = J["P_mix"]
X_exh      = J["X_exh"]
mdot_exh   = J["mdot_exh"]
h_exh      = J["h_exh"]
h_air2     = J["h_air2"]
TOTAL_AIR  = J["TOTAL_AIR"]
mdot_air1  = J["mdot_air1"]

X_AIR = "N2:0.78084, O2:0.20946"

# ---------- 2. Sweep the percentage ------------
fractions   = np.linspace(F_MIN, F_MAX, N_STEPS)      # fraction of TOTAL_AIR
phi2_list   = []
mdot_air2_list = []

for f in fractions:
    mdot_air2 = TOTAL_AIR - mdot_air1 - f*TOTAL_AIR
    if mdot_air2 <= 0.0:
        # store NaN so the array lengths stay consistent
        phi2_list.append(np.nan)
        mdot_air2_list.append(np.nan)
        continue

    X_mix, _, _ = mixture_properties(
        X1=X_exh, X2=X_AIR,
        mdot1=mdot_exh, mdot2=mdot_air2
    )

    h_mix = (h_exh*mdot_exh + h_air2*mdot_air2) / (mdot_exh + mdot_air2)

    gas = ct.Solution(MECH)
    gas.X, gas.HP = X_mix, (h_mix, P_mix)

    phi2_list.append(gas.equivalence_ratio())
    mdot_air2_list.append(mdot_air2)

phi2_arr     = np.array(phi2_list)
mdot_air2_arr = np.array(mdot_air2_list)

# ---------- 3. Save the map --------------------
np.savez(MAPFILE,
         frac=fractions,
         phi2=phi2_arr,
         mdot_air2=mdot_air2_arr)
print(f"Map written to {MAPFILE}")

# ---------- 4. Quick visual --------------------
fig, ax1 = plt.subplots()
ax1.plot(fractions*100, phi2_arr, label='$\\varphi_2$',color="blue")

# --- added visual cues ----------------------------------------------------- #
ax1.axvline(5.25, color='red')                                  # <<< added
ax1.axvspan(0, 5.25, color='red', alpha=0.15)                   # <<< added
ax1.axhline(0.8, color='red')                                   # <<< added
ax1.axhspan(0.8,1.2 , color='red', alpha=0.15)   # <<< added
# -------------------------------------------------------------------------- #

# -------- axis limits (insert here) -----------------
ax1.set_xlim(0, 80)          # x-axis: 0 – 30 % bleed-air
ax1.set_ylim(0.0, 1)       # left y-axis: φ₂ from 0 to 1.2

ax1.set_xlabel("Ammount of air used for bleed/cooling  [% of total air]")
ax1.set_ylabel("Global equivalence ratio in Lean CC $\\varphi_2$",color="blue")
ax1.grid(ls=":")

ax2 = ax1.twinx()
ax2.plot(fractions*100, mdot_air2_arr, color="C2",
         label='$\\dot m_{air2}$')
ax2.set_ylabel("$\\dot m_{air2}$: Mass entering the Lean CC  [kg s$^{-1}$]",color="C2")
fig.tight_layout()
plt.show()
