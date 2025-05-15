
"""
Parametric peak-temperature map for the rich stage
==================================================
Sweeps:
    φ  : 1 → 5
    Tin: 300 K → 900 K
Outputs:
    3-D matplotlib surface
"""

import cantera as ct
ct.suppress_thermo_warnings()

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
import h5py, json
from tqdm import tqdm

from mixture_properties import (
    mixture_properties,
    kerosene_to_h2,
    air_mass_flow_for_phi,
    MECH,
)

P_IN_ATM      = 12.1
ATM2PA        = 101_325.0
P_INLET       = P_IN_ATM * ATM2PA

TOTAL_AIR     = 4.635714          # kg/s
MDOT_KEROSENE = 0.08202600621432085
MDOT_H2       = kerosene_to_h2(MDOT_KEROSENE)

X_AIR = "N2:0.78084, O2:0.20946"
X_H2  = "H2:1"

WIDTH_RICH = 0.0004   # m

phi_vec   = np.linspace(1, 5, 100)
Tin_vec   = np.linspace(300.0, 900.0, 100) # Kelvin

T_peak    = np.empty((len(Tin_vec), len(phi_vec)))


print(f"Running {len(phi_vec)*len(Tin_vec)} flame solves …")

for iT, Tin in enumerate(tqdm(Tin_vec, desc="Tin grid")):
    for jφ, phi in enumerate(phi_vec):

        # 2.1 derive air flow for this φ
        mdot_air1 = air_mass_flow_for_phi(MDOT_H2, phi, X_air=X_AIR)
        if mdot_air1 > TOTAL_AIR:
            T_peak[iT, jφ] = np.nan
            continue


        X_mix, mdot_mix, _ = mixture_properties(
            mdot1=MDOT_H2, X1=X_H2,
            mdot2=mdot_air1, X2=X_AIR,
            T=Tin, P=P_INLET
        )


        g = ct.Solution(MECH)
        g.transport_model = "multicomponent"
        g.TPX = Tin, P_INLET, X_mix
        flame = ct.BurnerFlame(g, width=WIDTH_RICH)
        flame.burner.mdot = mdot_mix
        flame.set_refine_criteria(ratio=6.0, slope=0.01, curve=0.02)
        flame.solve(loglevel=0, auto=True, refine_grid=True)

        T_peak[iT, jφ] = float(flame.T.max())


Φ, TIN = np.meshgrid(phi_vec, Tin_vec)   # grid for plotting

fig = plt.figure(figsize=(7, 5))
ax  = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(
    Φ, TIN, T_peak,
    linewidth=0, antialiased=True, cmap="viridis"
)
ax.set_xlabel("Equivalence ratio  $\\varphi$")
ax.set_ylabel("Inlet $T_{mix}$  [K]")
ax.set_zlabel("Peak flame $T_{max}$  [K]")
ax.set_title("Rich-stage peak temperature map")
fig.colorbar(surf, shrink=0.6, label="$T_{max}$ [K]")
plt.tight_layout()
plt.show()

OUTDIR = Path("runs")
OUTDIR.mkdir(exist_ok=True)
outfile = OUTDIR / "rich_Tmax_map.h5"

with h5py.File(outfile, "w") as h5f:
    h5f.create_dataset("phi", data=phi_vec)
    h5f.create_dataset("Tin", data=Tin_vec)
    h5f.create_dataset("Tpeak", data=T_peak)
    meta = {
        "mechanism": MECH,
        "P_atm": P_IN_ATM,
        "fuel_mdot": MDOT_H2,
        "total_air": TOTAL_AIR,
        "timestamp": ct.now()
    }
    h5f.attrs["meta"] = json.dumps(meta, indent=2)

print(f"\n►► Data & metadata written to {outfile.resolve()}\n")

import numpy as np
np.savez("runs/rich_Tmax_map.npz", phi=phi_vec, Tin=Tin_vec, Tpeak=T_peak)
