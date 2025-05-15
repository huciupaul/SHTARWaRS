import cantera as ct
ct.suppress_thermo_warnings()

import sys, json
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

# ---- helper functions -----------------------------------------------------
from mixture_properties import (
    mixture_properties,
    kerosene_to_h2,
    air_mass_flow_for_phi,
    _as_mole_str,
    MECH,
)

# ---------------------------------------------------------------------------
# 0 – CONSTANTS & INPUTS
# ---------------------------------------------------------------------------
ATM2PA = 101_325.0
C2K    = 273.15

TOTAL_AIR     = 4.635545784891179         # kg/s
MDOT_KEROSENE = 0.08202600621432085
T_IN_C        = 330.0            # °C
P_IN_ATM      = 12.1             # atm
PHI_RICH      = 5                # φ₁

T_inlet = T_IN_C + C2K
P_inlet = P_IN_ATM * ATM2PA

# quench & lean-stage params
T_AIR2   = T_inlet
mdot_H2O = 0.0
T_H2O    = 873.0

X_AIR = "N2:0.78084, O2:0.20946"
X_H2  = "H2:1"

width_rich = 0.001   # m – domain length for free rich flame
width_lean = 0.0001    # m – domain length for free lean flame

PLOT_RICH = True
PLOT_LEAN = True

# ---------------------------------------------------------------------------
# 1 – DERIVED FLOWS & STATE
# ---------------------------------------------------------------------------
mdot_H2   = kerosene_to_h2(MDOT_KEROSENE)
mdot_air1 = air_mass_flow_for_phi(mdot_H2, PHI_RICH, X_air=X_AIR)
mdot_air2 = TOTAL_AIR - mdot_air1

if mdot_air2 < 0:
    sys.exit("ERROR: TOTAL_AIR is smaller than air needed for rich φ₁")

print("\n===== INPUT SUMMARY =====")
print(f"Mechanism file        : {MECH}")
print(f"P₀, T₀                : {P_IN_ATM:.2f} atm  {T_IN_C:.1f} °C")
print(f"φ₁ (rich)             : {PHI_RICH:.2f}")
print(f"H₂ mass-flow          : {mdot_H2:8.4f} kg/s")
print(f"Air for Rich stage    : {mdot_air1:8.4f} kg/s")
print(f"Air for Lean stage    : {mdot_air2:8.4f} kg/s")
print(f"Water mass-flow       : {mdot_H2O:8.4f} kg/s")
print("===========================================================\n")

# ---------------------------------------------------------------------------
# 2 – RICH PREMIXED MIXTURE & FREE FLAME
# ---------------------------------------------------------------------------
X_mix, mdot_rich_burn, _ = mixture_properties(
    mdot1=mdot_H2, X1=X_H2,
    mdot2=mdot_air1, X2=X_AIR,
    T=T_inlet, P=P_inlet
)

print("\n===== RICH REACTANTS SUMMARY =====")
print(f"P₀, T₀                : {P_inlet:.0f} Pa  {T_inlet:.1f} K")
print(f"Molar composition     : {X_mix}")
print("===========================================================\n")

gas = ct.Solution(MECH)
gas.transport_model = "multicomponent"
gas.TPX = T_inlet, P_inlet, X_mix

flame_rich = ct.FreeFlame(gas, width=width_rich)
flame_rich.set_refine_criteria(ratio=6.0, slope=0.01, curve=0.02)

flame_rich.solve(loglevel=0, auto=True, refine_grid=True)

# ---- diagnostics ----------------------------------------------------------
idx_NO  = gas.species_index("NO")
#idx_NO2  = gas.species_index("NO2")
idx_N2O = gas.species_index("N2O")
X_NO_out  = flame_rich.X[idx_NO,  -1] if idx_NO  >= 0 else 0.0
#X_NO2_out  = flame_rich.X[idx_NO2,  -1] if idx_NO2  >= 0 else 0.0
X_N2O_out = flame_rich.X[idx_N2O, -1] if idx_N2O >= 0 else 0.0

print("----- RICH FREE-FLAME RESULTS -----")
print(f"grid points           : {len(flame_rich.grid)}")
print(f"peak temperature      : {flame_rich.T.max():10.2f} K")
print(f"burnt-gas temperature : {flame_rich.T[-1]:10.2f} K")
print(f"X_NO  (burnt)         : {X_NO_out*1e6:10.3f} ppm")
print(f"X_N₂O (burnt)         : {X_N2O_out*1e6:10.3f} ppm")
#print(f"X_NO₂  (burnt)         : {X_NO2_out*1e6:10.3f} ppm")
print("-----------------------------------------------------------\n")

flame_rich.save("rich_flame_free.csv", basis="mole", overwrite=True)


# optional plots
if PLOT_RICH:
    plt.figure()
    plt.plot(flame_rich.grid, flame_rich.T)
    plt.xlabel("distance [m]"); plt.ylabel("T [K]")
    plt.title("Rich free flame temperature")
    plt.tight_layout()

    major = ("O2", "H2", "H2O")
    states = flame_rich.to_array()
    plt.figure()
    for sp in major:
        plt.plot(states.grid*1e3, states(sp).X, label=sp)
    plt.xlabel("distance [mm]"); plt.ylabel("mole fraction")
    plt.legend(); plt.tight_layout()
    plt.show()

    major = ("N2O", "NO"   )#   , "NO2")
    states = flame_rich.to_array()
    plt.figure()
    for sp in major:
        plt.plot(states.grid*1e3, states(sp).X, label=sp)
    plt.xlabel("distance [mm]"); plt.ylabel("mole fraction")
    plt.legend(); plt.tight_layout()
    plt.show()

# ---------------------------------------------------------------------------
# 3 – QUENCH: MIX EXHAUST + REMAINING AIR + WATER
# ---------------------------------------------------------------------------
P_mix = flame_rich.P

#Exhaust gases
gas_exh = ct.Solution(MECH)
gas_exh.TPX = flame_rich.T[-1], P_mix, flame_rich.X[:, -1]
mdot_gas_exh = mdot_rich_burn

#new air for quenching
gas_air2 = ct.Solution(MECH)
gas_air2.TPX = T_AIR2, P_mix, X_AIR
# mdot_air2 is the mass flow

#Water flow
gas_H2O = ct.Solution(MECH)
gas_H2O.TPX = T_H2O, P_mix, "H2O:1"
# mdot_H2O is the mass flow of water

X_exh = _as_mole_str({sp: x for sp, x in zip(gas_exh.species_names, gas_exh.X) if x > 0})

######## MIXTURE PART

#MIx exhaust with new air
X_mix1, mdot_mix1, _ = mixture_properties(
    X1=X_exh, X2=X_AIR,
    mdot1=mdot_gas_exh, mdot2=mdot_air2
)

# PART1: THE MOLAR FRACTIONS
#Mix the prev. mixture with water
X_mix_final, _, _ = mixture_properties(
    X1=X_mix1, X2="H2O:1",
    mdot1=mdot_mix1, mdot2=mdot_H2O
)


# PART2:  CALCULATE ENTHALPY OF MIXTURE
h_exh  = gas_exh.enthalpy_mass
h_air2 = gas_air2.enthalpy_mass
h_H2O  = gas_H2O.enthalpy_mass


mdot_total = mdot_gas_exh + mdot_air2 + mdot_H2O

h_mix = (
    h_exh*mdot_gas_exh +
    mdot_air2 * h_air2 +
    mdot_H2O  * h_H2O) / mdot_total

#initiate new gas
mixed_gas = ct.Solution(MECH)
mixed_gas.X = X_mix_final
mixed_gas.HP = h_mix, P_mix

export = dict(
    MECH      = MECH,
    P_mix     = P_mix,
    X_exh     = X_exh,
    mdot_exh  = mdot_rich_burn,          # or mdot_gas_exh
    h_exh     = h_exh,
    h_air2    = h_air2,
    TOTAL_AIR = TOTAL_AIR,
    mdot_air1 = mdot_air1,
)
with open("BLEED_AIR_VS_EQ RATIO_2/rich_to_lean.json", "w") as fp:
    json.dump(export, fp)
quit()

print("----- MIXER RESULTS -----")
print(f"mixed-gas temperature : {mixed_gas.T:10.2f} K")
print(f"lean φ₂               : {mixed_gas.equivalence_ratio():6.3f}")
print("-----------------------------------------------------------\n")


# ---------------------------------------------------------------------------
# 4 – LEAN FREE FLAME
# ---------------------------------------------------------------------------
gas_lean = ct.Solution(MECH)
gas_lean.transport_model = "mixture-averaged"
gas_lean.TPX = mixed_gas.T, P_mix, mixed_gas.X



# Set up flame object
flame_lean = ct.FreeFlame(gas_lean, width=width_lean)
flame_lean.set_refine_criteria(ratio=5.0, slope=0.001, curve=0.002)

# Solve with mixture-averaged transport model
flame_lean.transport_model = "mixture-averaged"
flame_lean.set_max_grid_points(flame_lean.flame,10000)
flame_lean.solve(loglevel=0, auto=True, refine_grid=True)
flame_lean.show()
print(f"mixture-averaged flamespeed = {flame_lean.velocity[0]:7f} m/s")


# Solve with multi-component transport properties
gas_lean.transport_model = "multicomponent"
flame_lean.transport_model = 'multicomponent'
flame_lean.set_max_grid_points(flame_lean.flame,10000)
flame_lean.set_refine_criteria(ratio=5.0, slope=0.01, curve=0.02)
flame_lean.solve(1)  # don't use 'auto' on subsequent solves
flame_lean.show()
print(f"multicomponent flamespeed = {flame_lean.velocity[0]:7f} m/s")
flame_lean.show()

idx_NO  = gas_lean.species_index("NO")
idx_N2O = gas_lean.species_index("N2O")
X_NO_lean  = flame_lean.X[idx_NO,  -1] if idx_NO  >= 0 else 0.0
X_N2O_lean = flame_lean.X[idx_N2O, -1] if idx_N2O >= 0 else 0.0

print("===== LEAN FREE-FLAME RESULTS =====")
print(f"grid points           : {len(flame_lean.grid)}")
print(f"peak temperature      : {flame_lean.T.max():10.2f} K")
print(f"burnt-gas temperature : {flame_lean.T[-1]:10.2f} K")
print(f"X_NO  (burnt)         : {X_NO_lean*1e6:10.2f} ppm")
print(f"X_N₂O (burnt)         : {X_N2O_lean*1e6:10.2f} ppm")
print("===========================================================\n")

flame = flame_lean.to_array()
flame.save("lean_flame_free.csv", basis="mole", overwrite=True)

# optional lean plots
if PLOT_LEAN:
    plt.figure()
    plt.plot(flame_lean.grid, flame_lean.T)
    plt.xlabel("distance [m]"); plt.ylabel("T [K]")
    plt.title("Lean free flame temperature")
    plt.tight_layout()

    major = ("O2", "H2", "H2O")
    states = flame_lean.to_array()
    plt.figure()
    for sp in major:
        plt.plot(states.grid*1e3, states(sp).X, label=sp)
    plt.xlabel("distance [mm]"); plt.ylabel("mole fraction")
    plt.legend(); plt.tight_layout()
    plt.show()

    states = flame_lean.to_array()
    x_mm = states.grid * 1e3

    # 1. Temperature & heat-release
    fig, ax1 = plt.subplots()
    ax1.plot(x_mm, flame_lean.heat_release_rate / 1e6, lw=1.5)
    ax1.set_xlabel("distance [mm]")
    ax1.set_ylabel("q̇ [$\mathrm{MW\,m^{-3}}$]")
    ax2 = ax1.twinx()
    ax2.plot(x_mm, flame_lean.T, "C3", lw=1.5)
    ax2.set_ylabel("T [K]")
    plt.title("Lean flame: heat-release & temperature")
    plt.tight_layout()



    # 5. NO & N2O
    plt.figure()
    for sp in ("NO", "N2O"):
        plt.semilogy(x_mm, np.maximum(states(sp).X, 1e-12), label=sp)
    plt.xlabel("distance [mm]");
    plt.ylabel("mole fraction (log)")
    plt.legend();
    plt.title("NOₓ – lean free flame");
    plt.tight_layout()

    plt.show()

# ---------------------------------------------------------------------------
# 5 – JSON summary
# ---------------------------------------------------------------------------
summary = {
    "input": {
        "mechanism": MECH,
        "P0_atm": P_IN_ATM,
        "T0_C": T_IN_C,
        "phi_rich": PHI_RICH,
        "total_air_kg_s": TOTAL_AIR,
        "mdot_kerosene_kg_s": MDOT_KEROSENE,
        "mdot_H2_kg_s": mdot_H2,
    },
    "flows": {
        "air_rich_kg_s": mdot_air1,
        "air_lean_kg_s": mdot_air2,
        "water_kg_s": mdot_H2O,
    },
    "rich_flame": {
        "grid_points": len(flame_rich.grid),
        "T_peak_K": float(flame_rich.T.max()),
        "T_exit_K": float(flame_rich.T[-1]),
        "X_NO_exit": float(X_NO_out),
        "X_N2O_exit": float(X_N2O_out),
    },
    "mixer": {
        "phi_lean_theoretical": float(mixed_gas.equivalence_ratio()),
        "T_mixed_K": float(mixed_gas.T),
    },
    "lean_flame": {
        "grid_points": len(flame_lean.grid),
        "T_peak_K": float(flame_lean.T.max()),
        "T_exit_K": float(flame_lean.T[-1]),
        "phi_lean_actual": float(mixed_gas.equivalence_ratio()),
        "X_NO_exit": float(X_NO_lean),
        "X_N2O_exit": float(X_N2O_lean),
    },
    "timestamp": datetime.utcnow().isoformat(timespec="seconds") + "Z"
}

summary_file = Path("run_summary.json")
with summary_file.open("w") as fp:
    json.dump(summary, fp, indent=2)
print(f"►► Summary written to {summary_file.resolve()}")