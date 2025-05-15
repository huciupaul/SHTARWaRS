import cantera as ct
ct.suppress_thermo_warnings()
import sys
from mixture_properties import mixture_properties,kerosene_to_h2,air_mass_flow_for_phi,_as_mole_str,MECH
import matplotlib.pyplot as plt


from datetime import datetime
from pathlib import Path
import json

# --------------------------------------------------------------------------- #

# -----------------------------  CONSTANTS  --------------------------------- #
ATM2PA = 101_325.0
C2K    = 273.15

# --------------------------  INPUT PARAMETERS  ----------------------------- #
TOTAL_AIR     = 4.635714 #5.0847705                 # kg/s  – total air available
MDOT_KEROSENE =  0.08202600621432085 #0.0818986224      # # kg/s  – reference kerosene flow
T_IN_C        = 330.0             # °C    – inlet T of rich stage
P_IN_ATM      = 12.1              # atm   – inlet P
PHI_RICH      = 5              # φ₁

T_inlet = T_IN_C + C2K
P_inlet = P_IN_ATM*ATM2PA

# quench
T_AIR2        = T_inlet             # K
MDOT_H2O      = 0.0               # kg/s  – water vapour
T_H2O         = 873.0             # K

# mixtures
X_AIR = "N2:0.78084, O2:0.20946"
X_H2  = "H2:1"

width_rich = 0.001          # m
width_lean = 0.02          # m

PLOT_LEAN = True
PLOT_RICH = True
# ========================================================================== #
# 0 – DERIVED FLOWS & STATE
# ========================================================================== #


mdot_H2   = kerosene_to_h2(MDOT_KEROSENE)

mdot_air1 = air_mass_flow_for_phi(mdot_H2, PHI_RICH, X_air=X_AIR)
mdot_air2 = TOTAL_AIR - mdot_air1          # quench air

if mdot_air2 < 0:
    sys.exit("ERROR: TOTAL_AIR is smaller than air needed for rich φ₁")

print("\n===== INPUT SUMMARY =====")
print(f"Mechanism file        : {MECH}")
print(f"P₀, T₀                : {P_IN_ATM:.2f} atm  {T_IN_C:.1f} °C")
print(f"φ₁ (rich)             : {PHI_RICH:.2f}")
print(f"H₂ mass-flow          : {mdot_H2:8.4f} kg/s")
print(f"Air for Rich mass-flow       : {mdot_air1:8.4f} kg/s")
print(f"Air for Lean mass-flow       : {mdot_air2:8.4f} kg/s")
print(f"Water mass-flow       : {MDOT_H2O:8.4f} kg/s")
print("=========================\n")



# ========================================================================== #
# 1 – MIX H₂ + AIR-1  (rich mixture for BurnerFlame)
# ========================================================================== #
X_mix, mdot_mix, _ = mixture_properties(
    mdot1=mdot_H2, X1=X_H2,
    mdot2=mdot_air1, X2=X_AIR,T=T_inlet,P=P_inlet
)

print("\n===== RICH REACTANTS SUMMARY =====")
print(f"P₀, T₀                :                        {P_inlet:.2f} Pa,  {T_inlet:.1f} K")
print(f"The reacting mixtures Molar Ratio is:          {X_mix}.")
print(f"The reacting mixtures MASS FLOW RATE is:       {mdot_mix:.4f} Kg/s.")
print(f"The reacting mixtures Molar Mass is:       {_:.4f} Kg/mol.")

# ========================================================================== #
# 2 – BUILD & SOLVE RICH BURNER-STABILIZED FLAME
# ========================================================================== #
gas = ct.Solution(MECH)   # emptry gas object

gas.transport_model = "multicomponent"
gas.TPX = T_inlet, P_inlet, X_mix

print(gas())

flame = ct.BurnerFlame(gas, width=width_rich)
flame.burner.mdot = mdot_mix
flame.set_refine_criteria(ratio=6.0, slope=0.005, curve=0.01)


################################################################################################################################################
###################  WATER IS CREATED HERE ???  ##############################
########################################################################

# one solve
print(flame.max_grid_points)
flame.set_max_grid_points(flame.flame,1200)
print("Flame refinements increased to ",flame.max_grid_points," Points")
flame.solve(loglevel=0, auto=True, refine_grid=True)



###################################

# rich-stage diagnostics
idx_NO  = flame.gas.species_index("NO")
idx_N2O = flame.gas.species_index("N2O")

X_NO_out  = flame.X[idx_NO,  -1] if idx_NO  >= 0 else 0.0
X_N2O_out = flame.X[idx_N2O, -1] if idx_N2O >= 0 else 0.0

print("----- RICH FLAME RESULTS -----")
print(f"grid points           : {len(flame.grid)}")
print(f"peak temperature      : {flame.T.max():10.2f}  K")
print(f"exit temperature      : {flame.T[-1]:10.2f}  K")
print(f"exit X_NO  , ppm      : {X_NO_out :10.3e}  ({X_NO_out*1e6:8.3f})")
print(f"exit X_N₂O, ppm       : {X_N2O_out:10.3e}  ({X_N2O_out*1e6:8.3f})")
print("-------------------------------\n")

# optional CSV save
flame.save("rich_flame.csv", basis="mole", overwrite=True)


# ========================================================================== #
# 4 – PLOTS  (optional)
# ========================================================================== #
if PLOT_RICH:
    plt.figure()
    plt.plot(flame.grid, flame.T)
    plt.xlabel("distance from burner [m]")
    plt.ylabel("T [K]")
    plt.title("Rich flame temperature")
    plt.tight_layout()

    major = ("O2",  "H2", "H2O")
    states = flame.to_array()
    plt.figure()
    for sp in major:
        plt.plot(states.grid*1e3, states(sp).X, label=sp)
    plt.xlabel("distance from burner [mm]")
    plt.ylabel("mole fraction")
    plt.legend()
    plt.tight_layout()

    major = ("NO", "N2O")#,"NO2")
    states = flame.to_array()
    plt.figure()
    for sp in major:
        plt.plot(states.grid*1e3, states(sp).X, label=sp)
    plt.xlabel("distance from burner [mm]")
    plt.ylabel("mole fraction")
    plt.legend()
    plt.tight_layout()

    fig, ax1 = plt.subplots()
    ax1.plot(flame.grid*1e3, flame.heat_release_rate/1e6)
    ax1.set_xlabel("distance [mm]")
    ax1.set_ylabel("heat-release [MW/m³]")
    ax2 = ax1.twinx()
    ax2.plot(flame.grid*1e3, flame.T, color="C3")
    ax2.set_ylabel("T [K]")
    plt.tight_layout()
    plt.show()


# ========================================================================== #
# 3 – QUENCH: MIX EXHAUST + REMAINING AIR + WATER
# ========================================================================== #


P_mix = flame.P

gas_exh = ct.Solution(MECH)
gas_exh.TPX = flame.T[-1], P_mix, flame.X[:, -1]

gas_air2 = ct.Solution(MECH)
gas_air2.TPX = T_AIR2, P_mix, X_AIR

gas_H2O = ct.Solution(MECH)
gas_H2O.TPX = T_H2O, P_mix, "H2O:1"


# two-step mix via helper
X_exh = _as_mole_str({sp: x for sp, x in zip(gas_exh.species_names, gas_exh.X) if x > 0})


X_mix1, mdot_mix1, _ = mixture_properties(X1=X_exh, X2=X_AIR,
                                          mdot1=mdot_mix, mdot2=mdot_air2)

X_mix_final, mdot_mix_final, _ = mixture_properties(X1=X_mix1, X2="H2O:1",
                                                    mdot1=mdot_mix1, mdot2=MDOT_H2O)

# adiabatic mixed temperature
h_mix = (
    mdot_mix      * gas_exh.enthalpy_mass +
    mdot_air2     * gas_air2.enthalpy_mass +
    MDOT_H2O      * gas_H2O.enthalpy_mass
) / mdot_mix_final

mixed_gas = ct.Solution(MECH)
mixed_gas.X =  X_mix_final
mixed_gas.HP  = h_mix, P_mix

print("----- MIXER RESULTS -----")
print(f"air-2 mass-flow       : {mdot_air2:10.4f} kg/s")
print(f"air-2 Temperature      : {T_AIR2:10.4f} kg/s")
print(f"Exhaust mass-flow       : {mdot_mix:10.4f} kg/s")
print(f"Exhaust temperature       : {gas_exh.T:10.4f} kg/s")
print(f"water mass-flow       : {MDOT_H2O:10.4f} kg/s")
print(f"water temperature       : {T_H2O:10.4f} kg/s")
print(f"mixed gas temperature : {mixed_gas.T:10.2f}  K")
print(f"lean φ₂     : {mixed_gas.equivalence_ratio():6.3f}")
print("-------------------------\n")
print(mixed_gas())





gas_lean = ct.Solution(MECH)
gas_lean.transport_model = "multicomponent"
gas_lean.TPX = mixed_gas.T, P_mix, mixed_gas.X

lean_flame = ct.BurnerFlame(gas_lean, width=width_lean)
lean_flame.burner.mdot = mdot_mix_final          # kg s⁻¹
lean_flame.set_refine_criteria(ratio=6.0, slope=0.01, curve=0.03)

# first solve (transport = mixture-averaged is fine here)
lean_flame.solve(loglevel=0, auto=True, refine_grid=True)

# ---------------------------------------------------------------------------
# 7.1  Post-process results
# ---------------------------------------------------------------------------
idx_NO  = gas_lean.species_index("NO")
idx_N2O = gas_lean.species_index("N2O")

X_NO_lean  = lean_flame.X[idx_NO,  -1] if idx_NO  >= 0 else 0.0
X_N2O_lean = lean_flame.X[idx_N2O, -1] if idx_N2O >= 0 else 0.0

print("===== LEAN FLAME RESULTS =====")
print(f"grid points           : {len(lean_flame.grid)}")
print(f"peak temperature      : {lean_flame.T.max():10.2f}  K")
print(f"exit temperature      : {lean_flame.T[-1]:10.2f}  K")
print(f"lean φ₂ (air only)    : {mixed_gas.equivalence_ratio():6.3f}")
print(f"exit X_NO  , ppm      : {X_NO_lean :10.3e}  ({X_NO_lean*1e6:8.2f})")
print(f"exit X_N₂O, ppm       : {X_N2O_lean:10.3e}  ({X_N2O_lean*1e6:8.2f})")
print("==============================\n")

# optional: save lean profile
lean_flame.save("lean_flame.csv", basis="mole", overwrite=True)



if PLOT_LEAN:
    # 4.1 Temperature profile
    plt.figure()
    plt.plot(lean_flame.grid, lean_flame.T, lw=2)
    plt.xlabel("distance from burner [m]")
    plt.ylabel("T [K]")
    plt.title("Lean flame temperature")
    plt.tight_layout()

    # 4.2 Major species
    major = ("O2", "H2", "H2O")
    states = lean_flame.to_array()
    plt.figure()
    for sp in major:
        plt.plot(states.grid*1e3, states(sp).X, label=sp)
    plt.xlabel("distance from burner [mm]")
    plt.ylabel("mole fraction")
    plt.legend()
    plt.tight_layout()


    minor = ("NO", "N2O")
    plt.figure()
    for sp in minor:
        plt.plot(states.grid*1e3, states(sp).X, label=sp)
    plt.xlabel("distance from burner [mm]")
    plt.ylabel("mole fraction")
    plt.legend()
    plt.tight_layout()

    # 4.4 Local heat-release overlayed with temperature
    fig, ax1 = plt.subplots()
    ax1.plot(lean_flame.grid*1e3,
             lean_flame.heat_release_rate/1e6,
             label="q̇", lw=2)
    ax1.set_xlabel("distance [mm]")
    ax1.set_ylabel("heat-release [MW m⁻³]")
    ax2 = ax1.twinx()
    ax2.plot(lean_flame.grid*1e3, lean_flame.T,
             label="T", lw=2, color="C3")
    ax2.set_ylabel("T [K]")
    plt.title("Lean flame: heat-release & temperature")
    plt.tight_layout()
    plt.show()


# do a final JSON storing the results
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
        "water_kg_s": MDOT_H2O,
    },
    "rich_flame": {
        "grid_points": len(flame.grid),
        "T_peak_K": float(flame.T.max()),
        "T_exit_K": float(flame.T[-1]),
        "X_NO_exit": float(X_NO_out),
        "X_N2O_exit": float(X_N2O_out),
    },
    "mixer": {
        "phi_lean_theoretical": float(mixed_gas.equivalence_ratio()),
        "T_mixed_K": float(mixed_gas.T),
    },
    "lean_flame": {
        "grid_points": len(lean_flame.grid),
        "T_peak_K": float(lean_flame.T.max()),
        "T_exit_K": float(lean_flame.T[-1]),
        "phi_lean_actual": float(mixed_gas.equivalence_ratio()),
        "X_NO_exit": float(X_NO_lean),
        "X_N2O_exit": float(X_N2O_lean),
    },
    "timestamp": datetime.utcnow().isoformat(timespec="seconds") + "Z"
}

summary_file = Path("../run_summary.json")
with summary_file.open("w") as fp:
    json.dump(summary, fp, indent=2)
print(f"\n►► Summary written to {summary_file.resolve()}\n")