import cantera as ct
from mixture_properties import mixture_properties
from typing import Union, Mapping
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np


"""
This script is an H2–air burner‐flame. 
 A first, artificially pre-heated solution is
converged and then used as the initial guess for the 330 °C inlet. 
"""

################################################################################
# 0 – CONSTANTS & SUPPORT FUNCTIONS                                            #
################################################################################

MECH = "CapursoMechanism.yaml"   # renamed for clarity but still the same file
LHV_KEROSENE = 43.15e6   # J kg⁻¹
LHV_H2       = 120.0e6   # J kg⁻¹
ATM2PA       = 101_325.0
C2K          = 273.15     # °C → K

TOTAL_Air = 6

def kerosene_to_h2(mdot_kerosene: float,
                   *,
                   lhv_kerosene: float = LHV_KEROSENE,
                   lhv_h2: float = LHV_H2) -> float:
    """
        This function inputs the amount of kerosene burnt and outputs the equivalent ammount
        of hydrogen needed for combustion using LHV.
    """
    if mdot_kerosene < 0:
        raise ValueError("Mass-flow must be non-negative.")
    return mdot_kerosene * lhv_kerosene / lhv_h2


def _as_mole_str(x: Union[str, Mapping[str, float]]) -> str:
    """
    Returns a mapping as a format of Cantera readable string library.
    """
    if isinstance(x, str):
        return x
    return ", ".join(f"{k}:{v}" for k, v in x.items())


def air_mass_flow_for_phi(
    mdot_H2: float,
    phi: float,
    *,
    X_air: Union[str, Mapping[str, float]],
    mech: str = MECH,
    T: float = 298.15,
    P: float = 101_325.0,
) -> float:
    """
    Given an equivalence ratio, and a mass flow rate of hydrogen, (and a dictionary of the mole fraction of air desired)
    this function delivers the ammount of mass flow of air necesary.
    """
    if mdot_H2 < 0:
        raise ValueError("mdot_H2 must be non-negative.")
    if phi <= 0:
        raise ValueError("phi must be positive.")

    gas = ct.Solution(mech)
    i_H2 = gas.species_index("H2")
    MW_H2 = gas.molecular_weights[i_H2] / 1000.0  # kg mol⁻¹

    n_H2 = mdot_H2 / MW_H2           # mol s⁻¹ H2
    n_O2_stoich = 0.5 * n_H2         # mol s⁻¹ O2

    MW_O2 = gas.molecular_weights[gas.species_index("O2")] / 1000.0
    mdot_O2_req = n_O2_stoich * MW_O2

    gas.TPX = T, P, _as_mole_str(X_air)
    Y_O2_air = gas.mass_fraction_dict()["O2"]

    mdot_air_stoich = mdot_O2_req / Y_O2_air
    return mdot_air_stoich / phi






# -----------------------------------------------------------------------------
# 1 – INLET CONDITIONS  (raise these two to trigger ignition)
# -----------------------------------------------------------------------------
P_0_atm    = 12.1           # atm   ← increase if desired
T_0_C      = 330            # °C    ← increase if desired
EQUIV_R    = 2.5            # φ


## -----------------------------------------------------------------------------
# 2 – COMPOSITION & MASS FLOWS
# -----------------------------------------------------------------------------
X_air = "N2:0.78084, O2:0.20946"
X_H2  = "H2:1"

mdot_kerosene = 0.0818986224               # kg s⁻¹
mdot_H2  = kerosene_to_h2(mdot_kerosene)
mdot_air = air_mass_flow_for_phi(mdot_H2, EQUIV_R, X_air=X_air)

X_mix, mdot_mix, _ = mixture_properties(X1=X_H2, X2=X_air,
                                        mdot1=mdot_H2, mdot2=mdot_air)

# -----------------------------------------------------------------------------
# # 3 – GRID
# # -----------------------------------------------------------------------------
width        = 0.01                        # m
initial_grid = np.linspace(0.0, width, 800)

# -----------------------------------------------------------------------------
# 4 – BUILD & SOLVE THE FLAME
# -----------------------------------------------------------------------------
gas = ct.Solution(MECH)
gas.TPX = T_0_C + C2K, P_0_atm * ATM2PA, X_mix

gas.transport_model = "multicomponent"
print(gas())


flame = ct.BurnerFlame(gas, grid=initial_grid)
flame.burner.mdot = mdot_mix


# enable reacting solution
flame.energy_enabled = True
flame.transport_model = "multicomponent"
flame.solve(loglevel=1, auto=False)

print(gas())
# -----------------------------------------------------------------------------
# 5 – SAVE & PLOT
# -----------------------------------------------------------------------------
flame.save("burner_flame_clean.csv",
           basis="mole",       # or "mass" or "siu"
           overwrite=True)


plt.figure()
plt.plot(flame.grid, flame.T)
plt.xlabel("Distance from burner [m]")
plt.ylabel("Temperature [K]")
plt.title("Temperature profile")
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot species mole fractions
plt.figure()
for species in ['H2', 'O2', 'H2O', 'OH']:
    plt.plot(f.grid, data(species).X, label=f'{species} Mole Fraction')
plt.xlabel('Distance (m)')
plt.ylabel('Mole Fraction')
plt.title('Species Mole Fractions')
plt.legend()
plt.grid()

# Show the plots
fig, ax = plt.subplots()
major = ('O2', 'N2O', 'NO','H2','H2O')
states = flame.to_array()
ax.plot(states.grid * 1000, states(*major).X, label=major)
ax.set(xlabel='distance from burner [mm]', ylabel='mole fractions')

ax.legend()
plt.show()




fig, ax1 = plt.subplots()

ax1.plot(flame.grid * 1000, flame.heat_release_rate / 1e6, color='C4')
ax1.set_ylabel('heat release rate [MW/m³]', color='C4')

ax1.set(xlabel='distance from burner [mm]')

ax2 = ax1.twinx()
ax2.plot(flame.grid * 1000, flame.T, color='C3')
ax2.set_ylabel('temperature [K]', color='C3')
plt.show()






###############################################################################
# 6 – MIX EXHAUST  +  QUENCH AIR  +  WATER   (with mixture_properties)
###############################################################################
# ===== USER INPUTS for the quench stage ======================================
mdot_air2 = TOTAL_Air - mdot_air   # kg/s  (remaining or quench air)
T_air2    = 400   # K     (temperature of that air)
mdot_H2O  = 0.5   # kg/s  (water vapour; if liquid, pre-convert to vapour mass)
T_H2O     = 400  # K     (≥ 373 K so it is gas phase)
P_mix     = flame.P        # keep the same pressure for the mixer
# ============================================================================

# ---------------------------------------------------------------------------
# 6.1  Build Solution objects for each feed stream
# ---------------------------------------------------------------------------
# --- exhaust from the rich burner-flame -------------------------------------
gas_exh = ct.Solution(MECH)
gas_exh.TPX = flame.T[-1], P_mix, flame.X[:, -1]   # exit state
print(gas_exh())

# --- quench air stream ------------------------------------------------------
gas_air2 = ct.Solution(MECH)
print(T_air2,"....", P_mix,"....", X_air)
gas_air2.TPX = T_air2, P_mix, X_air
quit()
# --- water vapour stream ----------------------------------------------------
gas_H2O = ct.Solution(MECH)
gas_H2O.TPX = T_H2O, P_mix, "H2O:1"

# ---------------------------------------------------------------------------
# 6.2  Use mixture_properties  (two streams at a time)
# ---------------------------------------------------------------------------
# First mix: (exhaust + quench-air)  →  X_mix1 , mdot_mix1
X_exh = _as_mole_str(
    {sp: x for sp, x in zip(gas_exh.species_names, gas_exh.X) if x > 0.0}
)
X_mix1, mdot_mix1, _ = mixture_properties(
    X1=X_exh,
    X2=X_air,
    mdot1=mdot_mix,      # from the rich flame
    mdot2=mdot_air2
)

# Second mix: (previous mixture + water)  →  X_mix_final
X_mix_final, mdot_mix_final, _ = mixture_properties(
    X1=X_mix1,
    X2="H2O:1",
    mdot1=mdot_mix1,
    mdot2=mdot_H2O
)

# ---------------------------------------------------------------------------
# 6.3  Compute adiabatic-mix temperature (mass-weighted enthalpy)
# ---------------------------------------------------------------------------
h_mix = (
    mdot_mix      * gas_exh.enthalpy_mass +
    mdot_air2     * gas_air2.enthalpy_mass +
    mdot_H2O      * gas_H2O.enthalpy_mass
) / mdot_mix_final

# ---------------------------------------------------------------------------
# 6.4  Create the final mixed gas Solution object
# ---------------------------------------------------------------------------
mixed_gas = ct.Solution(MECH)
mixed_gas.TPX = 500.0, P_mix, X_mix_final      # rough initial T
mixed_gas.HP  = h_mix, P_mix                   # enforce adiabatic mix

print("\n===== MIXER RESULTS =====")
print(f"Mass-flow total      : {mdot_mix_final:10.4f}  kg/s")
print(f"Adiabatic mix T      : {mixed_gas.T:10.2f}  K")
print(f"Lean-stage φ (air only): {mixed_gas.equivalence_ratio():6.3f}")
print("=========================\n")


print(mixed_gas())



