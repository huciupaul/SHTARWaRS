import cantera as ct
from mixture_properties import mixture_properties
from typing import Union, Mapping
import matplotlib.pyplot as plt
from pathlib import Path


"""
This script models an RQL combustor with H2 and air.
"""

################################################################################
# 0 – CONSTANTS & SUPPORT FUNCTIONS                                            #
################################################################################

#Mechanism
MECH = "CapursoMechanism.yaml"   # renamed for clarity but still the same file

# Constants
ATM2PA       = 101_325.0
C2K          = 273.15     # °C → K



TOTAL_Air = 6








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


<<<<<<< HEAD
<<<<<<< HEAD






"""
This script simulates a burner flame using Cantera. It sets up the initial conditions, 
calculates the required mass flow rates for hydrogen and air, and then solves the burner flame problem. 
The results are saved to a file and plotted.

0: At CC inlet. Inlet conditions
1: At Rich Combustion exhaust
2: At mixer inlet
3: At mixer outlet.
4: At Lean Burn inlet
5: At Lean Burn Outlet
"""

############################ Inlet conditions ###################################################

X_air = 'N2:0.78084,  O2:0.20946,  AR:0.00934,  CO2:0.000407 '
X_H2 = 'H2:1'

mdot_Kerosene = 0.0818986224 #Kg/s at take-off


# Quantities @ CC entrance
P_0 = 12.1 #Atmospheres
T_0 = 330 # CELSIUS
TOTAL_mdot_air = 12



P_0_pascals = P_0 * 101325
T_0_kelvin = T_0 + 273.15


RICH_EQUIVALENCE_RATIO = 2.5





mdot_H2 = kerosene_to_h2(mdot_Kerosene)
print(mdot_H2)
mdot_air = air_mass_flow_for_phi(mdot_H2=mdot_H2,phi=RICH_EQUIVALENCE_RATIO,X_air=X_air)

print(mdot_air)

(X_mix,mdot_mixture,M_mix) = mixture_properties(X1=X_H2,X2=X_air,mdot1=mdot_H2,mdot2=mdot_air)





gas = ct.Solution('SanDiego.yaml')
gas.TPX = T_0, P_0_pascals, X_mix


width = 0.2  # m

f = ct.BurnerFlame(gas, width=width)
f.burner.mdot = mdot_mixture
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.05)
f.show()

loglevel = 1

################# First Solve ######################
f.transport_model = 'mixture-averaged'
f.solve(loglevel, auto=True)

if "native" in ct.hdf_support():
    output = Path() / "burner_flame.h5"
else:
    output = Path() / "burner_flame.yaml"
output.unlink(missing_ok=True)

f.save(output, name="mix", description="solution with mixture-averaged transport")



################# Subsequent Solves #########################
f.transport_model = 'multicomponent'
f.solve(loglevel)  # don't use 'auto' on subsequent solves
f.show()
f.save(output, name="multi", description="solution with multicomponent transport")

f.save('burner_flame.csv', basis="mole", overwrite=True)

# Load the saved data
data = ct.SolutionArray(gas)
data.read_csv('burner_flame.csv')

# Plot temperature profile
=======
=======
>>>>>>> af1e39737ec0a5ccd90e2399aff2223894492a21
gas.transport_model = "multicomponent"
print(gas())


<<<<<<< HEAD
flame = ct.BurnerFlame(gas, grid=initial_grid)
flame.burner.mdot = mdot_mix

=======
flame = ct.BurnerFlame(gas,width=width)
flame.burner.mdot = mdot_mix

# Set refinement criteria (you can tweak these tolerances)
flame.set_refine_criteria(ratio=6.0, slope=0.005, curve=0.01)

# Solve with automatic grid refinement
flame.solve(loglevel=1, auto=True, refine_grid=True)
>>>>>>> af1e39737ec0a5ccd90e2399aff2223894492a21

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

<<<<<<< HEAD
>>>>>>> b3dfde696c6dceaa02ce7c3ce14dea7522bec72f
=======
>>>>>>> af1e39737ec0a5ccd90e2399aff2223894492a21
plt.figure()
plt.plot(flame.grid, flame.T)
plt.xlabel("Distance from burner [m]")
plt.ylabel("Temperature [K]")
plt.title("Temperature profile")
plt.grid(True)
plt.tight_layout()
plt.show()

<<<<<<< HEAD
<<<<<<< HEAD
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
=======
=======
>>>>>>> af1e39737ec0a5ccd90e2399aff2223894492a21
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
<<<<<<< HEAD
>>>>>>> b3dfde696c6dceaa02ce7c3ce14dea7522bec72f
=======
>>>>>>> af1e39737ec0a5ccd90e2399aff2223894492a21
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



