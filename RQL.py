import cantera as ct
from mixture_properties import mixture_properties
from typing import Union, Mapping
import matplotlib.pyplot as plt
from pathlib import Path


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
width        = 0.001                        # m
initial_grid = np.linspace(0.0, width, 800)

# -----------------------------------------------------------------------------
# 4 – BUILD & SOLVE THE FLAME
# -----------------------------------------------------------------------------
gas = ct.Solution(MECH)
gas.TPX = T_0_C + C2K, P_0_atm * ATM2PA, X_mix


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

>>>>>>> b3dfde696c6dceaa02ce7c3ce14dea7522bec72f
plt.figure()
plt.plot(flame.grid, flame.T)
plt.xlabel("Distance from burner [m]")
plt.ylabel("Temperature [K]")
plt.title("Temperature profile")
plt.grid(True)
plt.tight_layout()
plt.show()

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
>>>>>>> b3dfde696c6dceaa02ce7c3ce14dea7522bec72f
plt.show()



quit()








































##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################







import pandas as pd
import cantera as ct

# ────────────────────────────────────────────────────────────────────────────────
# 1) Read your combustion results from CSV
# ────────────────────────────────────────────────────────────────────────────────
# Assume 'combustion_results.csv' has columns: time, T (K), P (Pa), and species mass
df = pd.read_csv('burner_flame.csv')

# Take the last row as the "combustion product" state
prod = df.iloc[-1]
T_prod = prod['T']

# build a dict of species mass fractions
species_cols = [c for c in df.columns if c not in ('time','T','P')]
Y_prod = {s: prod[s] for s in species_cols}

# Create a Cantera Solution for the product
gas_prod = ct.Solution('SanDiegoNitrogen.yaml')
gas_prod.TPX = T_prod, P_0_pascals, Y_prod

# ────────────────────────────────────────────────────────────────────────────────
# 2) Define the “air” and “water” streams
# ────────────────────────────────────────────────────────────────────────────────
# Air: 21% O2, 79% N2 by mole, at T_air
T_air = 350.0  # K, your choice
gas_air = ct.Solution('gri30.yaml')
gas_air.TPX = T_air, P_0_pascals, 'O2:0.21, N2:0.79'

# Water vapor: pure H2O at T_water
T_water = 300.0  # K, your choice
gas_water = ct.Solution('gri30.yaml')
gas_water.TPX = T_water, P_0_pascals, 'H2O:1.0'

# ────────────────────────────────────────────────────────────────────────────────
# 3) Set up the 0D mixer reactor network (no reaction)
# ────────────────────────────────────────────────────────────────────────────────
# Create reservoirs for each stream
res_prod  = ct.Reservoir(gas_prod)
res_air   = ct.Reservoir(gas_air)
res_water = ct.Reservoir(gas_water)

# Create the “mixer” reactor: IdealGasReactor at constant volume
gas_mix = ct.Solution('gri30.yaml')
reactor = ct.IdealGasReactor(gas_mix, energy='on', volume=1.0)

# Hook up three MassFlowControllers to feed the reactor
# (choose mass flow rates in kg/s to get the mix ratio you want)
mdot_prod  = 1.0    # e.g. 1 kg/s of product
mdot_air   = 0.5    # 0.5 kg/s of air
mdot_water = 0.1    # 0.1 kg/s of water

mfc1 = ct.MassFlowController(res_prod,  reactor, mdot=mdot_prod)
mfc2 = ct.MassFlowController(res_air,   reactor, mdot=mdot_air)
mfc3 = ct.MassFlowController(res_water, reactor, mdot=mdot_water)

# Use a PressureController (valve) to let the mixed gas escape
# and hold the reactor near P_prod
outlet = ct.Reservoir(gas_mix)
# “master” can be any upstream controller; here we link it to mfc1
pc = ct.PressureController(reservoir=outlet, reactor=reactor, master=mfc1, K=1.0)

# Build the network and advance to steady‐state (e.g. 2 s)
net = ct.ReactorNet([reactor])
net.advance(2.0)

# ────────────────────────────────────────────────────────────────────────────────
# 4) Report the mixed‐stream state
# ────────────────────────────────────────────────────────────────────────────────
print(f"Mixed reactor temperature: {reactor.T:.2f} K")
print(f"Mixed reactor pressure:    {reactor.thermo.P/ct.one_atm:.3f} atm")
print("Mixed composition (mole fractions):")
for i, sp in enumerate(gas_mix.species_names):
    x = reactor.thermo.X[i]
    if x > 1e-6:
        print(f"  {sp:6s}: {x:.4f}")


