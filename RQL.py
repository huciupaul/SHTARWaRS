import cantera as ct
from mixture_properties import mixture_properties
from typing import Union, Mapping
import matplotlib.pyplot as plt
from pathlib import Path




######################## SUPPORT FUNCTIONS ######################################################

LHV_KEROSENE = 43.15e6   # J kg⁻¹   (≈ 43.15 MJ kg⁻¹)
LHV_H2       = 120.0e6   # J kg⁻¹   (≈ 120 MJ kg⁻¹)

def kerosene_to_h2(mdot_kerosene: float,
                   *,
                   lhv_kerosene: float = LHV_KEROSENE,
                   lhv_h2: float = LHV_H2) -> float:
    """
    Parameters
    ----------
    mdot_kerosene : float
        Mass-flow of Jet-A / kerosene [kg s⁻¹].

    Returns
    -------
    float
        Hydrogen mass-flow [kg s⁻¹] that releases the same heat at its LHV.
    """
    if mdot_kerosene < 0:
        raise ValueError("Mass-flow must be non-negative.")
    return mdot_kerosene * lhv_kerosene / lhv_h2

def _as_mole_str(x: Union[str, Mapping[str, float]]) -> str:
    """Return *x* as a Cantera mole-fraction string."""
    if isinstance(x, str):
        return x
    return ", ".join(f"{k}:{v}" for k, v in x.items())


def air_mass_flow_for_phi(
    mdot_H2: float,
    phi: float,
    *,
    X_air: Union[str, Mapping[str, float]] ,
    mech: str = "SanDiego.yaml",
    T: float = 298.15,
    P: float = 101_325.0,
) -> float:
    """
    Mass-flow of **air** (kg/s) needed to burn *mdot_H2* (kg/s) at the
    requested equivalence ratio *phi*.

    φ = (F/A)_actual / (F/A)_stoich.
    Therefore  mdot_air = mdot_air_stoich / φ.

    Parameters
    ----------
    mdot_H2 : float
        Hydrogen mass-flow [kg s⁻¹].
    phi : float
        Equivalence ratio (φ > 1 rich, φ < 1 lean).
    X_air : str | dict
        Air composition (mole basis).  Can be the default string or a dict.
    mech : str
        Mechanism file containing H2, O2, N2, etc.  Default: SanDiego.yaml.
    T, P : float
        State at which to create the `Solution` object (only for property look-up).

    Returns
    -------
    float
        Required air mass-flow [kg s⁻¹].

    Raises
    ------
    ValueError
        If *phi* ≤ 0 or *mdot_H2* < 0.
    """
    if mdot_H2 < 0:
        raise ValueError("mdot_H2 must be non-negative.")
    if phi <= 0:
        raise ValueError("phi must be positive.")

    gas = ct.Solution(mech)
    i_H2 = gas.species_index("H2")
    MW_H2 = gas.molecular_weights[i_H2] / 1000.0  # kg mol⁻¹

    # Stoichiometric O2 requirement: ½ mol O2 per mol H2  (2 H₂ + O₂ → 2 H₂O)
    n_H2 = mdot_H2 / MW_H2        # mol s⁻¹ H2
    n_O2_stoich = 0.5 * n_H2      # mol s⁻¹ O2
    MW_O2 = gas.molecular_weights[gas.species_index("O2")] / 1000.0
    mdot_O2_req = n_O2_stoich * MW_O2

    # Look up O2 mass-fraction in the *specified* air mixture
    gas.TPX = T, P, _as_mole_str(X_air)
    Y_O2_air = gas.mass_fraction_dict()["O2"]  # dimensionless

    # Stoichiometric air mass-flow (kg/s)
    mdot_air_stoich = mdot_O2_req / Y_O2_air

    # Desired air mass-flow for the given φ
    return mdot_air_stoich / phi
















"""


0: At CC inlet. Inlet conditions
1: At Rich Combustion exhaust
2: At mixer inlet
3:At mixer outlet.
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
T_0_kelvin = T_0 +273.15


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
plt.figure()
plt.plot(f.grid, data.T, label='Temperature (K)')
plt.xlabel('Distance (m)')
plt.ylabel('Temperature (K)')
plt.title('Temperature Profile')
plt.legend()
plt.grid()

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
plt.show()




















































