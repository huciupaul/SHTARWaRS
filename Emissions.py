# A burner-stabilized lean premixed hydrogen-oxygen flame

from pathlib import Path
import matplotlib.pyplot as plt
import cantera as ct




# Default dry‑air composition, mole fractions (ISO 2533 “Standard Atmosphere”
# with the current best CO₂ abundance 407ppm)
X_air = 'N2:0.78084,  O2:0.20946,  AR:0.00934,  CO2:0.000407 '


# Species name used in the mechanism file
_H2_SPECIES = "H2"

p = 0.05 * ct.one_atm
tburner = 373.0
mdot = 0.06
reactants = 'H2:1.5, O2:1, AR:7'  # premixed gas composition
width = 0.5  # m
loglevel = 1  # amount of diagnostic output (0 to 5)

gas = ct.Solution('SanDiego.yaml')
gas.TPX = tburner, p, reactants

f = ct.BurnerFlame(gas, width=width)
f.burner.mdot = mdot
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.05)
f.show()


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