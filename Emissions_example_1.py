# A burner-stabilized lean premixed hydrogen-oxygen flame at low pressure Cantera example

from pathlib import Path
import matplotlib.pyplot as plt
import cantera as ct

p = 0.05 * ct.one_atm
tburner = 373.0
mdot = 0.06
reactants = 'H2:1.5, O2:1, AR:7'  # premixed gas composition
width = 0.5  # m
loglevel = 1  # amount of diagnostic output (0 to 5)

gas = ct.Solution('h2o2.yaml')
gas.TPX = tburner, p, reactants

f = ct.BurnerFlame(gas, width=width)
f.burner.mdot = mdot
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show()

f.transport_model = 'mixture-averaged'
f.solve(loglevel, auto=True)

if "native" in ct.hdf_support():
    output = Path() / "burner_flame.h5"
else:
    output = Path() / "burner_flame.yaml"
output.unlink(missing_ok=True)

f.save(output, name="mix", description="solution with mixture-averaged transport")

f.transport_model = 'multicomponent'
f.solve(loglevel)  # don't use 'auto' on subsequent solves
f.show()
f.save(output, name="multi", description="solution with multicomponent transport")

f.save('burner_flame.csv', basis="mole", overwrite=True)

# Plot temperature and species profiles
z = f.flame.grid  # spatial grid
temperature = f.T  # temperature profile

# Plot temperature
plt.figure()
plt.plot(z, temperature, label='Temperature (K)')
plt.xlabel('Distance (m)')
plt.ylabel('Temperature (K)')
plt.title('Temperature Profile')
plt.legend()
plt.grid()

# Plot species profiles
plt.figure()
for species in ['H2', 'O2', 'H2O', 'AR']:
    plt.plot(z, f.X[gas.species_index(species)], label=species)
plt.xlabel('Distance (m)')
plt.ylabel('Mole Fraction')
plt.title('Species Profiles')
plt.legend()
plt.grid()

plt.show()