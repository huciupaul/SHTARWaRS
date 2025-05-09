import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# USER-DEFINED INPUTS (can be modified for different cases):
hydrogen_mass_flow = 0.1    # kg/s, total hydrogen fuel mass flow rate
phi_rich = 2.0             # Equivalence ratio in the rich combustor stage (fuel/air ratio relative to stoichiometric)
water_mass = 0.02          # kg/s (or kg per batch) of liquid water injected in the quench stage

# MECHANISM SELECTION:
# Use a mechanism with H2/O2 combustion and NOx chemistry. The user may specify a different mechanism file if desired.
# Examples: "gri30.yaml" (GRI-Mech 3.0 for methane, includes NOx),
#           "Mechanism.xml" (a custom detailed mechanism like UCSD or Glarborg NOx mechanism).
mechanism_file = "SanDiegoNitrogen.yaml"  # default mechanism; change as needed

# Initialize gas object for the chosen mechanism
gas = ct.Solution(mechanism_file)

# Simulation parameters:
pressure = ct.one_atm  # operating pressure (constant throughout stages). Use 1 atm as default.
# (Note: Pressure can be adjusted if needed, e.g., for gas turbine conditions. Make sure mechanism is valid at that pressure.)

# Calculate stoichiometric oxygen required for the given hydrogen flow:
# Reaction: 2 H2 + O2 -> 2 H2O. So 1 mole H2 needs 0.5 mole O2.
# First, find moles of H2 per second from the mass flow:
H2_mw = gas.molecular_weights[gas.species_index("H2")]  # molecular weight of H2 (≈2.016 g/mol)
mol_H2_per_s = hydrogen_mass_flow / H2_mw              # [kmol/s] of H2 (since mass in kg, mw in kg/kmol)
# Stoichiometric O2 (kmol/s) needed to fully burn that H2:
mol_O2_stoich = 0.5 * mol_H2_per_s

# Determine the air flow to achieve the rich equivalence ratio:
# Equivalence ratio phi = (fuel/O2)_actual / (fuel/O2)_stoich.
# Rearranged: actual O2 = stoich O2 / phi.
mol_O2_rich = mol_O2_stoich / phi_rich   # O2 moles/s in the rich stage (less than stoich if phi > 1)
# Corresponding N2 moles/s that come with this O2 (assuming dry air with O2:N2 ≈ 1:3.76 by mole):
mol_N2_rich = mol_O2_rich * (3.76)       # N2 moles/s in the rich stage air (3.76 N2 per O2)

# Set up initial state for the rich combustor:
# We create a gas mixture with H2, O2, N2 at the rich equivalence ratio.
# Assume the reactants enter at ambient temperature (e.g. 300 K). We'll then force ignition by raising the temperature.
initial_T = 300.0  # K, initial reactants temperature before ignition (can be ambient)
# Define the composition for the rich stage (as mole fractions):
gas.TPX = initial_T, pressure, {
    "H2": mol_H2_per_s,        # use moles in proportion to mass flow (for 1 s of fuel)
    "O2": mol_O2_rich,
    "N2": mol_N2_rich
}
# The above sets the state for a mixture representing 1 second of flow (for batch simulation).
# We now raise the temperature to force ignition (since the mixture may not autoignite at 300 K).
gas.TP = 1000.0, pressure  # increase T to ~1000 K to initiate combustion in the rich stage

# RICH COMBUSTOR REACTOR SIMULATION:
rich_reactor = ct.IdealGasConstPressureReactor(gas)  # constant-pressure adiabatic reactor with the rich mixture
sim = ct.ReactorNet([rich_reactor])

# Integrate reactor in time until reaction reaches steady state (no more temperature change).
# We'll integrate until a certain time or until dT/dt is very small.
time = 0.0
end_time = 0.5  # seconds (this is ample time for completion given H2 fast kinetics; adjust if needed)
# We will take small time steps and monitor temperature change for steady-state.
prev_temp = rich_reactor.T
for t in np.linspace(0, end_time, 200):
    sim.advance(t)
    # Check if temperature change is negligible:
    if abs(rich_reactor.T - prev_temp) < 1e-7:
        break
    prev_temp = rich_reactor.T

# Record rich stage results:
T_rich = rich_reactor.T
rich_gas = rich_reactor.thermo  # gas after rich combustor
# Mole fractions of NO, NO2, N2O after rich stage:
X_NO_rich = rich_gas[X := "NO"] if "NO" in rich_gas.species_names else 0.0
X_NO2_rich = rich_gas[X := "NO2"] if "NO2" in rich_gas.species_names else 0.0
X_N2O_rich = rich_gas[X := "N2O"] if "N2O" in rich_gas.species_names else 0.0

# QUENCH (MIXING) STAGE:
# Perfectly mix the rich combustor effluent with additional air and water (at lower temperature).
# Additional air to add (to go from phi_rich to lean):
# We assume we want to fully burn the remaining H2 (and have some excess O2). For simplicity, add enough air to make overall phi < 1.
# One approach: add the remaining stoichiometric O2 for leftover H2 (so total O2 = stoich), plus some excess (for lean burn).
# Here, we'll add O2 equal to the remaining stoichiometric requirement. This will yield overall phi ≈ 1 (slightly lean if water is added).
mol_O2_quench = mol_O2_stoich - mol_O2_rich   # O2 moles needed to reach stoichiometric for total fuel
# (If an explicitly lean final phi is desired, multiply mol_O2_quench by (>1) factor to supply excess O2.)
# Add corresponding N2 with this O2:
mol_N2_quench = mol_O2_quench * 3.76

# Create a gas object for the quench additives (air and water vapor).
air_water_gas = ct.Solution(mechanism_file)
# We will include water as H2O vapor in the composition for energy balance purposes (water will be vaporized).
# Assume air and water are at ambient conditions (e.g., 300 K) before mixing.
T_air = 300.0  # K
# Convert injected water mass to moles of H2O:
H2O_mw = gas.molecular_weights[gas.species_index("H2O")]
mol_H2O_quench = water_mass / H2O_mw  # moles of water added
# Set state of the air+water mixture (before mixing) at T_air:
air_water_gas.TPX = T_air, pressure, {
    "O2": mol_O2_quench,
    "N2": mol_N2_quench,
    "H2O": mol_H2O_quench
}

# Now perform adiabatic mixing of rich gas and quench gas at constant pressure.
# We assume no reaction during mixing, so total enthalpy is conserved.
# Calculate total moles of each stream for mixing:
# For rich_reactor, we can get moles via mass and mean molecular weight:
moles_rich = rich_gas.mass / rich_gas.mean_molecular_weight  # total moles in rich gas
moles_quench = mol_O2_quench + mol_N2_quench + mol_H2O_quench  # total moles in added air+water
# Get total enthalpy (J) of each stream:
H_rich = rich_gas.enthalpy_mass * rich_gas.mass   # (J) total enthalpy of rich gas
H_quench = air_water_gas.enthalpy_mass * air_water_gas.mass  # (J) total enthalpy of added air+water
H_total = H_rich + H_quench
# Total moles after mixing:
total_moles_mix = moles_rich + moles_quench
# Set the mixed gas state by specifying pressure, total moles composition, and finding temperature by energy balance:
# We combine the species moles from rich_gas and quench gas:
mixed_species_moles = {}
for i, sp in enumerate(gas.species_names):
    # moles from rich gas:
    X_rich = rich_gas.mole_fraction(i)
    mixed_species_moles[sp] = X_rich * moles_rich
    # moles from quench gas:
    if sp in air_water_gas.species_names:
        X_quench = air_water_gas[sp].X  # mole fraction of sp in quench gas
        mixed_species_moles[sp] += X_quench * moles_quench

# Set the composition of gas for the mixture (with an arbitrary T initially),
gas.TPX = 500.0, pressure, mixed_species_moles  # initial guess T=500 K
# Solve for the adiabatic mixing temperature: find T such that enthalpy matches H_total at constant P.
gas.HP = H_total / total_moles_mix, pressure  # set state by total specific enthalpy and pressure
T_quench = gas.T  # this is the temperature after mixing (before lean combustion)
# Update composition after quench (gas already has the mixture composition)
mixed_gas = ct.Solution(mechanism_file)
mixed_gas.TPX = T_quench, pressure, gas.X  # use the computed mixture state
# (At this point, mixed_gas contains the composition and state entering the lean combustor.)

# Capture NOx mole fractions after quench (no new formation, just diluted from rich stage):
X_NO_quench = mixed_gas[X := "NO"] if "NO" in mixed_gas.species_names else 0.0
X_NO2_quench = mixed_gas[X := "NO2"] if "NO2" in mixed_gas.species_names else 0.0
X_N2O_quench = mixed_gas[X := "N2O"] if "N2O" in mixed_gas.species_names else 0.0

# LEAN COMBUSTOR REACTOR SIMULATION:
# Set initial state for lean combustor with the mixed gas composition.
# If the post-quench temperature is not high enough to ignite the lean mixture, we force ignition similarly by raising T.
lean_initial_T = T_quench
if lean_initial_T < 1000.0:
    lean_initial_T = 1000.0  # ensure a high enough start temperature for ignition (adjust as needed)
gas.TPX = lean_initial_T, pressure, mixed_gas.X  # set gas to the mixed composition with elevated T if needed

# Create constant-pressure reactor for the lean stage and integrate to steady state
lean_reactor = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([lean_reactor])
time = 0.0
end_time = 0.5  # integrate up to 0.5 s (should reach steady state well before this)
prev_temp = lean_reactor.T
for t in np.linspace(0, end_time, 200):
    sim.advance(t)
    if abs(lean_reactor.T - prev_temp) < 1e-7:
        break
    prev_temp = lean_reactor.T

# Record lean stage results:
T_lean = lean_reactor.T
lean_gas = lean_reactor.thermo
X_NO_lean = lean_gas[X := "NO"] if "NO" in lean_gas.species_names else 0.0
X_NO2_lean = lean_gas[X := "NO2"] if "NO2" in lean_gas.species_names else 0.0
X_N2O_lean = lean_gas[X := "N2O"] if "N2O" in lean_gas.species_names else 0.0

# OUTPUT RESULTS:
print(f"Rich stage temperature: {T_rich:.1f} K")
print(f"Rich stage NO mole fraction: {X_NO_rich:.3e}")
print(f"Rich stage NO2 mole fraction: {X_NO2_rich:.3e}")
print(f"Rich stage N2O mole fraction: {X_N2O_rich:.3e}\n")

print(f"After quench (mix) temperature: {T_quench:.1f} K")
print(f"After quench NO mole fraction: {X_NO_quench:.3e}")
print(f"After quench NO2 mole fraction: {X_NO2_quench:.3e}")
print(f"After quench N2O mole fraction: {X_N2O_quench:.3e}\n")

print(f"Lean stage temperature: {T_lean:.1f} K")
print(f"Lean stage NO mole fraction: {X_NO_lean:.3e}")
print(f"Lean stage NO2 mole fraction: {X_NO2_lean:.3e}")
print(f"Lean stage N2O mole fraction: {X_N2O_lean:.3e}\n")

# Tabulate results for clarity:
stages = ["Rich Combustor", "After Quench", "Lean Combustor"]
temps = [T_rich, T_quench, T_lean]
NO_vals = [X_NO_rich, X_NO_quench, X_NO_lean]
NO2_vals = [X_NO2_rich, X_NO2_quench, X_NO2_lean]
N2O_vals = [X_N2O_rich, X_N2O_quench, X_N2O_lean]

print("Stage\t\tTemperature (K)\tX_NO\t\tX_NO2\t\tX_N2O")
for i, stage in enumerate(stages):
    print(f"{stage:16s}\t{temps[i]:.1f}\t\t{NO_vals[i]:.3e}\t{NO2_vals[i]:.3e}\t{N2O_vals[i]:.3e}")

# PLOTTING RESULTS:
# Plot temperature and NOx mole fractions at each stage.
plt.figure(figsize=(6,4))
plt.subplot(1,2,1)
plt.bar(stages, temps, color='tab:orange')
plt.ylabel("Temperature (K)")
plt.title("Temperature at Each Stage")
plt.xticks(rotation=20)

plt.subplot(1,2,2)
# Use logarithmic scale for NOx because values span orders of magnitude
plt.yscale('log')
plt.plot(stages, NO_vals, marker='o', label='NO')
plt.plot(stages, NO2_vals, marker='s', label='NO2')
plt.plot(stages, N2O_vals, marker='^', label='N2O')
plt.ylabel("Mole Fraction (log scale)")
plt.title("NOx Mole Fractions")
plt.xticks(rotation=20)
plt.legend()
plt.tight_layout()
plt.show()
