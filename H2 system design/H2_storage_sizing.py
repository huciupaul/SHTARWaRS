from CoolProp.CoolProp import PropsSI
import numpy as np

# Constants
LHV_H2    = 120 * 1e6    # J/kg
LHV_ja1   = 42.8 * 1e6  # J/kg

# Flight parameters:
m_ja1     = 2022  # kg , mass of conventional fuel

# Fuel consumption of the ja1 engine:
SFC_ja1   = np.array([0.509, 0.680])
SFC_h2    = SFC_ja1 * LHV_ja1 / LHV_H2  
H_ja1     = 18540
H_h2      = 60920
J         = 778.24

eta_tot_cc_ja1 = 550 * 3600 / (SFC_ja1 * H_ja1 * J)
eta_tot_cc_h2 = 550 * 3600 / (SFC_h2 * H_h2 * J)

# print(f"eta_cc_ja1: {eta_tot_cc_ja1}") 
# print(f"eta_cc_h2: {eta_tot_cc_h2}")


# FC efficiencies:
eta_stack = 0.70  # Efficiency of the fuel cell
eta_bop   = 0.98  # Efficiency of the balance of plant components
eta_ll    = 0.98  # Efficiency of the line loss
eta_em    = 0.99  # Efficiency of the electric motor
eta_s     = 0.95  # Efficiency of the shaft system
eta_fc = eta_stack * eta_bop * eta_ll * eta_em * eta_s  # Efficiency of the fuel cell system


# 1. Estimating total energy required for a flight:
#   1.1 From power profile:
V_hol = 54.0167 # m/s
E_tot_fp = 6347 * 1e6  # J, estimated energy required for flight (from power profile)
eta_prop = 0.85  # Propulsion efficiency
S = 28.79  # m^2, wing area

q_hol = 0.5 * 1.17295 * V_hol**2  # m^2/s^2, dynamic pressure
T = q_hol
P_hol = T * V_hol / eta_prop
cs_23_vel = np.array([[5*3600, V_hol], [30*3600, V_hol]])

# print(f"Estimated energy required for flight from flight profile:  MJ")
# print(E_tot_fp * 1e-6)



#   1.2 From current ja1 mass:
E_tot_current_mass = m_ja1 * eta_tot_cc_ja1 * LHV_ja1
# print(f"Estimated energy required for flight from current fuel mass (MJ): ")
# print(E_tot_current_mass * 1e-6)

E_tot = E_tot_fp

# 2. Estimating mass of hydrogen required for flight:
#   2.1 Using only combustion engine:
m_H2_cc = E_tot / (eta_tot_cc_h2 * LHV_H2)  # kg
# print(f"Mass of H2 required for combustion engine (kg):")
# print(m_H2_cc)

#   2.2 Using only fuel cell:
m_H2_fc = E_tot / (eta_fc * LHV_H2)  # kg
# print(f"Mass of H2 required for fuel cell (kg):")
# print(m_H2_fc)


h2 = 'ParaHydrogen'
# Preliminary volume estimation for H2 storage tank
# Assuming mass of H2 is m_H2 and it is stored as liquid hydrogen:
m_h2 = 0.5 * m_H2_cc[0] + 0.5 * m_H2_fc  # kg, mass of hydrogen

pressure_final = 11 * 1e5  # Pa, pressure of hydrogen storage tank
pressure_start = 1 * 1e5  # Pa, initial pressure of hydrogen storage tank

temp_final = PropsSI('T', 'P', pressure_final, 'Q', 0, h2)  # K, temperature of hydrogen at 11 atm and 20 K	
temp_init = PropsSI('T', 'P', pressure_start, 'Q', 0, h2)  # K, temperature of hydrogen at 1 atm and 20 K

rho_gas_final = PropsSI('D', 'P', pressure_final, 'Q', 1, h2)  # kg/m^3, density of hydrogen at 11 atm and 20 K

# Q_final = rho_gas * (1-fl)

print(f'final temperature: {temp_final} K')

H_final = PropsSI('H', 'P', pressure_final, 'Q', 1, h2)  # J/kg,

print(H_final)

# Uncompressed volume of hydrogen gas at 1 atm and 287.15 K
m_h2 = 0.5 * m_H2_cc[0] + 0.5 * m_H2_fc  # kg, mass of hydrogen
print(f"Mass of H2: {m_h2} kg")
rho_uncomp = PropsSI('D', 'P', 101325, 'T', 287.15, 'Hydrogen')  # kg/m^3, density of hydrogen at 1 atm and 287.15 K
V_uncomp = m_h2 / rho_uncomp  # m^3, volume of hydrogen at 1 atm and 287.15 K
print(f"Uncompressed H2 Volume: {V_uncomp} m^3")

# Compressed volume of hydrogen gas at 950 bar and 287.15 K
pressure_final = 700 * 1e5  # Pa, pressure of hydrogen storage tank
rho_comp = PropsSI('D', 'P', pressure_final, 'T', 287.15, 'Hydrogen')  # kg/m^3, density of hydrogen at 950 bar and 287.15 K
V_comp = m_h2 / rho_comp  # m^3, volume of hydrogen at 950 bar and 287.15 K
print(f"Compressed H2 Volume: {V_comp} m^3")

# Cryo-compressed volume of hydrogen gas at 950 bar and 25-110K temperature range
rho_ccomp_high = PropsSI('D', 'P', 50e5, 'T', 35, h2)
rho_ccomp_low = PropsSI('D', 'P', 700e5, 'T', 35, h2)
V_ccomp_low, V_ccomp_high = m_h2/rho_ccomp_high, m_h2/rho_ccomp_low
print(f"Cryo-Compressed H2 Volume: [{V_ccomp_low}, {V_ccomp_high}] m^3")


# Liuqid hydrogen volume at 1 atm and 20 K
pressure_final = 5 * 1e5  # Pa, pressure of hydrogen storage tank
rho_liquid = PropsSI('D', 'P', pressure_final, 'T', 20, h2)  # kg/m^3, density of hydrogen at 1 atm and 20 K
V_liquid = m_h2 / rho_liquid  # m^3, volume of hydrogen at 1 atm and 20 K
print(f"Liquid H2 Volume: {V_liquid} m^3")