from CoolProp.CoolProp import PropsSI

'''
Class which defines the flight conditions for the cooling system

'''

class FlightCondition:
    def __init__(self, name, T_amb, P_amb, RH_amb, V, power_required, power_split, thermal_efficiency, propulsive_efficiency, P_cc, T_cc):
        self.name = name
        self.T_amb = T_amb  # Ambient temperature in K
        self.P_amb = P_amb  # Ambient pressure in Pa
        self.RH_amb = RH_amb  # Relative humidity of the ambient air
        self.gamma = PropsSI('CPMASS', 'T', self.T_amb, 'P', self.P_amb, 'Air')/PropsSI('CVMASS', 'T', self.T_amb, 'P', self.P_amb, 'Air')
        self.V = V  # Velocity in m/s
        self.Mach = V / PropsSI('A', 'T', self.T_amb, 'P', self.P_amb, 'Air')
        self.power_required = power_required  # Power required in W
        self.power_split = power_split   # Proportion of power outputted from the fuel cell
        self.thermal_efficiency = thermal_efficiency
        self.propulsive_efficiency = propulsive_efficiency
        self.P_cc = P_cc
        self.T_cc = T_cc
        self.T_tot = self.T_amb * (1 + (self.gamma - 1)/2 * (self.Mach) ** 2)  # Total pressure in Pa
        self.P_tot = self.P_amb * (1 + (self.gamma - 1)/2 * (self.Mach) ** 2) ** (self.gamma/(self.gamma - 1)) # Total pressure in Pa

