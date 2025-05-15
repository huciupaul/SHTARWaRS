'''
Class which defines the flight conditions for the cooling system

'''

class FlightCondition:
    def __init__(self, name, T_amb, P_amb, RH_amb, V, power_required, power_split, thermal_efficiency, propulsive_efficiency, P_cc):
        self.name = name
        self.T_amb = T_amb  # Ambient temperature in K
        self.P_amb = P_amb  # Ambient pressure in Pa
        self.RH_amb = RH_amb  # Relative humidity of the ambient air
        self.V = V  # Velocity in m/s
        self.power_required = power_required  # Power required in W
        self.power_split = power_split   # Proportion of power outputted from the fuel cell
        self.thermal_efficiency = thermal_efficiency
        self.propulsive_efficiency = propulsive_efficiency
        self.P_cc = P_cc

