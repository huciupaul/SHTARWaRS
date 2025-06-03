'''
A class which contains the information of the fuel cell


'''

class FuelCell:
    def __init__(self, name, stack_efficiency, T, P_A, P_C, RH_A, RH_C, stoic_ratio_A, stoic_ratio_C, spec_power, spec_vol_power):
        self.name = name
        self.stack_efficiency = stack_efficiency
        self.T = T # Temperature in K
        self.P_A = P_A # Pressure of the anode in Pa
        self.P_C = P_C # Pressure of the cathode in Pa
        self.RH_A = RH_A # Relative humidity of the anode
        self.RH_C = RH_C # Relative humidity of the cathode
        self.stoic_ratio_A = stoic_ratio_A # Stoichiometric ratio of the anode
        self.stoic_ratio_C = stoic_ratio_C # Stoichiometric ratio of the cathode
        self.spec_power = spec_power # Specific power of the fuel cell in W/kg
        self.spec_vol_power = spec_vol_power # Specific volume power of the fuel cell in W/m^3

