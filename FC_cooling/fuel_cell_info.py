'''
A class which contains the information of the fuel cell


'''

class FuelCell:
    def __init__(self, name, stack_efficiency, T, P_A, P_C, RH_A, RH_C, stoic_ratio_A, stoic_ratio_C):
        self.name = name
        self.stack_efficiency = stack_efficiency
        self.T = T
        self.P_A = P_A
        self.P_C = P_C
        self.RH_A = RH_A
        self.RH_C = RH_C
        self.stoic_ratio_A = stoic_ratio_A
        self.stoic_ratio_C = stoic_ratio_C

