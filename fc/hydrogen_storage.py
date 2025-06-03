from CoolProp.CoolProp import PropsSI


'''
Class which contains the information of the state of hydrogen stored

'''

class HydrogenStorage:
    def __init__(self, name, T, P):
        self.name = name
        self.T = T  # Temperature in K
        self.P = P  # Pressure in Pa
        self.H_init = PropsSI('H', 'T', self.T, 'P', self.P, 'Hydrogen')
