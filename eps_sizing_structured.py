import numpy as np
from common.constants import P_TOGA # update once there
import matplotlib.pyplot as plt

class EPS:
    def __init__(self,
                 P_TOGA: float = P_TOGA, # TOGA power [W]
                 fc_split: float = 0.0,
                 Cable_length: float = 20, # Cable length [m]
                 Spec_W_cable: float = 0.0015, # Cables specific weight [kg/kW/m]
                 lambd: float = 0.5,
                 ): 
        self.P_TOGA = P_TOGA
        self.fc_split = fc_split
        self.P_req = P_TOGA * fc_split
        self.L = Cable_length
        self.Spec_W_cable = Spec_W_cable
        self.lambd = lambd

    def motor_weight(self):
        """Calculate the weight of the EPS based on the required power."""
        return self.P_req / 10 / 1000

    def inverter_weight(self):
        return 0.97 * self.P_req / 20 / 1000

    def cable_weight(self):
        return self.L * self.Spec_W_cable * self.P_req / 1000

    def heat_dissipation_motor(self):
        """Calculate the heat dissipation of the EPS."""
        return self.P_req * (1/0.97 - 1)

    def heat_dissipation_inverter(self):
        """Calculate the heat dissipation of the EPS."""
        return self.P_req * (1/0.99 - 1)

    

if __name__ == "__main__":
    eps = EPS(fc_split = fc_split, Cable_length = 20, Spec_W_cable = 0.0015, lambd = 0.5)
    print("Motor weight:", eps.motor_weight())
    print("Inverter weight:", eps.inverter_weight())
    print("Cable weight:", eps.cable_weight())
    print("Total EPS weight:", eps.motor_weight() + eps.inverter_weight() + eps.cable_weight())
    print("Heat dissipation motor:", eps.heat_dissipation_motor())
    print("Heat dissipation inverter:", eps.heat_dissipation_inverter())
    print("Total EPS heat dissipation:", eps.heat_dissipation_motor() + eps.heat_dissipation_inverter())


'''
    def motor_sizing(self):
        """Calculate the EPS sizing parameters."""
        c1 = 0.000805
        c2 = 0.00161
        c3 = 0.1934
        hs = 0.025
        P_req2 = self.P_req / 2
        r = (self.lambd * P_req2)**(1/3) * c1
        l = c2 * 1/self.lambd * (self.lambd * P_req2)**(1/3)
        hy = r * c3
        r_out = r + hs + hy
        r_tot = ((0.737 * (r_out ** 2)) - 0.580 * r_out + 1.1599) * r_out
        return r_tot, l

    def motor_volume(radius, length):
        """Calculate the volume of the EPS."""
        return np.pi * radius**2 * length

    def inverter_volume(self):
        return 0.97 * self.P_req / 18.7
    
    def fc_power_needed(self):
        """Calculate the total power needed from the fuel cell."""
        return self.P_req / (0.97*0.99)
'''