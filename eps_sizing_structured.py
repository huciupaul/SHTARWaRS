import numpy as np
from common.constants import P_TOGA # update once there
import matplotlib.pyplot as plt

class EPS:
    def __init__(self,
                 P_TOGA: float = P_TOGA, # TOGA power [W]
                 fc_split: float = 0.0,
                 Cable_length: float = 20, # Cable length [m]   TODO: estimate the length better
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
        return self.P_req / 10 / 1000

    def inverter_weight(self):
        return 0.97 * self.P_req / 20 / 1000

    def cable_weight(self):
        return self.L * self.Spec_W_cable * self.P_req / 1000

    def heat_dissipation_motor(self):
        return self.P_req * (1/0.97 - 1)

    def heat_dissipation_inverter(self):
        return self.P_req * (1/0.99 - 1)
    
    def __motor_sizing(self):
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

    def motor_volume(self):
        r_tot, l = self.__motor_sizing()
        return np.pi * r_tot**2 * l

    def inverter_volume(self):
        return 0.97 * self.P_req / 18.7

    
def main(fc_split: float = 0.0):
    """
    Main function to calculate the EPS components weight, heat dissipation and volume.
    Args:
        fc_split (float): Fraction of TOGA power allocated to the fuel cell. 
    Returns:
        comp (np.ndarray): Array containing the weights, heat dissipations and volumes of the EPS components.
        tot (np.ndarray): Array containing the total weight and heat dissipation of the EPS (volume excluded as not relevant).
    """
    eps = EPS(fc_split=fc_split)

    comp = np.zeros(7)
    tot = np.zeros(2)

    comp = eps.motor_weight(), eps.inverter_weight(), eps.cable_weight(), eps.heat_dissipation_motor(), eps.heat_dissipation_inverter(), eps.motor_volume(), eps.inverter_volume()
    tot = np.sum(comp[:3]), np.sum(comp[3:5]) # volume not summed as not relevant

    return comp, tot


if __name__ == "__main__":
    comp, tot = main(fc_split=0.5)

    print("Motor weight (kg):", comp[0])
    print("Inverter weight (kg):", comp[1])
    print("Cable weight (kg):", comp[2])
    print("Total EPS weight (kg):", tot[0])
    print("Motor heat dissipation (W):", comp[3])
    print("Inverter heat dissipation (W):", comp[4])
    print("Total EPS heat dissipation (W):", tot[1])
    print("Motor volume (m^3):", comp[5])
    print("Inverter volume (m^3):", comp[6])