import numpy as np
import matplotlib.pyplot as plt

P_req = 1000  # Required power in watts

def eps_weight(P):
    """Calculate the weight of the EPS based on the required power."""
    return P_req / 10 + 0.97 * P_req /25

def heat_dissipation(P):
    """Calculate the heat dissipation of the EPS."""
    return P_req * (1 - 0.97) + P_req * (1 - 0.97)**2

def fc_power_needed(P_req):
    """Calculate the total power needed from the fuel cell."""
    return P_req / (0.97**2)

def motor_sizing(lambd, P_req):
    """Calculate the EPS sizing parameters."""
    c1 = 0.0004
    c2 = 0.0008
    c3 = 1.0744
    hs = 0.0366
    P_req = P_req / 2
    r = (lambd * P_req)**(1/3) * c1
    l = (lambd * P_req)**(1/3) * c2 * 1/lambd
    hy = r * c3
    r_tot = r + hs + hy
    return r_tot, l

def motor_volume(radius, length):
    """Calculate the volume of the EPS."""
    return np.pi * radius**2 * length