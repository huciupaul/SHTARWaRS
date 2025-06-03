import numpy as np
import matplotlib.pyplot as plt
from PIL.WmfImagePlugin import WmfHandler

P_req = 3000000#[W] Required power in watts

def motor_weight(P):
    """Calculate the weight of the EPS based on the required power."""
    return P_req / 10 / 1000

def inverter_weight(P):
    return 0.97 * P_req /20 / 1000

def heat_dissipation(P):
    """Calculate the heat dissipation of the EPS."""
    return P_req * (1 - 0.97) + P_req * (1 - 0.97)**2

def fc_power_needed(P_req):
    """Calculate the total power needed from the fuel cell."""
    return P_req / (0.97**2)

def motor_sizing(lambd, P_req):
    """Calculate the EPS sizing parameters."""
    c1 = 0.000805
    c2 = 0.00161
    c3 = 0.1934
    hs = 0.025
    P_req = P_req / 2
    r = (lambd * P_req)**(1/3) * c1
    l = c2 * 1/lambd * (lambd * P_req)**(1/3)
    hy = r * c3
    r_tot = r + hs + hy
    return r_tot, l

def motor_volume(radius, length):
    """Calculate the volume of the EPS."""
    return np.pi * radius**2 * length

def inverter_volume(P):
    return 0.97 * P / 18.7

x = motor_sizing(0.5, P_req)
V = motor_volume(x[0], x[1])
Wm = motor_weight(P_req)
Wi = inverter_weight(P_req)
print(Wm, Wi)
print(x)
print(V)