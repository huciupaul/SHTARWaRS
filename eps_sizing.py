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

def eps_sizing(P_req):
    """Calculate the EPS sizing parameters."""