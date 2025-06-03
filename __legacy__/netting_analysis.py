import numpy as np

def netting_thickness(p_load, r_in, strength, FVF):
    
    #Calc membrane loads
    N_theta_p = p_load * r_in
    N_phi_p = p_load * r_in / 2
    strength = strength * .38 * .57 * .75 #may add safety factors
    angle = 30 #helical winding angle in degrees
    
    #Calc netting thickness in hoop and helical direction
    t_heli = N_phi_p / (strength * np.cos(angle)**2)
    t_hoop = (N_theta_p - N_phi_p * np.tan(angle)**2) / (strength)

    #Calc minimum thickness based on FVF
    t_min = (t_heli + t_hoop) / FVF

    return t_min