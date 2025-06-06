import numpy as np
from typing import Tuple

from DSE_1.global_constants import Beechcraft_1900D, rho_cargo

def main(L_tank: float, r_tank: float, V_tank: float):
    # write a function of the position of the back of the H2 tank w.r.t. the nose 
    X_tank_back = ...
    X_tank_front = X_tank_back - L_tank

    if space_in_front_of_tank > 0:
        push_tank_fwd
        
    V_cargo_aft = Beechcraft_1900D['V_cargo_aft'] + volume_from_removing_seats - volume_tank_with_space_underneath_and_above
    M_cargo_aft = V_cargo_aft * rho_cargo

    return X_tank, M_cargo_aft, X_cargo_aft, num_PAX