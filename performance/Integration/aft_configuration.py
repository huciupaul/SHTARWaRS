import numpy as np
import math as m
from typing import Tuple

from DSE_1.global_constants import Beechcraft_1900D, rho_cargo, seat_pitch

def determine_position_in_cone(L_tank: float, d_tank: float) -> Tuple[float, float, float]:
    a = (Beechcraft_1900D['X_aft_cone_beg'] - Beechcraft_1900D['X_aft_cone_end']) / (Beechcraft_1900D['d_aft_cone_beg'] - Beechcraft_1900D['d_aft_cone_end'])
    b = Beechcraft_1900D['X_aft_cone_beg'] - a * Beechcraft_1900D['d_aft_cone_beg']
    
    X_tank_back = a * d_tank + b
    X_tank_front = X_tank_back - L_tank
    X_aft_cyl_beg = Beechcraft_1900D['X_aft_cyl_beg']

    return X_tank_front, X_tank_back, X_aft_cyl_beg

def determine_diameter_in_cone(X_position: float) -> float:
    a = (Beechcraft_1900D['X_aft_cone_beg'] - Beechcraft_1900D['X_aft_cone_end']) / (Beechcraft_1900D['d_aft_cone_beg'] - Beechcraft_1900D['d_aft_cone_end'])
    b = Beechcraft_1900D['X_aft_cone_beg'] - a * Beechcraft_1900D['d_aft_cone_beg']
    
    d =  (X_position - b) / a

    return d

def check_for_seating_interference(X_tank_front: float) -> Tuple[float, int]:
    if Beechcraft_1900D['X_aft_cyl_beg'] - X_tank_front > 0:
        num_rows_removed = m.ceil((Beechcraft_1900D['X_aft_cyl_beg'] - X_tank_front) / Beechcraft_1900D['seat_pitch'])
        num_PAX = Beechcraft_1900D['num_PAX'] - 2 * num_rows_removed
        X_aft_cyl_beg = Beechcraft_1900D['X_aft_cyl_beg'] - num_rows_removed * seat_pitch
    else:
        num_PAX = Beechcraft_1900D['num_PAX']
        X_aft_cyl_beg = Beechcraft_1900D['X_aft_cyl_beg']

    return X_aft_cyl_beg, num_PAX

def shift_tank_fwd(X_tank_front: float, X_tank_back: float, X_aft_cyl_beg: float, L_tank: float, d_tank_TMS: float) -> Tuple[float, float]:
    if X_aft_cyl_beg - X_tank_front < 0:
        X_tank_front = X_aft_cyl_beg
        X_tank_back = X_tank_front + L_tank
        d_tank_TMS = determine_diameter_in_cone(X_tank_back)        

    return X_tank_front, X_tank_back, d_tank_TMS

def tank_and_TMS_positioning_specs(L_tank: float, d_tank_TMS: float) -> Tuple[float, float]:
    L_tank_cyl = Beechcraft_1900D['l_aft_cyl'] if Beechcraft_1900D['l_aft_cyl'] >= L_tank else L_tank
    L_tank_cone = L_tank - L_tank_cyl

    # Cylindical part
    V_tank_TMS_cyl = L_tank_cyl * Beechcraft_1900D['h_aft_cyl_ave']/2 * Beechcraft_1900D['w_aft_cyl']/2 * m.pi
    X_tank_TMS_cyl = Beechcraft_1900D['X_aft_cyl_beg'] + L_tank_cyl / 2
    
    # Conical part
    V_tank_TMS_cone = 1/3 * m.pi * L_tank_cone * (Beechcraft_1900D['d_aft_cone_beg']**2 + d_tank_TMS**2 + Beechcraft_1900D['d_aft_cone_beg'] * d_tank_TMS)/4
    X_tank_TMS_cone = L_tank_cone/4 * (Beechcraft_1900D['X_aft_cone_beg']**2 + 2 * Beechcraft_1900D['X_aft_cone_beg'] * d_tank_TMS + 3 * d_tank_TMS**2) / (Beechcraft_1900D['X_aft_cone_beg']**2 + Beechcraft_1900D['X_aft_cone_beg'] * d_tank_TMS + d_tank_TMS**2)

    # Total tank and TMS volume and cg position
    V_tank_TMS = V_tank_TMS_cyl + V_tank_TMS_cone
    X_tank_TMS = (V_tank_TMS_cyl * X_tank_TMS_cyl + V_tank_TMS_cone * X_tank_TMS_cone) / V_tank_TMS

    return L_tank_cone, d_tank_TMS, V_tank_TMS, X_tank_TMS

def calculate_cargo_specs(L_tank_cone: float, d_cargo_front: float) -> Tuple[float, float, float]:
    L_aft_cargo = Beechcraft_1900D['l_aft_cone'] - L_tank_cone

    V_aft_cargo = 1/3 * m.pi * L_aft_cargo * (d_cargo_front**2 + Beechcraft_1900D['d_aft_cone_end']**2 + d_cargo_front * Beechcraft_1900D['d_aft_cone_end'])/4
    M_aft_cargo = rho_cargo * V_aft_cargo
    X_aft_cargo = L_aft_cargo/4 * (d_cargo_front**2  + 2 * d_cargo_front * Beechcraft_1900D['X_aft_cone_end'] + 3 * Beechcraft_1900D['X_aft_cone_end']**2) / (d_cargo_front**2 + d_cargo_front * Beechcraft_1900D['X_aft_cone_end'] + Beechcraft_1900D['X_aft_cone_end']**2)

    return V_aft_cargo, M_aft_cargo, X_aft_cargo




def main(L_tank: float, d_tank: float):
    # Place the tank as aft as possible
    X_tank_front, X_tank_back, X_aft_cyl_beg = determine_position_in_cone(L_tank, d_tank)

    # Check if the tank interferes with the seating space, remove seats if necessary
    X_aft_cyl_beg, num_PAX = check_for_seating_interference(X_tank_front)
    
    # Move tank forward if possible
    X_tank_front, X_tank_back = shift_tank_fwd(X_tank_front, X_tank_back, X_aft_cyl_beg, L_tank)

    # Caluclate tank+TMS volume and cg position
    L_tank_cone, d_tank_TMS, V_tank_TMS, X_tank_TMS = tank_and_TMS_positioning_specs(L_tank, d_tank)
    
    # Calculate cargo volume, mass and cg position
    V_aft_cargo, M_aft_cargo, X_aft_cargo = calculate_cargo_specs(L_tank_cone, d_tank_TMS)

    aft_integration_data = {
        'X_tank_front': X_tank_front,
        'X_tank_back': X_tank_back,
        'X_tank_TMS': X_tank_TMS,
        'V_tank_TMS': V_tank_TMS,
        
        'X_aft_cargo': X_aft_cargo,
        'M_aft_cargo': M_aft_cargo,
        'V_aft_cargo': V_aft_cargo,

        'num_PAX': num_PAX
    }
    
    return aft_integration_data