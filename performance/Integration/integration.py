import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import matplotlib.pyplot as plt

from SHTARWaRS.global_constants import Beechcraft_1900D, seat_pitch, rho_cargo, M_PAX, \
    X_most_fwd, X_most_aft, X_cargo_fwd, X_first_seat, \
        X_wing_end, V_cargo_fwd, V_wing,\
            X_wing
from SHTARWaRS.performance.Integration import aft_configuration as acfg

def X_PAX(num_PAX) -> float:
    """Calculate the center of gravity position of the passenger seats."""
    return ((num_PAX-3)*M_PAX * (X_first_seat + (num_PAX-5)/2 * seat_pitch / 2) + 3*M_PAX * (X_first_seat + (num_PAX-3)/2 * seat_pitch)) / (num_PAX*M_PAX)

def X_OEW(num_PAX, M_FC, M_TMS_fwd, M_TMS_aft, M_TMS_mid, M_EPS, M_tank, X_tank_TMS, X_tank_front) -> float:
    """Calculate the center of gravity position of the original aircraft's empty weight."""
    """Original aircraft OEW"""
    M_cargo_og = V_cargo_fwd * rho_cargo + Beechcraft_1900D['M_cargo_aft']
    X_cargo_og = (V_cargo_fwd * rho_cargo * X_cargo_fwd + Beechcraft_1900D['M_cargo_aft'] * Beechcraft_1900D['X_cargo_aft']) / M_cargo_og

    M_PAX_og = num_PAX*M_PAX
    X_PAX_og = X_PAX(num_PAX)

    M_payload_og = M_cargo_og + M_PAX_og
    X_payload_og = (M_PAX_og * X_PAX_og + M_cargo_og*X_cargo_og) / M_payload_og

    # CG position of the OEW of Beechcraft 1900D based on its most forward MTOW CG position
    # TODO: rethink this choice
    X_OEW_og = (Beechcraft_1900D['MTOW'] * X_most_aft - Beechcraft_1900D['M_fuel'] * X_wing - M_payload_og * X_payload_og) / Beechcraft_1900D['OEW']

    OEW_H2D2 = Beechcraft_1900D['OEW'] + M_FC + M_EPS + M_TMS_fwd + M_TMS_aft + M_tank + M_TMS_mid
    X_OEW_H2D2 = (Beechcraft_1900D['OEW'] * X_OEW_og + (M_FC + M_EPS + M_TMS_fwd) * X_wing + (M_TMS_aft + M_tank) * X_tank_TMS + M_TMS_mid * (X_wing_end + X_tank_front) / 2) / OEW_H2D2

    return OEW_H2D2, X_OEW_H2D2

def __cargo(X_cargo_fwd, X_cargo_aft, M_cargo_fwd, M_cargo_aft, X_OEW, OEW):
    """Calculate the center of gravity position and weight of OEW+cargo."""
    W_cargo_front = np.zeros(3)
    X_cargo_front = np.zeros(3)
    W_cargo_back = np.zeros(3)
    X_cargo_back = np.zeros(3)

    """Initial point: OEW"""
    W_cargo_front[0] = OEW
    X_cargo_front[0] = X_OEW
    W_cargo_back[0] = W_cargo_front[0]
    X_cargo_back[0] = X_cargo_front[0]

    """Cargo loaded: front or aft first"""
    W_cargo_front[1] = W_cargo_front[0] + M_cargo_fwd
    X_cargo_front[1] = ((X_cargo_fwd * M_cargo_fwd) + (X_OEW * OEW)) / W_cargo_front[1]
    W_cargo_back[1] = OEW + M_cargo_aft
    X_cargo_back[1] = ((X_cargo_aft * M_cargo_aft) + (X_OEW * OEW)) / W_cargo_back[1]

    """Cargo loaded: the rest"""
    W_cargo_front[2] = W_cargo_front[1] + M_cargo_aft
    X_cargo_front[2] = ((X_cargo_aft * M_cargo_aft) + (X_cargo_front[1] * W_cargo_front[1])) / W_cargo_front[2]
    W_cargo_back[2] = W_cargo_front[2]
    X_cargo_back[2] = X_cargo_front[2]

    return X_cargo_front, W_cargo_front, X_cargo_back, W_cargo_back

def __passengers(X_cargo_front, W_cargo_front, num_PAX):
    """Calculate the center of gravity position and weight of OEW+cargo+PAX."""
    X_start = X_cargo_front[-1]
    W_start = W_cargo_front[-1]
    num_2PAX_rows = (num_PAX - 3) // 2
    num_3PAX_rows = 1

    """Front to back"""
    X_seat_front = np.zeros(num_2PAX_rows + num_3PAX_rows)
    W_seat_front = np.zeros(num_2PAX_rows + num_3PAX_rows)
    X_seat_front[0] = X_start
    W_seat_front[0] = W_start

    for i in range(1, num_2PAX_rows):
        W_seat_front[i] = W_seat_front[i - 1] + (2 * M_PAX)
        X_seat_front[i] = (((X_first_seat + (i * seat_pitch)) * 2 * M_PAX) + 
                    (X_seat_front[i - 1] * W_seat_front[i - 1])) / W_seat_front[i]
    
    W_seat_front[-1] = W_seat_front[-2] + (3 * M_PAX)
    X_seat_front[-1] = ((X_first_seat + (num_2PAX_rows * seat_pitch)) * 3 * M_PAX +
                        (X_seat_front[-2] * W_seat_front[-2])) / W_seat_front[-1]

    """Back to front"""
    X_seat_back = np.zeros(num_2PAX_rows + num_3PAX_rows)
    W_seat_back = np.zeros(num_2PAX_rows + num_3PAX_rows)
    X_seat_back[0] = X_start
    W_seat_back[0] = W_start
    X_last_seat = X_first_seat + (num_2PAX_rows * seat_pitch)

    W_seat_back[1] = W_seat_back[0] + (3 * M_PAX)
    X_seat_back[1] = (X_last_seat * 3 * M_PAX +
                    (X_seat_back[0] * W_seat_back[0])) / W_seat_back[1]

    for i in range(2, num_2PAX_rows + 1):
        W_seat_back[i] = W_seat_back[i - 1] + (2 * M_PAX)
        X_seat_back[i] = (((X_first_seat - (i * seat_pitch)) * 2 * M_PAX) + 
                        (X_seat_back[i - 1] * W_seat_back[i - 1])) / W_seat_back[i]

    return X_seat_front, W_seat_front, X_seat_back, W_seat_back

def __fuel(X_seat_front, W_seat_front, X_wing, M_fuel):
    """Calculate the center of gravity position and weight of OEW+cargo+PAX+fuel."""
    X_fuel = np.zeros(2)
    W_fuel = np.zeros(2)
    X_fuel[0] = X_seat_front[-1]
    W_fuel[0] = W_seat_front[-1]

    W_fuel[1] = W_fuel[0] + M_fuel
    X_fuel[1] = ((X_wing * W_fuel) + (X_fuel[0] * W_fuel[0])) / W_fuel[1]

    return X_fuel, W_fuel

def min_max_X_cg_positions(
    X_cargo_aft, M_cargo_aft,
    num_PAX, M_fuel, M_FC,
    M_TMS_fwd, M_TMS_aft,
    M_TMS_mid, M_EPS, M_tank,
    X_tank_TMS, X_tank_front
    ):
    """Find the minimum and maximum center of gravity positions."""
    OEW_H2D2, X_OEW_H2D2 = X_OEW(num_PAX, M_FC, M_TMS_fwd, M_TMS_aft, M_TMS_mid, M_EPS, M_tank, X_tank_TMS, X_tank_front) 
    X_cargo_front, W_cargo_front, X_cargo_back, W_cargo_back = __cargo(Beechcraft_1900D['X_cargo_fwd'], X_cargo_aft, Beechcraft_1900D['M_cargo_fwd'], M_cargo_aft, X_OEW_H2D2, OEW_H2D2)
    X_seat_front, W_seat_front, X_seat_back, W_seat_back = __passengers(X_cargo_front, W_cargo_front, num_PAX)
    X_fuel, W_fuel = __fuel(X_seat_front, W_seat_front, X_wing, M_fuel)

    arrays = [X_cargo_front, X_cargo_back, X_seat_front, X_seat_back, X_fuel]
    array_names = ["X_cargo", "X_seat_front", "X_seat_back", "X_fuel"]

    min_cg = np.min(np.hstack(arrays))
    max_cg = np.max(np.hstack(arrays))

    min_index = np.argmin(np.hstack(arrays))
    max_index = np.argmax(np.hstack(arrays))

    cumulative_lengths = np.cumsum([len(arr) for arr in arrays])

    def __find_source_array_and_index(flat_index):
        """Finds the original array and the index within that array."""
        for i, length in enumerate(cumulative_lengths):
            if flat_index < length:
                original_index = flat_index - (cumulative_lengths[i - 1] if i > 0 else 0)
                return array_names[i], original_index
        return None, None

    # Get source array names and indices
    min_source, min_array_index = __find_source_array_and_index(min_index)
    max_source, max_array_index = __find_source_array_and_index(max_index)

    # Apply 2% margin
    min_cg_with_margin = min_cg*0.98
    max_cg_with_margin = max_cg*1.02

    print(f"Min value with 2% margin: {min_cg*0.98}, found in {min_source} at index {min_array_index}")
    print(f"Max value with 2% margin: {max_cg*1.02}, found in {max_source} at index {max_array_index}")

    return min_cg, max_cg, min_cg_with_margin, max_cg_with_margin

# '''Plotting'''
# def plot_loading_diagram():
#     X_cargo, W_cargo = __cargo(X_c, W_c, X_OEW, OEW)
#     X_seat_front, W_seat_front, X_seat_back, W_seat_back = __passengers(X_cargo, W_cargo, W_s, num_seats, X_most_forward_seat, seat_spacing)
#     X_fuel, W_fuel = __fuel(X_seat_back, W_seat_back, X_f, W_f)
#     min_cg, max_cg, min_cg_with_margin, max_cg_with_margin = __min_max_X_cg_positions(X_cargo, X_seat_front, X_seat_back, X_fuel)

#     print(f'OEW CG: ', X_OEW)
#     print(f'Cargo CG: ', X_cargo)
#     print(f'Fuel CG: ', X_fuel[-1])

#     plt.figure(figsize=(8, 6))
#     plt.scatter(X_cargo, W_cargo, color='blue', label="Cargo")
#     plt.scatter(X_seat_front, W_seat_front, color='red', label="Passengers Front to Back")
#     plt.scatter(X_seat_back, W_seat_back, color='green', label="Passengers Back to Front")
#     plt.scatter(X_fuel, W_fuel, color='black', label="Fuel")

#     plt.plot(X_seat_front, W_seat_front, color='red')
#     plt.plot(X_seat_back, W_seat_back, color='green')
#     plt.plot(X_fuel, W_fuel, color='black')
#     plt.plot([X_OEW, X_cargo], [OEW, W_cargo], 'blue')

#     plt.xlabel("X_cg [m]")
#     plt.ylabel("Weight [kg]")
#     plt.title("Aircraft CG Position vs. Weight")
#     plt.legend()
#     plt.grid(True)
#     plt.show()   

# class Retrofitted_aircraft:
#     MTOW: float
#     X_MTOW: float
#     M_fuel: float
#     X_fuel: float
#     X_cargo_fwd: float
#     X_cargo_aft: float
#     num_PAX: float
#     OEW: float
#     M_EPS: float
#     X_EPS: float

#     M_FC: float
#     X_FC: float
#     M_storage: float
#     X_storage: float
#     M_TMS: float
#     X_TMS: float

# class CG_calculation:
#     def __init__(self,
#                  original_aircraft: Original_aircraft,
#                  retrofitted_aircraft: Retrofitted_aircraft
#                  ):
#         self.og = original_aircraft
#         self.retro = retrofitted_aircraft
        

def main(design: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Main integration function to return the constrained design tensor, array of N_PAX, and array of m_cargo.
    Args:
        design (np.ndarray): Design tensor with shape (N, M, P, Q) where:
            N: Power splits (FC use vs Turboprop use)
            M: FC TOGA throttle setting
            P: FC cruise throttle setting
            Q: Design variables: 
                < m_EPS, m_FC, m_H2, m_storage, m_TMS, V_FC, V_storage, V_ELMO, MTOW, L_storage, D_storage >
    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
        Constrained design tensor, array of N_PAX, and array of m_cargo.    
    """
    # Create np.array of relevant design variables
    var_idx = np.array([0, 1, 2, 3, 4, 8, 9, 10, 11, 12, 13])
    m_EPS, m_FC, m_H2, m_storage, \
        V_FC, L_storage, D_storage, \
            m_cargo, m_TMS_front, m_TMS_aft, m_TMS_mid = design[:, :, :, var_idx]
            
    # Apply volumetric constraints
    # TODO: Revise margin
    margin = 1.05
    valid_design = V_FC*margin/2 <= V_wing  # Fuel cell volume constraint
    
    # Get aircraft configuration
    ac_config = acfg.main(
        L_tank=L_storage,
        d_tank=D_storage
    )
    
    N_PAX = ac_config['num_PAX']
    m_cargo = ac_config['M_aft_cargo']
    X_aft_cargo = ac_config['X_aft_cargo']
    X_tank_front = ac_config['X_tank_front']
    X_tank_TMS = ac_config['X_tank_TMS']
    
    _, _, min_cg_with_margin, max_cg_with_margin = min_max_X_cg_positions(
        X_cargo_aft=X_aft_cargo,
        M_cargo_aft=m_cargo,
        num_PAX=N_PAX,
        M_fuel=m_H2,
        M_FC=m_FC,
        M_TMS_fwd=m_TMS_front,
        M_TMS_aft=m_TMS_aft,
        M_TMS_mid=m_TMS_mid,
        M_EPS=m_EPS,
        M_tank=m_storage,
        X_tank_TMS=X_tank_TMS,
        X_tank_front=X_tank_front
    )
    
    valid_design &= (
        (X_most_fwd <= min_cg_with_margin <= X_most_aft) & 
        (X_most_fwd <= max_cg_with_margin <= X_most_aft)
    )
    
    return valid_design, N_PAX, m_cargo