import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import matplotlib.pyplot as plt

from DSE_1.global_constants import M_PAX, X_cargo_fwd, V_cargo_fwd, V_cargo_aft, V_cargo, X_first_seat,  seat_pitch


@dataclass
class Aircraft_data:
    MTOW: float
    X_MTOW: float
    M_fuel: float
    X_fuel: float
    X_cargo_fwd: float
    X_cargo_aft: float
    num_PAX: float
    OEW: float
    M_EPS: float
    X_EPS: float

    M_FC: float=None
    X_FC: float=None
    M_storage: float=None
    X_storage: float=None
    M_TMS: float=None
    X_TMS: float=None


    @property
    def X_PAX(self) -> float:
        return ((self.num_PAX-3)*M_PAX * (X_first_seat + (self.num_PAX-5)/2 * seat_pitch / 2) + 3*M_PAX * (X_first_seat + (self.num_PAX-3)/2 * seat_pitch)) / (self.num_PAX*M_PAX)

    @property
    def M_cargo(self, V_cargo_specific) -> float:
        return 939 * V_cargo_specific/V_cargo
    
    @property
    def X_OEW(self) -> float:
        M_cargo_tot = self.M_cargo(V_cargo_fwd) + self.M_cargo(V_cargo_aft)
        M_payload = self.num_PAX*M_PAX + M_cargo_tot

        X_cargo_tot = (self.M_cargo(V_cargo_fwd)*self.X_cargo_fwd + self.M_cargo(V_cargo_aft)*self.X_cargo_aft) / M_cargo_tot
        X_payload = (self.num_PAX*M_PAX*self.X_PAX + M_cargo_tot*X_cargo_tot) / M_payload

        return (self.MTOW*self.X_MTOW - self.M_fuel*self.X_fuel - M_payload*X_payload) / self.OEW



class CG_calculation:
    def __init__(self,
                 original_aircraft: Original_aircraft,
                 retrofitted_aircraft: Retrofitted_aircraft
                 ):
        self.og = original_aircraft
        self.retro = retrofitted_aircraft
        

def main():
    Beechcraft_1900D = Aircraft_data(
        MTOW=7766
        X_MTOW=7.51
        M_fuel=337
        X_fuel=7.29
        X_cargo_fwd=X_cargo_fwd
        X_cargo_aft=12.64
        num_PAX=19
        X_PAX=7.74
        OEW=4894
    )

    H2D2 = Aircraft_data(

    )






















"""
    '''Function Definitions'''
    def cargo(def_X_c, def_W_c, def_X_OEW, def_OEW):

        W_cargo = def_OEW + def_W_c
        X_cargo = ((def_X_OEW * def_OEW) + (def_X_c * def_W_c)) / W_cargo

        return X_cargo, W_cargo


    def passengers(def_X_cargo, def_W_cargo, def_W_seat, def_num_seats, def_X_most_forward_seat, def_seat_spacing):

        X_start = def_X_cargo
        W_start = def_W_cargo
        num_rows = (def_num_seats - 1) // 2

        # Front to back
        X_seat_front = np.zeros(num_rows)
        W_seat_front = np.zeros(num_rows)
        X_seat_front[0] = X_start
        W_seat_front[0] = W_start

        for i in range(1, num_rows + 1):
            W_seat_front[i] = W_seat_front[i - 1] + (2 * def_W_seat)
            X_seat_front[i] = (((def_X_most_forward_seat + (i * def_seat_spacing)) * 2 * def_W_seat) + 
                        (X_seat_front[i - 1] * W_seat_front[i - 1])) / W_seat_front[i]

        # Back to front
        X_seat_back = np.zeros(num_rows)
        W_seat_back = np.zeros(num_rows)
        X_seat_back[0] = X_start
        W_seat_back[0] = W_start

        for i in range(1, num_rows + 1):
            W_seat_back[i] = W_seat_back[i - 1] + (2 * def_W_seat)
            X_seat_back[i] = (((def_X_most_forward_seat - (i * seat_spacing)) * 2 * def_W_seat) + 
                            (X_seat_back[i - 1] * W_seat_back[i - 1])) / W_seat_back[i]

        # returning: front to back, back to front
        return X_seat_front, W_seat_front, X_seat_back, W_seat_back


    def fuel(def_X_seat_back, def_W_seat_back, def_X_fuel, def_W_fuel):
        
        X_start = def_X_seat_back[-1]
        W_start = def_W_seat_back[-1]

        X_fuel = np.zeros(2)
        W_fuel = np.zeros(2)
        X_fuel[0] = X_start
        W_fuel[0] = W_start

        W_fuel[1] = W_start + def_W_fuel
        X_fuel[1] = ((def_X_fuel * def_W_fuel) + (X_start * W_start)) / W_fuel[1]

        return X_fuel, W_fuel


    def min_max_X_cg_positions(def_X_cargo, def_X_seat_front, def_X_seat_back, def_X_fuel):
        # Combine all x_cg positions
        arrays = [def_X_cargo, def_X_seat_front, def_X_seat_back, def_X_fuel]
        array_names = ["X_cargo", "X_seat_front", "X_seat_back", "X_fuel"]

        min_cg = np.min(np.hstack([def_X_cargo, def_X_seat_front, def_X_seat_back, def_X_fuel]))
        max_cg = np.max(np.hstack([def_X_cargo, def_X_seat_front, def_X_seat_back, def_X_fuel]))

        min_index = np.argmin(np.hstack([def_X_cargo, def_X_seat_front, def_X_seat_back, def_X_fuel]))
        max_index = np.argmax(np.hstack([def_X_cargo, def_X_seat_front, def_X_seat_back, def_X_fuel]))

        cumulative_lengths = np.cumsum([len(arr) for arr in arrays])

        def find_source_array_and_index(flat_index):
            '''Finds the original array and the index within that array.'''
            for i, length in enumerate(cumulative_lengths):
                if flat_index < length:
                    original_index = flat_index - (cumulative_lengths[i - 1] if i > 0 else 0)
                    return array_names[i], original_index
            return None, None

        # Get source array names and indices
        min_source, min_array_index = find_source_array_and_index(min_index)
        max_source, max_array_index = find_source_array_and_index(max_index)

        print(f"Min value with 2% margin: {min_cg*0.98}, found in {min_source} at index {min_array_index}")
        print(f"Max value with 2% margin: {max_cg*1.02}, found in {max_source} at index {max_array_index}")

        return min_cg, max_cg


'''Plotting'''
# Loading diagram of the reference aircraft
X_cargo, W_cargo = cargo(X_c, W_c, X_OEW, OEW)
X_seat_front, W_seat_front, X_seat_back, W_seat_back = passengers(X_cargo, W_cargo, W_s, num_seats, X_most_forward_seat, seat_spacing)
X_fuel, W_fuel = fuel(X_seat_back, W_seat_back, X_f, W_f)
min_cg, max_cg = min_max_X_cg_positions(X_cargo, X_seat_front, X_seat_back, X_fuel)

print(f'OEW CG: ', X_OEW)
print(f'Cargo CG: ', X_cargo)
print(f'Fuel CG: ', X_fuel[-1])

plt.figure(figsize=(8, 6))
plt.scatter(X_cargo, W_cargo, color='blue', label="Cargo")
plt.scatter(X_seat_front, W_seat_front, color='red', label="Passengers Front to Back")
plt.scatter(X_seat_back, W_seat_back, color='green', label="Passengers Back to Front")
plt.scatter(X_fuel, W_fuel, color='black', label="Fuel")

plt.plot(X_seat_front, W_seat_front, color='red')
plt.plot(X_seat_back, W_seat_back, color='green')
plt.plot(X_fuel, W_fuel, color='black')
plt.plot([X_OEW, X_cargo], [OEW, W_cargo], 'blue')

plt.xlabel("X_cg [m]")
plt.ylabel("Weight [kg]")
plt.title("Aircraft CG Position vs. Weight")
plt.legend()
plt.grid(True)
plt.show()
"""