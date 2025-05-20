import numpy as np
import matplotlib.pyplot as plt
import os

'''Reference Data'''
OEW =              # kg
X_OEW =            # m

# Cargo
X_c =             # m
W_c =               # kg

# Passengers
W_s = 84            # kg
num_seats = 19               # seats
X_most_forward_seat = 4.8375 # m
seat_spacing = 0.7525  # m

# Fuel
W_f = 368.3           # kg
X_f =        # m


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
        """Finds the original array and the index within that array."""
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
plt.scatter(X_seat_front, W_seat_front, color='red', label="Window Passengers Front to Back")
plt.scatter(X_seat_back, W_seat_back, color='green', label="Window Passengers Back to Front")
plt.scatter(X_seat_aisle, W_seat_aisle, color='orange', label="Aisle Passengers Front to Back")
plt.scatter(X_seat_back_aisle, W_seat_back_aisle, color='purple', label="Aisle Passengers Back to Front")
plt.scatter(X_fuel, W_fuel_value, color='black', label="Fuel")

plt.plot(X_seat, W_seat, color='red')
plt.plot(X_seat_back, W_seat_back, color='green')
plt.plot(X_seat_aisle, W_seat_aisle, color='orange')
plt.plot(X_seat_back_aisle, W_seat_back_aisle, color='purple')
plt.plot(X_fuel, W_fuel_value, color='black')

plt.plot([X_cargo[0], X_cargo[1]], [OEW, W_OEW_fc], 'b-')  # OEW to Front Cargo
plt.plot([X_cargo[1], X_cargo[3]], [W_OEW_fc, W_OEW_c], 'b-')  # Front Cargo to Both Cargo
plt.plot([X_cargo[0], X_cargo[2]], [OEW, W_OEW_rc], 'b-')  # OEW to Rear Cargo
plt.plot([X_cargo[2], X_cargo[3]], [W_OEW_rc, W_OEW_c], 'b-')  # Rear Cargo to Both Cargo

plt.xlabel("X_cg/MAC [-]")
plt.ylabel("Weight [kg]")
plt.title("Aircraft CG Position vs. Weight")
plt.legend()
plt.grid(True)


# Loading diagram with comparison
X_cargo_new, W_cargo_new, XCG_CW_new = cargo(XCG_fc_new, W_fc_new, XCG_rc_new, W_rc_new, X_OEW_with_Batt_new, OEW_with_Batt_new, X_LEMAC_new, MAC_new)
W_OEW_c_new, W_OEW_fc_new, W_OEW_rc_new = W_cargo_new[3], W_cargo_new[1], W_cargo_new[2]
X_seat_new, W_seat_new, X_seat_back_new, W_seat_back_new, X_seat_aisle_new, W_seat_aisle_new, X_seat_back_aisle_new, W_seat_back_aisle_new = passengers(XCG_CW_new, W_OEW_c_new, W_seat_value_new, num_seats_new, X_most_forward_seat_new, X_most_aft_seat_new, X_LEMAC_new, MAC_new)
X_fuel_new, W_fuel_value_new = fuel(X_seat_aisle_new, W_seat_aisle_new, XCG_fuel_new, W_fuel_new, X_LEMAC_new, MAC_new)
min_cg_new, max_cg_new = min_max_X_cg_positions(X_cargo_new, X_seat_new, X_seat_back_new, X_seat_aisle_new, X_seat_back_aisle_new, X_fuel_new)


print(f'OEW CG: ', X_OEW_with_Batt_new, convert_to_MAC_reference_inverse(X_OEW_with_Batt_new, X_LEMAC_new, MAC_new))
print(f'Front cargo CG: ', X_cargo_new[1], convert_to_MAC_reference_inverse(X_cargo_new[1], X_LEMAC, MAC))
print(f'Rear cargo CG: ', X_cargo_new[2], convert_to_MAC_reference_inverse(X_cargo_new[2], X_LEMAC, MAC))
print(f'Fuel CG: ', X_fuel_new[-1], convert_to_MAC_reference_inverse(X_fuel_new[-1], X_LEMAC_new, MAC_new))

plt.figure(figsize=(8, 6))
# Reference
plt.scatter(X_cargo, W_cargo, color='gray', label="Original Aircraft")
plt.scatter(X_seat, W_seat, color='gray')
plt.scatter(X_seat_back, W_seat_back, color='gray')
plt.scatter(X_seat_aisle, W_seat_aisle, color='gray')
plt.scatter(X_seat_back_aisle, W_seat_back_aisle, color='gray')
plt.scatter(X_fuel, W_fuel_value, color='gray')

plt.plot(X_seat, W_seat, color='gray')
plt.plot(X_seat_back, W_seat_back, color='gray')
plt.plot(X_seat_aisle, W_seat_aisle, color='gray')
plt.plot(X_seat_back_aisle, W_seat_back_aisle, color='gray')
plt.plot(X_fuel, W_fuel_value, color='gray')

plt.plot([X_cargo[0], X_cargo[1]], [OEW, W_OEW_fc], 'gray')  # OEW to Front Cargo
plt.plot([X_cargo[1], X_cargo[3]], [W_OEW_fc, W_OEW_c], 'gray')  # Front Cargo to Both Cargo
plt.plot([X_cargo[0], X_cargo[2]], [OEW, W_OEW_rc], 'gray')  # OEW to Rear Cargo
plt.plot([X_cargo[2], X_cargo[3]], [W_OEW_rc, W_OEW_c], 'gray')  # Rear Cargo to Both Cargo

# Modified
plt.scatter(X_cargo_new, W_cargo_new, color='blue', label="Cargo")
plt.scatter(X_seat_new, W_seat_new, color='red', label="Window Passengers Front to Back")
plt.scatter(X_seat_back_new, W_seat_back_new, color='green', label="Window Passengers Back to Front")
plt.scatter(X_seat_aisle_new, W_seat_aisle_new, color='orange', label="Aisle Passengers Front to Back")
plt.scatter(X_seat_back_aisle_new, W_seat_back_aisle_new, color='purple', label="Aisle Passengers Back to Front")
plt.scatter(X_fuel_new, W_fuel_value_new, color='black', label="Fuel")

plt.plot(X_seat_new, W_seat_new, color='red')
plt.plot(X_seat_back_new, W_seat_back_new, color='green')
plt.plot(X_seat_aisle_new, W_seat_aisle_new, color='orange')
plt.plot(X_seat_back_aisle_new, W_seat_back_aisle_new, color='purple')
plt.plot(X_fuel_new, W_fuel_value_new, color='black')

plt.plot([X_cargo_new[0], X_cargo_new[1]], [OEW_with_Batt_new, W_OEW_fc_new], 'b-')  # OEW to Front Cargo
plt.plot([X_cargo_new[1], X_cargo_new[3]], [W_OEW_fc_new, W_OEW_c_new], 'b-')  # Front Cargo to Both Cargo
plt.plot([X_cargo_new[0], X_cargo_new[2]], [OEW_with_Batt_new, W_OEW_rc_new], 'b-')  # OEW to Rear Cargo
plt.plot([X_cargo_new[2], X_cargo_new[3]], [W_OEW_rc_new, W_OEW_c_new], 'b-')  # Rear Cargo to Both Cargo

plt.xlabel("X_cg/MAC [-]")
plt.ylabel("Weight [kg]")
plt.title("Aircraft CG Position vs. Weight")
plt.legend()
plt.grid(True)

# save_path = os.path.join('Output', f"Potato_diagram.pdf")
# plt.savefig(save_path, format='pdf')

plt.show()