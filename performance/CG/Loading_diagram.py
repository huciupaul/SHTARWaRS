import numpy as np
import matplotlib.pyplot as plt
import os

'''Reference Data'''
OEW = 13600             # kg
X_OEW = 12.25           # m
X_LEMAC = 11.242        # m
MAC = 2.303             # m

# Front cargo
XCG_fc = 4.29           # m
XCG_rc = 21.49          # m

# Rear cargo
W_fc = 560 * (4.6/9.4)
W_rc = 560 * (4.8/9.4)

# Passengers
W_seat_value = 95       # kg
num_seats = 72          # seats
X_most_forward_seat = 6.71 # m
X_most_aft_seat = 18.71 # m

# Fuel
W_fuel = 2000           # kg
XCG_fuel = 12.068       # m

'''Modified Aircraft Data'''
OEW_with_Batt_new = 14219.1       # kg
X_OEW_with_Batt_new = 12.71 # 13.56       # m
X_LEMAC_new = 11.242 # 12.580    # m
MAC_new = 2.303        # m

# Front cargo
XCG_fc_new = 4.29       # m
W_fc_new = 436 * (4.6/9.4) # kg

# Rear cargo
XCG_rc_new = 21.49      # m
W_rc_new = 436 * (4.8/9.4) # kg

# Passengers
W_seat_value_new = 95   # kg
num_seats_new = 56      # seats
# seats spread evenly (change if needed)
X_most_forward_seat_new = 6.71 # m
X_most_aft_seat_new = 16.71 # m

# Fuel
W_fuel_new = 3025       # kg
XCG_fuel_new = 12.068      # m


'''Function Definitions'''
def convert_to_MAC_reference(def_X, def_X_LEMAC, def_MAC):
    return (def_X - def_X_LEMAC)/def_MAC


def convert_to_MAC_reference_inverse(def_X, def_X_LEMAC, def_MAC):
    return def_X_LEMAC + (def_X * def_MAC)


def cargo(def_XCG_fc, def_W_fc, def_XCG_rc, def_W_rc, def_X_OEW, def_OEW, def_X_LEMAC, def_MAC):

    # Only front
    XCG_FW = ((def_X_OEW * def_OEW) + (def_XCG_fc * def_W_fc)) / (def_OEW + def_W_fc)
    W_OEW_fc = def_OEW + def_W_fc

    # Only rear
    XCG_RW = ((def_X_OEW * def_OEW) + (def_XCG_rc * def_W_rc)) / (def_OEW + def_W_rc)
    W_OEW_rc = def_OEW + def_W_rc

    # Combined
    XCG_CW = ((def_X_OEW * def_OEW) + (def_XCG_fc * def_W_fc) + (def_XCG_rc * def_W_rc)) / (def_OEW + def_W_fc + def_W_rc)
    W_OEW_c = def_OEW + def_W_fc + def_W_rc

    X_cargo = [def_X_OEW, XCG_FW, XCG_RW, XCG_CW]
    for i in range(len(X_cargo)):
        X_cargo[i] = convert_to_MAC_reference(X_cargo[i], def_X_LEMAC, def_MAC)
    W_cargo = [def_OEW, W_OEW_fc, W_OEW_rc, W_OEW_c]
    
    return X_cargo, W_cargo, XCG_CW


def passengers(def_XCG_CW, def_W_OEW_c, def_W_seat_value, def_num_seats, def_X_most_forward_seat, def_X_most_aft_seat, def_X_LEMAC, def_MAC):

    # Window seats
    X_start = def_XCG_CW 
    W_start = def_W_OEW_c
    num_rows = (def_num_seats // 4) + 1
    seat_spacing = (def_X_most_aft_seat - def_X_most_forward_seat) / (num_rows) # m

    # Front to back window seats
    X_seat = np.zeros(num_rows)
    W_seat = np.zeros(num_rows)
    X_seat[0] = X_start
    W_seat[0] = W_start
    for i in range(1, num_rows):
        W_seat[i] = W_seat[i - 1] + (2 * def_W_seat_value)
        X_seat[i] = (((X_most_forward_seat + (i * seat_spacing)) * 2 * def_W_seat_value) + 
                    (X_seat[i - 1] * W_seat[i - 1])) / (2 * def_W_seat_value + W_seat[i - 1])

    for i in range(len(X_seat)):
        X_seat[i] = convert_to_MAC_reference(X_seat[i], def_X_LEMAC, def_MAC)

    # Back to front window seats
    X_seat_back = np.zeros(num_rows)
    W_seat_back = np.zeros(num_rows)
    X_seat_back[0] = X_start
    W_seat_back[0] = W_start
    for i in range(1, num_rows):
        W_seat_back[i] = W_seat_back[i - 1] + (2 * def_W_seat_value)
        X_seat_back[i] = (((def_X_most_aft_seat - (i * seat_spacing)) * 2 * def_W_seat_value) + 
                        (X_seat_back[i - 1] * W_seat_back[i - 1])) / (2 * def_W_seat_value + W_seat_back[i - 1])

    for i in range(len(X_seat_back)):
        X_seat_back[i] = convert_to_MAC_reference(X_seat_back[i], def_X_LEMAC, def_MAC)

    # Aisle seats
    X_start = convert_to_MAC_reference_inverse(X_seat[-1], def_X_LEMAC, def_MAC)
    W_start = W_seat[-1]

    # Front to back aisle seats
    X_seat_aisle = np.zeros(num_rows)
    W_seat_aisle = np.zeros(num_rows)
    X_seat_aisle[0] = X_start
    W_seat_aisle[0] = W_start
    for i in range(1, num_rows):
        W_seat_aisle[i] = W_seat_aisle[i - 1] + (2 * def_W_seat_value)
        X_seat_aisle[i] = (((def_X_most_forward_seat + (i * seat_spacing)) * 2 * def_W_seat_value) + 
                        (X_seat_aisle[i - 1] * W_seat_aisle[i - 1])) / (2 * def_W_seat_value + W_seat_aisle[i - 1])

    for i in range(len(X_seat_aisle)):
        X_seat_aisle[i] = convert_to_MAC_reference(X_seat_aisle[i], def_X_LEMAC, def_MAC)

    X_start = convert_to_MAC_reference_inverse(X_seat_back[-1], def_X_LEMAC, def_MAC)

    # Back to front aisle seats
    X_seat_back_aisle = np.zeros(num_rows)
    W_seat_back_aisle = np.zeros(num_rows)
    X_seat_back_aisle[0] = X_start
    W_seat_back_aisle[0] = W_start
    for i in range(1, num_rows):
        W_seat_back_aisle[i] = W_seat_back_aisle[i - 1] + (2 * def_W_seat_value)
        X_seat_back_aisle[i] = (((def_X_most_aft_seat - (i * seat_spacing)) * 2 * def_W_seat_value) + 
                                (X_seat_back_aisle[i - 1] * W_seat_back_aisle[i - 1])) / (2 * def_W_seat_value + W_seat_back_aisle[i - 1])

    for i in range(len(X_seat_back)):
        X_seat_back_aisle[i] = convert_to_MAC_reference(X_seat_back_aisle[i], def_X_LEMAC, def_MAC)

    # returning: front to back window, back to front window, front to back aisle, back to front aisle
    return X_seat, W_seat, X_seat_back, W_seat_back, X_seat_aisle, W_seat_aisle, X_seat_back_aisle, W_seat_back_aisle


def fuel(def_X_seat_aisle, def_W_seat_aisle, def_XCG_fuel, def_W_fuel, def_X_LEMAC, def_MAC):
    X_start = convert_to_MAC_reference_inverse(def_X_seat_aisle[-1], def_X_LEMAC, def_MAC)
    W_start = def_W_seat_aisle[-1]

    X_fuel = np.zeros(2)
    W_fuel_value = np.zeros(2)
    X_fuel[0] = X_start
    W_fuel_value[0] = W_start

    W_fuel_value[1] = W_fuel_value[0] + def_W_fuel
    X_fuel[1] = (((def_XCG_fuel * def_W_fuel) + (X_start * W_start)) / (def_W_fuel + W_start))
    for i in range(len(X_fuel)):
        X_fuel[i] = convert_to_MAC_reference(X_fuel[i], def_X_LEMAC, def_MAC)

    return X_fuel, W_fuel_value


def min_max_X_cg_positions(def_X_cargo, def_X_seat, def_X_seat_back, def_X_seat_aisle, def_X_seat_back_aisle, def_X_fuel):
    # Combine all x_cg positions
    arrays = [def_X_cargo, def_X_seat, def_X_seat_back, def_X_seat_aisle, def_X_seat_back_aisle, def_X_fuel]
    array_names = ["X_cargo", "X_seat", "X_seat_back", "X_seat_aisle", "X_seat_back_aisle", "X_fuel"]

    min_cg = np.min(np.hstack([def_X_cargo, def_X_seat, def_X_seat_back, def_X_seat_aisle, def_X_seat_back_aisle, def_X_fuel]))
    max_cg = np.max(np.hstack([def_X_cargo, def_X_seat, def_X_seat_back, def_X_seat_aisle, def_X_seat_back_aisle, def_X_fuel]))

    min_index = np.argmin(np.hstack([def_X_cargo, def_X_seat, def_X_seat_back, def_X_seat_aisle, def_X_seat_back_aisle, def_X_fuel]))
    max_index = np.argmax(np.hstack([def_X_cargo, def_X_seat, def_X_seat_back, def_X_seat_aisle, def_X_seat_back_aisle, def_X_fuel]))

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

    print(f"Min value: {min_cg*0.98}, found in {min_source} at index {min_array_index}")
    print(f"Max value: {max_cg*1.02}, found in {max_source} at index {max_array_index}")

    min_cg_conv = convert_to_MAC_reference_inverse(min_cg*0.98, 11.242, 2.303)
    max_cg_conv = convert_to_MAC_reference_inverse(max_cg*1.02, 11.242, 2.303)
    print(f'min val: ', min_cg_conv)
    print(f'max val: ', max_cg_conv)

    return min_cg, max_cg


'''Plotting'''
# Loading diagram of the reference aircraft
X_cargo, W_cargo, XCG_CW = cargo(XCG_fc, W_fc, XCG_rc, W_rc, X_OEW, OEW, X_LEMAC, MAC)
W_OEW_c, W_OEW_fc, W_OEW_rc = W_cargo[3], W_cargo[1], W_cargo[2]
X_seat, W_seat, X_seat_back, W_seat_back, X_seat_aisle, W_seat_aisle, X_seat_back_aisle, W_seat_back_aisle = passengers(XCG_CW, W_OEW_c, W_seat_value, num_seats, X_most_forward_seat, X_most_aft_seat, X_LEMAC, MAC)
X_fuel, W_fuel_value = fuel(X_seat_aisle, W_seat_aisle, XCG_fuel, W_fuel, X_LEMAC, MAC)
min_cg, max_cg = min_max_X_cg_positions(X_cargo, X_seat, X_seat_back, X_seat_aisle, X_seat_back_aisle, X_fuel)

print(f'OEW CG: ', X_OEW, convert_to_MAC_reference_inverse(X_OEW, X_LEMAC, MAC))
print(f'Front cargo CG: ', X_cargo[1], convert_to_MAC_reference_inverse(X_cargo[1], X_LEMAC, MAC))
print(f'Rear cargo CG: ', X_cargo[2], convert_to_MAC_reference_inverse(X_cargo[2], X_LEMAC, MAC))
print(f'Fuel CG: ', X_fuel[-1], convert_to_MAC_reference_inverse(X_fuel[-1], X_LEMAC, MAC))

plt.figure(figsize=(8, 6))
plt.scatter(X_cargo, W_cargo, color='blue', label="Cargo")
plt.scatter(X_seat, W_seat, color='red', label="Window Passengers Front to Back")
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