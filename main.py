# file: main.py
# desc: This script is the main entry point for the SHTARWaRS design tool.

# External imports
import os
import csv
import argparse
import warnings
import matplotlib.pyplot as plt
import numpy as np
from global_constants import * # Import global constants
import pickle

# Local imports
from eps.eps_sizing import eps_main  # Import EPS sizing function
from fpp.flight import fpp_main  # Import FPP sizing function
from storage.tank import main_storage  # Import storage sizing function



def main(minimum, maximum, no_of_splits, max_iter):
    # Check if the output directory exists, if not create it
    output_dir = "data/logs"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Define the power splits and fuel cell percentages
    power_splits = np.linspace(minimum, maximum, no_of_splits)  # 10 splits from 0 to 1
    fc_toga_percentages = np.linspace(minimum, maximum, no_of_splits)  # Example percentages of TOGA power for the fuel cell
    fc_cruise_percentages = np.linspace(minimum, maximum, no_of_splits)  # Example percentages of CRUISE power for the fuel cell

    # Initialize the result tensor with zeros
    result_tensor = np.zeros((
    len(power_splits),
    len(fc_toga_percentages),
    len(fc_cruise_percentages),
    14  # Length of the tensor
    ))
    
    loading_tensor = np.zeros((
        len(power_splits),
        len(fc_toga_percentages),
        len(fc_cruise_percentages),
        4, 4
    ))

    for i_split, split in enumerate(power_splits):
        
        # EPS
        m_eps, Qdot_eps, V_elmo = eps_main(split)

        for i_toga, fc_toga_percentage in enumerate(fc_toga_percentages):

            for i_cruise, fc_cruise_percentage in enumerate(fc_cruise_percentages):

                # Initialize variables for the first iteration of convergence
                MTOW = Beechcraft_1900D['OEW'] # Original OEW
                delta_ap, c_d_rad = 0.0, 0.0
                m_eps_prev, m_fc_tms_prev, m_h2_prev, m_sto_prev, m_cargo_prev = 0.0, 326.63, 315.39, 141.38, 714.05 # from midterm
                MTOW += (m_eps + m_fc_tms_prev + m_h2_prev + m_sto_prev)  # Initial MTOW
                MTOW_prev = 0.0

                # Convergence of MTOW, stop if the change is less than 1% of the previous MTOW or reached max number of iterations
                for i in range(max_iter):
                    if np.abs((MTOW-MTOW_prev)/MTOW_prev) < 0.01 and i > 0:
                        print(f"Converged after {i} iterations with MTOW: {MTOW:.2f} kg")
                        print(f"Split: {split:.2f}, TOGA: {fc_toga_percentage:.2f}, Cruise: {fc_cruise_percentage:.2f}")
                        break
                    
                    # FPP
                    TMS_inputs, m_h2, FC_outputs, mission_profile, loading_vector = fpp_main(split, fc_toga_percentage, fc_cruise_percentage, MTOW, c_d_rad, delta_ap, 10)
                    m_fc = FC_outputs['m_fc']
                    V_fc = FC_outputs['V_fc']
                    
                    # Storage
                    m_sto, V_sto, t1, dv, t2, length_sto, diameter_sto = main_storage(m_h2)
                    
                    # Update delta_AP, c_D_rad based on the TMS code and get:
                    m_tms_front = 0
                    m_tms_aft = 0
                    m_tms_mid = 0
                    m_cargo = 0

                    # --------- ADD INTEGRATION HERE ---------
                    
                    #MTOW update
                    MTOW_prev = MTOW
                    MTOW += ((m_eps - m_eps_prev) + (m_h2 - m_h2_prev) + (m_sto - m_sto_prev) + \
                        (m_tms_front + m_tms_aft + m_tms_mid + m_fc - m_fc_tms_prev))
                    # Update the values based on the previous iteration
                    m_eps_prev = m_eps
                    m_fc_tms_prev = m_tms_aft + m_tms_front + m_tms_mid + m_fc
                    m_h2_prev = m_h2
                    m_sto_prev = m_sto
                    m_cargo_prev = m_cargo

                # After convergence, store the results in the tensor, then in the 4D tensor
                tensor = np.array([m_eps, m_fc, m_h2, m_sto,V_fc, V_sto, V_elmo, MTOW, length_sto, diameter_sto, m_cargo, m_tms_front, m_tms_aft, m_tms_mid])
                result_tensor[i_split, i_toga, i_cruise, :] = tensor
                loading_tensor[i_split, i_toga, i_cruise, :, :] = loading_vector

    # Save the result tensor to a pickle file
    with open("data/logs/result_tensor.pkl", "wb") as f:
        pickle.dump(result_tensor, f)
    
    with open("data/logs/loading_tensor.pkl", "wb") as f:
        pickle.dump(loading_tensor, f)
    
    print("Result tensor saved to data/logs/result_tensor.pkl")
    print("Loading tensor saved to data/logs/loading_tensor.pkl")

    return None       
    

#main()

# with open("data/logs/result_tensor.pkl", "rb") as f:
#     loaded_tensor = pickle.load(f)

# print(loaded_tensor[0,0,1,:])
# print(loaded_tensor.shape)  # Print the shape of the tensor for verification


if __name__=="__main__":
    main(
        minimum=0.1,
        maximum=1.0,
        no_of_splits=20,
        max_iter=100
    )