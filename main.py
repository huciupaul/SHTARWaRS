# file: main.py
# desc: This script is the main entry point for the SHTARWaRS design tool.

# External imports
import os
import csv
import argparse
import warnings
import matplotlib.pyplot as plt
import numpy as np
from global_constants import TOGA, MTOW_orig # Import global constants
import pickle

from eps.eps_sizing import eps_main  # Import EPS sizing function
from fpp.flight import main_fpp  # Import FPP sizing function
from storage.tank import main_storage  # Import storage sizing function

# Internal imports
#...

def main():
    # Check if the output directory exists, if not create it
    output_dir = "data/logs"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Define power split array
    # power_splits = np.linspace(args.min, args.max, args.n_splits+1)
    # fc_toga_percentages = np.linspace(0, 1, 200)  # Percentage of TOGA power for the fuel cell

    power_splits = np.array([0.2])  # Example: 10 splits from 0 to 1
    fc_toga_percentages = np.array([0.2])  # Example percentages of TOGA power for the fuel cell


    for split in power_splits:

        m_eps, Qdot_eps, V_elmo = eps_main(split)

        for fc_toga_percentage in fc_toga_percentages:

            for fc_cruise_percentage in fc_toga_percentages:

                
                MTOW = MTOW_orig  # Original Maximum Take-Off Weight 
                delta_ap, c_d_rad = 0, 0
                m_eps_prev, m_fc_prev, m_h2_prev, m_sto_prev, m_tms_prev = 0, 0, 0, 0, 0

                for i in range(10):
                    
                    if i > 0:
                        # Update the values based on the previous iteration
                        m_eps_prev = m_eps
                        m_fc_prev = m_fc
                        m_h2_prev = m_h2
                        m_sto_prev = m_sto
                        m_tms_prev = m_tms
                        

                    # Perform optimization steps here
                    TMS_inputs, m_h2, FC_outputs, mission_profile = main_fpp(split, fc_toga_percentage, fc_cruise_percentage, MTOW, c_d_rad, delta_ap, 0.1)
                    print("Test", m_h2, FC_outputs, split, fc_toga_percentage, fc_cruise_percentage)
                    #m_sto, V_sto, t1, dv, t2, length_sto, diameter_sto = main_storage(m_h2)

                    # Update delta_AP, c_D_rad based on the TMS code
                    m_tms = 0
                    m_sto = 0
                    m_fc = 0
                    #MTOW update
                    MTOW += (m_eps - m_eps_prev + m_fc - m_fc_prev + m_h2 - m_h2_prev + m_sto - m_sto_prev + m_tms - m_tms_prev)
    return MTOW


print(main())

                
                ### TENSOR SHAPE DESCRIPTION ###
                # The tensor is expected to have the shape (N, M, P, Q) where:
                # N: split
                # M: fc_toga_percentage
                # P: fc_cruise_percentage
                # Q: Design variables: 
                #   < m_eps, m_fc, m_h2, m_sto, m_tms, V_FC, V_sto, V_elmo, MTOW, length_Sto, diameter_sto >

                # tensor = np.array([
                #         m_eps, m_fc, m_h2, m_sto, m_tms,
                #         V_fc, V_sto, V_elmo, MTOW, length_sto, diameter_sto
                #     ])

                # # === Save tensor to a pickle file ===
                # file_name = f"tensor_split{split}_toga{fc_toga_percentage}_cruise{fc_cruise_percentage}.pkl"
                # file_path = os.path.join(output_dir, file_name)

                # with open(file_path, 'wb') as f:
                #     pickle.dump(tensor, f)

                


# if __name__=="__main__":
#     args = argparse.ArgumentParser(description="SHTARWaRS Design Tool")
#     args.add_argument("--min", type=float, default=0.0, help="Minimum power split value")
#     args.add_argument("--max", type=float, default=1.0, help="Maximum power split value")
#     args.add_argument("--n_splits", type=int, default=10, help="Number of power splits")
#     args.add_argument("--max_iter", type=int, default=100, help="Maximum number of iterations for optimization")
#     args.add_argument("--overwrite", action='store_true', help="Overwrite existing results")
#     args.add_argument("--plot", action='store_true', help="Plot the results")
    
#     args = args.parse_args()
    
#     main(args)