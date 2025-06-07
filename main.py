# file: main.py
# desc: This script is the main entry point for the SHTARWaRS design tool.

# External imports
import os
import csv
import argparse
import warnings
import matplotlib.pyplot as plt
import numpy as np
from global_constants import TOGA, MTOW_mod # Import global constants

# Internal imports
#...

def main(args):
    # Check if the output directory exists, if not create it
    output_dir = "data/logs"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Define power split array
    power_splits = np.linspace(args.min, args.max, args.n_splits+1)
    fc_toga_percentages = np.linspace(0, 1, 200)  # Percentage of TOGA power for the fuel cell

    for split in power_splits:

        #m_eps, Qdot_eps, V_elmo = eps_main(split, TOGA)

        for fc_toga_percentage in fc_toga_percentages:

            for fc_cruise_percentage in range(fc_toga_percentage):

                
                #MTOW = MTOW_mod_init 
                #MTOW_mod = m_eps + m_fc + m_tms + m_h2 + m_sto
                #delta_ap, c_d_rad = 0, 0

                for i in range(args.max_iter):
                    
                    # Perform optimization steps here
                    #TMS_inputs, m_h2, m_fc, mission_profile = main_fpp(split, fc_toga_percentage, fc_cruise_percentage, MTOW, c_d_rad, delta_ap, dt: float=0.1)
                    #m_sto, V_sto, t1, dv, t2, length_sto, radius_sto = main_storage(m_h2)

                    # Update MTOW_mod, delta_AP, c_D_rad based on the optimization logic
                    #MTOW = OEW - () + (m_eps + m_fc + m_tms + m_h2 + m_sto)

                




    
    
    

if __name__=="__main__":
    args = argparse.ArgumentParser(description="SHTARWaRS Design Tool")
    args.add_argument("--min", type=float, default=0.0, help="Minimum power split value")
    args.add_argument("--max", type=float, default=1.0, help="Maximum power split value")
    args.add_argument("--n_splits", type=int, default=10, help="Number of power splits")
    args.add_argument("--max_iter", type=int, default=100, help="Maximum number of iterations for optimization")
    args.add_argument("--overwrite", action='store_true', help="Overwrite existing results")
    args.add_argument("--plot", action='store_true', help="Plot the results")
    
    args = args.parse_args()
    
    main(args)