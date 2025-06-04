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
    
    # Warn for logs to be overwritten
    # for split in power_splits:
    #     log_file = os.path.join(output_dir, f"log_{split:.5f}.csv")
    #     if os.path.exists(log_file) and not args.overwrite:
    #         warnings.warn(f"Log file {log_file} already exists. Use --overwrite to overwrite it.")
    #     elif os.path.exists(log_file) and args.overwrite:
    #         warnings.warn(f"Overwriting log file {log_file}.")
    #     else:
    #         pass

    for split in power_splits:
        #m_eps, Qdot_eps = eps_main(split, TOGA)

        #P_cc, airspeed, T_cc = fpp_simple(TOGA, split)
        #mdot_h2_cc, T_cc, pressure_cc = Turboprop.function(P_cc)
        for fc_toga_percentage in fc_toga_percentages:

            #m_fc, m_tms, c_d_tms, delta_ap = fuelcell_main(split, TOGA, fc_toga_percentage, mdot_h2_cc, T_cc, pressure_cc)

            #MTOW_mod = MTOW_mod_init
            #MTOW_mod = m_eps + m_fc + m_tms + m_h2 + m_sto
            for fc_cruise_percentage in range(fc_toga_percentage):
                #m_h2 = fpp(split, fc_toga_percentage, fc_cruise_percentage, MTOW_mod, c_d_tms)
                #m_sto, vol_sto, l_sto = sto_main(m_h2)

                




    
    
    

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