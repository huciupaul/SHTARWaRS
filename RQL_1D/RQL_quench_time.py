import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.join(current_dir, '..')
sys.path.insert(0, os.path.abspath(parent_dir))

import cantera as ct
ct.suppress_thermo_warnings()

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from typing import List, Tuple

from MECH import MECH
from mixture_properties import mixture_properties_with_temp

'''
This script is used to estimate the NOx formation in the quench zone of a RQL combustor.
'''

# Input gases
class gas_obj:
    mdot: float # kg/s
    X: str      # Cantera mole fraction string
    T: float    # K
    P: float    # Pa
    def __init__(self, mdot: float, X: str, T: float, P: float):
        self.mdot = mdot
        self.X = X
        self.T = T
        self.P = P

        self.gas = ct.Solution(MECH)
        self.gas.TPX = T, P, X


def main(hdot_h2o_step: float = 1e-3) -> Tuple[List[float], np.ndarray]:
    #----------Inputs----------
    # Exhaust of the rich zone
    rich_ex = gas_obj(
            mdot=0.23071218802678672,
            X="N2:0.27396853378378067, H2:0.5784925158354742, H:0.0005133907896591396, O2:7.393652373124469e-10, O:4.591482831748451e-09, H2O:0.14701811528905126, OH:7.430992453322465e-06, H2O2:1.1742925935728515e-10, HO2:3.6285427485688574e-12, NO:5.4364866743951985e-09, N2O:9.901157119223157e-12, N:4.953608429224339e-12, NH:1.4251848100599631e-11, NNH:1.977256219409285e-09, NH2:4.148258601540983e-10",
            T=734.22,
            P=12.1 * ct.one_atm
            )

    # Air bleeded, to be mixed with rich exhaust to create stoichiometric mixture
    air = gas_obj(
            mdot=1.0,
            X="N2:0.78084, O2:0.20946",
            T=603.0,
            P=rich_ex.P
            )

    n_h2 = rich_ex.mdot * rich_ex.gas.X[1]
    n_o2 = 0.5 * n_h2  # stoichiometric ratio for H2 combustion, assuming no O2 in exhaust
    air.mdot = n_o2 / air.gas.X[3]

    mdot_mix1, T_mix1, P_mix1, X_mix1 = mixture_properties_with_temp(
        stream1=rich_ex,
        stream2=air,
        MECH=MECH
    )

    stoich_mix = gas_obj(
            mdot=mdot_mix1,
            X=X_mix1,
            T=T_mix1,
            P=P_mix1
            )

    #----------Outputs----------
    NOx = []
    water_content = []
    temps_before_reaction = []
    temps_after_reaction = []

    mdot_h2o = np.arange(0, 0.1 + 1e-2, hdot_h2o_step)  # kg/s



    lh2o = gas_obj(
            mdot=0.0,
            X="H2O:1",
            T=353.0,
            P=rich_ex.P
            )
    
    mdot_mix, T_mix, P_mix, X_mix = mixture_properties_with_temp(
        stream1=stoich_mix,
        stream2=lh2o,
        MECH=MECH
    )

    stoich_with_water = gas_obj(
        mdot=mdot_mix,
        X=X_mix,
        T=T_mix,
        P=P_mix
    )

    # Set up the Cantera gas object
    reactor = ct.IdealGasConstPressureReactor(stoich_with_water.gas)
    sim = ct.ReactorNet([reactor])
    temps_before_reaction.append(reactor.thermo.T)

    times = []
    NOx_formed = []
    temp_evolution = []

    O2 = []
    N2 = []
    
    for t in tqdm(np.arange(1e-4, 5.0, 1e-4)):
        sim.advance(t)
        times.append(t)
        NOx_formed.append(reactor.thermo['NO'].X[0] * 1e6)
        # O2.append(reactor.thermo['O2'].X[0])
        # N2.append(reactor.thermo['N2'].X[0])
        temp_evolution.append(reactor.thermo.T)

    print(np.array(NOx_formed))
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    # Ensure water_content matches the length of NOx and temperature arrays

    ax1.plot(np.array(times), np.array(NOx_formed), color='r')
    # ax1.plot(np.array(times), np.array(O2), color='b')
    # ax1.plot(np.array(times), np.array(N2), color='g')
    ax1.set_title('NOx emissions')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('NOx [ppm]')

    ax2.plot(np.array(times), np.array(temp_evolution))
    ax2.set_title('Temp evolution')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Temp [K]')

    plt.tight_layout()
    plt.show()

        

    # #----------Iteration----------
    # for m_dot_h2o_i in tqdm(mdot_h2o):
    #     lh2o = gas_obj(
    #             mdot=m_dot_h2o_i,
    #             X="H2O:1",
    #             T=353.0,
    #             P=rich_ex.P
    #             )
        
    #     mdot_mix, T_mix, P_mix, X_mix = mixture_properties_with_temp(
    #         stream1=stoich_mix,
    #         stream2=lh2o,
    #         MECH=MECH
    #     )

    #     stoich_with_water = gas_obj(
    #         mdot=mdot_mix,
    #         X=X_mix,
    #         T=T_mix,
    #         P=P_mix
    #     )
    #     water_content.append(lh2o.mdot/stoich_with_water.mdot * 100)
        
    #     # Set up the Cantera gas object
    #     reactor = ct.IdealGasConstPressureReactor(stoich_with_water.gas)
    #     sim = ct.ReactorNet([reactor])
    #     temps_before_reaction.append(reactor.thermo.T)

    #     try:
    #         res = sim.advance_to_steady_state(max_steps=1e3, return_residuals=True)
    #         # advance_to_steady_state(self, max_steps=10000, residual_threshold=0.0, atol=0.0, return_residuals=False)
    #         # print(res)
    #         print(f"Temp after reaction: ", reactor.thermo.T)

    #     except Exception as e:
    #         print("Error at m_dot_h2o_i =", m_dot_h2o_i)
    #         print("T =", stoich_with_water.gas.T)
    #         print("P =", stoich_with_water.gas.P)
    #         print("X =", stoich_with_water.gas.X)
        
    #         if reactor.thermo.T > 20000 or reactor.thermo.T < 200:
    #             print(f"Unphysical T = {reactor.thermo.T}, skipping this iteration.")
    #             print(f"Reactor state before failure:\nT = {reactor.thermo.T}, Y = {reactor.thermo.Y}")
    #             print(f"Sum(Y) = {sum(reactor.thermo.Y)}")

                
    #     temps_after_reaction.append(reactor.thermo.T)
    #     NOx.append(reactor.thermo['NO'].X[0] * 1e4) # ppm

    # NOx = np.array(NOx)
    # water_content = np.array(water_content)
    # temps_before_reaction = np.array(temps_before_reaction)
    # temps_after_reaction = np.array(temps_after_reaction)

    # return NOx, water_content, temps_before_reaction, temps_after_reaction


if __name__ == "__main__":
    # NOx, water_content, temps_before_reaction, temps_after_reaction = main(hdot_h2o_step=1e-3)
    main()
    NOx_limit = 25 # ppm

    # print(f'NOx produced: ', NOx)
    # print(f'Temps before: ', temps_before_reaction, 'Temps after: ', temps_after_reaction)

    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))  # adjust figsize as needed

    # # Ensure water_content matches the length of NOx and temperature arrays
    # water_content_plot = water_content[:len(NOx)]

    # ax1.plot(water_content_plot, NOx)
    # ax1.axhline(y=NOx_limit, color='r')
    # ax1.set_title('NOx emissions')
    # ax1.set_xlabel('Water content in quench [%]')
    # ax1.set_ylabel('NOx [ppm]')

    # ax2.plot(water_content, temps_before_reaction, color='b')
    # ax2.plot(water_content_plot, temps_after_reaction, color='r')
    # ax2.set_title('Temperature')
    # ax2.set_xlabel('Water content in quench [%]')
    # ax2.set_ylabel('Temp [K]')

    # plt.tight_layout()
    # plt.show()