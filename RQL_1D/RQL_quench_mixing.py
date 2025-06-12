import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.join(current_dir, '..')
sys.path.insert(0, os.path.abspath(parent_dir))

import cantera as ct
ct.suppress_thermo_warnings()

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple

from MECH import MECH
from mixture_properties import mix_streams_const_HP

'''
This script is used to estimate the maximum mixing time in the quench zone of a RQL combustor 
to keep the NOx emissions under the specified limit.
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


def main(time_limit: float = 10.0, dt: float = 0.01, hdot_h2o_step: float = 1e-3) -> Tuple[List[float], np.ndarray]:
    #----------Inputs----------
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

    X_mix1, T_mix1, mdot_mix1 = mix_streams_const_HP(
        stream1=rich_ex,
        stream2=air,
        P=rich_ex.P
        )

    stoich_mix = gas_obj(
            mdot=mdot_mix1,
            X=X_mix1,
            T=T_mix1,
            P=rich_ex.P
            )

    #----------Outputs----------
    times_to_reach_limit_NOx = []
    NOx = []
    water_content = []

    mdot_h2o = np.arange(0, 0.2 + 1e-2, hdot_h2o_step)  # kg/s
    NOx_limit = 2.5 * 1e-5 # mole fraction equivalent to 25 ppm

    for m_dot_h2o_i in tqdm(mdot_h2o):
        lh2o = gas_obj(
                mdot=m_dot_h2o_i,
                X="H2O:1",
                T=333.0,
                P=rich_ex.P
                )
        
        # # Liquid water properties
        # water = ct.Water()
        # # Set to saturated liquid at 12.1 atm
        # water.set_sat_p(lh2o.P, x=0.0)
        # h_liq = water.enthalpy_mass
        # # Set to saturated vapor at 12.1 atm
        # water.set_sat_p(lh2o.P, x=1.0)
        # h_vap = water.enthalpy_mass
        # h_fg = h_vap - h_liq
        # lh2o.gas.enthalpy_mass -= h_fg
        
        X_mix, T_mix, mdot_mix = mix_streams_const_HP(
            stream1=stoich_mix,
            stream2=lh2o,
            P=rich_ex.P
        )

        water_content.append(lh2o.mdot/mdot_mix)

        stoich_with_water = gas_obj(
            mdot=mdot_mix,
            X=X_mix,
            T=T_mix,
            P=rich_ex.P
        )

        # Create the reactor
        reactor = ct.IdealGasConstPressureReactor(stoich_with_water.gas)
        sim = ct.ReactorNet([reactor])

        # Advance until steady state (change in NOx < tol or max time reached)
        NO_index = stoich_with_water.gas.species_index('NO')
        NOx_prev = 0.0
        t = 0.0
        tol = 1e-9

        while t < time_limit:
            sim.advance(t + dt)
            NOx_now = reactor.thermo.X[NO_index]
            if abs(NOx_now - NOx_prev) < tol or NOx_now >= NOx_limit:
                break
            NOx_prev = NOx_now
            t += dt

        NOx.append(reactor.thermo.X[NO_index]) # This is the NOx produced (mole fraction)
        times_to_reach_limit_NOx.append(t)
    return np.array(times_to_reach_limit_NOx), np.array(mdot_h2o), np.array(NOx), np.array(water_content)


if __name__ == "__main__":
    times, mdot_h2o, NOx, water_content = main(time_limit=15.0, dt=1e-20, hdot_h2o_step=1e-2)
    print("Times to reach NOx limit:", times)
    print("NOx emissions:", NOx)

    # plt.figure()
    # plt.plot(water_content, NOx)
    # plt.show()








        # #----------Iteration----------
    # for m_dot_h2o_i in tqdm(mdot_h2o):
    #     NOx_emitted = 0.0
    #     t = 0.0

    #     lh2o = gas_obj(
    #             mdot=m_dot_h2o_i,
    #             X="H2O:1",
    #             T=333.0,
    #             P=rich_ex.P
    #             )
        
    #     # # Liquid water properties
    #     # water = ct.Water()
    #     # # Set to saturated liquid at 12.1 atm
    #     # water.set_sat_p(lh2o.P, x=0.0)
    #     # h_liq = water.enthalpy_mass
    #     # # Set to saturated vapor at 12.1 atm
    #     # water.set_sat_p(lh2o.P, x=1.0)
    #     # h_vap = water.enthalpy_mass
    #     # h_fg = h_vap - h_liq
    #     # lh2o.gas.enthalpy_mass -= h_fg
        
    #     X_mix, T_mix, mdot_mix = mix_streams_const_HP(
    #         stream1=stoich_mix,
    #         stream2=lh2o,
    #         P=rich_ex.P
    #     )

    #     water_content.append(lh2o.mdot/mdot_mix)

    #     stoich_with_water = gas_obj(
    #         mdot=mdot_mix,
    #         X=X_mix,
    #         T=T_mix,
    #         P=rich_ex.P
    #     )

    #     l = 0.1
    #     d = 0.001

    #     reactor = ct.NewReactor("IdealGasConstPressureReactor",
    #                          gas=stoich_with_water.gas)
    #     reactor.volume = np.pi * (d / 2) ** 2 * l
    #     reactor.initialize()
    #     NOx_prev = 0.0  # mole fraction of NO in the initial gas mixture


    #     # # Set up the Cantera gas object
    #     # reactor = ct.IdealGasConstPressureReactor(stoich_with_water.gas)
    #     # l = 0.1
    #     # d = 0.001
    #     # reactor.volume = np.pi * (d / 2) ** 2 * l
    #     # sim = ct.ReactorNet([reactor])
    #     # sim.initialize()

    #     while t < time_limit and reactor.thermo['NO'].X[0] < NOx_limit and reactor.thermo['NO'].X[0] != NOx_prev:
    #         try:
    #             # sim.advance(t + dt)

    #             # n_NO = reactor.thermo.X[NO_index] * reactor.thermo.density_mole * reactor.volume
    #             # NOx_emitted += max(n_NO - n_NO_prev, 0.0)
    #             # n_NO_prev = n_NO
    #             reactor.advance(t + dt)
    #             NOx_prev = reactor.thermo['NO'].X[0]  # mole fraction of NO in the gas mixture
    #             t += dt
    #         except Exception as e:
    #             print("Error at m_dot_h2o_i =", m_dot_h2o_i)
    #             print("T =", stoich_with_water.gas.T)
    #             print("P =", stoich_with_water.gas.P)
    #             print("X =", stoich_with_water.gas.X)
    #             raise
        
    #     times_to_reach_limit_NOx.append(t)
    #     NOx.append(reactor.thermo['NO'].X[0])  # mole fraction of NO in the final gas mixture
    #     # NOx.append(NOx_emitted)

    # times_to_reach_limit_NOx = np.array(times_to_reach_limit_NOx)
    # NOx = np.array(NOx)
    # water_content = np.array(water_content)

    # return times_to_reach_limit_NOx, mdot_h2o, NOx, water_content