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


def main(dt: float = 0.01) -> Tuple[List[float], np.ndarray]:
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

    mdot_h2o = np.arange(0, 0.2 + 1e-2, 1e-2)  # kg/s
    NOx_limit = 2 * 1e-9 # mole fraction
    time_limit = 10.0

    #----------Iteration----------
    for m_dot_h2o_i in tqdm(mdot_h2o):
        NOx_emitted = 0.0
        t = 0.0

        # TODO: Liquid water injected
        lh2o = gas_obj(
                mdot=m_dot_h2o_i,
                X="H2O:1",
                T=333.0,
                P=rich_ex.P
                )

        X_mix, T_mix, mdot_mix = mix_streams_const_HP(
            stream1=stoich_mix,
            stream2=lh2o,
            P=rich_ex.P
        )

        stoich_with_water = gas_obj(
            mdot=mdot_mix,
            X=X_mix,
            T=T_mix,
            P=rich_ex.P
        )
        print(f"m_dot_h2o_i = {m_dot_h2o_i:.3f} kg/s, T = {stoich_with_water.T:.2f} K, P = {stoich_with_water.P/ct.one_atm:.2f} atm, X = {stoich_with_water.X}")

        while NOx_emitted < NOx_limit and t < time_limit:
            try:
                # Set up 0D reactor
                reactor = ct.IdealGasConstPressureReactor(stoich_with_water.gas)
                sim = ct.ReactorNet([reactor])
                sim.advance(t)
                NOx_emitted += reactor.thermo['NO'].X[0]
                t += dt
            except Exception as e:
                print("Error at m_dot_h2o_i =", m_dot_h2o_i)
                print("T =", stoich_with_water.gas.T)
                print("P =", stoich_with_water.gas.P)
                print("X =", stoich_with_water.gas.X)
                raise

        times_to_reach_limit_NOx.append(t)
    times_to_reach_limit_NOx = np.array(times_to_reach_limit_NOx)

    return times_to_reach_limit_NOx, mdot_h2o

def plot(mdot_h2o: np.ndarray,
    times: np.ndarray) -> None:

    plt.figure(figsize=(8, 5))
    plt.plot(mdot_h2o, times, marker='o')
    plt.xlabel('Injected H2O Mass Flow Rate [kg/s]')
    plt.ylabel('Time to Reach NOx Limit [s]')
    plt.title('Effect of H2O Injection on Time to Reach NOx Limit')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    times, mdot_h2o = main(dt=0.001)
    plot(mdot_h2o, times)