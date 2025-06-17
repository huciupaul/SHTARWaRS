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
import pickle
from scipy.interpolate import RegularGridInterpolator, InterpolatedUnivariateSpline

from MECH import MECH
from mixture_properties import mixture_properties_with_temp
from global_constants import Power, T_rich, P_rich

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


def main(
        T_rich_ex: np.ndarray,
        mdot_h2o: np.ndarray,
        water_liquid: bool = True
        ) -> Tuple[List[float], np.ndarray]:

    NOx_normalized = np.full((len(T_rich_ex), len(mdot_h2o)), np.nan)
    NOx_at_autoignition = np.full((len(T_rich_ex), len(mdot_h2o)), np.nan)
    P_matrix = np.zeros((len(T_rich_ex), len(mdot_h2o)))
    time_max_dT_dlogt = np.zeros((len(T_rich_ex), len(mdot_h2o)))

    P_min, P_max, T_min, T_max = 482017, 1465306, 1579, 1730

    for i, T_rich_ex_i in tqdm(enumerate(T_rich_ex)):
        P_rich_ex = (P_max - P_min) / (T_max - T_min) * T_rich_ex_i + (P_min - (P_max - P_min) / (T_max - T_min) * T_min)

        # Exhaust of the rich zone
        rich_ex = gas_obj(
                mdot=0.23071218802678672,
                X="N2:0.27396853378378067, H2:0.5784925158354742, H:0.0005133907896591396, O2:7.393652373124469e-10, O:4.591482831748451e-09, H2O:0.14701811528905126, OH:7.430992453322465e-06, H2O2:1.1742925935728515e-10, HO2:3.6285427485688574e-12, NO:5.4364866743951985e-09, N2O:9.901157119223157e-12, N:4.953608429224339e-12, NH:1.4251848100599631e-11, NNH:1.977256219409285e-09, NH2:4.148258601540983e-10",
                T= T_rich_ex_i,
                P= P_rich_ex
                )

        # Air bled, to be mixed with rich exhaust to create stoichiometric mixture
        air = gas_obj(
                mdot=1.0,
                X="N2:0.78084, O2:0.20946",
                T=603.0,
                P=rich_ex.P
                )

        mdot_h2 = rich_ex.mdot * rich_ex.gas.X[1]
        mdot_o2 = 0.5 * mdot_h2  # stoichiometric ratio for H2 combustion, assuming no O2 in exhaust
        air.mdot = mdot_o2 / air.gas.X[3]

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
        
        for j, m_dot_h2o_j in enumerate(mdot_h2o):
            # Liquid water injection
            h2o = gas_obj(
                    mdot=m_dot_h2o_j,
                    X="H2O:1",
                    T=353.0,
                    P=rich_ex.P
                    )
            
            if water_liquid:
                vapor_h = h2o.gas.enthalpy_mass
                del_h = 2580 * 1e3 # J/kg, simple estimation
                h2o.HP = (vapor_h - del_h), h2o.gas.P
            

            mdot_mix, T_mix, P_mix, X_mix = mixture_properties_with_temp(
                stream1=stoich_mix,
                stream2=h2o,
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

            Temp = []
            X_nox = []
            time = []

            V = stoich_with_water.mdot / (2*np.pi*0.15*0.02 * stoich_with_water.gas.density)
            tau_max = 0.1/V
            # print(V, tau_max)

            for t in np.arange(1e-9, tau_max, 1e-6):
                try:
                    # sim.advance_to_steady_state()
                    sim.advance(t)
                    Temp.append(reactor.thermo.T)
                    X_nox.append(reactor.thermo['NO'].X[0])
                    time.append(t)
                    # times_to_steady[i, j] = sim.time

                except Exception as e:
                    print("Error at m_dot_h2o_i =", m_dot_h2o_j)
                    print("T =", stoich_with_water.gas.T)
                    print("P =", stoich_with_water.gas.P)
                    print("X =", stoich_with_water.gas.X)
                
                    if reactor.thermo.T > 20000 or reactor.thermo.T < 200:
                        print(f"Unphysical T = {reactor.thermo.T}, skipping this iteration.")
                        print(f"Reactor state before failure:\nT = {reactor.thermo.T}, Y = {reactor.thermo.Y}")
                        print(f"Sum(Y) = {sum(reactor.thermo.Y)}")
                        raise

            log_time = np.log10(time)
            dT_dlogt = np.diff(Temp) / np.diff(log_time)
            time_mid = 10**((log_time[:-1] + log_time[1:]) / 2)
            idx_max = np.argmax(dT_dlogt)

            # plt.figure()
            # plt.plot(time_mid, dT_dlogt)
            # plt.xscale('log')
            # plt.xlabel('Time [s] (log scale)')
            # plt.ylabel('dT/d(log t)')
            # plt.title('Derivative of Temperature vs Log(Time)')
            # plt.show()
            
            # plt.figure()
            # plt.plot(np.arange(1e-9, tau_max, 1e-6), Temp)
            # plt.xscale('log')
            # plt.show()            

            # X_NOx[i, j] = reactor.thermo['NO'].X[0] # molar fraction # * stoich_with_water.mdot # kg/s
            # X_O2[i, j] = reactor.thermo['O2'].X[0]
            # X_H2O[i, j] = reactor.thermo['H2O'].X[0]

            P_matrix[i, j] = reactor.thermo.P
            time_max_dT_dlogt[i, j] = time_mid[idx_max]

            O2_dry = reactor.thermo['O2'].X[0] / (1.0 - reactor.thermo['H2O'].X[0])
            F15 = (20.9 - 15.0) / (20.9 - 100.0 * O2_dry)
            
            NOx_normalized[i, j] = reactor.thermo['NO'].X[0] * F15 #* 1e6
            NOx_at_autoignition[i, j] = X_nox[idx_max] * F15 #* 1e6

    return NOx_normalized, NOx_at_autoignition, P_matrix, time_max_dT_dlogt


if __name__ == "__main__":
    T_rich_ex = np.arange(1570, 1740, 1)
    mdot_h2o = np.arange(0, 0.1 + 1e-2, 1e-2)

    NOx_normalized, NOx_at_autoignition, P_matrix, time_max_dT_dlogt = main(T_rich_ex=T_rich_ex, mdot_h2o=mdot_h2o, water_liquid=True)

    T_grid, W_grid = np.meshgrid(T_rich_ex, mdot_h2o, indexing='ij')


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(W_grid, T_grid, NOx_normalized, cmap='viridis')

    # Add 25 ppm contour (2.5e-5)
    contour = ax.contour(
        W_grid, T_grid, NOx_normalized,
        levels=[2.5e-5],
        colors='red',
        linewidths=2,
        offset=2.5e-5  # Place contour at the correct z-level
    )
    ax.clabel(contour, fmt={2.5e-5: '25 ppm'}, colors='red')

    ax.set_xlabel('Water injection [kg/s]')
    ax.set_ylabel('Rich zone exhaust temperature [K]')
    ax.set_zlabel('NOx formed in the quench zone [ppm]')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    plt.show()


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(W_grid, T_grid, NOx_at_autoignition, cmap='viridis')

    # Add 25 ppm contour (2.5e-5)
    contour = ax.contour(
        W_grid, T_grid, NOx_at_autoignition,
        levels=[2.5e-5],
        colors='red',
        linewidths=2,
        offset=2.5e-5  # Place contour at the correct z-level
    )
    ax.clabel(contour, fmt={2.5e-5: '25 ppm'}, colors='red')

    ax.set_xlabel('Water injection [kg/s]')
    ax.set_ylabel('Rich zone exhaust temperature [K]')
    ax.set_zlabel('NOx formed till the autoignition time [ppm]')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    plt.show()



    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    time_log = np.log10(time_max_dT_dlogt)
    surf = ax.plot_surface(W_grid, T_grid, time_log, cmap='viridis')

    # Choose tick positions (in log10 space)
    zticks = np.arange(np.floor(time_log.min()), np.ceil(time_log.max())+1)
    ax.set_zticks(zticks)
    ax.set_zticklabels([f"{10**z:.1e}" for z in zticks])

    ax.set_xlabel('Water injection [kg/s]')
    ax.set_ylabel('Rich zone exhaust temperature [K]')
    ax.set_zlabel('Autoignition time [s]')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    plt.show()


    NOx_interpolator = RegularGridInterpolator(
        (T_rich_ex, mdot_h2o), NOx_normalized, bounds_error=False, fill_value=None
    )

    os.makedirs('data/interpolants', exist_ok=True)

    # Pickleeeeeeeeeeeee the interpolator
    with open('data/interpolants/NOx_interpolator.pkl', 'wb') as f:
        pickle.dump(NOx_interpolator, f, protocol=pickle.HIGHEST_PROTOCOL)


    # fig, ax = plt.subplots(figsize=(6, 5))

    # # Time to autoignite heatmap
    # im = ax.imshow(
    #     time_max_dT_dlogt,
    #     aspect='auto',
    #     origin='lower',
    #     extent=[mdot_h2o[0], mdot_h2o[-1], T_rich_ex[0], T_rich_ex[-1]],
    #     cmap='plasma'
    # )
    # ax.set_xlabel('mdot_h2o [kg/s]')
    # ax.set_ylabel('T_rich_ex [K]')
    # ax.set_title('Temperature [K]')
    # fig.colorbar(im, ax=ax, label='NOx [molar ppm]')

    # # Overlay pressure contours
    # levels = 5  # or a list of pressure values
    # contours = ax.contour(
    #     mdot_h2o, T_rich_ex, P_matrix,
    #     levels=levels, colors='white', linewidths=1
    # )
    # ax.clabel(contours, inline=True, fontsize=8, fmt='%.0f Pa')

    # Add contour for 25 ppm (2.5e-5 mole fraction)
    # contour = ax.contour(
    #     mdot_h2o, T_rich_ex, NOx_at_15_O2,
    #     levels=[2.5e-5], colors='red', linewidths=2
    # )
    # ax.clabel(contour, fmt={2.5e-5: '25 ppm'}, colors='red')


    # plt.tight_layout()
    # plt.show()
    




'''

# Looping over powers, not inlet temps
def main(
        Power: np.ndarray,
        mdot_h2o: np.ndarray,
        water_liquid: bool = True
        ) -> Tuple[List[float], np.ndarray]:

    NOx_normalized = np.full((len(T_rich_ex), len(mdot_h2o)), np.nan)
    NOx_at_autoignition = np.full((len(T_rich_ex), len(mdot_h2o)), np.nan)
    P_matrix = np.zeros((len(T_rich_ex), len(mdot_h2o)))
    time_max_dT_dlogt = np.zeros((len(T_rich_ex), len(mdot_h2o)))

    temp_spline = InterpolatedUnivariateSpline(Power, T_rich, k=3)
    pres_spline = InterpolatedUnivariateSpline(Power, P_rich, k=3)

    for i, Power_i in tqdm(enumerate(Power)):
        # Exhaust of the rich zone
        rich_ex = gas_obj(
                mdot=0.23071218802678672,
                X="N2:0.27396853378378067, H2:0.5784925158354742, H:0.0005133907896591396, O2:7.393652373124469e-10, O:4.591482831748451e-09, H2O:0.14701811528905126, OH:7.430992453322465e-06, H2O2:1.1742925935728515e-10, HO2:3.6285427485688574e-12, NO:5.4364866743951985e-09, N2O:9.901157119223157e-12, N:4.953608429224339e-12, NH:1.4251848100599631e-11, NNH:1.977256219409285e-09, NH2:4.148258601540983e-10",
                T=temp_spline(Power_i),
                P=pres_spline(Power_i)
                )

        # Air bled, to be mixed with rich exhaust to create stoichiometric mixture
        air = gas_obj(
                mdot=1.0,
                X="N2:0.78084, O2:0.20946",
                T=603.0,
                P=rich_ex.P
                )

        mdot_h2 = rich_ex.mdot * rich_ex.gas.X[1]
        mdot_o2 = 0.5 * mdot_h2  # stoichiometric ratio for H2 combustion, assuming no O2 in exhaust
        air.mdot = mdot_o2 / air.gas.X[3]

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
        
        for j, m_dot_h2o_j in enumerate(mdot_h2o):
            # Liquid water injection
            h2o = gas_obj(
                    mdot=m_dot_h2o_j,
                    X="H2O:1",
                    T=353.0,
                    P=rich_ex.P
                    )
            
            if water_liquid:
                vapor_h = h2o.gas.enthalpy_mass
                del_h = 2580 * 1e3 # J/kg, simple estimation
                h2o.HP = (vapor_h - del_h), h2o.gas.P
            

            mdot_mix, T_mix, P_mix, X_mix = mixture_properties_with_temp(
                stream1=stoich_mix,
                stream2=h2o,
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

            Temp = []
            X_nox = []
            time = []

            V = stoich_with_water.mdot / (2*np.pi*0.15*0.02 * stoich_with_water.gas.density)
            tau_max = 0.1/V
            # print(V, tau_max)

            for t in np.arange(1e-9, tau_max, 1e-6):
                try:
                    # sim.advance_to_steady_state()
                    sim.advance(t)
                    Temp.append(reactor.thermo.T)
                    X_nox.append(reactor.thermo['NO'].X[0])
                    time.append(t)
                    # times_to_steady[i, j] = sim.time

                except Exception as e:
                    print("Error at m_dot_h2o_i =", m_dot_h2o_j)
                    print("T =", stoich_with_water.gas.T)
                    print("P =", stoich_with_water.gas.P)
                    print("X =", stoich_with_water.gas.X)
                
                    if reactor.thermo.T > 20000 or reactor.thermo.T < 200:
                        print(f"Unphysical T = {reactor.thermo.T}, skipping this iteration.")
                        print(f"Reactor state before failure:\nT = {reactor.thermo.T}, Y = {reactor.thermo.Y}")
                        print(f"Sum(Y) = {sum(reactor.thermo.Y)}")
                        raise

            log_time = np.log10(time)
            dT_dlogt = np.diff(Temp) / np.diff(log_time)
            time_mid = 10**((log_time[:-1] + log_time[1:]) / 2)
            idx_max = np.argmax(dT_dlogt)

            # plt.figure()
            # plt.plot(time_mid, dT_dlogt)
            # plt.xscale('log')
            # plt.xlabel('Time [s] (log scale)')
            # plt.ylabel('dT/d(log t)')
            # plt.title('Derivative of Temperature vs Log(Time)')
            # plt.show()
            
            # plt.figure()
            # plt.plot(np.arange(1e-9, tau_max, 1e-6), Temp)
            # plt.xscale('log')
            # plt.show()            

            # X_NOx[i, j] = reactor.thermo['NO'].X[0] # molar fraction # * stoich_with_water.mdot # kg/s
            # X_O2[i, j] = reactor.thermo['O2'].X[0]
            # X_H2O[i, j] = reactor.thermo['H2O'].X[0]

            P_matrix[i, j] = reactor.thermo.P
            time_max_dT_dlogt[i, j] = time_mid[idx_max]

            O2_dry = reactor.thermo['O2'].X[0] / (1.0 - reactor.thermo['H2O'].X[0])
            F15 = (20.9 - 15.0) / (20.9 - 100.0 * O2_dry)
            
            NOx_normalized[i, j] = reactor.thermo['NO'].X[0] * F15 * 1e6
            NOx_at_autoignition[i, j] = X_nox[idx_max] * F15 * 1e6

    return NOx_normalized, NOx_at_autoignition, P_matrix, time_max_dT_dlogt

    '''