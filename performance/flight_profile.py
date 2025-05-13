import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple


def f_eta_prop(eta_prop: float, V_inf: float, rho_inf: float, P: float, A_prop: float) -> float:
    """
    Propulsive efficiency equation across a range of flight conditions.
    :param eta_prop: Propulsive efficiency [-]
    :param V_inf: Aircraft airspeed [m/s]
    :param rho_inf: Air density [kg/m^3]
    :param P: Propulsive power [W]
    :param A_prop: Total propeller disk area (x N_prop) [m^2]
    :return: f(eta_prop) - eta_prop [-]
    """
    return 2/(1 + (1/V_inf)*np.sqrt(V_inf**2 + 2*P*eta_prop/(rho_inf*A_prop*V_inf))) - eta_prop
    
def get_P_cr(V_cruise: float, rho_inf: float, S: float, C_D0: float, e: float, AR: float, W: np.ndarray, A_prop: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the power required to maintain cruise speed.
    :param V_cruise: Cruise speed [m/s]
    :param rho_inf: Air density [kg/m^3]
    :param S: Wing area [m^2]
    :param C_D0: Zero-lift drag coefficient [-]
    :param e: Oswald efficiency factor [-]
    :param AR: Aspect ratio [-]
    :param W: Aircraft weight [N]
    :param A_prop: Total propeller disk area (x N_prop) [m^2]
    :return: P -> Power required [W], eta_prop -> Propulsive efficiency [-]
    """
    # Constant altitude: Lift = Weight
    C_L = W/(0.5*rho_inf*V_cruise**2*S)
    # Induced drag coefficient
    C_Di = C_L**2/(np.pi*e*AR)
    
    # Total drag coefficient
    C_D = C_D0 + C_Di
    
    # Constant velocity: Thrust = Drag
    T = 0.5*rho_inf*V_cruise**2*S*C_D
    
    # Propulsive efficiency
    eta_prop = 2/(1 + (1/V_cruise)*np.sqrt(V_cruise**2 + 2*T/(rho_inf*A_prop)))
    
    # Power required
    P = T*V_cruise/eta_prop

    return P, eta_prop

def get_T_02(T_0: np.ndarray, V_0: np.ndarray, c_pa: float) -> np.ndarray:
    """
    Calculate the temperature at the inlet of the compressor.
    :param T_0: Ambient temperature [K]
    :param V_0: Inlet velocity [m/s]
    :param c_pa: Specific heat of air at constant pressure [J/(kg*K)]
    :return: T_02 -> Temperature at the inlet of the compressor [K]
    """
    return T_0 + (V_0**2)/(2*c_pa)

def get_T_03(T_02: np.ndarray, PI_comp: float, eta_iso: float, k: float) -> np.ndarray:
    """
    Calculate the temperature at the exit of the compressor.
    :param T_02: Temperature at the inlet of the compressor [K]
    :param PI_comp: Pressure ratio across the compressor [-]
    :param eta_iso: Isentropic efficiency of the compressor [-]
    :param k: Specific heat ratio [-]
    :return: T_03 -> Temperature at the exit of the compressor [K]
    """
    return T_02 + (T_02/eta_iso)*(PI_comp**((k-1)/k) - 1)

def get_mdot_f(mdot_a: np.ndarray, T_04: np.ndarray, T_03: np.ndarray, c_pg: float, LHV_f: float, eta_cc: float) -> np.ndarray:
    """
    Calculate the fuel mass flow rate.
    :param mdot_a: Air mass flow rate [kg/s]
    :param T_04: Turbine exit temperature [K]
    :param T_03: Compressor exit temperature [K]
    :param c_pg: Specific heat of gas at constant pressure [J/(kg*K)]
    :param LHV_f: Lower heating value of fuel [J/kg]
    :param eta_cc: Compressor efficiency [-]
    :return: mdot_f -> Fuel mass flow rate [kg/s]
    """
    return mdot_a*c_pg*(T_04 - T_03)/(LHV_f*eta_cc)

def ISA(T_0: float, P_0: float, h: float) -> Tuple[float, float]:
    """
    Calculate the temperature and pressure at a given altitude using the International Standard Atmosphere (ISA) model.
    :param T_0: Sea level temperature [K]
    :param P_0: Sea level pressure [Pa]
    :param h: Altitude [m]
    :return: T -> Temperature at altitude [K], P -> Pressure at altitude [Pa], rho -> Density at altitude [kg/m^3]
    """
    # Constants
    R = 287.05  # Specific gas constant for air [J/(kg*K)]
    g = 9.81  # Gravitational acceleration [m/s^2]
    
    # Temperature lapse rate
    L = 0.0065  # Temperature lapse rate [K/m]
    
    # Calculate temperature and pressure at altitude
    T = T_0 - L*h
    P = P_0 * (T/T_0)**(g/(L*R))
    rho = P/(R*T)  # Density at altitude [kg/m^3]
    
    return T, P, rho

if __name__=="__main__":
    # Input params
    T_0 = 288.15  # Sea level temperature [K]
    P_0 = 101325  # Sea level pressure [Pa]
    Pw_TOGA = 1.908e6 # TOGA power [W]
    S = 28.79 # Wing area [m^2]
    C_D0 = 0.0215 # Zero-lift drag coefficient [-]
    b = 17.64 # Wing span [m]
    AR = b**2/S # Aspect ratio [-]
    e = 1.78*(1-0.045*AR**0.68) - 0.64 # Oswald efficiency factor [-]
    D_prop = 2.78 # Propeller diameter [m]
    A_prop = 2*(np.pi*(0.5*D_prop)**2) # Total propeller disk area (x N_prop) [m^2]
    
    h = np.array([
        0,      # Takeoff 
        2438.4, # Climb station 0 
        4876.8, # Climb station 1 
        7620.0, # Cruise
        304.8,  # Approach
        0       # Landing
    ])
    V = np.array([
        54.0167, # Takeoff
        139.929, # Climb station 0
        146.102, # Climb station 1
        142.501, # Cruise
        60.2,    # Approach
        60.2     # Landing
    ])
    ROC = np.array([
        797/60,  # Takeoff to climb station 0
        797/60,  # Climb station 0 to climb station 1
        797/60,  # Climb station 1 to cruise
        -6.4,    # Cruise to approach
        -3.1506  # Approach to landing
    ])
    Pw = np.array([
        1.0, # Takeoff (5 mins)
        0.95, # Rest of climb
        0.19  # Descent
    ])*Pw_TOGA # Power setting
    
    # Find time taken to get to each station
    t = np.diff(h)/ROC
    dt = 1 # Time step [s]
    
    # Extrapolate velocity profile across climbs to stations (start to end split across time)
    V_segments = [
        np.linspace(V[i], V[i+1], int(np.ceil(abs((h[i+1]-h[i])/(ROC[i] or 1))*1))) 
        for i in range(len(h)-1)
    ]

    alt_segments = [
        np.linspace(h[i], h[i+1], int(np.ceil(abs((h[i+1]-h[i])/(ROC[i] or 1))*1))) 
        for i in range(len(h)-1)
    ]
    
    #––––––––––––––––––––––––––––––––––––––––––––––––
    # 1) separate climb (legs 0–2) vs descent (legs 3–4)
    climb_legs = 3
    V_climb = np.hstack(V_segments[:climb_legs])
    alt_climb = np.hstack(alt_segments[:climb_legs])
    ROC_climb = np.hstack([ROC[i]*np.ones_like(V_segments[i]) for i in range(climb_legs)])
    gamma_climb = np.arcsin(ROC_climb / V_climb)
    V_ground_climb = V_climb * np.cos(gamma_climb)
    range_climb = np.cumsum(V_ground_climb) * dt
    total_climb_range = range_climb[-1]

    V_descent = np.hstack(V_segments[climb_legs:])
    alt_descent = np.hstack(alt_segments[climb_legs:])
    ROC_descent = np.hstack([ROC[i]*np.ones_like(V_segments[i]) for i in range(climb_legs, len(V_segments))])
    gamma_descent = np.arcsin(ROC_descent / V_descent)
    V_ground_descent = V_descent * np.cos(gamma_descent)
    range_descent_local = np.cumsum(V_ground_descent) * dt
    total_descent_range = range_descent_local[-1]

    #––––––––––––––––––––––––––––––––––––––––––––––––
    # 2) cruise leg to make up 707 km total
    R_tot = 707e3
    R_cruise = R_tot - total_climb_range - total_descent_range
    V_cr = V[3]
    n_cr = int(np.ceil((R_cruise / V_cr) / dt))

    V_cr_arr = V_cr * np.ones(n_cr)
    ROC_cr   = np.zeros(n_cr)
    alt_cr   = h[3]  * np.ones(n_cr)
    V_ground_cr = V_cr
    range_cr_local = np.arange(1, n_cr+1) * V_ground_cr * dt

    #––––––––––––––––––––––––––––––––––––––––––––––––
    # 3) assemble full profiles in the correct order
    V_full    = np.hstack([V_climb,    V_cr_arr,    V_descent])
    ROC_full  = np.hstack([ROC_climb,  ROC_cr,      ROC_descent])
    alt_full  = np.hstack([alt_climb,  alt_cr,      alt_descent])
    time_full = np.arange(len(V_full)) * dt
    range_full= np.hstack([
        range_climb,
        total_climb_range + range_cr_local,
        total_climb_range + R_cruise + range_descent_local
    ])

    # #––––––––––––––––––––––––––––––––––––––––––––––––
    # # 4) quick plots
    # plt.figure(figsize=(10,6))
    # plt.subplot(2,2,1)
    # plt.plot(time_full, V_full);       plt.title("Airspeed vs Time");    plt.xlabel("t [s]"); plt.ylabel("V [m/s]")
    # plt.subplot(2,2,2)
    # plt.plot(time_full, range_full);   plt.title("Ground Range vs Time");plt.xlabel("t [s]"); plt.ylabel("R [m]")
    # plt.subplot(2,2,3)
    # plt.plot(time_full, ROC_full);     plt.title("ROC vs Time");        plt.xlabel("t [s]"); plt.ylabel("ROC [m/s]")
    # plt.subplot(2,2,4)
    # plt.plot(time_full, alt_full);     plt.title("Altitude vs Time");    plt.xlabel("t [s]"); plt.ylabel("h [m]")
    # plt.tight_layout()
    # plt.show()
    
    #––––––––––––––––––––––––––––––––––––––––––––––––
    # 5) calculate power required for cruise
    
    