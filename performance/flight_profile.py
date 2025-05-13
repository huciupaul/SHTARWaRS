import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Dict, List, Optional
from dataclasses import dataclass


def f_eta_prop(eta_prop: float, V_inf: float, rho_inf: float, P: float, A_prop: float) -> float:
    """
    Propulsive efficiency equation across a range of flight conditions.
    
    Parameters:
        eta_prop: Propulsive efficiency [-]
        V_inf: Aircraft airspeed [m/s]
        rho_inf: Air density [kg/m^3]
        P: Propulsive power [W]
        A_prop: Total propeller disk area (x N_prop) [m^2]
        
    Returns:
        f(eta_prop) - eta_prop [-]
    """
    return 2/(1 + (1/V_inf)*np.sqrt(V_inf**2 + 2*P*eta_prop/(rho_inf*A_prop*V_inf))) - eta_prop
    

def get_P_cr(V_cruise: float, rho_inf: float, S: float, C_D0: float, e: float, AR: float, 
            W: np.ndarray, A_prop: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the power required to maintain cruise speed.
    
    Parameters:
        V_cruise: Cruise speed [m/s]
        rho_inf: Air density [kg/m^3]
        S: Wing area [m^2]
        C_D0: Zero-lift drag coefficient [-]
        e: Oswald efficiency factor [-]
        AR: Aspect ratio [-]
        W: Aircraft weight [N]
        A_prop: Total propeller disk area (x N_prop) [m^2]
        
    Returns:
        P: Power required [W]
        eta_prop: Propulsive efficiency [-]
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
    
    Parameters:
        T_0: Ambient temperature [K]
        V_0: Inlet velocity [m/s]
        c_pa: Specific heat of air at constant pressure [J/(kg*K)]
        
    Returns:
        T_02: Temperature at the inlet of the compressor [K]
    """
    return T_0 + (V_0**2)/(2*c_pa)


def get_T_03(T_02: np.ndarray, PI_comp: float, eta_iso: float, k: float) -> np.ndarray:
    """
    Calculate the temperature at the exit of the compressor.
    
    Parameters:
        T_02: Temperature at the inlet of the compressor [K]
        PI_comp: Pressure ratio across the compressor [-]
        eta_iso: Isentropic efficiency of the compressor [-]
        k: Specific heat ratio [-]
        
    Returns:
        T_03: Temperature at the exit of the compressor [K]
    """
    return T_02 + (T_02/eta_iso)*(PI_comp**((k-1)/k) - 1)


def get_mdot_f(mdot_a: np.ndarray, T_04: np.ndarray, T_03: np.ndarray, 
              c_pg: float, LHV_f: float, eta_cc: float) -> np.ndarray:
    """
    Calculate the fuel mass flow rate.
    
    Parameters:
        mdot_a: Air mass flow rate [kg/s]
        T_04: Turbine exit temperature [K]
        T_03: Compressor exit temperature [K]
        c_pg: Specific heat of gas at constant pressure [J/(kg*K)]
        LHV_f: Lower heating value of fuel [J/kg]
        eta_cc: Compressor efficiency [-]
        
    Returns:
        mdot_f: Fuel mass flow rate [kg/s]
    """
    return mdot_a*c_pg*(T_04 - T_03)/(LHV_f*eta_cc)


def ISA(T_0: float, P_0: float, h: float) -> Tuple[float, float, float]:
    """
    Calculate atmospheric properties at a given altitude using the International Standard Atmosphere (ISA) model.
    
    Parameters:
        T_0: Sea level temperature [K]
        P_0: Sea level pressure [Pa]
        h: Altitude [m]
        
    Returns:
        T: Temperature at altitude [K]
        P: Pressure at altitude [Pa]
        rho: Density at altitude [kg/m^3]
    """
    # Constants
    R = 287.05  # Specific gas constant for air [J/(kg*K)]
    g = 9.81    # Gravitational acceleration [m/s^2]
    L = 0.0065  # Temperature lapse rate [K/m]
    
    # Calculate temperature and pressure at altitude
    T = T_0 - L*h
    P = P_0 * (T/T_0)**(g/(L*R))
    rho = P/(R*T)
    
    return T, P, rho


@dataclass
class AircraftParameters:
    """Data class to store aircraft performance parameters."""
    T_0: float              # Sea level temperature [K]
    P_0: float              # Sea level pressure [Pa]
    power_TOGA: float       # Take-off/go-around power [W]
    wing_area: float        # Wing area [m^2]
    CD0: float              # Zero-lift drag coefficient [-]
    wingspan: float         # Wing span [m]
    propeller_diameter: float  # Propeller diameter [m]
    
    @property
    def aspect_ratio(self) -> float:
        """Calculate aircraft aspect ratio."""
        return self.wingspan**2/self.wing_area
    
    @property
    def oswald_efficiency(self) -> float:
        """Calculate Oswald efficiency factor."""
        return 1.78*(1-0.045*self.aspect_ratio**0.68) - 0.64
    
    @property
    def propeller_area(self) -> float:
        """Calculate total propeller disk area."""
        return 2*(np.pi*(0.5*self.propeller_diameter)**2)


class FlightProfileCalculator:
    """
    Calculates and manages flight profile parameters for a hybrid-electric aircraft.
    
    This class handles the computation of flight phases (climb, cruise, descent),
    and associated parameters such as airspeed, altitude, rate of climb, and range.
    """
    
    def __init__(self, aircraft: AircraftParameters):
        """
        Initialize the flight profile calculator with aircraft parameters.
        
        Parameters:
            aircraft: AircraftParameters object containing aircraft specifications
        """
        self.aircraft = aircraft
        self.waypoints = {
            'altitude': np.array([
                0,      # Takeoff 
                2438.4, # Climb station 0 
                4876.8, # Climb station 1 
                7620.0, # Cruise
                304.8,  # Approach
                0       # Landing
            ]),
            'airspeed': np.array([
                54.0167, # Takeoff
                139.929, # Climb station 0
                146.102, # Climb station 1
                142.501, # Cruise
                60.2,    # Approach
                60.2     # Landing
            ]),
            'ROC': np.array([
                797/60,  # Takeoff to climb station 0
                797/60,  # Climb station 0 to climb station 1
                797/60,  # Climb station 1 to cruise
                -6.4,    # Cruise to approach
                -3.1506  # Approach to landing
            ]),
            'power_setting': np.array([
                1.0,  # Takeoff (5 mins)
                0.95, # Rest of climb
                0.19  # Descent
            ]) * aircraft.power_TOGA  # Power setting
        }
        
        self.dt = 1  # Time step [s]
        self.total_range = 707e3  # Total range [m]
        self.climb_legs = 3  # Number of climb segments
        
        # Store calculated profile data
        self.profile_data = {}
    
    def calculate_flight_segments(self) -> None:
        """
        Calculate discretized flight segments from waypoints.
        """
        h = self.waypoints['altitude']
        v = self.waypoints['airspeed']
        roc = self.waypoints['ROC']
        
        # Calculate time-discretized segments
        self.velocity_segments = []
        self.altitude_segments = []
        
        for i in range(len(h)-1):
            # Calculate number of time steps for this segment
            segment_time = abs((h[i+1]-h[i])/(roc[i] if roc[i] != 0 else 1))
            num_steps = int(np.ceil(segment_time/self.dt))
            
            # Create linearly interpolated arrays for velocity and altitude
            v_segment = np.linspace(v[i], v[i+1], num_steps)
            h_segment = np.linspace(h[i], h[i+1], num_steps)
            
            self.velocity_segments.append(v_segment)
            self.altitude_segments.append(h_segment)
    
    def calculate_climb_profile(self) -> None:
        """
        Calculate climb profile parameters.
        """
        # Combine climb segments
        self.profile_data['V_climb'] = np.hstack(self.velocity_segments[:self.climb_legs])
        self.profile_data['alt_climb'] = np.hstack(self.altitude_segments[:self.climb_legs])
        
        # Compute rate of climb array
        roc_segments = []
        for i in range(self.climb_legs):
            roc_segments.append(
                self.waypoints['ROC'][i] * np.ones_like(self.velocity_segments[i])
            )
        self.profile_data['ROC_climb'] = np.hstack(roc_segments)
        
        # Calculate flight path angle and ground speed
        gamma_climb = np.arcsin(self.profile_data['ROC_climb'] / self.profile_data['V_climb'])
        self.profile_data['V_ground_climb'] = self.profile_data['V_climb'] * np.cos(gamma_climb)
        
        # Calculate range
        self.profile_data['range_climb'] = np.cumsum(self.profile_data['V_ground_climb']) * self.dt
        self.profile_data['total_climb_range'] = self.profile_data['range_climb'][-1]
    
    def calculate_descent_profile(self) -> None:
        """
        Calculate descent profile parameters.
        """
        # Combine descent segments
        self.profile_data['V_descent'] = np.hstack(
            self.velocity_segments[self.climb_legs:]
        )
        self.profile_data['alt_descent'] = np.hstack(
            self.altitude_segments[self.climb_legs:]
        )
        
        # Compute rate of climb array for descent
        roc_segments = []
        for i in range(self.climb_legs, len(self.velocity_segments)):
            roc_segments.append(
                self.waypoints['ROC'][i] * np.ones_like(self.velocity_segments[i])
            )
        self.profile_data['ROC_descent'] = np.hstack(roc_segments)
        
        # Calculate flight path angle and ground speed
        gamma_descent = np.arcsin(self.profile_data['ROC_descent'] / self.profile_data['V_descent'])
        self.profile_data['V_ground_descent'] = self.profile_data['V_descent'] * np.cos(gamma_descent)
        
        # Calculate range
        self.profile_data['range_descent_local'] = np.cumsum(self.profile_data['V_ground_descent']) * self.dt
        self.profile_data['total_descent_range'] = self.profile_data['range_descent_local'][-1]
    
    def calculate_cruise_profile(self) -> None:
        """
        Calculate cruise profile parameters.
        """
        # Calculate cruise range
        R_cruise = self.total_range - self.profile_data['total_climb_range'] - self.profile_data['total_descent_range']
        self.profile_data['R_cruise'] = R_cruise
        
        # Cruise airspeed
        V_cr = self.waypoints['airspeed'][3]
        
        # Calculate number of time steps for cruise
        n_cr = int(np.ceil((R_cruise / V_cr) / self.dt))
        
        # Create arrays for cruise segment
        self.profile_data['V_cr_arr'] = V_cr * np.ones(n_cr)
        self.profile_data['ROC_cr'] = np.zeros(n_cr)
        self.profile_data['alt_cr'] = self.waypoints['altitude'][3] * np.ones(n_cr)
        self.profile_data['V_ground_cr'] = V_cr
        self.profile_data['range_cr_local'] = np.arange(1, n_cr+1) * V_cr * self.dt
    
    def assemble_complete_profile(self) -> None:
        """
        Assemble the complete flight profile by combining climb, cruise, and descent segments.
        """
        # Combine all profile arrays
        self.profile_data['V_full'] = np.hstack([
            self.profile_data['V_climb'],
            self.profile_data['V_cr_arr'],
            self.profile_data['V_descent']
        ])
        
        self.profile_data['ROC_full'] = np.hstack([
            self.profile_data['ROC_climb'],
            self.profile_data['ROC_cr'],
            self.profile_data['ROC_descent']
        ])
        
        self.profile_data['alt_full'] = np.hstack([
            self.profile_data['alt_climb'],
            self.profile_data['alt_cr'],
            self.profile_data['alt_descent']
        ])
        
        # Create time and range arrays
        self.profile_data['time_full'] = np.arange(len(self.profile_data['V_full'])) * self.dt
        
        self.profile_data['range_full'] = np.hstack([
            self.profile_data['range_climb'],
            self.profile_data['total_climb_range'] + self.profile_data['range_cr_local'],
            self.profile_data['total_climb_range'] + self.profile_data['R_cruise'] + 
                self.profile_data['range_descent_local']
        ])
    
    def calculate_complete_profile(self) -> Dict:
        """
        Execute the full flight profile calculation process.
        
        Returns:
            Dict containing all calculated flight profile data
        """
        self.calculate_flight_segments()
        self.calculate_climb_profile()
        self.calculate_descent_profile()
        self.calculate_cruise_profile()
        self.assemble_complete_profile()
        return self.profile_data
    
    def plot_flight_profile(self) -> None:
        """
        Generate plots of the flight profile data.
        """
        # Create figure with subplots
        plt.figure(figsize=(12, 8))
        
        # Airspeed vs Time
        plt.subplot(2, 2, 1)
        plt.plot(self.profile_data['time_full'], self.profile_data['V_full'])
        plt.title("Airspeed vs Time")
        plt.xlabel("Time [s]")
        plt.ylabel("Airspeed [m/s]")
        plt.grid(True, alpha=0.3)
        
        # Ground Range vs Time
        plt.subplot(2, 2, 2)
        plt.plot(self.profile_data['time_full'], self.profile_data['range_full'] / 1000)  # Convert to km
        plt.title("Ground Range vs Time")
        plt.xlabel("Time [s]")
        plt.ylabel("Range [km]")
        plt.grid(True, alpha=0.3)
        
        # Rate of Climb vs Time
        plt.subplot(2, 2, 3)
        plt.plot(self.profile_data['time_full'], self.profile_data['ROC_full'])
        plt.title("Rate of Climb vs Time")
        plt.xlabel("Time [s]")
        plt.ylabel("ROC [m/s]")
        plt.grid(True, alpha=0.3)
        
        # Altitude vs Time
        plt.subplot(2, 2, 4)
        plt.plot(self.profile_data['time_full'], self.profile_data['alt_full'])
        plt.title("Altitude vs Time")
        plt.xlabel("Time [s]")
        plt.ylabel("Altitude [m]")
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    # Initialize aircraft parameters
    aircraft_params = AircraftParameters(
        T_0=288.15,                # Sea level temperature [K]
        P_0=101325,                # Sea level pressure [Pa]
        power_TOGA=1.908e6,        # TOGA power [W]
        wing_area=28.79,           # Wing area [m^2]
        CD0=0.0215,                # Zero-lift drag coefficient [-]
        wingspan=17.64,            # Wing span [m]
        propeller_diameter=2.78    # Propeller diameter [m]
    )
    
    # Create flight profile calculator and compute profile
    flight_calculator = FlightProfileCalculator(aircraft_params)
    profile = flight_calculator.calculate_complete_profile()
    
    # Generate plots
    flight_calculator.plot_flight_profile()
    
    # Future calculation of power required for cruise can be added here
    # flight_calculator.calculate_cruise_power_requirements()