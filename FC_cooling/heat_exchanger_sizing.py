import heat_removal_required

from CoolProp.CoolProp import PropsSI, Props
from pyfluids import Fluid, FluidsList

import numpy as np
from scipy.optimize import fsolve

class HeatExchanger:
    """
    Class which defines the heat exchanger for the cooling system
    """
    design_point: heat_removal_required.Design_point


    def __init__(self, design_point: heat_removal_required.Design_point):

        self.design_point = design_point

        # self.design_point.O2_cooling_required()
        # self.design_point.mass_flow_calculation()
        # self.design_point.heat_removal_available()

    def mass_flow(self, cp_coolant:float):
        """
        Calculate the mass flow rate of the coolant
        """
        self.Q_dot = self.design_point.fuel_cell.heat_exchange_fc
        self.T_in = self.design_point.fuel_cell.T
        self.T_amb = self.design_point.flight_condition.T_amb
        self.T_out = self.T_in - (self.T_in - self.T_amb)/2 # K Placeholder value because we don't know

        # Calculate the mass flow rate of the coolant
        self.m_dot_cool = self.Q_dot / (cp_coolant * (self.T_in - self.T_out))
        return self.m_dot_cool

    #def u_radiator(self, A_norm: float, delta_T_air_coolant: float):
        """
        Compute U_rad: overall heat transfer coefficient (W/m^2/K)
        """
        U_rad = self.Q_dot / (A_norm * delta_T_air_coolant)
        return U_rad

    def area_U(self, delta_T_air_coolant: float):
        self.A_U = self.Q_dot/delta_T_air_coolant 
        return self.A_U
    

    def D_int(self):
        """
        Solve for internal diameter D_i, assuming fixed wall thickness t
        """


    def T_coolant_avg(self, T_in: float, T_out: float):
        """
        Compute average temperature of coolant
        """
        return 0.5 * (T_in + T_out)

    def delta_T_air_coolant_arithmetic(self, T_coolant_avg: float, T_air_in: float):
        """"
        Compute temperature difference of air and coolant: arithmetic method
        """
        return T_coolant_avg - T_air_in

    def delta_T_air_coolant_LMTD(self, T_coolant_in: float, T_coolant_out: float,
                                  T_air_in: float, T_air_out: float, F: float = 1.0):
        """
        Compute temperature difference of air and coolant: LMTD method

        IMPORTANT: define F: most likely not equal to 1.0

        """
        delta_T1 = T_coolant_in - T_air_out
        delta_T2 = T_coolant_out - T_air_in
        delta_T_LMTD = F * (delta_T1 - delta_T2) / np.log(delta_T1 / delta_T2)

        return delta_T_LMTD

    def r_ratio(self, T_coolant_in: float, T_coolant_out: float, T_air_in: float, T_air_out: float):
        """"
        LMTD preferred to the arithmetic mean when r > 1.5
        """
        r_ratio = (T_coolant_in - T_air_in) / (T_coolant_out - T_air_out)
        if r_ratio > 1.5:
            print("r-ratio is higher than 1.5, use LMTD method")
        return r_ratio

    def R(self, T_coolant_in: float, T_coolant_out: float, T_air_in: float, T_air_out: float):
        R = (T_coolant_in - T_coolant_out) / (T_air_out - T_air_in)
        return R

    def P(self, T_coolant_in: float, T_air_in: float, T_air_out: float):    
        P = (T_air_out - T_air_in) / (T_coolant_in - T_air_in)
        return P



