import heat_removal_required
from CoolProp.CoolProp import PropsSI, Props
from pyfluids import Fluid, FluidsList



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
        


    def mass_flow(self):
        """
        Calculate the mass flow rate of the coolant
        """
        self.Q_dot = self.design_point.fuel_cell.heat_exchange_fc
        self.T_in = self.design_point.fuel_cell.T
        self.T_amb = self.design_point.flight_condition.T_amb
        self.T_out = self.T_in - (self.T_in - self.T_amb)/2 # K Placeholder value because we don't know
        print(f"T_out: {self.T_out} K")
        print(f"T_in: {self.T_in} K")
        print(f"T_amb: {self.T_amb} K")


        self.cp_coolant = 2300 # PropsSI('C', 'T', self.T_in, 'P', 101325, "TherminolD12") # J/kg/K Coolant subject to change
        print(f"cp_coolant: {self.cp_coolant} J/kg/K")
        # Calculate the mass flow rate of the coolant
        self.m_dot_cool = self.Q_dot / (self.cp_coolant * (self.T_in - self.T_out))

