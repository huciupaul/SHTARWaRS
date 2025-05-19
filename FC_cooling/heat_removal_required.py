from CoolProp.CoolProp import PropsSI
import numpy as np
import fuel_cell_info
import flight_condition
import hydrogen_storage

'''
Class which defines the design point of the fuel cell & flight condition
This class is used to store all the values found for the design point
'''

LHV_H2 = 120000000  # Lower heating value of hydrogen in J/kg

class Design_point:
    def __init__(self, fuel_cell, flight_condition, hydrogen_storage, P_C):
        self.fuel_cell = fuel_cell
        self.flight_condition = flight_condition
        self.hydrogen_storage = hydrogen_storage
        p_diff = self.fuel_cell.P_C - self.fuel_cell.P_A
        self.fuel_cell.P_C = P_C
        self.fuel_cell.P_A = P_C - p_diff

            
    def O2_cooling_required(self):
        # Calculate the heat removal required from the oxygen for the fuel cell
        H_air_initial = PropsSI('H', 'T', self.flight_condition.T_tot, 'P', self.flight_condition.P_tot, 'Air')
        gamma_air = PropsSI('CPMASS', 'T', self.flight_condition.T_tot, 'P', self.flight_condition.P_tot, 'Air') / PropsSI('CVMASS', 'T', self.flight_condition.T_tot, 'P', self.flight_condition.P_tot, 'Air')
        P_air_comp = self.fuel_cell.P_A
        T_air_comp = self.flight_condition.T_tot * (P_air_comp / self.flight_condition.P_tot)**((gamma_air - 1) / gamma_air) # Adiabatic compression
        H_air_comp = PropsSI('H', 'T', T_air_comp, 'P', P_air_comp, 'Air')
        H_air_fc = PropsSI('H', 'T', self.fuel_cell.T, 'P', self.fuel_cell.P_C, 'Air')
        compressor_efficiency = 0.85 ### REVIEW THIS VALUE

        # Calculate the energy required for air compression
        self.fuel_cell.air_compression_energy = (H_air_comp - H_air_initial) / compressor_efficiency  # J/kg
        
        # Calculate the heat removal required from the air
        self.fuel_cell.heat_change_air = (H_air_fc - H_air_comp)  # J/kg

        print(f"Power required for air compression: {self.fuel_cell.air_compression_energy * 10**(-6)} MJ/kg")
        print(f"Heat added to air after compression: {self.fuel_cell.heat_change_air * 10**(-6)} MJ/kg")


    def mass_flow_calculation(self):

        # Calculate air/hydrogen mass flow ratio
        self.fuel_cell.air_hydrogen_ratio = 15.999 / 1.00784 * self.fuel_cell.stoic_ratio_C / self.fuel_cell.stoic_ratio_A / 0.2314 # mass fraction of total air per total hydrogen
        

        # Calculate the required mass flow rate of gasses for the fuel cell
        self.fuel_cell.shaft_power_required = self.flight_condition.power_required * self.flight_condition.power_split # W
        self.fuel_cell.m_H2_tot = self.fuel_cell.shaft_power_required / (self.fuel_cell.stack_efficiency  * LHV_H2 - self.fuel_cell.air_compression_energy * self.fuel_cell.air_hydrogen_ratio)
        self.fuel_cell.m_H2_used = self.fuel_cell.m_H2_tot / self.fuel_cell.stoic_ratio_A
        self.fuel_cell.m_O2_used = self.fuel_cell.m_H2_used * 15.999/1.00784
        self.fuel_cell.m_O2_tot = self.fuel_cell.m_O2_used * self.fuel_cell.stoic_ratio_C
        self.fuel_cell.m_air = self.fuel_cell.m_O2_tot / 0.2314  # Divided by mass fraction of O2 in air

        
        self.fuel_cell.total_power_O2_comp = self.fuel_cell.air_compression_energy * self.fuel_cell.m_air # W
        print(f"Power required for air compression (total): {self.fuel_cell.total_power_O2_comp * 10**(-6)} MW")

        self.fuel_cell.m_tot = self.fuel_cell.m_air + self.fuel_cell.m_H2_tot
        self.fuel_cell.m_H2O_out = self.fuel_cell.m_H2_used + self.fuel_cell.m_O2_used
        self.fuel_cell.m_H2_out = self.fuel_cell.m_H2_tot - self.fuel_cell.m_H2_used
        self.fuel_cell.m_O2_out = self.fuel_cell.m_O2_tot - self.fuel_cell.m_O2_used
        self.fuel_cell.m_rest_out = self.fuel_cell.m_tot - self.fuel_cell.m_H2O_out - self.fuel_cell.m_H2_out - self.fuel_cell.m_O2_out

        print(f"Mass flow rate of hydrogen used in fuel cell: {self.fuel_cell.m_H2_tot} kg/s")
        print(f"Mass flow rate of used oxygen for fuel cell: {self.fuel_cell.m_O2_used} kg/s")
        print(f"Mass flow rate of air for fuel cell: {self.fuel_cell.m_air} kg/s")
        print(f"Mass flow rate of hydrogen for fuel cell: {self.fuel_cell.m_H2_tot} kg/s")
        print(f"Mass flow rate of water produced by fuel cell: {self.fuel_cell.m_H2O_out} kg/s")
        
        self.fuel_cell.reaction_efficiency = self.fuel_cell.stack_efficiency * self.fuel_cell.stoic_ratio_A # Efficiency of the reaction (to determine how much heat is generated)
        self.fuel_cell.total_efficiency = self.fuel_cell.stack_efficiency - self.fuel_cell.air_compression_energy * self.fuel_cell.air_hydrogen_ratio / (LHV_H2) # Efficiency of the fuel cell
        self.fuel_cell.heat_power = LHV_H2 * self.fuel_cell.m_H2_used * (1 - self.fuel_cell.reaction_efficiency)  # Heat power produced by the fuel cell
        self.fuel_cell.total_electrical_power = LHV_H2 * self.fuel_cell.m_H2_tot * self.fuel_cell.stack_efficiency  # Total electrical power produced by the fuel cell


        print(f"Total efficiency of the fuel cell: {self.fuel_cell.total_efficiency * 100} %")
        print(f"Shaft power produced by the fuel cell: {self.fuel_cell.shaft_power_required * 10**(-6)} MW")
        print(f"Total power produced by the fuel cell (+BOP): {self.fuel_cell.total_electrical_power * 10**(-6)} MW")
        print(f"Heat power produced by the fuel cell: {self.fuel_cell.heat_power * 10**(-6)} MW")

        # Calculate the required mass flow rate of hydrogen for the combustion chamber
        self.cc_power = self.flight_condition.power_required * (1 - self.flight_condition.power_split)
        self.m_H2_cc = self.cc_power / (self.flight_condition.thermal_efficiency * self.flight_condition.propulsive_efficiency * LHV_H2)

        # Calculate the total mass flow rate
        self.m_H2_tot = self.m_H2_cc + self.fuel_cell.m_H2_tot

        print(f"Mass flow rate of hydrogen used in fuel cell: {self.fuel_cell.m_H2_tot} kg/s")
        print(f"Mass flow rate of hydrogen for combustion chamber: {self.m_H2_cc} kg/s")
        print(f"Total mass flow rate of hydrogen: {self.m_H2_tot} kg/s")




    def heat_removal_available(self):
        # Calculate the heat removal available from the hydrogen for the fuel cell (combustion chamber treated separately)

        H_fin_fc = PropsSI('H', 'T', self.fuel_cell.T, 'P', self.fuel_cell.P_A, 'Hydrogen')
        self.HSP_fc_H2_sp = H_fin_fc - self.hydrogen_storage.H_init # J/kg
        self.HSP_fc_H2 = self.HSP_fc_H2_sp * self.fuel_cell.m_H2_tot # W
        print(f"Heat added to hydrogen for fuel cell: {self.HSP_fc_H2 * 10**(-6)} MW")

        # Calculate the heat removal available from the hydrogen for the combustion chamber
        H_fin_cc = PropsSI('H', 'T', self.flight_condition.T_cc, 'P', self.flight_condition.P_cc, 'Hydrogen')
        self.HSP_cc_H2_sp = H_fin_cc - self.hydrogen_storage.H_init # J/kg
        self.HSP_cc_H2 = self.HSP_cc_H2_sp * self.m_H2_cc
        
        # Calculate the heat removal available from the hydrogen for the combustion chamber until the fuel cell temp
        H_fc_cc = PropsSI('H', 'T', self.fuel_cell.T, 'P', self.flight_condition.P_cc, 'Hydrogen')
        self.HSP_cc_fc_H2_sp = H_fc_cc - self.hydrogen_storage.H_init # J/kg
        self.HSP_cc_fc_H2 = self.HSP_cc_fc_H2_sp * self.m_H2_cc # W
        print(f"Heat added to hydrogen for combustion chamber until fuel cell temperature: {self.HSP_cc_fc_H2 * 10**(-6)} MW")

        print(f"Heat added to hydrogen for combustion chamber: {self.HSP_cc_H2 * 10**(-6)} MW")

        # Calculate heat exchange required for the fuel cell
        self.fuel_cell.heat_exchange_fc = self.fuel_cell.heat_power - (self.HSP_fc_H2 + self.HSP_cc_fc_H2) # W
        print(f"Heat exchange required for the fuel cell: {self.fuel_cell.heat_exchange_fc * 10**(-6)} MW")

        self.fuel_cell.mass =  self.fuel_cell.shaft_power_required / self.fuel_cell.spec_power # kg
        self.fuel_cell.volume = self.fuel_cell.total_electrical_power / self.fuel_cell.spec_vol_power # m^3

        print(f"Mass of the fuel cell: {self.fuel_cell.mass} kg")
        print(f"Volume of the fuel cell: {self.fuel_cell.volume} m^3")


