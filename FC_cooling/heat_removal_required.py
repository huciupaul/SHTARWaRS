from CoolProp.CoolProp import PropsSI
import numpy as np
import fuel_cell_info
import flight_condition
import hydrogen_storage

'''
Class which defines the design point of the fuel cell & flight condition
This class is used to store all the values found for the design point
'''
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
        print(f"Heat added to from air after compression: {self.fuel_cell.heat_change_air * 10**(-6)} MJ/kg")


    def mass_flow_calculation(self):

        # Calculate air/hydrogen mass flow ratio
        self.fuel_cell.air_hydrogen_ratio = 15.999 / 1.00784 * self.fuel_cell.stoic_ratio_A / 0.2314 # mass fraction of air per hydrogen
        

        # Calculate the required mass flow rate of gasses for the fuel cell
        self.fuel_cell.power_required = self.flight_condition.power_required * self.flight_condition.power_split
        self.fuel_cell.m_H2_tot = self.fuel_cell.power_required / ((self.fuel_cell.stack_efficiency  * LHV_H2 - self.fuel_cell.air_compression_energy * self.fuel_cell.air_hydrogen_ratio) * self.flight_condition.propulsive_efficiency)
        self.fuel_cell.m_H2_used = self.fuel_cell.m_H2_tot / self.fuel_cell.stoic_ratio_A
        self.fuel_cell.m_O2_used = self.fuel_cell.m_H2_used * 15.999/1.00784
        self.fuel_cell.m_O2_tot = self.fuel_cell.m_O2_used * self.fuel_cell.stoic_ratio_C
        self.fuel_cell.m_air = self.fuel_cell.m_O2_tot / 0.2314  # Divided by mass fraction of O2 in air

        self.fuel_cell.m_tot = self.fuel_cell.m_air + self.fuel_cell.m_H2_tot
        self.fuel_cell.m_H2O_out = self.fuel_cell.m_H2_used + self.fuel_cell.m_O2_used
        self.fuel_cell.m_H2_out = self.fuel_cell.m_H2_tot - self.fuel_cell.m_H2_used
        self.fuel_cell.m_O2_out = self.fuel_cell.m_O2_tot - self.fuel_cell.m_O2_used
        self.fuel_cell.m_rest_out = self.fuel_cell.m_tot - self.fuel_cell.m_H2O_out - self.fuel_cell.m_H2_out - self.fuel_cell.m_O2_out


        print(f"Mass flow rate of used oxygen for fuel cell: {self.fuel_cell.m_O2_used} kg/s")
        print(f"Mass flow rate of air for fuel cell: {self.fuel_cell.m_air} kg/s")
        
        self.fuel_cell.reaction_efficiency = self.fuel_cell.stack_efficiency * self.fuel_cell.stoic_ratio_A # Efficiency of the reaction (to determine how much heat is generated)
        self.fuel_cell.total_efficiency = self.fuel_cell.stack_efficiency - self.fuel_cell.air_compression_energy * self.fuel_cell.air_hydrogen_ratio / (LHV_H2) # Efficiency of the fuel cell
        self.fuel_cell.heat_power = LHV_H2 * self.fuel_cell.m_H2_used * (1 - self.fuel_cell.reaction_efficiency)  # Heat power produced by the fuel cell


        print(f"Mass flow rate of hydrogen used in fuel cell: {self.fuel_cell.m_H2_tot} kg/s")
        print(f"Total efficiency of the fuel cell: {self.fuel_cell.total_efficiency * 100} %")
        print(f"Electrical power produced by the fuel cell: {self.fuel_cell.power_required * 10**(-6)} MW")
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



# Create three fuel cell types
LTPEM = fuel_cell_info.FuelCell(
    name = "LTPEM",
    stack_efficiency=0.6,
    T=273.15 + 80,  # K
    P_A=1.8 * 101325 + 0.06 * 10**5,  # Pa https://pdf.sciencedirectassets.com/271472/1-s2.0-S0360319921X00614/1-s2.0-S0360319921027634/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEFwaCXVzLWVhc3QtMSJIMEYCIQDn6w4UOTzRqoEPtaJbiJTNcGIZB9kDUYb2ecZ%2BmBXN9wIhAMOwbI4XLuHoL5%2Bsp27d4gD7Tv66lv%2Bm2yXTJiurmwktKrMFCBUQBRoMMDU5MDAzNTQ2ODY1IgzoT0L%2Boa5nr7Qm%2BbQqkAViGsDYeaEzBi%2BfWU4I%2BGGl1oz5YeF1BhCa3EI%2BsdekiMWno3GyGuSawP6UhKhrMmxI2IiRRuMFNbQODhnC9aXp%2BNApcCE6YnLxInhezcfKHTenH67lBh2Caol23dSWqQ7g6xNDYQSvWmt2o1GO%2B9thieoAoMwWWsIbcYh7JzjS9%2BWZossrFcPvxYWLhx978ah9UStWo3wfqn5VRnbAcGqwB29JmIoGWdeOxH99kXZLH8ErtK7rtZeXsovb27UXBkYUTtaVE0%2BTxzWsTiHUJUcknbgkdDvdZ0Tc5CTdMeepnDQSx2ad1teq43xKqIaaHEJ%2BZsFjcPt%2FaTn4AqvGYoAAEyYtrgBPJjD8nrFwpN%2B1VGYDNWUwBrL0N%2F67aoI29n%2BeeIj%2F48FjeRoINWoCtcMRFFJAFHJIqjvthQqxgicN%2F%2FwwlgbZplrrVJN2rq2%2Bj7bCV8IJKxu%2Bhm%2FjEWjBaU5%2FDLcf9gOjTE9Rcf1WBwemFImvlAMduz%2B%2BZnOnKRrNX95GoEJdkIBWHWOXcEA7uEhArudTaVqldwCVjs0VBRo95bUWgpkm78YoOF9%2BK5ETevqXTvceQWtAmVoWZzL5O1SxFYV9CeG1i2uX4y9%2B3%2BxqZlgbIo755dSJJD0w%2B9Oc8MPhSTYbUnXdF%2FiHJfIMDODertE8sKXKbKrFfB1c7ajcO7ptzCP56TxzF9qa8cyAshzlLj%2F%2Brz7FxxyPmV4toDhBSPS3gtg%2F9OPjk6AcCoHRuLjXXi6YKk323Rwpz%2Fy%2F7onraITAqwXvzFfb9GedKAeSuPz%2Bk4etVDmQS87AAMAMBN4v7MehJ4YnI5MCm%2BSaKKWfa0grEMX8z2YFexgk69lVJts0TfbnWKwt9bpgvZL0nDDa%2F5HBBjqwAaSfaaYmRRuwBONTr8vEiVqhE2xjVkEBLlVgXw3XtSO5cUHiAEh1a1f9IVwrIZSepVSyNyboY3yZ5%2B8dhkkkALCSJlGbudhItI4Z0Zc5rFIBIUBiTIbnkiknIz5PJRMRWBdBwKmzxWacUGvSJRD7SP%2BR8csCcd4nyTFeBPrx9J%2Fr%2Fsbsc2KbpghupePp7nREX5kxCbpIAxjlhUc%2FhpRqae4kQqzMcqNbL5NgGmnxwjr2&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250514T123737Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYTNM6Z4RK%2F20250514%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=218bd5bcc7fa27dde6283d4177a9db49807cbe6645b3470448d3e7297e4fa515&hash=23dba463727fe3e5de483bc687ec6fd5a87197a2d15fb33588c30bd19bbd42a8&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0360319921027634&tid=spdf-b7d767c9-1df7-4ff6-9758-ed1bc6e98c49&sid=8b7897ae44a1344a266b0d809b3b3978ce3agxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=140a5b5604065551020d01&rr=93fa72abae34312a&cc=nl
    P_C=1.8 * 101325,  # Pa
    RH_A=0.5,  # Relative humidity of the anode air https://pdf.sciencedirectassets.com/271472/1-s2.0-S0360319925X00203/1-s2.0-S036031992501448X/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEF0aCXVzLWVhc3QtMSJGMEQCIGFSKC6cUV%2FTlte1lJMGc60Z%2BUx0HrZZrzgCHQTgussLAiAxQgV%2BrFoV5RIZr8LAqFLn4apLvVIVefwEMBY6gWwoXiqzBQgWEAUaDDA1OTAwMzU0Njg2NSIMV2RLWboACEgN%2BRh%2BKpAFuplJmh%2FNY3JUamb%2BOqeyn5%2FoQNFDT1l4p6XWwwqy14tBbrJcydxGxgZyDvY594XUDl9E8%2FufOzVID%2BZ9vLBBYOlWtafzLF1qT4T%2Btr6UQJFHxH0NLRGUYfnl4xDBqA6m3VaR%2BJThYEx%2FKWv%2Fta1DOdnN6JQT3wViUMdkFXEJH%2F6TYh98QGwMvTfPAjhULZJROeOxnfi%2FaZR9ik5J%2F4i4WPiDX93YTSVg1glIv9ZkisyyW8gsMDiJDKgFfjYyETHqvwthoJQ85T2T4dDa8h4EilDODwoqVIarsl%2FCLutpBKMbpKcrm9xvXIjPWs24N1uEi7JYa%2BNnDfE8Wqv5nhN2MiAJqAQwsbz86%2B1VT1IsUyVJKKD%2FM6XwYOmXbnzuO3B8n6yllO6M6jbOh4g85%2FoiEHkp8iED3IN%2Bn1bCZDFWU0kUc5DLxaB8Xl3adQnvLssg5gBAx8xLiq1AFOtL2OlWCU1ghGL4s6dZvNr2sDUF2oiG5KV4njKAtzf4LgXO5ksxPmTN%2FQ%2B7byA3qfrvse14XOpB08gJVT3jRVSWuuxitQLh1MBl9n1jLHMDGnHjKeyWPySdefZIi5XnM3CGCdoBlBPQEWoDUpiDTz0wu%2BrikS7z8b0CfCiC3p%2FpuRDdnqFozW4lH%2Ba4Hw%2BandKw2sfp2irrBqINxOHgp2VNvdTuFBxW1%2B8A6oFaAtqlzebbISbRv4D1MoCmQ1tmJm09NQo7E38zvoHeGpzQq0ytwxRQ%2FVMNZfwovic8GebetFaePwQm47oThz9FNAmz4gkFaOi7fM4qJBCVuLVEK%2BE0DTS2FdEcn2N5opDTZ6wR86nrpGHqJ5Nmhhu%2B2KorHjUxnLaDymkXSylUGcZRnp7VdlxX6b8w%2FqiSwQY6sgHxTRT3%2F%2FIK3Gwb8Ft0vqxrFHuh4CR59Mt3o4kYs0nqCSpLKSii310q0YA4I5fn%2FO2n9BNJU0aXnUAGEIJVZu9LM9rBrqO6wZlMeZviNZhmXoRhrFKnJF3lSROZEQVzMGVDsk2sKCSA0B%2Fu03d8ZGrLSYq0Vl5tI%2FcHLgwwhor3yc4LQJZ93xFrV2hlZcXN1ev3wQih149PdA%2BaVWAIYcL%2BvE8heNyqs340O380jBIenDgi&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250514T135049Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY74AGHFSZ%2F20250514%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=a14d911c3d9cd53464a9ef0b436e42f5c3092f7e1a616e09785a83869d6d732e&hash=8aeb960286ea350c176c237d2ded6cf83e5cda255e1304799b323e34ddec6791&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S036031992501448X&tid=spdf-b31244d5-d4e7-42bc-a100-b06567e324ee&sid=8b7897ae44a1344a266b0d809b3b3978ce3agxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=140a5b5604060607065a55&rr=93fadde5585e1c7d&cc=nl
    RH_C=0.38,  # Relative humidity of the cathode air
    stoic_ratio_C=1.8,  # Stoichiometric ratio of the cathode air
    stoic_ratio_A=1.05  # Stoichiometric ratio of the anode air
)

HTPEM = fuel_cell_info.FuelCell(
    name = "HTPEM",
    stack_efficiency=0.6,
    T=273.15 + 160,  # K
    P_A=1.85 * 101325 + 0.06 * 10**5,  # Pa
    P_C=1.85 * 101325,  # Pa
    RH_A=0.5,  # Relative humidity of the anode air
    RH_C=0.38,  # Relative humidity of the cathode air
    stoic_ratio_C=1.8,  # Stoichiometric ratio of the cathode air
    stoic_ratio_A=1.05  # Stoichiometric ratio of the anode air
)

SOFC = fuel_cell_info.FuelCell(
    name = "SOFC",
    stack_efficiency=0.65,
    T=273.15 + 800,  # K
    P_A=2 * 101325 + 0.06 * 10**5,  # Pa
    P_C=2 * 101325,  # Pa
    RH_A=0.5,  # Relative humidity of the anode air
    RH_C=0.38,  # Relative humidity of the cathode air
    stoic_ratio_C=1.8,  # Stoichiometric ratio of the cathode air
    stoic_ratio_A=1.05  # Stoichiometric ratio of the anode air
)


# Create a flight condition object for takeoff
takeoff = flight_condition.FlightCondition(
    name = "Takeoff",
    T_amb=273.15 + 55,  # K
    P_amb=101325,  # Pa
    RH_amb=0,  # Relative humidity of the ambient air
    V=0,  # m/s
    power_required = 1908000, # Power required in W
    power_split = 0.5,  # proportion of power outputted from the fuel cell
    thermal_efficiency = 0.4,  # thermal efficiency of the combustion chamber
    propulsive_efficiency = 0.85, # propulsive efficiency of the propeller
    P_cc = 12.1 * 101325,  # Pa
    T_cc = 300 + 273.15 #K
)

# Create a hydrogen storage object
LH2 = hydrogen_storage.HydrogenStorage(
    name = "LH2",
    T = 20,  # K
    P = 6 * 101325  # Pa
)

CcH2 = hydrogen_storage.HydrogenStorage(
    name = "CcH2",
    T = 273.15 - 207,  # K
    P = 350 * 101325  # Pa
)

GCH2 = hydrogen_storage.HydrogenStorage(
    name = "GCH2",
    T = 273.15 + 20,  # K
    P = 700 * 101325  # Pa
)





    
# Create an instance of the Output class
HTPEM_TO = Design_point(fuel_cell=HTPEM, flight_condition=takeoff, hydrogen_storage=LH2, P_C = 1.85 * 101325)

LHV_H2 = 120000000  # Lower heating value of hydrogen in J/kg

HTPEM_TO.O2_cooling_required()
HTPEM_TO.mass_flow_calculation()
HTPEM_TO.heat_removal_available()
