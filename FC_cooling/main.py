import heat_exchanger_sizing
import fuel_cell_info
import flight_condition
import hydrogen_storage
import heat_removal_required



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
    power_split = 1,  # proportion of power outputted from the fuel cell
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
HTPEM_TO = heat_removal_required.Design_point(fuel_cell=LTPEM, flight_condition=takeoff, hydrogen_storage=LH2, P_C = 1.85 * 101325)

LHV_H2 = 120000000  # Lower heating value of hydrogen in J/kg
print(HTPEM_TO.fuel_cell.name)
print(HTPEM_TO.flight_condition.name)
HTPEM_TO.O2_cooling_required()
HTPEM_TO.mass_flow_calculation()
HTPEM_TO.heat_removal_available()

Heat_exchanger = heat_exchanger_sizing.HeatExchanger(design_point=HTPEM_TO)
Heat_exchanger.mass_flow()
print(f"Mass flow rate of coolant: {Heat_exchanger.m_dot_cool} kg/s")




T_outlet = (80 + 273.15) * (12/1.8)**((1.4 - 1)/1.4) # K
print(f"Outlet temperature of the coolant: {T_outlet} K")