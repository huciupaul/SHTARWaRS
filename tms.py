import math
import abc
import numpy as np
from CoolProp.CoolProp import PropsSI

# ------------------------------ PARAMETERS -----------------------------------

# Ambient Conditions --------------------------------
ambient_conditions= {
    'T': 283.5,  # Ambient temperature in Kelvin
    'Mach': 0.4,  # Flight Mach number
    'V': 100,  # Flight speed in m/s
    'rho': 1.225  # Air density at sea level in kg/m^3
}

# Power and Efficiency Parameters --------------------
total_power = 1.812 *1e6
power_split_fc = 0.45
fc_power = total_power * power_split_fc
cc_power = total_power * (1 - power_split_fc)
fc_eff = 0.6

# Battery
bat_p = 0.0
battery_heat_fraction = 0.05

#EM
em_power = 400000  # Watts
em_eff = 0.99 * 0.95

# TEG
efficiency_teg = 0.05 

# H2 ---------------------------------------------------
t_init = 20 # K
# FC Loop
h2_mf_fc = 0.0156  # kg/s per MW
p_init_fc = 6e5 
p_fin_fc = 3e5
t_fin_fc = 180+273.15 # K

# CC Loop
h2_mf_cc = 0.0294  # kg/s 
p_init_cc = 6e5
p_fin_cc = 12.1e5
t_fin_cc = 180+273.15 # K

# Ram Air HX -------------------------------
# Air properties
h_air = 200 
cp_air = 1005.0
T_air_in = ambient_conditions['T']
T_air_out = ambient_conditions['T'] + 20.0

# Coolant and HX Properties
h_cool = 1500 # average coolant heat transfer coefficient for double phase
U_ra = 1 / (1/h_air + 1/h_cool) # overall heat transfer coefficient
ra_coolant_temp = 350.0 # K

# Sizing Parameters
porosity = 0.10              
hole_diameter_mm = 1.0 
meredith_recovery = 0.50
inlet_cd = 0.02 

# Fan 
fan_eff = 0.7
ra_density = ambient_conditions['rho']  # kg/m³, air density
delta_pressure = 1000  # Pa, pressure drop across the fan


#Skin HX -----------------------------------------------	
area_wing = 4.6

#Coolant
k_in = 0.258 # thermal conductivity of coolant [W/m.K]
prandtl_in = 0.71 # fluid
Reynolds = 1e4 # fluid
f = (0.79 * math.log(Reynolds) - 1.64) ** -2
nu = ((f/8)*(Reynolds-1000)*prandtl_in)/(1+12.7*(f/8)**0.5 * (prandtl_in**0.66-1)) # Nusselt number [m²/s]
dh = 0.01 #tube diameter [m]
h_in = k_in * nu / dh
wing_coolant_temp = 320  # K

# Air
prandtl_air = 0.71
reynolds_air = 1e7
h_ext = ambient_conditions['rho'] * ambient_conditions['V'] * cp_air * 0.185 * (math.log10(reynolds_air))**-2.584 * prandtl_air**-0.66
recovery_factor = prandtl_air ** 0.33  

U_wing = 1/(1/h_in + 1/h_ext)


# HX 1 (FC-FC) -----------------------------------------------
# Coolant properties
coolant_cp = 3672.348 # (3672.347969811853, 4222.1975576689)  
k_coolant = 0.485  #  (0.48537754234800523, 0.6806615139526146)
mu_coolant = 0.0002636  # (0.0001982859557074688, 0.0002636263335097737) 

# H2 properties
#h2_cp =  # (9384.657038075966, 16485.931431151545)
mu_h2 = 1.4246e-05  # (1.494409455847204e-06, 1.4246629589288973e-05)
k_h2 = 0.0272 # (0.02721713884147831, 0.24471920308558517)


# Main Classes for Thermal Management System (TMS) -------------------------------	

class HeatSource(abc.ABC):
    def __init__(self, name, power_output_w):
        self.name = name
        self.power_output = power_output_w  # Watts
    @abc.abstractmethod
    def waste_heat(self):
        pass

class HeatSink(abc.ABC):
    def __init__(self, name):
        self.name = name
    @abc.abstractmethod
    def absorb_heat(self, heat_w):
        pass
 

# Main Subclasses for Heat Sources (TMS) -------------------------------

class FuelCell(HeatSource):
    def __init__(self, power_output_w, efficiency):
        super().__init__("FuelCell", power_output_w)
        self.efficiency = efficiency 
    def waste_heat(self):
        q = self.power_output * (1.0/self.efficiency - 1.0)
        return q  # Watts

class Battery(HeatSource):
    def __init__(self, power_output_w, heat_fraction):
        super().__init__("Battery", power_output_w)
        self.heat_fraction = heat_fraction
    def waste_heat(self):
        # Simplified battery model (constant for maximum)
        return self.power_output * self.heat_fraction

class ElectricMotor(HeatSource):
    def __init__(self, power_output_w, efficiency):
        super().__init__("Motor", power_output_w)
        self.efficiency = efficiency  
    def waste_heat(self):
        return self.power_output * (1.0/self.efficiency - 1.0)
    

# Main Subclasses for Heat Sinks (TMS) -------------------------------
    
class LH2Tank(HeatSink):
    def __init__(self, mass_flow_kg_s, T_final_K=273+180, P_final=3e5, P_init = 600000):
        super().__init__("LH2Tank")
        self.m_dot = mass_flow_kg_s  
        self.T_final = T_final_K         # K (before fuel cell)
        self.P_final = P_final
        self.P_init = P_init
    def absorb_heat(self,heat_w):
        H_fin = PropsSI('H', 'T', self.T_final, 'P', self.P_final, 'ParaHydrogen')  # Final enthalpy
        H_init = PropsSI('H', 'P', self.P_init, 'Q', 0, 'ParaHydrogen')  # Initial enthalpy
        Q_absorbed = self.m_dot * (H_fin - H_init)  # Heat absorbed in Joules
        return Q_absorbed
    
class SkinHeatExchanger(HeatSink):
    def __init__(self, area_m2, U_W_per_m2K, coolant_temp_K):
        super().__init__("SkinHX")
        self.area = area_m2
        self.U = U_W_per_m2K              # overall HT coefficient (W/m^2.K)
        self.coolant_temp = coolant_temp_K
        self.ambient_temp = None          # to be set based on flight condition

    def set_ambient(self, T_ambient_K, recovery_factor=0.9, mach=0.3):
        # Set effective ambient (adiabatic wall) temperature for convection calculations
        T_aw = T_ambient_K * (1 + recovery_factor*0.5*(1.4-1)*mach**2)
        self.ambient_temp = T_aw 
    def absorb_heat(self, heat_w):
        if self.ambient_temp is None:
            raise RuntimeError("Ambient conditions not set for SkinHeatExchanger")

        Q_capacity = self.U * self.area * (self.coolant_temp - self.ambient_temp)
        Q_absorbed = min(heat_w, Q_capacity)
        if Q_capacity > Q_absorbed:
            print(f"SkinHX: Absorbed {Q_absorbed:.2f} W, remaining capacity {Q_capacity - Q_absorbed:.2f} W")
        return Q_absorbed
    
class RamAirHeatExchanger(HeatSink):
    def __init__(self, U_W_per_m2K, coolant_temp_K):
        super().__init__("RamAirHX")
        self.U = U_W_per_m2K
        self.coolant_temp = coolant_temp_K

        # To be determined during sizing:
        self.required_area = 0.0    # m² of radiator face
        self.inlet_area = 0.0       # m² of total open hole area
        self.skin_area = 0.0        # m² of skin needed (porous area)
        self.hole_count = 0         # total number of small holes
        self.drag_N = 0.0           # N, raw drag before Meredith effect
        self.net_drag_N = 0.0       # N, after Meredith recovery

    def absorb_heat(self, heat_w):
        return heat_w

    def size_exchanger(self, heat_w, T_air_in_K, T_air_out_K):
        dT1 = self.coolant_temp - T_air_in_K
        dT2 = self.coolant_temp - T_air_out_K

        if abs(dT2 - dT1) < 1e-6:
            dT_lm = 0.5 * (dT1 + dT2)
        else:
            dT_lm = (dT2 - dT1) / math.log(dT2 / dT1)

        self.required_area = heat_w / (self.U * dT_lm)
        return self.required_area

    def compute_hole_and_skin_area(self,
                                   heat_w,
                                   T_air_in_K,
                                   T_air_out_K,
                                   flight_speed_m_s,
                                   air_density,
                                   porosity=0.10,
                                   hole_diameter_mm=1.0,
                                   inlet_CD=0.02,
                                   cp_air=1005.0,
                                   meredith_recovery=0.50):
        
        delta_T_air = T_air_out_K - T_air_in_K
        if delta_T_air <= 0:
            raise ValueError("T_air_out must exceed T_air_in to carry heat away.")

        m_dot_air = heat_w / (cp_air * delta_T_air)  # kg/s
        A_open = m_dot_air / (air_density * flight_speed_m_s)  # m²
        skin_area = A_open / porosity  # m² of total panel with holes
        if skin_area >= 0.10:
            d_m = hole_diameter_mm / 1000.0
            area_per_hole = math.pi * (d_m / 2) ** 2 # m² per hole
            hole_count = int(math.ceil(A_open / area_per_hole))
        else:
            hole_count = 0
            d_m = 0.0
            area_per_hole = 0.0

        # Base drag
        drag_N = 0.5 * air_density * (flight_speed_m_s ** 2) * inlet_CD * A_open

        # Net drag 
        net_drag_N = drag_N * (1.0 - meredith_recovery)

        # Store results
        self.inlet_area = A_open
        self.skin_area = skin_area
        self.hole_count = hole_count
        self.drag_N = drag_N
        self.net_drag_N = net_drag_N

        return {
            'm_dot_air': m_dot_air,
            'A_open': A_open,
            'skin_area': skin_area,
            'hole_count': hole_count,
            'drag_N': drag_N,
            'net_drag_N': net_drag_N
        }

class ThermoElectricGenerator(HeatSink):
    def __init__(self, hot_temp_K, cold_temp_K, efficiency):
        super().__init__("TEG")
        self.efficiency = efficiency
        self.hot_temp = hot_temp_K
        self.cold_temp = cold_temp_K
    def absorb_heat(self, heat_w):
        Q_absorbed = heat_w 
        P_out = self.efficiency * Q_absorbed
        
        Q_rejected = Q_absorbed - P_out
        print(f"TEG: Absorbed {Q_absorbed:.2f} W, Power output {P_out:.2f} W, Heat rejected {Q_rejected:.2f} W")
        return P_out  # return the heat that still needs disposal


# -------------------------------- Additional Components -------------------------------
class Fan():
    def __init__(self, fan_eff, ra_mf, ra_density, delta_pressure):
        self.fan_eff = fan_eff
        self.ra_mf = ra_mf
        self.ra_density = ra_density
        self.delta_pressure = delta_pressure
    def power(self):
        return (self.ra_mf * self.delta_pressure) / (self.fan_eff * self.ra_density)
    
class Pump(): # Only for liquid coolant systems
    def __init__(self, pump_eff, coolant_mf, coolant_density, delta_pressure):
        self.pump_eff = pump_eff
        self.coolant_mf = coolant_mf
        self.coolant_density = coolant_density
        self.delta_pressure = delta_pressure
    def power(self):
        return (self.coolant_mf * self.delta_pressure) / (self.pump_eff * self.coolant_density)

class HX():
    def __init__(self, name, fluid_cold, fluid_hot, Q_req):
        self.name = name
        self.fluid_cold = fluid_cold
        self.fluid_hot = fluid_hot
        self.Q_req = Q_req

    def size(self):
        # Plate propeties
        plate_thickness = 0.6 * 1e-3 #m
        plate_thermal_conductivity = 17.5 # W/(m·K), SS
        size_factor = 1.15 #(1.15-1.25)
        gap_bt_plates = 3e-3
        N_plates = 8


        # Pipe Properties
        N_passes = 1
        dh_coolant = 2*gap_bt_plates / size_factor
        dh_h2 = dh_coolant
        


        #initial guess for mass flow rate of coolant
        self.fluid_hot.mf_calculated = self.Q_req / (self.fluid_hot.cp * (self.fluid_hot.T - self.fluid_cold.T))
        
        area = 2.3 # assumed area of heat exchanger plate in m²
        H_hx = 500 # Overall heat exchange coefficient for HX [W/m².K]
        self.fluid_cold.cp = 9384
        # Calculation
        Q = 1000
        while abs(Q-self.Q_req) > 1e-5:
        #for i in range(5):

            
            t_hot_out = self.Q_req / (self.fluid_hot.mf_calculated * self.fluid_hot.cp) + self.fluid_hot.T


            t_hot_in = self.fluid_hot.T
            t_cold_in = self.fluid_cold.T
            #t_cold_out = self.fluid_hot.T - 5
            t_cold_out = t_cold_in + 445/9.384
            delta_t1 = t_hot_out - t_cold_in
            delta_t2 = t_hot_in - t_cold_out
            
            lmtd = (delta_t2 - delta_t1) / math.log(delta_t2 / delta_t1)
            S_eff = self.Q_req / (lmtd * H_hx)
            #s_eff = area * size_factor #assumed area of heat exchanger plate
            #N_plates = math.ceil(S_eff / s_eff)
            s_eff = (area - 2* math.pi * (dh_coolant/2)**2 ) * size_factor
            N_plates = min(math.ceil(S_eff / s_eff),8)
            n_channels = math.ceil((N_plates-1) / 2)
            volume = N_plates * area * (plate_thickness+gap_bt_plates)  # m³
            
            # Coolant
            C_h_coolant = 0.348
            y_coolant = 0.663
            R_fh = 0.0003526
            Pr_coolant = self.fluid_hot.mu * self.fluid_hot.cp / self.fluid_hot.k
            Re_coolant = self.fluid_hot.mf_calculated * 4/  (math.pi * dh_coolant * self.fluid_hot.mu) 
            Nu_coolant = C_h_coolant * Re_coolant**y_coolant * Pr_coolant**0.33 * 1 # Assuming (mu/mu_w)**0.17 = 1, given the tubes are very narrow
            hc_hot = Nu_coolant * self.fluid_hot.k / dh_coolant

            # H2
            C_h_h2 = 0.348  # Placeholder
            y_h2 = 0.663 # Placeholder
            R_fc = 8.815 * 1e-5 # for vapour
            Pr_h2 = self.fluid_cold.mu * self.fluid_cold.cp / self.fluid_cold.k  # Prandtl number for H2
            Re_h2 = self.fluid_cold.mf_calculated * 4/  (math.pi * dh_h2 * self.fluid_cold.mu) 
            Nu_h2 = C_h_h2 * Re_h2**y_h2 * Pr_h2**0.33 * 1 # Assuming (mu/mu_w)**0.17 = 1, given the tubes are very narrow
            hc_cold = Nu_h2 * self.fluid_cold.k / dh_h2


            # Recalculate overall heat transfer coefficient
            H_hx = 1/ (1/hc_cold + 1/hc_hot + R_fc + R_fh + plate_thickness/plate_thermal_conductivity)

            Q = H_hx * s_eff * N_plates * lmtd  * N_passes
            print("Try Q:{H_hx*}")
            factor = self.Q_req / Q 
            print(f"Iteration: Q = {Q:.2f} W, Required Q = {self.Q_req:.2f} W, Factor = {factor:.2f}")
            print(f"Iterarion: Area = {area:.2f} m², Volume = {volume:.2f} m³, N_plates = {N_plates}, n_channels = {n_channels}")
            hole_area = math.pi * (dh_coolant/2)**2
            ratio =  hole_area / area
            print(f"Pipe diameter: {dh_coolant:.3f} m")
            print(f"Coolant mass flow rate: {self.fluid_hot.mf_calculated:.2f} kg/s")
            print(f"Total Surface Required: {S_eff:.2f} m²") 

            area *= factor            

        return self.fluid_hot.mf_calculated, area, volume
        
class Compressor():
    def __init__(self, comp_efficiency, pressure_ratio):
        self.efficiency = comp_efficiency
        self.pressure_ratio = pressure_ratio

    def power(self, mass_flow_rate, inlet_temperature, gamma):
        # Calculate the power required for the compressor
        H2_molar_mass = 2.016  # kg/kmol for hydrogen
        R = 8314 / H2_molar_mass  # Specific gas constant for hydrogen in J/(kg*K)
        T_out = inlet_temperature * (self.pressure_ratio ** ((gamma - 1) / gamma))  # Isentropic relation

        '''SOMEBODY CHANGE PRESSURE OF HYDROGEN FROM AMBIENT TO "BEFORE COMPRESSOR C1" TO GET PROPER cp ''' 
        
        cp_hydrogen = PropsSI('C', 'T', inlet_temperature, 'P', 101325, 'Hydrogen')  # Specific heat capacity at constant pressure
        # power = (mass_flow_rate * R * inlet_temperature * (T_out - inlet_temperature)) / self.efficiency
        power = mass_flow_rate * cp_hydrogen * (T_out - inlet_temperature) / self.efficiency  # Power in Watts
        return power

class Fluid():
    def __init__(self, name, T, P, C,cp,mf,k, mu):
        self.name = name
        self.T = T  # Temperature in Kelvin
        self.P = P  # Pressure in Pascals
        self.C = C  # Heat capacity in J/(kg*K)
        self.cp = cp
        self.mf_calculated = C/cp
        self.mf_given = mf
        self.k = k
        self.mu = mu

# ------------------------------ FUNCTIONS -----------------------------------

def size_thermal_management(heat_sources, sinks):
    total_heat = sum(src.waste_heat() for src in heat_sources)
    print(f"Total waste heat = {total_heat/1000:.1f} kW")

    remaining_heat = total_heat
    for sink in sinks:
        if remaining_heat <= 0:
            break
        if isinstance(sink, SkinHeatExchanger):
            sink.set_ambient(T_ambient_K=ambient_conditions['T'], mach=ambient_conditions['Mach'], recovery_factor=recovery_factor)
        if not isinstance(sink, RamAirHeatExchanger):
            absorbed = sink.absorb_heat(remaining_heat)
            print(f"{sink.name} absorbed {absorbed/1000:.1f} kW")
            remaining_heat -= absorbed

    if remaining_heat > 0:
        # assume last sink is RamAirHX
        for sink in sinks:
            if isinstance(sink, RamAirHeatExchanger):
                
                area = sink.size_exchanger(heat_w=remaining_heat,T_air_in_K=T_air_in,T_air_out_K=T_air_out)
                print(f"Radiator face area required: {area:.2f} m²")
                sizing = sink.compute_hole_and_skin_area(
                    heat_w=remaining_heat,
                    T_air_in_K=T_air_in,
                    T_air_out_K=T_air_out,
                    flight_speed_m_s=ambient_conditions['V'],
                    air_density=ambient_conditions['rho'],
                    porosity=porosity,
                    hole_diameter_mm=hole_diameter_mm,
                    inlet_CD=inlet_cd,       
                    cp_air=cp_air,
                    meredith_recovery=meredith_recovery
                )
                print(f"Air mass flow required: {sizing['m_dot_air']:.2f} kg/s")
                print(f"Total open hole area: {sizing['A_open']:.3f} m²")
                print(f"Total skin area (10% porosity): {sizing['skin_area']:.2f} m²")
                print(f"Number of {hole_diameter_mm:.1f} mm holes: {sizing['hole_count']}")
                print(f"Raw inlet drag: {sizing['drag_N']:.1f} N")
                print(f"Net drag after Meredith (50% recovery): {sizing['net_drag_N']:.1f} N")

                fan = Fan(
                    fan_eff=fan_eff,
                    ra_mf=sizing['m_dot_air'],
                    ra_density=ambient_conditions['rho'],
                    delta_pressure=delta_pressure
                )
                fan_power = fan.power()
                print(f"Fan power required during Take-off: {fan_power/1000:.2f} kW")

                remaining_heat = 0
    if remaining_heat > 0:
        print(f"WARNING: {remaining_heat/1000:.1f} kW of heat could not be dissipated!")



# ------------------------------ MAIN PROGRAM -----------------------------------


# Heat Sources
fuel_cell = FuelCell(power_output_w=fc_power, efficiency=fc_eff)
battery = Battery(power_output_w=bat_p, heat_fraction=battery_heat_fraction)
motor = ElectricMotor(power_output_w=em_power, efficiency=em_eff)

heat_sources = [fuel_cell, battery, motor]

# Heat Sinks
lh2_tank = LH2Tank(mass_flow_kg_s=(2*h2_mf_fc*fc_power/1e6),T_final_K=t_fin_fc,P_final=p_fin_fc,P_init=p_init_fc) # for both wings (FC)
lh2_tank_1 = LH2Tank(mass_flow_kg_s=2*h2_mf_cc,T_final_K=t_fin_cc, P_final = p_fin_cc, P_init = p_init_cc) # for 2 comb. (CC)
skin_hx = SkinHeatExchanger(area_m2=area_wing, U_W_per_m2K=U_wing, coolant_temp_K=wing_coolant_temp) # for 2 wings
ram_air_hx = RamAirHeatExchanger(U_W_per_m2K=U_ra, coolant_temp_K=ra_coolant_temp)
teg = ThermoElectricGenerator(hot_temp_K=wing_coolant_temp, cold_temp_K=ambient_conditions['T'], efficiency=efficiency_teg)


sinks = [lh2_tank, lh2_tank_1,teg, skin_hx, ram_air_hx]


# RUN -----------------------------------------------------------
if __name__ == "__main__":
    size_thermal_management(heat_sources, sinks)



    # Sizing
    C_h2 = lh2_tank.absorb_heat(fuel_cell.waste_heat())/(t_fin_fc - t_init)
    print(f"Capacity rate of LH2: {C_h2:.2f})") 
    fluid_cold = Fluid(name="LH2", T=t_init, P=p_init_fc, C=C_h2, cp = C_h2/h2_mf_fc, k = k_h2, mu = mu_h2, mf = h2_mf_fc)
    fluid_hot = Fluid(name="FuelCellCoolant", T=t_fin_fc, P=p_fin_fc, C=0, cp = h_cool, mf = 0, k=k_coolant, mu=mu_coolant) 

    hx = HX(name="FuelCellHX", fluid_cold=fluid_cold, fluid_hot=fluid_hot, Q_req=lh2_tank.absorb_heat(fuel_cell.waste_heat()))
    hx_cool_flow_mass, hx_area, hx_volume = hx.size()
    print(f"Coolant Mass Flow: {hx_cool_flow_mass:.2f} kg/s")
    print(f"Heat Exchanger Area: {hx_area:.2f} m²")
    print(f"Heat Exchanger Volume: {hx_volume:.2f} m³")
