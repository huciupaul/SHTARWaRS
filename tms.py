import abc
import numpy as np
from CoolProp.CoolProp import get_global_param_string
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
import tms_plotting 
import math

# ------------------------------ PARAMETERS -----------------------------------
# TODO: 
# - Delete all parameters that are not used in the code
# - List all constants clearly and cite from literature


# Ambient Conditions --------------------------------
ambient_conditions= {
    'T': 238.651,  # Ambient temperature in K
    'M': 0.4,  # Flight Mach number
    'V': 123.8644,  # Flight speed in m/s
    'rho': 0.55  # Air density at cruise in kg/m^3
}

# Parameters
gamma = 1.4

# TEG
efficiency_teg = 0.05 

# Ram Air HX -------------------------------
# Air properties
h_air = 250         # REF: Shah 
cp_air = 1005.0
T_air_in = ambient_conditions['T']
T_air_out = ambient_conditions['T'] + 20.0

# Coolant and HX Properties 
ra_coolant_temp = 165 + 20 + 273.15 # K
        

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
area_wing = 2.3

# Air
prandtl_air = 0.71
reynolds_air = 1e7  # depends on temperature?
h_ext_w = ambient_conditions['rho'] * ambient_conditions['V'] * cp_air * 0.185 * (np.log10(reynolds_air))**-2.584 * prandtl_air**-0.66
recovery_factor = prandtl_air ** 0.33  

# Overall heat coefficients WING and RADIATOR
#U_wing = 1/(1/h_cool + 1/h_ext_w)
#U_ra = 1 / (1/h_air + 1/h_cool)

deltaT_fc = 20
HEX_1_deltaT = 20 # Assumed temperature increase during evaporation


# Main Classes for Heat Sinks (TMS) -------------------------------
    
class SkinHeatExchanger():
    def __init__(self, area, U, fluid):
        super().__init__()
        self.area = area
        self.U = U             # overall HT coefficient (W/m^2.K)
        self.coolant_temp = fluid.T
        self.ambient_temp = None          # to be set based on flight condition
        self.fluid = fluid                

    def set_ambient(self, T_ambient_K, recovery_factor=0.9, mach=0.3):
        # Set effective ambient (adiabatic wall) temperature for convection calculations
        T_aw = T_ambient_K * (1 + recovery_factor*0.5*(1.4-1)*mach**2)
        self.ambient_temp = T_aw 
    def absorb_heat(self, heat_w):
        if self.ambient_temp is None:
            raise RuntimeError("Ambient conditions not set for SkinHeatExchanger")

        Q_capacity = self.U * self.area * (self.coolant_temp - self.ambient_temp)
        Q_absorbed = min(heat_w, Q_capacity)
        t_coolant_out = self.coolant_temp - Q_absorbed / (self.fluid.cp * self.fluid.mf_given)
        if Q_absorbed < heat_w:
            return Q_absorbed, t_coolant_out
        else:
            print("No need for Radiator")
            return Q_absorbed, t_coolant_out
    
class RamAirHeatExchanger():
    def __init__(self, coolant_temp_K, fluid):
        #self.U = U_W_per_m2K
        self.coolant_temp = coolant_temp_K
        self.Pr = fluid.Pr
        self.dyn_visc = fluid.dyn_visc
        self.mf_coolant = fluid.mf_given
        self.k_cool = fluid.k

        # To be determined during sizing:
        self.U_ra = None
        self.required_area = 0.0    # m² of radiator face
        self.inlet_area = 0.0       # m² of total open hole area
        self.skin_area = 0.0        # m² of skin needed (porous area)
        self.hole_count = 0         # total number of small holes
        self.TR = 0.0
        #self.drag_N = 0.0           # N, raw drag before Meredith effect
        self.net_drag_N = 0.0       # N, after Meredith recovery
    
    def U_rad(self):

        # PRINT fluid parameters
        #print("k_cool:", self.k_cool)
        #print("Pr:", self.Pr)
        #print("mf_coolant:", self.mf_coolant)
        #print("dyn_visc:", self.dyn_visc)

        # Coolant constants
        dh = 0.002    # 4*A/P  hydraulic diameter

        # Reynolds calc
        Re = 4*self.mf_coolant/(np.pi * dh * self.dyn_visc)
        if Re < 2300:
            raise ValueError("Reynolds number in laminar regime, use a turbulent/combined relation correlation.")
        
        # Nusselt calc
        f = (0.79 * np.log(Re) - 1.64) ** -2
        nu = 0.023*Re**(4/5)*self.Pr**0.3        #nu = ((f/8)*(Re-1000)*self.Pr)/(1+12.7*(f/8)**0.5 * (self.Pr**0.66-1)) # Nusselt number for turbulent flow in tubes [m²/s]
        #print("Nusselt number of the coolant:", nu)
        h_cool = 1500#self.k_cool * nu / dh
        wing_coolant_temp = 320  # K
        #print("Heat transfer coefficient of the coolant:", h_cool)

        self.U_ra = 1 / (1/h_air + 1/h_cool)
        #print("Overall heat transfer coefficient", self.U_ra)

        return self.U_ra


    def absorb_heat(self, heat_w):
        return heat_w

    def size_exchanger(self, heat_w, T_air_in_K):
        self.U_rad()

        # Iteration for delta T
        deltaT_grid = np.linspace(1, 150, 100 )        
        tol = 0.1                                       # “almost 1” → ±1 %

        for T in deltaT_grid:
            print(T)

            dT1 = self.coolant_temp - T_air_in_K                      
            dT2 = self.coolant_temp - (T_air_in_K + T)                

            dT_lm = (dT2 - dT1) / np.log(dT2 / dT1)

            self.required_area = heat_w / (self.U_ra * dT_lm)
            beta  = 800
            vol_rad = self.required_area/beta
            front_area_rad = vol_rad/0.05
            length_rad = np.sqrt(front_area_rad)
            mf_air = front_area_rad * ambient_conditions['V'] * ambient_conditions['rho']

            # --- ratio we want to drive to 1 -----------------------------------------
            ratio = (mf_air * cp_air * T) / (self.U_ra * dT_lm)

            if abs(ratio - 1.0) < tol:           # close enough
                print(f"Converged: ΔT_air = {T:.3f} K  (ratio = {ratio:.4f})")
                print("Heat area radiator", self.required_area, "W/m²K")
                print("Length of the radiator", length_rad, "m")
                print("Volume of the radiator", vol_rad, "m³")
                break
            

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
            area_per_hole = np.pi * (d_m / 2) ** 2 # m² per hole
            hole_count = int(np.ceil(A_open / area_per_hole))
        else:
            hole_count = 0
            d_m = 0.0
            area_per_hole = 0.0

        # Total gross drag
        mdot_air = ambient_conditions['rho'] * ambient_conditions['V'] * A_open
        D_g = mdot_air * ambient_conditions['V']

        #drag_N = 0.5 * air_density * (flight_speed_m_s ** 2) * inlet_CD * A_open

        # Net drag 
        #net_drag_N = TR * (1.0 - meredith_recovery)

        # Store results
        self.inlet_area = A_open
        self.skin_area = skin_area
        self.hole_count = hole_count
        self.D_g = D_g
        #self.drag_N = drag_N
        #self.net_drag_N = net_drag_N

        return {
            'm_dot_air': m_dot_air,
            'A_open': A_open,
            'skin_area': skin_area,
            'hole_count': hole_count
            #'drag_N': drag_N,
            #'net_drag_N': net_drag_N
        }

    


    # THRUST RECOVERY
    def thrust_ratio(
        deltaT_HX0,                 # can be a scalar or a NumPy array
        gamma: float,
        R: float,
        M0: float,
        T0: float,
        eta_p07: float ,
        CD_d_CD_sp: float
    ) -> np.ndarray:
        """
        Vectorised thrust-ratio calculation.

        Parameters
        ----------
        deltaT_HX0 : float or array-like
            Heat-exchanger temperature rise(s) [K].
        gamma      : float
            Ratio of specific heats.
        R          : float
            Gas constant [J kg⁻¹ K⁻¹].
        M0         : float
            Flight Mach number.
        T0         : float
            Ambient static temperature [K].
        eta_p07    : float
            Propulsive efficiency term.
        CD_d_CD_sp : float
            Sum (C_D,d + C_D,sp).

        Returns
        -------
        TR : ndarray
            Thrust ratio(s) corresponding to `deltaT_HX0`.
        """
        deltaT_HX0 = np.asarray(deltaT_HX0, dtype=float)

        comp_ratio = (1.0 + 0.5 * (gamma - 1.0) * M0**2) * eta_p07 ** ((gamma - 1.0) / gamma)

        # Guard against unphysical D <= 1
        term2 = 1.0 - 1.0 / comp_ratio

        numerator = (2.0 * gamma / (gamma - 1.0) * R) * term2 * ( T0 * (1.0 + 0.5 * (gamma - 1.0) * M0**2) + deltaT_HX0)

        denominator = M0 * np.sqrt(gamma * R * T0) * (1.0 + CD_d_CD_sp / 2.0)
        
        TR = np.sqrt(numerator) / denominator - 1.0
        #TR = np.where(value > 0, np.sqrt(value) - 1.0, np.nan)  # nan for impossible points
        return TR
    

    # -- sweep ΔT_HX0 from 0 to 110 K --------------------------------------
    deltaT = np.linspace(0, 100, 100)
    eta_values = [0.995, 0.92, 0.90, 0.85]

    for eta in eta_values:
        TR = thrust_ratio(
            deltaT,
            gamma = 1.4,
            R     = 287.058,
            M0    = ambient_conditions["M"],
            T0    = ambient_conditions["T"],
            eta_p07       = eta,   # optional override
            CD_d_CD_sp    = 0.11    # optional override
        )
        # plt.plot(deltaT, TR, label=f"ηₚ,07 = {eta:.3f}")

    # print(TR[-1])
    # -- plot ---------------------------------------------------------------
    #plt.plot(deltaT, TR, marker="o")
    '''
    plt.xlabel("ΔT_HX⁰ (K)")
    plt.ylabel("Thrust Ratio (TR)")
    plt.title("Temperature difference of Radiator vs TR (for different pressure ratios)")
    plt.grid(True)
    plt.legend()
    plt.show()
    '''


class ThermoElectricGenerator():
    def __init__(self, hot_temp_K, cold_temp_K, efficiency):
        super().__init__("TEG")
        self.efficiency = efficiency
        self.hot_temp = hot_temp_K
        self.cold_temp = cold_temp_K
    def absorb_heat(self, heat_w):
        Q_absorbed = heat_w 
        P_out = self.efficiency * Q_absorbed
        
        Q_rejected = Q_absorbed - P_out
        # print(f"TEG: Absorbed {Q_absorbed:.2f} W, Power output {P_out:.2f} W, Heat rejected {Q_rejected:.2f} W")
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
    
class Pump(): 
    def __init__(self, fluid, pump_eff, delta_pressure):
        self.pump_eff = pump_eff
        self.coolant_mf = fluid.mf_given
        self.coolant_density = fluid.rho
        self.delta_pressure = delta_pressure
    def power(self):
        return (self.coolant_mf * self.delta_pressure) / (self.pump_eff * self.coolant_density)

class HEX():
    def __init__(self, name, fluid_cold, fluid_hot, Q_req):
        self.name = name
        self.fluid_cold = fluid_cold
        self.fluid_hot = fluid_hot
        self.Q_req = Q_req

    def size(self):
        # Plate properties
        plate_thickness = 3.6 * 1e-3 #m 
        plate_thermal_conductivity = 17.5 # W/(m·K), SS
        size_factor = 1.15 #(1.15-1.25)
        gap_bt_plates = 3e-3  # from excel tool
        N_plates = 10          # from Michelle 

        # Pipe Properties
        N_passes = 1
        dh_coolant = 2*gap_bt_plates / size_factor  # diameter of coolant pipes
        dh_h2 = dh_coolant                          # diameter of hydrogen pipes
        

        #initial guess for mass flow rate of coolant
        self.fluid_hot.mf_calculated = 10 # self.Q_req / (self.fluid_hot.cp * (self.fluid_hot.T - self.fluid_cold.T))
        
        L_h2 = 447000 # J/kg
        area = 1.3 # assumed area of heat exchanger plate in m² (total)
        H_hx = 500 # Overall heat exchange coefficient for HX [W/m².K]
        self.fluid_cold.cp = 9384

        # Calculation
        Q = 1000
        iteration = 0
        while abs(Q-self.Q_req) > 1e-5:
        #for i in range(5):

            t_hot_out = self.Q_req / (self.fluid_hot.mf_calculated * self.fluid_hot.cp) + self.fluid_hot.T


            t_hot_in = self.fluid_hot.T
            t_cold_in = self.fluid_cold.T
            #t_cold_out = self.fluid_hot.T - 5
            t_cold_out = t_cold_in + (self.fluid_cold.mf_given/self.fluid_hot.mf_given )*(L_h2/self.fluid_hot.cp)
            delta_t1 = t_hot_out - t_cold_in
            delta_t2 = t_hot_in - t_cold_out
            
            lmtd = (delta_t2 - delta_t1) / np.log(delta_t2 / delta_t1)
            S_eff = self.Q_req / (lmtd * H_hx)
            #s_eff = area * size_factor #assumed area of heat exchanger plate
            #N_plates = np.ceil(S_eff / s_eff)
            s_eff = (area - 4* np.pi * (dh_coolant/2)**2 ) * size_factor  # TOTAL EFF AREA
            #N_plates = min(np.ceil(S_eff / s_eff),15)
            n_channels = np.ceil((N_plates-1) / 2)
            volume = N_plates * area * (plate_thickness+gap_bt_plates)  # m³
            
            # Coolant
            C_h_coolant = 0.348 # Chevron angle of 30deg
            y_coolant = 0.663   # Chevron angle of 30deg
            R_fh = 0.0003526
            Pr_coolant = self.fluid_hot.mu * self.fluid_hot.cp / self.fluid_hot.k
            Re_coolant = self.fluid_hot.mf_calculated * 4/  (np.pi * dh_coolant * self.fluid_hot.mu) 
            Nu_coolant = C_h_coolant * Re_coolant**y_coolant * Pr_coolant**0.33 * 1 # Assuming (mu/mu_w)**0.17 = 1, given the tubes are very narrow
            hc_hot = Nu_coolant * self.fluid_hot.k / dh_coolant

            # H2
            C_h_h2 = 0.348  # Chevron angle of 30deg
            y_h2 = 0.663    # Chevron angle of 30deg
            R_fc = 8.815 * 1e-5 # for vapour
            Pr_h2 = self.fluid_cold.Pr  # Prandtl number for H2
            Re_h2 = self.fluid_cold.mf_calculated * 4 / (np.pi * dh_h2 * self.fluid_cold.mu) 
            Nu_h2 = C_h_h2 * Re_h2**y_h2 * Pr_h2**0.33 * 1 # Assuming (mu/mu_w)**0.17 = 1, given the tubes are very narrow
            hc_cold = Nu_h2 * self.fluid_cold.k / dh_h2


            # Recalculate overall heat transfer coefficient
            H_hx = 1/ (1/hc_cold + 1/hc_hot + R_fc + R_fh + plate_thickness/plate_thermal_conductivity)

            # Recalculate the heat transferred
            Q = H_hx * s_eff * lmtd * N_plates  * N_passes
            # print("Try Q:{H_hx*}")
            factor = self.Q_req / Q 
            # print(f"Iteration: Q = {Q:.2f} W, Required Q = {self.Q_req:.2f} W, Factor = {factor:.2f}")
            # print(f"Iteration: Area = {area:.2f} m², Volume = {volume:.2f} m³, N_plates = {N_plates}, n_channels = {n_channels}")
            hole_area = np.pi * (dh_coolant/2)**2
            ratio =  hole_area / area
            # print(f"Pipe diameter: {dh_coolant:.3f} m")
            # print(f"Coolant mass flow rate: {self.fluid_hot.mf_calculated:.2f} kg/s")
            # print(f"Total Surface Required: {S_eff:.2f} m²") 

            area *= factor    
            iteration += 1

        # TODO: Add pressure drop -------------------- REDO FUNCTION---------------------
        

        self.fluid_hot.T = t_hot_in
        cool_in = self.fluid_hot
        cool_in.mf_given *= 0.5
        self.fluid_hot.T = t_hot_out
        cool_out = cool_in
        self.fluid_cold.T = t_cold_out
        cold_out = self.fluid_cold

        return cool_in, cool_out, cold_out, area, volume
        
        
        

class Compressor():
    def __init__(self, comp_efficiency, pressure_ratio, fluid):
        self.efficiency = comp_efficiency
        self.pressure_ratio = pressure_ratio
        self.fluid = fluid  

    def power(self):
        inlet_temperature = self.fluid.T  
        gamma = self.fluid.gamma 
        T_out = inlet_temperature * (1 + 1 / self.efficiency * (self.pressure_ratio ** ((gamma - 1) / gamma) - 1))  # Isentropic relation

        cp = self.fluid.cp  
        mdot = self.fluid.mf_given
        power_req = mdot * cp * (T_out - inlet_temperature) / self.efficiency  # Power in Watts
        return power_req
    
    def mass(self, power):
        return power * 0.0400683 + 5.17242

class Turbine():
    def __init__(self, turbine_efficiency, fluid):
        self.efficiency = turbine_efficiency
        self.fluid = fluid

    def power(self, T_in, T_out):
        # Calculate the power produced by the turbine
        cp_gas = self.fluid.cp
        mdot = self.fluid.mf_given
        power_provide = mdot * cp_gas * (T_in - T_out) / self.efficiency  # Power in Watts\
        return power_provide
    
    def mass(self, power):
        # MASSIVE assumption of turbine mass being calculated in the same way as compressor
        return power * 0.0400683 + 5.17242

class Valve():
    def __init__(self, fluid, valve_efficiency = 0.9 ):
        self.efficiency = valve_efficiency
        self.fluid = fluid

    def valve_mass(self):
        mdot_coolant = self.fluid.mf_given
        return (0.568 * mdot_coolant ** 0.5541)

class Fluid():
    def __init__(self, name, T, P, mf, fluid_type, C = None):
        sat_temp = 0
        gas = True
        if P < PropsSI('PCRIT', fluid_type):
            sat_temp = PropsSI('T', 'P', P, 'Q', 0, fluid_type)  # Saturation temperature at given pressure
        if abs(T - sat_temp) < 1e-3:
            gas = False
        self.name = name
        self.T = T  # Temperature in Kelvin
        self.P = P  # Pressure in Pascals
        self.cp = PropsSI('C', 'P', P, 'T', T, fluid_type) if gas else PropsSI('C', 'P', P, 'Q', 0, fluid_type)  # Specific heat capacity at constant pressure
        self.C = C if C is not None else self.cp * mf  # Heat capacity rate
        self.mf_calculated = self.C/self.cp
        self.mf_given = mf
        self.k = PropsSI('CONDUCTIVITY', 'P', P, 'T', T, fluid_type) if gas else PropsSI('CONDUCTIVITY', 'P', P, 'Q', 0, fluid_type)  # Thermal conductivity of LH2 at storage pressure
        self.dyn_visc = PropsSI('VISCOSITY', 'P', P, 'T', T, fluid_type) if gas else PropsSI('VISCOSITY', 'P', P, 'Q', 0, fluid_type) # Viscosity of LH2 at storage pressure
        self.gamma = PropsSI('CP0MASS', 'P', P, 'T', T, fluid_type) / PropsSI('CVMASS', 'P', P, 'T', T, fluid_type) if gas else PropsSI('CP0MASS', 'P', P, 'Q', 0, fluid_type) / PropsSI('CVMASS', 'P', P, 'Q', 0, fluid_type) # Heat capacity ratio of LH2 at storage pressure
        self.Pr = PropsSI('PRANDTL', 'P', P, 'T', T, fluid_type) if gas else PropsSI('PRANDTL', 'P', P, 'Q', 0, fluid_type) # Prandtl number of LH2 at storage pressure
        self.rho = PropsSI('D', 'P', P, 'T', T, fluid_type) if gas else PropsSI('D', 'P', P, 'Q', 0, fluid_type) # Density of LH2 at storage pressure
        self.vdot = mf / self.rho  # Volumetric flow rate in m³/s
        self.kin_visc = self.dyn_visc / self.rho  # Kinematic viscosity of LH2 at storage pressure
        self.fluid_type = fluid_type
        self.mu = self.dyn_visc  # Dynamic viscosity of the fluid
    
class Pipe():
    def __init__(self, length, d_in, fluid, type = 'normal'):
        """Initialize the Heat Pipe Analyzer"""
        self.mat_density = 7900
        self.length = length
        self.d_in = d_in
        self.d_out = self.d_in * 1.1
        self.d_is = self.d_out * 1.1
        self.temp_fluid = fluid.T
        self.temp_air = 273.15 + 20
        self.therm_conduc_fluid = fluid.k
        self.prandtl_number = fluid.Pr
        self.mass_flow = fluid.mf_given
        self.spec_heat = fluid.cp
        self.heat_trans_air = 10
        self.therm_conduc_pipe = 15
        self.therm_conduc_insulation = 15
        self.density_fluid = fluid.rho
        self.kinematic_visc = fluid.kin_visc
        self.material_roughness = 0.004e-3
        self.fluid = fluid  # Store the fluid object for later use

    def analyze_heat_pipe(self):
        """
        Combined analysis of heat pipe: calculate both pressure drop and temperature distribution

        Parameters:
        length: pipe length
        d_in: inner diameter
        d_out: outer diameter
        d_is: insulation outer diameter
        temp_fluid: initial fluid temperature
        temp_air: ambient air temperature
        therm_conduc_fluid: thermal conductivity of fluid
        prandtl_number: Prandtl number
        mass_flow: mass flow rate
        spec_heat: specific heat of fluid
        heat_trans_air: heat transfer coefficient air
        therm_conduc_pipe: thermal conductivity of pipe material
        therm_conduc_insulation: thermal conductivity of insulation material
        density_fluid: fluid density
        kinematic_visc: kinematic viscosity of fluid
        material_roughness: absolute roughness of pipe material

        Returns:
        dict: {'pressure_drop': float, 'temperature_profile': np.array, 'reynolds': float}
        """
        # Calculate fluid velocity
        velocity_fluid = self.mass_flow / (self.density_fluid * ((self.d_in / 2) ** 2) * np.pi)

        # Calculate Reynolds number
        reynolds = velocity_fluid * self.d_in / self.kinematic_visc
        

        # ==================== PRESSURE DROP CALCULATION ====================
        # Check for laminar or turbulent flow and calculate Darcy-Weisbach coefficient
        if reynolds < 2300:
            darcy_weisbach = 64 / reynolds
        else:
            x_darcy = (-2.457 * np.log((7 / reynolds) ** 0.9 + 0.27 * self.material_roughness / self.d_in)) ** 16
            y_darcy = (37530 / reynolds) ** 16
            darcy_weisbach = 8 * ((8 / reynolds) ** 12 + (x_darcy + y_darcy) ** -1.5) ** (1 / 12)

        # Calculate pressure drop at the end of the pipe
        pressure_drop = self.length * (darcy_weisbach * self.density_fluid * velocity_fluid ** 2) / (self.d_in * 2)

        # ==================== TEMPERATURE PROFILE CALCULATION ====================
        # Initialize temperature array
        temperature = [self.temp_fluid]
        energies = 0
        steps = 100

        # Calculate heat transfer coefficient
        nu_turb = 0.0223 * (reynolds ** 0.8) * self.prandtl_number ** 0.4
        heat_trans_co = nu_turb * self.therm_conduc_fluid / self.d_in

        # Calculate surface areas
        surface_in = np.pi * self.d_in * (self.length)
        surface_out = np.pi * self.d_is * (self.length / steps)

        # Calculate thermal resistances
        r_fluid = (1 / (self.mass_flow * self.spec_heat * (1 - np.exp(-1 * heat_trans_co * surface_in / (self.mass_flow * self.spec_heat))))) * steps
        r_air = 1 / (self.heat_trans_air * surface_out)
        r_wall = np.log(self.d_out / self.d_in) / (2 * np.pi * self.therm_conduc_pipe * (self.length / steps))
        r_insulation = np.log(self.d_is / self.d_out) / (2 * np.pi * self.therm_conduc_insulation * (self.length / steps))
        r_tot = r_fluid + r_air + r_wall + r_insulation


        # Calculate time in pipe
        time_in_pipe = self.length / velocity_fluid


        temps_fluid = self.temp_fluid
        mass_section = (self.length / steps) * (self.d_in / 2) ** 2 * np.pi * self.density_fluid

        # Calculate temperature at each step
        for i in range(steps):
            delta_q = (temps_fluid - self.temp_air) / r_tot
            delta_t = time_in_pipe * delta_q / (mass_section * self.spec_heat) / steps
            temps_fluid = temps_fluid - delta_t
            temperature.append(temps_fluid)
            energies = energies + delta_q
        
        energy_lost = energies
        if type == 'normal':
            self.fluid.P -= pressure_drop
            self.fluid.T = temperature[-1]  # Update fluid temperature to final value
        if type == 'inv':
            self.fluid.P += pressure_drop
            self.fluid.T += (temperature[0]- temperature[-1]) 

        return self.fluid

    def mass(self):
        """
        Calculate the mass of the pipe
        Returns:
        float: mass of the pipe in kg
        """
        mass = self.mat_density * self.length * np.pi * ((self.d_out / 2) ** 2 - (self.d_in / 2) ** 2)
        return mass


def size_pipes_h2(h2_mf_fc, h2_mf_cc, p_sto,fluid,diam_est):
    sf_delta_p = 2e5
    pressure_drop = p_sto -4.511e5 - sf_delta_p 
    total_pressure_drop = 0

    # First half
    length1 = 6.3425
    diam1 = diam_est
    mf1 = (h2_mf_fc + h2_mf_cc)
    v_fluid1 = mf1 / (fluid.rho * ((diam1 / 2) ** 2) * np.pi)
    reynolds1 = v_fluid1 * diam1 / fluid.kin_visc

    if reynolds1 < 2300:
            darcy_weisbach = 64 / reynolds1
    else:
        x_darcy = (-2.457 * np.log((7 / reynolds1) ** 0.9 + 0.27 * 0.004e-3 / diam1)) ** 16
        y_darcy = (37530 / reynolds1) ** 16
        darcy_weisbach = 8 * ((8 / reynolds1) ** 12 + (x_darcy + y_darcy) ** -1.5) ** (1 / 12)
    pressure_drop1 = length1 * (darcy_weisbach * fluid.rho * v_fluid1 ** 2) / (diam1 * 2)

    # Second half
    diam2 = diam_est * h2_mf_cc / (h2_mf_cc + h2_mf_fc) 
    length2 = 3.225
    mf2 = h2_mf_cc
    v_fluid2 = mf2 / (fluid.rho * ((diam2 / 2) ** 2) * np.pi)
    reynolds2 = v_fluid2 * diam2 / fluid.kin_visc

    if reynolds2 < 2300:
            darcy_weisbach = 64 / reynolds2
    else:
        x_darcy = (-2.457 * np.log((7 / reynolds2) ** 0.9 + 0.27 * 0.004e-3 / diam2)) ** 16
        y_darcy = (37530 / reynolds2) ** 16
        darcy_weisbach = 8 * ((8 / reynolds2) ** 12 + (x_darcy + y_darcy) ** -1.5) ** (1 / 12)
    pressure_drop2 = length2 * (darcy_weisbach * fluid.rho * v_fluid2 ** 2) / (diam2 * 2)

    # Total pressure drop
    total_pressure_drop = pressure_drop1 + pressure_drop2
    press_err = total_pressure_drop - pressure_drop
    if press_err < 0:
        print(f"Selected pipe diameter {diam_est:.3f} m is adequate, pressure drop is {total_pressure_drop:.2f} Pa, with an estimated maximum allowed of {pressure_drop:.2f} Pa")
    else: 
        print(f"Increase pipe diameter, pressure drop is {total_pressure_drop:.2f} Pa, with an estimated maximum allowed of {pressure_drop:.2f} Pa")


    
# ------------------------------ MAIN PROGRAM -----------------------------------


def main(Q_dot_fc, Q_dot_eps, p_fc, p_cc, h2_mf_fc, h2_mf_cc, T_fc, T_cc, air_mf_fc, T_amb, rho_amb, V_amb, p_amb, h2_mf_rec, air_out_fc, p_sto):
    p_cool = 5.7e5
    fc_press_drop_cool = 100000
    cool_0 = Fluid(name="FuelCellCoolantGeneric", T=T_fc, P=p_cool, C=None, mf=10, fluid_type='Water')  # Coolant to FC
    cool_mf_per_fc = Q_dot_fc / (cool_0.cp * deltaT_fc * 2)
    cool_0.mf_given = cool_mf_per_fc  # Mass flow rate of coolant to FC


    sources = np.array([Q_dot_fc, Q_dot_eps])
    Q_dot_rem = np.sum(sources)

    h2_mf_fc = h2_mf_fc * 2  # Total H2 mass flow rate to fuel cell (both wings)
    h2_mf_cc = h2_mf_cc * 2  # Total H2 mass flow rate to combustion chamber (both wings)

    # Calculate pipe diameter for H2
    h2_test = Fluid(
        name="H2_Test",
        T=PropsSI('T', 'P', p_sto, 'Q', 1, 'ParaHydrogen'),  # Q=1 for saturated vapor (gaseous H2)
        P=p_sto,
        mf=h2_mf_fc + h2_mf_cc,
        fluid_type='ParaHydrogen'
    )
    diam_est = 0.02
    size_pipes_h2(h2_mf_fc, h2_mf_cc, p_sto,h2_test, diam_est)

    # Initialize -----------------------------
    T_tank = PropsSI('T', 'P', p_sto, 'Q', 0, 'ParaHydrogen')  # Initial temperature of LH2 tank
    h2_1 = Fluid(name="H2_1", T=T_tank, P=p_sto, C = None, mf = h2_mf_fc + h2_mf_cc, fluid_type='ParaHydrogen')

    # Pipe h2 12
    pipe_h2_12 = Pipe(1, 0.05, h2_1)  # 1 m length, 50 mm diameter
    h2_2 = pipe_h2_12.analyze_heat_pipe()

    # Pipe h2 34
    T_sat_h2 = PropsSI('T', 'P', h2_2.P, 'Q', 0, h2_2.fluid_type)  # Saturation temperature of H2 at pressure P
    h2_3 = Fluid(name="H2_3", T=T_sat_h2 + HEX_1_deltaT, P=p_sto, C = None, mf = h2_mf_fc + h2_mf_cc, fluid_type='ParaHydrogen')
    pipe_h2_34 = Pipe(1, 0.05, h2_3)  # 1 m length, 50 mm diameter
    h2_4 = pipe_h2_34.analyze_heat_pipe()

    # Superheater HEX
    cool_19 = Fluid(name="FuelCellCoolant", T=T_fc, P=p_cool-fc_press_drop_cool, C=None, mf = cool_mf_per_fc * 2, fluid_type='Water') 
    Q_h2_heat = (PropsSI('H', 'P', h2_2.P, 'T', T_fc, h2_2.fluid_type) - PropsSI('H', 'P', h2_2.P, 'T', T_sat_h2 + HEX_1_deltaT, h2_2.fluid_type)) * h2_2.mf_given   # Heat of vaporization
    hx_heat = HEX(name="Vaporizer", fluid_cold=h2_4, fluid_hot=cool_19, Q_req=Q_h2_heat)
    cool_19,cool_20,h2_5, hx_heat_area, hx_heat_volume = hx_heat.size()
    Q_dot_rem -= Q_h2_heat 

    # Pipe cool 20-1
    pipe_cool_20_1 = Pipe(1, 0.05, cool_20)  # 1 m length, 50 mm diameter
    cool_1 = pipe_cool_20_1.analyze_heat_pipe()

    # Vaporizer HEX
    Q_vap = (PropsSI('H', 'P', h2_2.P, 'T', T_sat_h2 + HEX_1_deltaT, h2_2.fluid_type) - PropsSI('H', 'P', h2_2.P, 'Q', 0, h2_2.fluid_type)) * h2_2.mf_given   # Heat of vaporization
    hx_vap = HEX(name="Vaporizer", fluid_cold=h2_2, fluid_hot=cool_1, Q_req=Q_vap)
    cool_1,cool_2,h2_3, hx_vap_area, hx_vap_volume = hx_vap.size()
    Q_dot_rem -= Q_vap

    # Pipe h2 5-6
    pipe_h2_5_6 = Pipe(1, 0.05, h2_5)  # 1 m length, 50 mm diameter
    h2_6 = pipe_h2_5_6.analyze_heat_pipe()

    # Intersection 6-11-7
    h2_7 = Fluid(name="H2_7", T=h2_6.T, P=h2_6.P, C=0, mf=h2_mf_fc, fluid_type='ParaHydrogen')  # H2 to fuel cell
    h2_11 = Fluid(name="H2_7", T=h2_6.T, P=h2_6.P, C=0, mf=h2_mf_cc, fluid_type='ParaHydrogen')  # H2 to CC
    valve_6711 = Valve(fluid=h2_6, valve_efficiency=0.9)  # Valve for H2 to fuel cell
    valve_6711_mass = valve_6711.valve_mass()  # Mass of the valve

    # PATH to CC ----------------------------------------------
    # Pipe h2 11-12
    pipe_h2_11_12 = Pipe(1, 0.05, h2_11)  # 1 m length, 50 mm diameter
    h2_12 = pipe_h2_11_12.analyze_heat_pipe()

    # Intersection 12-13-26
    comp_14_15_PI = 5
    split_13_26 = 1 - ((p_cc+1)/h2_12.P - comp_14_15_PI)/(h2_12.P - comp_14_15_PI)
    h2_13 = Fluid(name="H2_13", T=h2_12.T, P=h2_12.P, C=0, mf=h2_mf_cc * split_13_26, fluid_type='ParaHydrogen')  # H2 to Compressor
    h2_26 = Fluid(name="H2_26", T=h2_12.T, P=h2_12.P, C=0, mf=h2_mf_cc * (1-split_13_26), fluid_type='ParaHydrogen')  # H2 around Compressor
    valve_121326 = Valve(fluid=h2_12, valve_efficiency=0.9)  # Valve for H2 to CC
    valve_121326_mass = valve_121326.valve_mass()  # Mass of the valve

    # Pipe h2 13-14
    pipe_h2_13_14 = Pipe(1, 0.05, h2_13)  # 1 m length, 50 mm diameter
    h2_14 = pipe_h2_13_14.analyze_heat_pipe()

    # Pipe h2 26-18
    pipe_h2_26_18 = Pipe(1, 0.05, h2_26)  # 1 m length, 50 mm diameter
    h2_18 = pipe_h2_26_18.analyze_heat_pipe()

    # Compressor 14-15
    comp_14_15 = Compressor(comp_efficiency=0.9, pressure_ratio=comp_14_15_PI, fluid=h2_14)
    comp14_15_power = comp_14_15.power()
    comp_14_15_mass = comp_14_15.mass(comp14_15_power)
    h2_15 = Fluid(name="H2_15", T=h2_14.T, P=h2_14.P * comp_14_15_PI, C=0, mf=h2_14.mf_given, fluid_type='ParaHydrogen')

    # Pipe h2 15-16
    pipe_h2_15_16 = Pipe(1, 0.05, h2_15)  # 1 m length, 50 mm diameter
    h2_16 = pipe_h2_15_16.analyze_heat_pipe()

    # Intersection 16-17-18
    h2_17 = Fluid(name="H2_17", T=h2_16.T, P=p_cc + 1e5, C=0, mf=h2_mf_cc, fluid_type='ParaHydrogen')  # H2 to CC

    # Pipe h2 17-19
    pipe_h2_17_19 = Pipe(1, 0.05, h2_17)  # 1 m length, 50 mm diameter
    h2_19 = pipe_h2_17_19.analyze_heat_pipe()

    # Valve 19-20
    valve_1920 = Valve(fluid=h2_19, valve_efficiency=0.9) 
    valve_1920_mass = valve_1920.valve_mass()  # Mass of the valve
    h2_20 = Fluid(name="H2_20", T=h2_19.T, P=p_cc, C=0, mf=h2_mf_cc, fluid_type='ParaHydrogen')  # H2 to CC

    # Pipe h2 20-21
    h2_20_copy = Fluid(name="H2_20_copy", T=h2_20.T, P=h2_20.P, C=0, mf=h2_mf_cc/2, fluid_type='ParaHydrogen')
    pipe_h2_20_21 = Pipe(1, 0.05, h2_20_copy)  # 1 m length, 50 mm diameter
    h2_21 = pipe_h2_20_21.analyze_heat_pipe()

    # PATH to FC -------------------------------------------------
    pipe_h2_7_8 = Pipe(1, 0.05, h2_7)  # 1 m length, 50 mm diameter
    h2_8 = pipe_h2_7_8.analyze_heat_pipe()
    h2_8_copy = Fluid(name="H2_8_copy", T=h2_8.T, P=h2_8.P, C=0, mf=h2_mf_fc/2, fluid_type='ParaHydrogen')
    
    valve_8_9 = Valve(fluid=h2_8, valve_efficiency=0.9)  # Valve for H2 to fuel cell
    valve_8_9_mass = valve_8_9.valve_mass()  # Mass of the valve
    h2_9 = Fluid(name="H2_9", T=h2_8.T, P=p_fc, C=0, mf=h2_mf_fc, fluid_type='ParaHydrogen')  # H2 to fuel cell

    # Pipe h2 9-10
    h2_9_copy = Fluid(name="H2_9_copy", T=h2_9.T, P=h2_9.P, C=0, mf=h2_mf_fc/2, fluid_type='ParaHydrogen')
    pipe_h2_9_10 = Pipe(1, 0.05, h2_9_copy)  # 1 m length, 50 mm diameter
    h2_10 = pipe_h2_9_10.analyze_heat_pipe()

    # FC Circulation 
    fc_press_drop_h2 = 100000
    h2_22 = Fluid(name="H2_22", T=T_fc, P=p_fc - fc_press_drop_h2 , C=0, mf=h2_mf_rec, fluid_type='ParaHydrogen')  

    # Pipe h2 22-23
    pipe_h2_22_23 = Pipe(1, 0.05, h2_22)  # 1 m length, 50 mm diameter
    h2_23 = pipe_h2_22_23.analyze_heat_pipe()

    # Compressor 23-24
    comp_23_24_PI = p_fc / h2_23.P
    comp_23_24 = Compressor(comp_efficiency=0.9, pressure_ratio=comp_23_24_PI, fluid=h2_23)
    comp23_24_power = comp_23_24.power()
    comp_23_24_mass = comp_23_24.mass(comp23_24_power)
    h2_24 = Fluid(name="H2_24", T=h2_23.T, P=h2_23.P * comp_23_24_PI, C=0, mf=h2_mf_fc/2, fluid_type='ParaHydrogen')

    # Pipe h2 24-25
    pipe_h2_24_25 = Pipe(1, 0.05, h2_24)  # 1 m length, 50 mm diameter
    h2_25 = pipe_h2_24_25.analyze_heat_pipe()

    # FC Cooling ----------------------------------------

    # Pipe cool 2-3
    pipe_cool_2_3 = Pipe(1, 0.05, cool_2)  # 1 m length, 50 mm diameter
    cool_3 = pipe_cool_2_3.analyze_heat_pipe()

    # Intersection 3-4A-4B
    cool_4 = Fluid(name="Cool_4", T=cool_3.T, P=cool_3.P, C=0, mf=cool_3.mf_given/2, fluid_type='Water') 

    # Pipe cool 4-5
    pipe_cool_4_5 = Pipe(1, 0.05, cool_4)  # 1 m length, 50 mm diameter
    cool_5 = pipe_cool_4_5.analyze_heat_pipe()

    # Pipe 19-22
    pipe_cool_19_22 = Pipe(1, 0.05, cool_19, type ='inv')  # 1 m length, 50 mm diameter
    cool_22 = pipe_cool_19_22.analyze_heat_pipe()

    # Intersection 22-21A-21B
    cool_21 = Fluid(name="Cool_21", T=cool_22.T, P=cool_22.P, C=0, mf=cool_22.mf_given/2, fluid_type='Water')

    # Pipe 21-20
    pipe_cool_21_20 = Pipe(1, 0.05, cool_21, type='inv')  # 1 m length, 50 mm diameter
    cool_20 = pipe_cool_21_20.analyze_heat_pipe()

    # Pipe 12-13
    cool_12 = Fluid(name="Cool_12", T=T_fc, P=p_cool, C=0, mf=cool_mf_per_fc, fluid_type='Water')  # Coolant to FC
    pipe_cool_12_13 = Pipe(1, 0.05, cool_12)  # 1 m length, 50 mm diameter
    cool_13 = pipe_cool_12_13.analyze_heat_pipe()

    # Valve 13-20-14
    valve_132014 = Valve(fluid=cool_13, valve_efficiency=0.9)  # Valve for coolant to FC
    valve_132014_mass = valve_132014.valve_mass()  # Mass of the valve
    cool_14 = Fluid(name="Cool_14", T=cool_13.T, P=cool_13.P, C=None, mf=cool_13.mf_given - cool_20.mf_given, fluid_type='Water')  # Coolant to FC
    
    # Pipe 14-15
    pipe_cool_14_15 = Pipe(1, 0.05, cool_14)  # 1 m length, 50 mm diameter
    cool_15 = pipe_cool_14_15.analyze_heat_pipe()

    # Valve 15-16-17
    valve_151617 = Valve(fluid=cool_15, valve_efficiency=0.9)  # Valve for coolant to FC    
    valve_151617_mass = valve_151617.valve_mass()  # Mass of the valve
    
    # Calculate bypass flow rate:
    T_not_HEX = (cool_12.mf_given * (cool_12.T - deltaT_fc) - cool_20.mf_given * cool_5.T) / (cool_12.mf_given - cool_20.mf_given)  # Temperature of coolant to FC without HEX
    cool_17_mf_minimum = Q_dot_rem / (cool_15.cp * (cool_15.T - T_amb))  
    cool_mf_bypass = (T_not_HEX * cool_14.mf_given - cool_17_mf_minimum * T_amb) / cool_15.T
    cool_17_mf = np.maximum(cool_17_mf_minimum, cool_15.mf_given - cool_mf_bypass)  # Minimum mass flow rate of coolant to FC

    cool_17 = Fluid(name="Cool_17", T=cool_15.T, P=cool_15.P, C=0, mf=cool_17_mf, fluid_type='Water')  # Coolant to FC
    cool_16 = Fluid(name="Cool_16", T=cool_15.T, P=cool_15.P, C=0, mf=cool_15.mf_given - cool_17_mf, fluid_type='Water')  # Coolant to FC
    
    # Pipe 17-18
    pipe_cool_17_18 = Pipe(1, 0.05, cool_17)  # 1 m length, 50 mm diameter    
    cool_18 = pipe_cool_17_18.analyze_heat_pipe()

    # Skin Heat Exchanger
    # ---------------------------------- COMPLETE THE SKIN HX SIZING HERE ----------------------------------
    skin_hx = SkinHeatExchanger(area = area_wing, fluid = cool_18, U = 80)
    skin_hx.set_ambient(T_ambient_K = ambient_conditions['T'])
    Q_abs, t_cool_out = skin_hx.absorb_heat(Q_dot_rem)
    cool_23 = Fluid(name="Cool_23", T=t_cool_out, P=cool_18.P, C=0, mf=cool_18.mf_given, fluid_type='Water')  
    Q_dot_rem -= Q_abs  

    #TODO: Fix if lofic and add pipes, also add radiator
    if Q_dot_rem > 0:
        #TODO: Implement radiator
        pass
    else:
        # No need for radiator, all waste heat is absorbed
        pass
    
    cool_24 = cool_23
    cool_24.T -= Q_dot_rem / (cool_24.mf_given * cool_24.cp)

    # Intersection 24-5-25
    cool_25 = Fluid(name="Cool_25", T=(cool_24.T * cool_24.mf_given + cool_5.T * cool_5.mf_given) / (cool_5.mf_given + cool_24.mf_given), P=cool_24.P, C=0, mf=cool_24.mf_given + cool_5.mf_given, fluid_type='Water') 
    
    # Radiator 
    Rad_1 = RamAirHeatExchanger(coolant_temp_K=ra_coolant_temp, fluid=cool_23)
    radiator_area = Rad_1.size_exchanger(Q_dot_rem/2, T_amb )        # T_amb + 0.5 * (cool_23.T - T_amb)
    #print(f"Heat exchange area of the rad: {radiator_area}")

    # Electric heater:
    if cool_25.T <= T_fc - deltaT_fc:
        Q_dot_heater = (cool_25.mf_given * cool_25.cp * (T_fc - deltaT_fc - cool_25.T))  # Heat required to heat the coolant to FC temperature
    else:
        Q_dot_heater = 0
    cool_26 = Fluid(name="Cool_26", T= cool_25.T + Q_dot_heater / (cool_25.cp * cool_25.mf_given), P=cool_25.P, C=0, mf=cool_25.mf_given, fluid_type='Water') 

    # Coolant pump
    pump_26_27 = Pump(fluid=cool_26, pump_eff=0.9, delta_pressure=p_cool - cool_26.P)  # Pump for coolant to FC
    pump_26_27_power = pump_26_27.power()
    cool_27 = Fluid(name="Cool_27", T=cool_26.T, P=p_cool, C=0, mf=cool_26.mf_given, fluid_type='Water')  # Coolant to FC after pump

    # Store all Fluid instances in a dictionary with their name as key
    fluids_dict = {}
    for var_name, var_value in locals().items():
        if isinstance(var_value, Fluid):
            fluids_dict[var_value.name] = {
                'massflow': var_value.mf_given,
                'temperature': var_value.T,
                'pressure': var_value.P
            }
    return fluids_dict

# RUN -----------------------------------------------------------
if __name__ == "__main__":
    
    # Dummy values for main function inputs
    Q_dot_fc = 500000         # Waste heat from fuel cell [W]
    Q_dot_eps = 20000        # Waste heat from EPS [W]
    p_fc = 1.8e5               # Fuel cell pressure [Pa]
    p_cc = 12e5             # Combustion chamber pressure [Pa]
    h2_mf_fc = 0.01          # H2 mass flow to FC [kg/s]
    h2_mf_cc = 0.01          # H2 mass flow to CC [kg/s]
    T_fc = 273.15+160            # Fuel cell temperature [K]
    T_cc = T_fc           # Combustion chamber temperature [K]
    air_mf_fc = 1.0          # Air mass flow to FC [kg/s]
    T_amb = 273.15+15           # Ambient temperature [K]
    rho_amb = 1.225          # Ambient air density [kg/m^3]
    V_amb = 120.0            # Ambient velocity [m/s]
    p_amb = 101325           # Ambient pressure [Pa]
    h2_mf_rec = 0.005         # H2 recirculation mass flow [kg/s]
    air_out_fc = 1         # Air out of FC [kg/s]
    p_sto = 7e5              # Storage pressure [Pa]

    fluids = main(Q_dot_fc, Q_dot_eps, p_fc, p_cc, h2_mf_fc, h2_mf_cc, T_fc, T_cc, air_mf_fc, T_amb, rho_amb, V_amb, p_amb, h2_mf_rec, air_out_fc, p_sto)

    #tms_plotting.fluids_table(fluids)