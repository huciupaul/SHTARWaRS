import heat_exchanger_sizing
import fuel_cell_info
import flight_condition
import hydrogen_storage
import heat_removal_required

class MainSimulation:
    def __init__(self):
        # Define fuel cells
        self.LTPEM = fuel_cell_info.FuelCell(
            name="LTPEM",
            stack_efficiency=0.6,
            T=273.15 + 80,
            P_A=1.8 * 101325 + 0.06 * 1e5,
            P_C=1.8 * 101325,
            RH_A=0.5,
            RH_C=0.38,
            stoic_ratio_C=1.8,
            stoic_ratio_A=1.05
        )

        # Define flight condition
        self.takeoff = flight_condition.FlightCondition(
            name="Takeoff",
            T_amb=273.15 + 55,
            P_amb=101325,
            RH_amb=0,
            V=0,
            power_required=1_908_000,
            power_split=1,
            thermal_efficiency=0.4,
            propulsive_efficiency=0.8,
            P_cc=12.1 * 101325,
            T_cc=573.15
        )

        # Define hydrogen storage
        self.LH2 = hydrogen_storage.HydrogenStorage(name="LH2", T=20, P=6 * 101325)

        # Set design point
        self.design_point = heat_removal_required.Design_point(
            fuel_cell=self.LTPEM,
            flight_condition=self.takeoff,
            hydrogen_storage=self.LH2,
            P_C=1.85 * 101325
        )

        # Initialize HeatExchanger
        self.heat_exchanger = heat_exchanger_sizing.HeatExchanger(design_point=self.design_point)

    def run(self):
        print(self.design_point.fuel_cell.name)
        print(self.design_point.flight_condition.name)
        
        self.design_point.O2_cooling_required()
        self.design_point.mass_flow_calculation()
        self.design_point.heat_removal_available()

        # Coolant constants: Ethylene Glycol solution 60% (lowest freezing point)
        # https://www.engineeringtoolbox.com/ethylene-glycol-d_146.html
        cp_coolant =  3424                  # J/kg*deg 
        rho_coolant = 1068                  # at 60deg
        boiling_coolant = 111.1 + 273.15    # K

        self.heat_exchanger.mass_flow(cp_coolant)
        print(f"Mass flow rate of coolant: {self.heat_exchanger.m_dot_cool:.4f} kg/s")

        # Constants
        T_in = self.heat_exchanger.T_in
        T_out = self.heat_exchanger.T_out
        T_air_in = self.design_point.flight_condition.T_amb
        T_air_out = 60 + 273.15 # K Placeholder
        delta_T_air_coolant = T_in - T_out 
        F = 0.97                 # Using R and P chart: F is correction factor for single-pass, cross HE with unmixed fluids

        # Material constants
        k = 205           #  Aluminum (typically alloy 1100) -->  corrosion resistance in deionized water coolant loops + lightweight

        # Equations
        #U_rad = self.heat_exchanger.u_radiator(A_norm, delta_T_air_coolant)
        #print(f"U_rad: {U_rad:.2f} W/m^2/K")

        area_U = self.heat_exchanger.area_U(delta_T_air_coolant)
        D_int = self.heat_exchanger.D_int()
        print(f"Wall thickness and inner diameter: {D_int} m")

        T_coolant_avg = self.heat_exchanger.T_coolant_avg(T_in, T_out)

        delta_T_air_coolant_arithmetic = self.heat_exchanger.delta_T_air_coolant_arithmetic(T_coolant_avg, T_air_in)
        print(f"Delta T of air and coolant (arithmetic): {delta_T_air_coolant_arithmetic:.2f} K")

        delta_T_air_coolant_LMTD = self.heat_exchanger.delta_T_air_coolant_LMTD(T_in, T_out, T_air_in, T_air_out, F)
        print(f"Delta T of air and coolant (LMTD): {delta_T_air_coolant_LMTD:.2f} K")

        r_ratio = self.heat_exchanger.r_ratio(T_in, T_out, T_air_in, T_air_out)
        print(f"R ratio: {r_ratio:.2f} ")

        R = self.heat_exchanger.R(T_in, T_out, T_air_in, T_air_out)
        print(f"R: {R:.2f} ")

        P = self.heat_exchanger.P(T_in, T_air_in, T_air_out)
        print(f"P: {P:.2f} ")

if __name__ == "__main__":
    sim = MainSimulation()
    sim.run()