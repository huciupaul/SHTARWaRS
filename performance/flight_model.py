import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

#  CONSTANTS & ISA MODEL 
R_AIR = 287.05      # [J/(kg·K)]
G_0 = 9.80665       # [m/s²]
LAPSE = 0.0065      # [K/m]
k_air = 1.4         # [-] Specific heat ratio of air
A_inlet = 0.015446456799602548*np.pi*2.78**2/4

def isa_atmosphere(h: np.ndarray, T_sl: float = 288.15, P_sl: float = 101_325.0) -> Tuple[float, float, float]:
    """Thin-layer ISA (no tropopause).
    Returns T [K], P [Pa], rho [kg/m^3], and speed of sound."""
    T = T_sl - LAPSE * h
    P = P_sl * (T / T_sl) ** (G_0 / (LAPSE * R_AIR))
    rho = P / (R_AIR * T)
    a = np.sqrt(k_air * R_AIR * T)  # speed of sound
    return T, P, rho, a

#  BASIC AIRCRAFT DESCRIPTION 
@dataclass
class Aircraft:
    name: str
    wing_area: float
    wing_span: float
    CD0: float
    prop_diameter: float
    MTOW: float
    TOGA: float = 1.908e+6
    n_prop: int = 2

    @property
    def aspect_ratio(self) -> float:
        return self.wing_span ** 2 / self.wing_area

    @property
    def oswald(self) -> float:
        return 1.78 * (1 - 0.045 * self.aspect_ratio ** 0.68) - 0.64

    @property
    def prop_area(self) -> float:
        return self.n_prop * np.pi * (self.prop_diameter * 0.5) ** 2

#  WAYPOINT SPEC 
@dataclass
class Waypoint:
    name: str
    altitude: float   # [m]
    airspeed: float   # [m/s]
    roc: float        # [m/s]  (+ climb, – descent, 0 level)
    hold_time: float = 0.0  # [s]

#  POWER SPEC
@dataclass
class Powerpoint:
    name: str
    power: float=None # [W]
    time: float=None  # [s]
    until_phase: str=None   # e.g. "cruise"
    
#  TURBOPROP MODEL
class Turboprop:
    """Full model of a turboprop engine. The model is based on the following
    assumptions:
    - Brayton cycle
    - Constant inlet, compressor, combustion chamber, and turbine efficiencies
    - Constant specific heat ratio of air and combustion products
    - Constant specific heat of combustion products
    - Constant pressure ratio across the compressor and combustion chamber
    - Constant turbine inlet temperature
    - ISA atmosphere
    - No bleed air extraction
    - No variable geometry
    """
    
    def __init__(
        self,
        h0: np.ndarray, # [m] Altitude
        V0: np.ndarray, # [m/s] Airspeed
        A_inlet: float, # [m^2] Inlet area
        eta_in: float,  # [-] Inlet efficiency
        PI_comp: float, # [-] Compressor pressure ratio
        eta_comp: float, # [-] Compressor efficiency
        PI_cc: float,  # [-] Combustion chamber pressure ratio
        eta_cc: float, # [-] Combustion chamber efficiency
        LHV_fuel: float, # [J/kg] Lower heating value of fuel
        T04: float,  # [K] Turbine inlet temperature
        eta_turb: float, # [-] Turbine efficiency
        eta_mech: float, # [-] Mechanical efficiency
        c_pa: float, # [J/(kg·K)] Specific heat of air
        k_air: float,   # [-] Specific heat ratio of air
        c_pg: float, # [J/(kg·K)] Specific heat of combustion products
        k_gas: float  # [-] Specific heat ratio of combustion products
    ) -> None:
        
        self.h0 = h0
        self.V0 = V0
        self.T0, self.P0, self.rho0, self.a0 = isa_atmosphere(h0)
        self.M0 = V0 / self.a0
        self.mdot_air = self.rho0 * V0 * A_inlet
        self.eta_in = eta_in
        self.PI_comp = PI_comp
        self.eta_comp = eta_comp
        self.PI_cc = PI_cc
        self.eta_cc = eta_cc
        self.LHV_fuel = LHV_fuel
        self.T04 = T04
        self.eta_turb = eta_turb
        self.eta_mech = eta_mech
        self.c_pa = c_pa
        self.k_air = k_air
        self.c_pg = c_pg
        self.k_gas = k_gas
        
        self.mdot_fuel = None
        self.eta_th = None
    
    #  PUBLIC INTERFACE
    @property
    def profile(self) -> Dict[str, float]:
        """Returns the engine profile."""
        if self.mdot_fuel is None or self.eta_th is None:
            self.__build()
        return {
            "mdot_fuel": self.mdot_fuel,
            "eta_th": self.eta_th
        }
    
    #  INTERNAL HELPERS
    @staticmethod
    def __02(
        T0: np.ndarray,
        P0: np.ndarray,
        M0: np.ndarray,
        eta_in: float,
        k_air: float
        ) -> np.ndarray:
        """Station 2 (compressor inlet) conditions."""
        T02 = T0*(1 + (k_air-1)*M0**2/2)
        P02 = P0*(1 + eta_in*(k_air-1)*M0**2/2)**(k_air/(k_air-1))
        return T02, P02
    
    @staticmethod
    def __03(
        T02: np.ndarray,
        P02: np.ndarray,
        PI_comp: float,
        eta_comp: float,
        k_air: float
        ) -> np.ndarray:
        """Station 3 (compressor exit) conditions."""
        T03 = T02*(1 + (1/eta_comp)*(PI_comp**((k_air-1)/k_air)-1))
        P03 = P02*PI_comp
        return T03, P03
    
    @staticmethod
    def __04(
        T04: np.ndarray,
        P03: np.ndarray,
        PI_cc: float
    ):
        """Station 4 (combustion chamber exit) conditions."""
        T04 = T04
        P04 = P03*PI_cc
        return T04, P04
    
    @staticmethod
    def __mdot_fuel(
        mdot_air: np.ndarray,
        c_pg: float,
        T03: np.ndarray,
        T04: np.ndarray,
        LHV_fuel: float,
        eta_cc: float
    ):
        """Mass flow rate of fuel."""
        mdot_fuel = mdot_air*(c_pg*(T04-T03)/(LHV_fuel*eta_cc))
        return mdot_fuel
    
    @staticmethod
    def __05(
        mdot_air: np.ndarray,
        mdot_fuel: np.ndarray,
        c_pa: float,
        c_pg: float,
        T02: np.ndarray,
        T03: np.ndarray,
        T04: np.ndarray,
        P04: np.ndarray,
        eta_mech: float,
        eta_turb: float,
        k_gas: float
    ):
        """Station 5 (turbine exit) conditions."""
        mdot = mdot_air + mdot_fuel
        Wdot_comp = mdot_air*c_pa*(T03-T02)
        Wdot_turb = Wdot_comp/eta_mech
        T05 = T04 - Wdot_turb/(mdot*c_pg)
        P05 = P04*(1-(1/eta_turb)*(1-T05/T04))**(k_gas/(k_gas-1))
        return T05, P05
    
    @staticmethod
    def __eta_th(
        mdot_air: np.ndarray,
        mdot_fuel: np.ndarray,
        T03: np.ndarray,
        T04: np.ndarray,
        T05: np.ndarray,
        P0: np.ndarray,
        P05: np.ndarray,
        V0: np.ndarray,
        c_pg: float,
        k_gas: float
    ):
        """Thermal efficiency."""
        mdot = mdot_air + mdot_fuel
        Wdot_gg = mdot*c_pg*T05*(1-(P0/P05)**((k_gas-1)/k_gas))-0.5*mdot*V0**2
        eta_th = Wdot_gg/(mdot*c_pg*(T04-T03))
        return eta_th
    
    #  CORE BUILDER
    def __build(self) -> None:
        """Build the turboprop model."""
        # Station 0 (inlet)
        T0 = self.T0
        P0 = self.P0
        V0 = self.V0
        M0 = self.M0
        mdot_air = self.mdot_air

        # Station 2 (compressor inlet)
        T02, P02 = self.__02(T0, P0, M0, self.eta_in, self.k_air)

        # Station 3 (compressor exit)
        T03, P03 = self.__03(T02, P02, self.PI_comp, self.eta_comp, self.k_air)

        # Station 4 (combustion chamber exit)
        T04, P04 = self.__04(self.T04, P03, self.PI_cc)

        # Mass flow rate of fuel
        mdot_fuel = self.__mdot_fuel(mdot_air, self.c_pg, T03, T04, self.LHV_fuel, self.eta_cc)

        # Station 5 (turbine exit)
        T05, P05 = self.__05(mdot_air, mdot_fuel, self.c_pa, self.c_pg, T02, T03, T04, P04,
                             self.eta_mech, self.eta_turb, self.k_gas)

        # Thermal efficiency
        eta_th = self.__eta_th(mdot_air, mdot_fuel, T03, T04, T05,
                               P0, P05, V0,
                               self.c_pg,
                               self.k_gas)
        
        # Store results
        self.mdot_fuel = mdot_fuel
        self.eta_th = eta_th
        
    
#  FLIGHT MISSION (KINEMATICS ONLY) 
class FlightMission:
    """Builds a time-resolved kinematic mission profile *and* auto-sizes the
    cruise leg to satisfy a required **ground range**.

    - Specify *total_range_m* (straight-line + holds) if you want the class to
      work out how long to sit in cruise.
    - Waypoint labelled *cruise_wp_name* anchors altitude & speed for that leg.
      Any *hold_time* at that waypoint is inserted **before** the auto-sized
      cruise block.
    """

    def __init__(
        self,
        waypoints: List[Waypoint],
        dt: float = 1.0,
        total_range_m: Optional[float] = None,
        cruise_wp_name: str = "cruise",
    ) -> None:
        if len(waypoints) < 2:
            raise ValueError("Need at least two waypoints.")
        self.wps = waypoints
        self.dt = dt
        self.total_range = total_range_m
        self.cruise_wp_name = cruise_wp_name.lower()
        self._profile: Dict[str, np.ndarray] = {}

    #  PUBLIC
    @property
    def profile(self) -> Dict[str, np.ndarray]:
        if not self._profile:
            self.__build()
        return self._profile

    #  INTERNAL HELPERS 
    @staticmethod
    def __ground_speed(V: float, roc: float) -> float:
        """Horizontal speed neglecting wind."""
        if V == 0:
            return 0.0
        gamma = np.arcsin(np.clip(roc / V, -1.0, 1.0))
        return V * np.cos(gamma)

    #  CORE BUILDER 
    def __build(self) -> None:
        # Identify cruise waypoint 
        try:
            i_cruise = next(i for i, wp in enumerate(self.wps) if wp.name.lower() == self.cruise_wp_name)
        except StopIteration:
            if self.total_range is not None:
                raise ValueError("Cruise waypoint not found; needed for range sizing.")
            i_cruise = None  # maybe no cruise sizing requested

        # Containers 
        V_full, h_full, roc_full, time_full = [], [], [], []
        t = 0.0
        dist_accum = 0.0  # [m]

        # Utility to append one step 
        def __append(V_s: float, h_s: float, roc_s: float):
            nonlocal t, dist_accum
            V_full.append(V_s)
            h_full.append(h_s)
            roc_full.append(roc_s)
            time_full.append(t)
            dist_accum += self.__ground_speed(V_s, roc_s) * self.dt
            t += self.dt

        # ------------------------- BUILD CLIMB + CRUISE HOLD -------------------
        for idx, wp in enumerate(self.wps):
            # Insert hold at *current* waypoint (except if idx==0 where we start)
            if wp.hold_time > 0 and not (idx == 0):
                n_hold = int(np.ceil(wp.hold_time / self.dt))
                for _ in range(n_hold):
                    __append(wp.airspeed, wp.altitude, 0.0)

            # Last waypoint → nothing to transition to
            if idx == len(self.wps) - 1:
                break

            nxt = self.wps[idx + 1]

            # If we are at cruise wp AND range sizing requested → defer transition
            if idx == i_cruise and self.total_range is not None:
                break  # store state, add descent later

            # Normal transition (climb/descent/level)
            dh = nxt.altitude - wp.altitude
            roc_cmd = wp.roc if wp.roc != 0 else nxt.roc
            if abs(dh) < 1e-6:
                continue  # same altitude & roc zero ⇒ nothing to do
            if roc_cmd == 0:
                raise ValueError(f"No ROC defined for leg {wp.name} → {nxt.name}.")
            n_steps = int(np.ceil(abs(dh / roc_cmd) / self.dt))
            alt_seg = np.linspace(wp.altitude, nxt.altitude, n_steps, endpoint=False)
            V_seg = np.linspace(wp.airspeed, nxt.airspeed, n_steps, endpoint=False)
            for V_s, h_s in zip(V_seg, alt_seg):
                __append(V_s, h_s, roc_cmd)

        # ------------------------- AUTO-SIZE CRUISE -----------------------------
        if self.total_range is not None:
            # Build the **descent + remaining legs** separately to know distance
            V_post, h_post, roc_post = [], [], []
            dist_post = 0.0
            wp_start = self.wps[i_cruise]
            for idx in range(i_cruise, len(self.wps) - 1):
                wp = self.wps[idx]
                nxt = self.wps[idx + 1]

                # hold at wp (except the cruise wp – already handled above)
                if idx != i_cruise and wp.hold_time > 0:
                    n_hold = int(np.ceil(wp.hold_time / self.dt))
                    for _ in range(n_hold):
                        V_post.append(wp.airspeed)
                        h_post.append(wp.altitude)
                        roc_post.append(0.0)
                        dist_post += self.__ground_speed(wp.airspeed, 0.0) * self.dt

                # transition to nxt
                dh = nxt.altitude - wp.altitude
                roc_cmd = wp.roc if wp.roc != 0 else nxt.roc
                if abs(dh) < 1e-6:
                    continue
                if roc_cmd == 0:
                    raise ValueError(f"No ROC for leg {wp.name} → {nxt.name}.")
                n_steps = int(np.ceil(abs(dh / roc_cmd) / self.dt))
                alt_seg = np.linspace(wp.altitude, nxt.altitude, n_steps, endpoint=False)
                V_seg = np.linspace(wp.airspeed, nxt.airspeed, n_steps, endpoint=False)
                for V_s, h_s in zip(V_seg, alt_seg):
                    V_post.append(V_s)
                    h_post.append(h_s)
                    roc_post.append(roc_cmd)
                    dist_post += self.__ground_speed(V_s, roc_cmd) * self.dt

            # Distance left for cruise ------------------------------------
            dist_remaining = self.total_range - dist_accum - dist_post
            if dist_remaining < 0:
                raise ValueError(
                    f"Mission (climb/descent/holds) exceeds total range by {-dist_remaining/1000:.1f} km."
                )
            if dist_remaining > 0:
                V_cruise = wp_start.airspeed
                n_cruise = int(np.ceil(dist_remaining / (V_cruise * self.dt)))
                for _ in range(n_cruise):
                    __append(V_cruise, wp_start.altitude, 0.0)
            # Append the pre-built post-segments ---------------------------
            for V_s, h_s, roc_s in zip(V_post, h_post, roc_post):
                __append(V_s, h_s, roc_s)

        # If no total_range given & we broke at cruise wp above (unlikely) – resume build normally
        if self.total_range is None and i_cruise is not None and t == 0:
            raise RuntimeError("Builder logic error – mission truncated.")

        # Export 
        self._profile = {
            "time": np.asarray(time_full),
            "V": np.asarray(V_full),
            "alt": np.asarray(h_full),
            "ROC": np.asarray(roc_full),
        }

    #  QUICK PLOT 
    def quicklook(self) -> None:
        import matplotlib.pyplot as plt
        p = self.profile
        fig, axs = plt.subplots(1, 3, figsize=(15, 4))
        axs[0].plot(p["time"], p["alt"])
        axs[0].set_ylabel("Altitude [m]")
        axs[1].plot(p["time"], p["V"])
        axs[1].set_ylabel("TAS [m/s]")
        axs[2].plot(p["time"], p["ROC"])
        axs[2].set_ylabel("ROC [m/s]")
        for ax in axs:
            ax.set_xlabel("Time [s]")
            ax.grid(alpha=0.3)
        fig.tight_layout()
        plt.show()

if __name__ == "__main__":
    wps = [
        Waypoint("takeoff", 0.0, 54.0, 4.0),
        Waypoint("climb-1", 2_438.4, 140.0, 4.0),
        Waypoint("climb-2", 4_876.8, 146.0, 4.0),
        Waypoint("cruise", 7_620.0, 142.5, 0.0),   # anchor for auto-sized cruise
        Waypoint("IAF", 900.0, 70.0, -6.0, hold_time=240.),
        Waypoint("final", 0.0, 60.0, -3.0),
    ]
    pws = [
        Powerpoint("takeoff", 1.0, time=5*60),
        Powerpoint("climb", 0.95, until_phase="cruise"),
        Powerpoint("cruise"), # anchor for auto-sized cruise power
        Powerpoint("descent", 0.19, until_phase="final")
    ]
    ac_model = Aircraft(
        name="Beechcraft 1900D",
        wing_area=28.79,
        wing_span=17.64,
        CD0=0.0215,
        prop_diameter=2.78,
        MTOW=10_000.0
    )
    # Ask for a 707 km mission
    mission = FlightMission(wps, dt=1.0, total_range_m=707e3)
    mission.quicklook()
    # print("A_inlet", A_inlet)

