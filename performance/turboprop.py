import numpy as np
from typing import Tuple

# Local imports
from common.constants import A_inlet, R_AIR

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
        name: str,
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
        
        self.name = name
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
    
    #  COMPUTE
    def compute(self, T0: np.ndarray, P0: np.ndarray, rho0:np.ndarray, V0: np.ndarray, eta_prop: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compute the turboprop engine performance."""
        
        # Atmospheric conditions
        M0 = V0/np.sqrt(self.k_air*R_AIR*T0)
        mdot_air = rho0 * V0 * A_inlet/eta_prop
        
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
        
        return mdot_fuel, eta_th