import numpy as np
from typing import Tuple

# Local imports
from .constants import LAPSE, G_0, R_AIR, k_air

def isa_atmosphere(h: np.ndarray, T_sl: float = 288.15, P_sl: float = 101_325.0) -> Tuple[float, float, float]:
    """Thin-layer ISA (no tropopause).
    Returns T [K], P [Pa], rho [kg/m^3], and speed of sound."""
    T = T_sl - LAPSE * h
    P = P_sl * (T / T_sl) ** (G_0 / (LAPSE * R_AIR))
    rho = P / (R_AIR * T)
    a = np.sqrt(k_air * R_AIR * T)  # speed of sound
    return T, P, rho, a