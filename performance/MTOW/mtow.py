# file: performance/MTOW/mtow.py
# desc: Check MTOW against original aircraft constraints analysis and return a score.

# Global imports
import numpy as np
from typing import Tuple, Union

# Local imports
from SHTARWaRS.global_constants import OG_constraints, MTOW_orig, C_d_0_orig
from SHTARWaRS.performance.MTOW import ca
 
_cons      = ca.constraint_curves(C_d_0_orig, MTOW_orig) 
_WS_x      = _cons["WS_x"]
_PW_cl     = _cons["PW_cl"]
_PW_cr     = _cons["PW_cr"]
_WS_stall  = _cons["WS_stall"]

_PW_req_curve = np.maximum(_PW_cl, _PW_cr)


def mtow_score(
    loading_vector: np.ndarray
    ) -> np.ndarray:
    """Check MTOW against original aircraft constraints analysis and return a score.

    Args:
        MTOW (Union[float, np.ndarray]): Maximum Take-Off Weight.
        PW (Union[float, np.ndarray]): Power loading.
        SW (Union[float, np.ndarray]): Wing loading.

    Returns:
        np.ndarray: Score for the design based on MTOW, power loading, and wing loading.
    """
    
    PW_req_at_WS = np.interp(
        WS,
        _WS_x,
        _PW_req_curve,
        left=np.inf,
        right=np.inf
    )
    
    valid_design = (PW >= PW_req_at_SW) & (SW <= _WS_stall)
    
    # For valid designs, calculate the distance from the original design point
    d = np.where(
        valid_design,
        np.sqrt(
            ((PW - PW_OG_TO) / PW_OG_TO) ** 2 +
            ((SW - SW_OG_TO) / SW_OG_TO) ** 2
        )
    )
    
    score = np.where(
        valid_design,
        1/ (1 + d),
        0.0
    )
    
    return score
