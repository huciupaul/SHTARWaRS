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
# Original 

_PW_req_curve_TO = np.maximum(_PW_cl, _PW_cr)
_PW_req_curve_cr = _cons["PW_req_curve_cr"]


def power_and_wing_loading(
    design: np.ndarray,
    loading: np.ndarray
    ) -> np.ndarray:
    """Check loading against original aircraft constraints analysis and return a mask.

    Args:
        design (np.ndarray): Design tensor with shape (N, M, P, Q) where:
            N: Power splits (FC use vs Turboprop use)
            M: FC TOGA throttle setting
            P: FC cruise throttle setting
            Q: Design variables: 
                < m_EPS, m_FC, m_H2, m_storage, m_TMS, V_FC, V_storage, V_ELMO, MTOW, L_storage, D_storage >
                
        loading (np.ndarray): Mass loading tensor with shape (N, M, P, Q, 4) where:
            N: Power splits (FC use vs Turboprop use)
            M: FC TOGA throttle setting
            P: FC cruise throttle setting
            Q: Flight-phase dependent loading vector:
                < PW_1, WS_1, PW_2, WS_2 > x 4
                where:
                    PW_1: Power loading at phase start
                    WS_1: Wing loading at phase start
                    PW_2: Power loading at phase end
                    WS_2: Wing loading at phase end
    
    Returns:
        np.ndarray: Boolean mask indicating whether the tensor indices meet the constraints.
    """
    PW_1 = loading[:, :, :, :, 0]  # Power loading at phase start
    WS_1 = loading[:, :, :, :, 1]
    PW_2 = loading[:, :, :, :, 2]  # Power loading at phase end
    WS_2 = loading[:, :, :, :, 3]
    
    # Check the first two rows for Take-Off (TO) and Second-Stage Climb
    PW_req_at_WS_TO = np.interp(
        WS_1[:2],
        _WS_x,
        _PW_req_curve_TO,
        left=np.inf,
        right=np.inf
    )
    
    # Check cruise and hold
    PW_req_at_WS_cr = np.interp(
        WS_2[2:],
        _WS_x,
        _PW_req_curve_cr,
        left=np.inf,
        right=np.inf
    )
    
    valid_design = (
        (np.all(PW_1[:2] >= PW_req_at_WS_TO)) &  # TO and Second-Stage Climb
        (np.all(PW_2[2:] >= PW_req_at_WS_cr)) &  # Cruise and Hold
        (np.all(WS_1 >= _WS_stall)) &  # Stall speed check
        (np.all(WS_2 >= _WS_stall))    # Stall speed check
    )
    
    # Design vectors that do not meet the constraints are filled with NaN
    de = np.full_like(design, np.nan, dtype=np.float64)
    de[valid_design] = design[valid_design]
    
    
    return de
