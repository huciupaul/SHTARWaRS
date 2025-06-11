# file: performance/MTOW/mtow.py
# desc: Check MTOW against original aircraft constraints analysis and return a score.

# Global imports
import numpy as np
from typing import Tuple, Union

# Local imports
from global_constants import MTOW_orig, C_d_0_orig
from performance.MTOW import ca
 
_cons      = ca.constraint_curves(C_d_0_orig, MTOW_orig) 
_WS_x      = _cons["WS_x"]
_PW_cl     = _cons["PW_cl"]
_PW_cr     = _cons["PW_cr"]
_WS_stall  = _cons["WS_stall"]
# Original 

_PW_req_curve_TO = np.maximum(_PW_cl, _PW_cr)
_PW_req_curve_cr = _PW_cr


def power_and_wing_loading(
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
    # Extract the 4 load‐cases from the last axis (shape = N×M×P×Q×4)
    pw_ws = loading[..., :4]               # (N, M, P, Q, 4)
    pw_ws = np.moveaxis(pw_ws, -1, 0)      # (4, N, M, P, Q)
    PW_1, WS_1, PW_2, WS_2 = pw_ws         # unpack channels

    # Now PW_1, WS_1, … all have shape (N, M, P, Q)
    # and you can proceed exactly as before:
    PW_req_at_WS_TO = np.interp(
        WS_1[:2].reshape(-1),
        _WS_x, _PW_req_curve_TO,
        left=np.inf, right=np.inf
    ).reshape(2, *WS_1.shape[1:])

    PW_req_at_WS_cr = np.interp(
        WS_2[2:].reshape(-1),
        _WS_x, _PW_req_curve_cr,
        left=np.inf, right=np.inf
    ).reshape(WS_2.shape[0]-2, *WS_2.shape[1:])

    valid_design = (
        (np.all(PW_1[:2] >= PW_req_at_WS_TO)) &
        (np.all(PW_2[2:] >= PW_req_at_WS_cr)) &
        (np.all(WS_1 >= _WS_stall)) &
        (np.all(WS_2 >= _WS_stall))
    )

    return valid_design
