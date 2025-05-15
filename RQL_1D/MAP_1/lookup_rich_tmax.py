# lookup_rich_Tmax.py
import numpy as np
from scipy.interpolate import RegularGridInterpolator          # SciPy ≥1.2
import h5py                                                    # or numpy.load

# -------- 2.1  load arrays -----------------------------------
with h5py.File("runs/rich_Tmax_map.h5", "r") as f:
    phi_vec = f["phi"][...]       # shape (nφ,)
    Tin_vec = f["Tin"][...]       # shape (nT,)
    T_peak  = f["Tpeak"][...]     # shape (nT, nφ)

# construct interpolator

interp_obj = RegularGridInterpolator(
    (Tin_vec, phi_vec), T_peak,
    method="quadratic", bounds_error=False, fill_value=np.nan
)                                                   # :contentReference[oaicite:1]{index=1}

def Tmax(phi_in, Tin_in):
    """
    Return interpolated rich-zone peak temperature [K]
    for arbitrary equivalence ratio φ_in and inlet T Tin_in [K].

    Returns NaN if the query lies outside the pre-computed domain.
    """
    Tin_in = float(Tin_in)
    phi_in = float(phi_in)
    return float(interp_obj((Tin_in, phi_in)))