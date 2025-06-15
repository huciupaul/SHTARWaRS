# file: analysis.py
# desc: This script analyzes the 4D design space tensor of the SHTARWaRS design tool.

# Global imports
import numpy as np
import pickle
import copy
from typing import Tuple

# Local imports
from global_constants import Beechcraft_1900D, NOx_pPax_TO, NOx_pPax_cruise
from performance.MTOW.mtow import power_and_wing_loading as pws
from performance.Integration.integration import main as cg_excursion
from cost import calc_cost as f_cost
from sus import get_sus as f_gwp

# Cached constants
M_CARGO_AFT = Beechcraft_1900D["M_cargo_aft"]  # [kg] Cargo mass in the aft compartment

def load_tensor(file_path: str) -> np.ndarray:
    """
    Load the 4D design space tensor from a pickle file.
    
    Args:
        file_path (str): Path to the pickle file containing the design space tensor.
        
    Returns:
        np.ndarray: The loaded 4D design space tensor.
    """
    with open(file_path, 'rb') as f:
        tensor = pickle.load(f)
    return tensor

def objective_function(weight: float=0.5, design: np.ndarray=None, N_PAX: np.ndarray=None) -> Tuple[np.ndarray, np.ndarray]:
    """Objective function to collapse the 4D (N, M, P, Q) design
    space tensor to both a 3D scalar field and a 4D (N, M, P, 2) tensor.
    The 3D scalar field is the weighted sum of cost and GWP, while the 4D
    tensor contains the cost and GWP values for each design.

    Args:
        weight (float, optional): Cost weight from 0 to 1. Defaults to 0.5.
        design (np.ndarray, optional): Design (N, M, P, Q) tensor. Defaults to None.

    Returns:
        Tuple[np.ndarray, np.ndarray]: A tuple containing the 3D scalar field 
        and the 4D tensor.
    """
    # Initialize the selection tensor and scalar field
    selection_tensor = np.full(design.shape[:-1] + (2,), np.nan)
    scalar_field = np.full(design.shape[:-1], np.nan)
    
    # Arrays for GWP
    GWP_fc = design[..., 20]  # GWP of the fuel cell production/disposal
    GWP_sto = design[..., 21]  # GWP of the storage system production/disposal
    GWP_eps = design[..., 22]  # GWP of the EPS production/disposal
    m_h2_nom = design[..., 18]   # Nominal mass of hydrogen used
    m_nox = design[..., 14]  # Mass of NOx produced
    
    # Arrays for cost
    fc_cost = design[..., 17]  # Cost of the fuel cell
    P_eps = design[..., 19]  # Max power used by the electric propulsion system
    m_sto = design[..., 3]  # Mass of the storage system
    m_h2_nom = design[..., 18]  # Nominal mass of hydrogen used
    m_h2 = design[..., 2]  # Total mass of hydrogen stored
    
    # Calculate GWP and cost for each design (per passenger)
    GWP = f_gwp(GWP_fc, GWP_sto, GWP_eps, m_h2_nom, m_nox)/N_PAX
    cost = f_cost(fc_cost, P_eps, m_sto, m_h2_nom, m_h2)/N_PAX
    # Store GWP and cost in the selection tensor
    selection_tensor[..., 0] = cost
    selection_tensor[..., 1] = GWP
    
    # Normalize the GWP and cost values
    GWP_norm = GWP / np.nanmax(GWP)
    cost_norm = cost / np.nanmax(cost)
    
    # Calculate the scalar field as a weighted sum of cost and GWP
    scalar_field = weight * cost_norm + (1 - weight) * GWP_norm
    
    return scalar_field, selection_tensor

def main(weight: float=0.5,
         dim_bounds: np.ndarray=np.array([[0.1, 1.0],
                                          [0.1, 1.0],
                                          [0.0, 1.0]]),
         fpath_design: str='data/logs/result_tensor.pkl',
         fpath_loading: str='data/logs/loading_tensor.pkl') -> None:
    """Main function to collapse the 4D design space to a 3D scalar field using the provided weights.

    Args:
        weights (tuple, optional): Integration, MTOW, N_PAX, and m_cargo weights. Defaults to (0.5, 0.25, 0.125, 0.125).
        dim_bounds (np.ndarray, optional): Bounds for the dimensions of the design space tensor. Defaults to [[0.1, 1.0], [0.1, 1.0], [0.0, 1.0]].
        fpath (str, optional): Filepath to the design space tensor created from main.py. Defaults to 'data/design_space.pkl'.
    """
    
    # Load the 4D design space tensor
    try:
        tensor = load_tensor(fpath_design)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {fpath_design}")
    
    # Load the 4D loading space tensor
    try:
        loading = load_tensor(fpath_loading)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {fpath_loading}")
    
    print(f"Design space tensor shape: {tensor.shape}")
    print(f"Loading space tensor shape: {loading.shape}")
    
    # Create axes for the tensor dimensions
    N, M, P, Q = tensor.shape
    splits = np.linspace(dim_bounds[0, 0], dim_bounds[0, 1], N)
    toga_throttle = np.linspace(dim_bounds[1, 0], dim_bounds[1, 1], M)
    cruise_throttle = np.linspace(dim_bounds[2, 0], dim_bounds[2, 1], P)
        
    ### TENSOR SHAPE DESCRIPTION ###
    # The tensor is expected to have the shape (N, M, P, Q) where:
    # N: Power splits (FC use vs Turboprop use)
    # M: FC TOGA throttle setting
    # P: FC cruise throttle setting
    # Q: Design variables: 
    #   < m_EPS, m_FC, m_H2, m_storage, m_TMS, V_FC, V_storage, V_ELMO, MTOW, L_storage, D_storage >
    # Keep a 'design' tensor that maintains the original tensor structure
    design      = copy.deepcopy(tensor)
    # Initialize a 'tradeoff' tensor with the same N, M, P dimensions and a 4th dimension containing 4 np.nan values

    # Apply design space constraints
    valid_integration, N_PAX, m_cargo = cg_excursion(design)
    valid_mass = pws(loading)
    valid_NOx   = (design[..., 15]/N_PAX < NOx_pPax_TO) & \
        (design[..., 16]/N_PAX < NOx_pPax_cruise)  # NOx emissions constraints
    
    valid3d = valid_integration & valid_mass  & valid_NOx
    # Keep only valid designs, fill rest with np.nan
    # lift into a 4th dim so it matches (N, M, P, Q)
    valid4d = valid3d[..., None]        # shape (N, M, P, 1)
    valid4d = np.broadcast_to(valid4d, design.shape)

    # now prune
    design[~valid4d] = np.nan
    
    print(f"Percentage pruned by integration: {np.sum(~valid_integration) / np.prod(design.shape[:-1]) * 100:.2f}%")
    print(f"Percentage pruned by mass loading: {np.sum(~valid_mass) / np.prod(design.shape[:-1]) * 100:.2f}%")
    print(f"Percentage pruned by NOx emissions: {np.sum(~valid_NOx) / np.prod(design.shape[:-1]) * 100:.2f}%")
    print(f"Total percentage pruned from design space: {np.sum(~valid4d) / np.prod(design.shape) * 100:.2f}%")
    
    # Calculate the scalar field and selection tensor
    scalar_field, selection_tensor = objective_function(weight, design, N_PAX)
    
    # Save the scalar field and selection tensor to a pickle file
    with open('data/logs/CostnSus.pkl', 'wb') as f:
        pickle.dump(selection_tensor, f)
    
    with open('data/logs/SWAG_score.pkl', 'wb') as f:
        pickle.dump(scalar_field, f)
    
    print(f"Scalar field shape: {scalar_field.shape}")
    
    
    
if __name__ == "__main__":
    main(
        weight=0.5,
        dim_bounds=np.array([[0.1, 1.0],
                             [0.1, 1.0],
                             [0.1, 1.0]]),
        fpath_design='data/logs/result_tensor.pkl',
        fpath_loading='data/logs/loading_tensor.pkl'
    )