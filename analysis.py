# file: analysis.py
# desc: This script analyzes the 4D design space tensor of the SHTARWaRS design tool.

# Global imports
import numpy as np
import pickle
import copy

# Local imports
# from SHTARWaRS.performance import *

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

def integration(design: np.ndarray) -> np.ndarray:
    """Evaluate design integration into existing aircraft, score the results and normalize.

    Args:
        design (np.ndarray): design tensor with shape (N, M, P, Q) where:
            N: Power splits (FC use vs Turboprop use)
            M: FC TOGA throttle setting
            P: FC cruise throttle setting
            Q: Design variables: 
                < m_EPS, m_FC, m_H2, m_storage, m_TMS, V_FC, V_storage, V_ELMO, MTOW, A_radiator, L_storage, m_cargo >

    Returns:
        np.ndarray: Vector of integration scores
    """
    # Create np.array of relevant design variables
    var_idx = np.array([0, 1, 2, 3, 4])
    m_EPS, m_FC, m_H2, m_storage, m_TMS = design[:, :, :, var_idx]
     

def main(weights: tuple=(0.5, 0.25, 0.125, 0.125),
         dim_bounds: np.ndarray=np.array([[0.1, 1.0],
                                          [0.1, 1.0],
                                          [0.0, 1.0]]),
         fpath: str='data/design_space.pkl'):
    """Main function to collapse the 4D design space to a 3D scalar field using the provided weights.

    Args:
        weights (tuple, optional): Integration, MTOW, N_PAX, and m_cargo weights. Defaults to (0.5, 0.25, 0.125, 0.125).
        dim_bounds (np.ndarray, optional): Bounds for the dimensions of the design space tensor. Defaults to [[0.1, 1.0], [0.1, 1.0], [0.0, 1.0]].
        fpath (str, optional): Filepath to the design space tensor created from main.py. Defaults to 'data/design_space.pkl'.
    """
    # Quick validation of weights
    if sum(weights) != 1.0:
        raise ValueError("Weights must sum to 1.0")
    if len(weights) != 4:
        raise ValueError("Weights must be a tuple of length 4")
    
    # Load the 4D design space tensor
    try:
        tensor = load_tensor(fpath)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {fpath}")
    
        
    ### TENSOR SHAPE DESCRIPTION ###
    # The tensor is expected to have the shape (N, M, P, Q) where:
    # N: Power splits (FC use vs Turboprop use)
    # M: FC TOGA throttle setting
    # P: FC cruise throttle setting
    # Q: Design variables: 
    #   < m_EPS, m_FC, m_H2, m_storage, m_TMS, V_FC, V_storage, V_ELMO, MTOW, A_radiator, L_storage, m_cargo >
    
    # Keep a 'design' tensor that maintains the original tensor structure
    design      = copy.deepcopy(tensor)
    # Initialize a 'tradeoff' tensor with the same N, M, P dimensions and a 4th dimension containing 4 np.nan values
    tradeoff    = np.full(tensor.shape[:-1] + (4,), np.nan)
    tradeoff[:, :, :] = tensor[:, :, :]
    