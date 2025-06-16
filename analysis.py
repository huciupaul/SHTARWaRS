# file: analysis.py
# desc: This script analyzes the 4D design space tensor of the SHTARWaRS design tool.

# Global imports
import numpy as np
import pickle
import copy
from typing import Tuple
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull


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

def pareto_front(selection_tensor: np.ndarray, fpath: str="data/plots/pareto_front.png", dpi=600, 
                 remove_outliers: bool=True, outlier_threshold: float=2.0,
                 bounds: np.ndarray=np.array([
                     [0.1, 1.0],
                     [0.1, 1.0],
                     [0.1, 1.0]
                     ])) -> None:
    """Create a Pareto front plot from the selection tensor.
    Args:
        selection_tensor (np.ndarray): The 4D selection tensor containing cost and GWP values.
        fpath (str, optional): Filepath to save the Pareto front plot. Defaults to "data/plots/pareto_front.png".
        dpi (int, optional): DPI for the saved plot. Defaults to 600.
        remove_outliers (bool, optional): Whether to remove outliers. Defaults to True.
        outlier_threshold (float, optional): IQR multiplier for outlier detection. Defaults to 2.0.
        bounds (np.ndarray, optional): Bounds for the dimensions of the design space tensor. Defaults to [[0.1, 1.0], [0.1, 1.0], [0.0, 1.0]].
    """
    # Create parameter arrays based on tensor dimensions and bounds
    N, M, P, _ = selection_tensor.shape
    splits = np.linspace(bounds[0, 0], bounds[0, 1], N)
    toga_throttle = np.linspace(bounds[1, 0], bounds[1, 1], M)
    cruise_throttle = np.linspace(bounds[2, 0], bounds[2, 1], P)
    
    # Track original indices before flattening
    indices = np.zeros(selection_tensor[..., 0].shape + (3,), dtype=int)
    for i in range(N):
        for j in range(M):
            for k in range(P):
                indices[i, j, k] = [i, j, k]
    
    # Extract cost and GWP from the selection tensor
    cost = selection_tensor[..., 0].flatten()
    gwp = selection_tensor[..., 1].flatten()
    flat_indices = indices.reshape(-1, 3)
    
    # Remove NaN values and keep corresponding indices
    mask = ~np.isnan(cost) & ~np.isnan(gwp)
    cost = cost[mask]/1e6  # Convert cost to M€
    gwp = gwp[mask]
    flat_indices = flat_indices[mask]
    
    # Create data points array
    points = np.column_stack((cost, gwp))
    
    # Remove outliers if requested
    if remove_outliers:
        # Calculate IQR for cost and GWP
        q1_cost, q3_cost = np.percentile(cost, [25, 75])
        q1_gwp, q3_gwp = np.percentile(gwp, [25, 75])
        iqr_cost = q3_cost - q1_cost
        iqr_gwp = q3_gwp - q1_gwp
        
        # Define outlier bounds
        cost_lower = q1_cost - outlier_threshold * iqr_cost
        cost_upper = q3_cost + outlier_threshold * iqr_cost
        gwp_lower = q1_gwp - outlier_threshold * iqr_gwp
        gwp_upper = q3_gwp + outlier_threshold * iqr_gwp
        
        # Filter outliers
        outlier_mask = ((cost >= cost_lower) & (cost <= cost_upper) & 
                (gwp >= gwp_lower) & (gwp <= gwp_upper))
        
        # Save original points for later reference
        all_points = points.copy()
        original_count = len(points)
        
        # Apply outlier filter
        points = points[outlier_mask]
        flat_indices = flat_indices[outlier_mask]
        outlier_count = original_count - len(points)
        print(f"Removed {outlier_count} outliers ({outlier_count/original_count:.1%} of data)")
    
    # Find Pareto-optimal points and their indices
    pareto_points = []
    pareto_indices = []
    for i, point in enumerate(points):
        is_dominated = False
        for other_point in points:
            # Check if other_point dominates point (lower cost AND lower GWP)
            if (other_point[0] < point[0] and other_point[1] <= point[1]) or \
               (other_point[0] <= point[0] and other_point[1] < point[1]):
                is_dominated = True
                break
        if not is_dominated:
            pareto_points.append(point)
            pareto_indices.append(flat_indices[i])
    
    pareto_points = np.array(pareto_points)
    pareto_indices = np.array(pareto_indices)
    print(f"Found {len(pareto_points)} Pareto-optimal points")
    
    # Sort Pareto points by cost for drawing the front
    sorted_indices = np.argsort(pareto_points[:, 0])
    sorted_pareto = pareto_points[sorted_indices]
    
    # Create a scatter plot for all designs
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    
    # Plot all design points
    plt.scatter(points[:, 0], points[:, 1], c='blue', alpha=0.3, label='Designs')
    
    # Plot Pareto points with different style
    plt.scatter(pareto_points[:, 0], pareto_points[:, 1], 
                c='red', edgecolors='black', s=70, alpha=0.9, label='Pareto-optimal designs')
    
    # Draw the Pareto front line
    plt.plot(sorted_pareto[:, 0], sorted_pareto[:, 1], 'r-', linewidth=2, label='Pareto front')
    
    # Generate convex hull Pareto front (boundary front)
    if len(pareto_points) >= 3:  # ConvexHull requires at least 3 points
        hull = ConvexHull(pareto_points)
        # Get the vertices of the convex hull in order
        hull_vertices = pareto_points[hull.vertices]
        hull_indices = pareto_indices[hull.vertices]
        # Sort by x-coordinate to draw from left to right
        sort_order = np.argsort(hull_vertices[:, 0])
        hull_vertices = hull_vertices[sort_order]
        hull_indices = hull_indices[sort_order]
        
        # Keep only the lower part of the hull (the actual Pareto boundary)
        hull_y_max = hull_vertices[:, 1].max()
        lower_mask = hull_vertices[:, 1] < hull_y_max
        lower_hull = hull_vertices[lower_mask]
        lower_hull_indices = hull_indices[lower_mask]
        
        if len(lower_hull) >= 2:  # Need at least 2 points to draw a line
            # Plot convex hull boundary
            plt.plot(lower_hull[:, 0], lower_hull[:, 1], 'g--', linewidth=2, 
                    label='Convex Pareto boundary')
            
            # Label key points on the convex hull boundary
            for i, (point, idx) in enumerate(zip(lower_hull, lower_hull_indices)):
                # Convert indices to parameter values
                split_val = splits[idx[0]]
                toga_val = toga_throttle[idx[1]]
                cruise_val = cruise_throttle[idx[2]]
                
                plt.annotate(f"({split_val:.2f}, {toga_val:.2f}, {cruise_val:.2f})",
                           xy=(point[0], point[1]),
                           xytext=(10, 5 + (i % 3) * 15),  # Stagger labels to avoid overlap
                           textcoords='offset points',
                           fontsize=8,
                           arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2', color='black'))
    
    # Add additional information
    min_cost_idx = np.argmin(pareto_points[:, 0])
    min_gwp_idx = np.argmin(pareto_points[:, 1])
    
    # Convert min cost and min GWP indices to parameter values
    min_cost_params = (
        splits[pareto_indices[min_cost_idx][0]],
        toga_throttle[pareto_indices[min_cost_idx][1]],
        cruise_throttle[pareto_indices[min_cost_idx][2]]
    )
    
    min_gwp_params = (
        splits[pareto_indices[min_gwp_idx][0]],
        toga_throttle[pareto_indices[min_gwp_idx][1]],
        cruise_throttle[pareto_indices[min_gwp_idx][2]]
    )
    
    # plt.annotate(f'Min cost: {pareto_points[min_cost_idx, 0]:.3f} M€/PAX\n({min_cost_params[0]:.2f}, {min_cost_params[1]:.2f}, {min_cost_params[2]:.2f})',
    #             xy=(pareto_points[min_cost_idx, 0], pareto_points[min_cost_idx, 1]),
    #             xytext=(10, 15), textcoords='offset points',
    #             arrowprops=dict(arrowstyle='->'))
    
    # plt.annotate(f'Min GWP: {pareto_points[min_gwp_idx, 1]:.1f} kg CO₂e/PAX\n({min_gwp_params[0]:.2f}, {min_gwp_params[1]:.2f}, {min_gwp_params[2]:.2f})',
    #             xy=(pareto_points[min_gwp_idx, 0], pareto_points[min_gwp_idx, 1]),
    #             xytext=(10, -25), textcoords='offset points',
    #             arrowprops=dict(arrowstyle='->'))
    
    # # Add statistics about the design space
    # plt.text(0.02, 0.98, 
    #         f'Total designs: {len(points)}\nPareto-optimal designs: {len(pareto_points)}',
    #         transform=plt.gca().transAxes, 
    #         verticalalignment='top',
    #         bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    # Add labels and title
    plt.xlabel('Lifetime Cost [M€/PAX]')
    plt.ylabel(r'GWP [kg CO$_2$e/PAX]')
    plt.title('Pareto Front of Cost vs GWP')
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(fpath, dpi=dpi)
    plt.close()

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
    
    try:
        convergence = load_tensor('data/logs/convergence_tensor.pkl')
    except FileNotFoundError:
        raise FileNotFoundError("Convergence tensor not found. Please run the main.py script to generate it.")
    
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
    valid4d = np.broadcast_to(valid4d, design.shape[:-1] + (2,))  # shape (N, M, P, 2)
    # valid4d = np.broadcast_to(valid4d, design.shape)

    # # now prune
    # design[~valid4d] = np.nan
    
    print(f"Percentage pruned by integration: {np.sum(~valid_integration) / np.prod(design.shape[:-1]) * 100:.2f}%")
    print(f"Percentage pruned by mass loading: {np.sum(~valid_mass) / np.prod(design.shape[:-1]) * 100:.2f}%")
    print(f"Percentage pruned by NOx emissions: {np.sum(~valid_NOx) / np.prod(design.shape[:-1]) * 100:.2f}%")
    
    # Calculate the scalar field and selection tensor
    scalar_field, selection_tensor = objective_function(weight, design, N_PAX)
    
    # Create a Pareto front plot
    pareto_front(selection_tensor, fpath="data/plots/pareto_front_no_cons.png", dpi=600)
    constrained_selection_tensor = selection_tensor.copy()
    constrained_selection_tensor[~valid4d] = np.nan  # Apply the valid mask to the selection tensor
    pareto_front(constrained_selection_tensor, fpath="data/plots/pareto_front.png", dpi=600)
    
    # Save the scalar field and selection tensor to a pickle file
    with open('data/logs/CostnSus.pkl', 'wb') as f:
        pickle.dump(selection_tensor, f)
    
    with open('data/logs/SWAG_score.pkl', 'wb') as f:
        pickle.dump(scalar_field, f)
    
    # Determine optimal point index
    optimal_index = np.nanargmin(scalar_field)
    optimal_params = np.unravel_index(optimal_index, scalar_field.shape)
    optimal_split = splits[optimal_params[0]]
    optimal_toga = toga_throttle[optimal_params[1]]
    optimal_cruise = cruise_throttle[optimal_params[2]]
    
    print(f"Optimal point found at index {optimal_index} with parameters:"
          f"\n  Power split: {optimal_split:.2f}"
          f"\n  TOGA throttle: {optimal_toga:.2f}"
          f"\n  Cruise throttle: {optimal_cruise:.2f}")
    
    optimal_convergence = convergence[optimal_params]
    print(optimal_convergence)
    plt.plot(optimal_convergence)
    plt.show()
    
    print(f"Scalar field shape: {scalar_field.shape}")
    
    
    
if __name__ == "__main__":
    main(
        weight=0.5,
        dim_bounds=np.array([[0.1, 0.9],
                             [0.1, 1.0],
                             [0.1, 1.0]]),
        fpath_design='data/logs/result_tensor.pkl',
        fpath_loading='data/logs/loading_tensor.pkl'
    )