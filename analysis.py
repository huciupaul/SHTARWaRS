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
from scipy.stats import gaussian_kde
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
import plotly.graph_objects as go
import pandas as pd



# Local imports
from global_constants import Beechcraft_1900D, NOx_pPax_TO, NOx_pPax_cruise
from performance.MTOW.mtow import power_and_wing_loading as pws
from performance.Integration.integration import main as cg_excursion
from cost import calc_cost as f_cost
from sus import get_sus as f_gwp
from cost import calc_cost_breakdown

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
    # print(fc_cost[..., 0])
    P_eps = design[..., 19]  # Max power used by the electric propulsion system
    m_sto = design[..., 3]  # Mass of the storage system
    m_h2_nom = design[..., 18]  # Nominal mass of hydrogen used
    m_h2 = design[..., 2]  # Total mass of hydrogen stored
    
    # Calculate GWP and cost for each design (per passenger)
    GWP = f_gwp(GWP_fc, GWP_sto, GWP_eps, m_h2_nom, m_nox)/N_PAX
    cost = f_cost(fc_cost, P_eps, m_sto, m_h2, m_h2_nom)/N_PAX
    cost = cost / 1e6  # Convert cost to M€    

    # Normalize cost and GWP via L^2 norm
    # obj = np.stack([cost, GWP], axis=-1)
    cost_flat = cost
    gwp_flat = GWP
    
    c = cost_flat
    g = gwp_flat
    c_min, c_max = c.min(), c.max()
    g_min, g_max = g.min(), g.max()

    cost_norm = (c - c_min) / (c_max - c_min)
    gwp_norm  = (g - g_min) / (g_max - g_min)
    # Calculate the scalar field as a weighted sum of cost and GWP
    scalar_field = weight * cost_norm + (1 - weight) * gwp_norm
    
    selection_tensor[..., 0] = cost 
    selection_tensor[..., 1] = GWP
    
    return scalar_field, selection_tensor


def score_distribution(
        selection_tensor: np.ndarray,
        fpath: str = "data/plots/cost_gwp_distrib.png",
        dpi: int = 600,
        method: str | None = None,        # 'auto', 'hex', 'kde', 'scatter'
        bins: int = 60,                   # for hex / hist
    ) -> None:
    """
    Visualise the distribution of cost vs GWP with a smarter default than
    a blunt 2-D histogram.

    Parameters
    ----------
    method : str | None
        'hex', 'kde', 'scatter', or None => choose automatically.
    """

    cost = selection_tensor[..., 0].ravel()
    gwp  = selection_tensor[..., 1].ravel()
    mask = np.isfinite(cost) & np.isfinite(gwp)
    cost, gwp = cost[mask], gwp[mask]
    N = len(cost)

    if method is None:
        # choose automatically
        if N > 3e5:
            method = "hex"
        elif N > 5e4:
            method = "scatter"
        else:
            method = "kde"

    sns.set_theme(style="whitegrid")
    g = sns.JointGrid(x=cost, y=gwp, height=6, space=0)

    if method == "hex":
        hb = g.ax_joint.hexbin(cost, gwp, gridsize=bins, cmap="viridis",
                               norm=plt.LogNorm(), mincnt=1)
        cb = g.figure.colorbar(hb, ax=g.ax_joint, label="counts (log-scaled)")
    elif method == "kde":
        # joint KDE
        xy = np.vstack([cost, gwp])
        kde = gaussian_kde(xy)(xy)
        idx = kde.argsort()
        g.ax_joint.scatter(cost[idx], gwp[idx], c=kde[idx], s=8,
                           cmap="viridis", edgecolors="none")
        # filled contour lines (levels = 10%, 30%, 50%, 70%, 90%)
        levels = kde.max() * np.array([.1, .3, .5, .7, .9])
        g.ax_joint.tricontourf(cost, gwp, kde, levels=levels,
                               alpha=.25, cmap="viridis")
        cb = g.figure.colorbar(plt.cm.ScalarMappable(cmap="viridis"),
                            ax=g.ax_joint, label="relative density")
    elif method == "scatter":
        # compute local density for colour
        xy = np.vstack([cost, gwp])
        kde = gaussian_kde(xy)(xy)
        idx = kde.argsort()
        g.ax_joint.scatter(cost[idx], gwp[idx], c=kde[idx], s=4,
                           cmap="viridis", edgecolors="none")
        cb = g.figure.colorbar(plt.cm.ScalarMappable(cmap="viridis"),
                            ax=g.ax_joint, label="relative density")
    else:
        raise ValueError("method must be 'hex', 'kde', 'scatter', or None")

    # marginals
    g.ax_marg_x.hist(cost, bins=bins, color="grey", alpha=.7)
    g.ax_marg_y.hist(gwp,  bins=bins, color="grey", alpha=.7, orientation="horizontal")

    # labels
    g.ax_joint.set_xlabel("Lifetime cost [M€/PAX]")
    g.ax_joint.set_ylabel(r"GWP [kg CO$_2$e/PAX]")
    # g.figure.suptitle(f"Cost vs GWP distribution  (N = {N:,})", y=.98)
    g.figure.tight_layout()
    g.figure.show()
    g.figure.savefig(fpath, dpi=dpi)
    plt.close(g.figure)
    
    
def print_design_point(point_name, tensor_point, n_pax=None):
    """Print all design variables for a specific point with proper formatting."""
    print(f"\n=== {point_name.upper()} POINT DESIGN VARIABLES ===")
    
    # Masses
    print(f"Mass Variables:")
    print(f"  m_eps:          {tensor_point[0]} kg  (Electric Propulsion System)")
    print(f"  m_fc:           {tensor_point[1]} kg  (Fuel Cell)")
    print(f"  m_h2:           {tensor_point[2]} kg  (Hydrogen)")
    print(f"  m_sto:          {tensor_point[3]} kg  (Storage System)")
    print(f"  m_tms_front:    {tensor_point[11]} kg  (Thermal Management - Front)")
    print(f"  m_tms_aft:      {tensor_point[12]} kg  (Thermal Management - Aft)")
    print(f"  m_tms_mid:      {tensor_point[13]} kg  (Thermal Management - Mid)")
    print(f"  m_h2_nom:       {tensor_point[18]} kg  (Nominal Hydrogen Usage)")
    
    # Volumes
    print(f"\nVolume Variables:")
    print(f"  V_fc:           {tensor_point[4]} m³  (Fuel Cell Volume)")
    print(f"  V_sto:          {tensor_point[5]} m³  (Storage System Volume)")
    print(f"  V_elmo:         {tensor_point[6]} m³  (Electric Motor Volume)")
    
    # Aircraft properties
    print(f"\nAircraft Properties:")
    print(f"  MTOW:           {tensor_point[7]} kg  (Maximum Take-Off Weight)")
    print(f"  length_sto:     {tensor_point[8]} m   (Storage System Length)")
    print(f"  diameter_sto:   {tensor_point[9]} m   (Storage System Diameter)")
    print(f"  M_aft_cargo:    {tensor_point[10]} kg (Aft Cargo Mass)")
    print(f"  N_PAX:          {n_pax} pax   (Number of Passengers)") if n_pax is not None else print(f"  N_PAX:         {tensor_point[23]:.2f} pax  (Number of Passengers)")
    
    # Emissions & Power
    print(f"\nEmissions & Power:")
    print(f"  m_nox:          {tensor_point[14]} kg  (NOx Mass)")
    print(f"  NOx_TO:         {tensor_point[15]} kg/s  (NOx at Takeoff)")
    print(f"  NOx_cruise:     {tensor_point[16]} kg/s  (NOx at Cruise)")
    print(f"  P_elmo:         {tensor_point[19]} W  (Electric Motor Power)")
    
    # Cost & Environmental Impact
    print(f"\nCost & Environmental Impact:")
    print(f"  cost_fc:        {tensor_point[17]/1e6} M€  (Fuel Cell Cost)")
    print(f"  co2_fc:         {tensor_point[20]} kg CO₂e  (Fuel Cell GWP)")
    print(f"  co2_sto:        {tensor_point[21]} kg CO₂e  (Storage System GWP)")
    print(f"  co2_eps:        {tensor_point[22]} kg CO₂e  (Electric Propulsion System GWP)")

def _pareto_mask(selection_tensor: np.ndarray) -> np.ndarray:
    """Return a boolean mask of Pareto-optimal points (dominance on cost & gwp)."""
    # Get the original shape (excluding last dimension)
    original_shape = selection_tensor.shape[:-1]
    
    # Flatten the arrays
    cost = selection_tensor[..., 0].ravel()
    gwp = selection_tensor[..., 1].ravel()
    
    # Get indices of finite values
    valid_mask = np.isfinite(cost) & np.isfinite(gwp)
    valid_indices = np.where(valid_mask)[0]
    valid_cost = cost[valid_indices]
    valid_gwp = gwp[valid_indices]
    
    # Find Pareto front points
    pts = np.column_stack((valid_cost, valid_gwp))
    sort_idx = np.argsort(pts[:, 0])           # ascending cost
    pts_sorted = pts[sort_idx]
    gwp_min_run = np.minimum.accumulate(pts_sorted[:, 1])
    is_pareto_sorted = pts_sorted[:, 1] == gwp_min_run
    
    # Map back to original indices
    pareto_indices = valid_indices[sort_idx[is_pareto_sorted]]
    
    # Create full mask
    full_mask = np.zeros_like(cost, dtype=bool)
    full_mask[pareto_indices] = True
    
    # Reshape to original dimensions
    return full_mask.reshape(original_shape)

def interactive_pareto_front(
        selection_tensor: np.ndarray,
        design_tensor: np.ndarray,
        splits: np.ndarray,
        toga_throttle: np.ndarray,
        cruise_throttle: np.ndarray,
        N_PAX: np.ndarray,
        fpath: str = "data/plots/interactive_pareto.html",
        remove_outliers: bool = False,
        outlier_threshold: float = 1.5,
        star_boy: np.ndarray | None = None,
):
    """
    Create an interactive Plotly plot of the Pareto front where hovering over points
    shows design variables.
    """
    # Extract cost and GWP data
    cost = selection_tensor[..., 0].ravel()
    gwp = selection_tensor[..., 1].ravel()
    mask = np.isfinite(cost) & np.isfinite(gwp)
    cost, gwp = cost[mask], gwp[mask]
    
    # Get indices of valid points for mapping back to design tensor
    valid_indices = np.where(mask)[0]
    original_shape = selection_tensor.shape[:-1]
    point_indices = np.array([np.unravel_index(idx, original_shape) for idx in valid_indices])
    
    if remove_outliers:
        def _cut(x):
            q1, q3 = np.percentile(x, [25, 75])
            iqr = q3 - q1
            lo, hi = q1 - outlier_threshold*iqr, q3 + outlier_threshold*iqr
            return (x >= lo) & (x <= hi)
        keep = _cut(cost) & _cut(gwp)
        cost, gwp = cost[keep], gwp[keep]
        point_indices = point_indices[keep]
    
    # Prepare data for Pareto front
    pts = np.column_stack((cost, gwp))
    sort_idx = np.argsort(pts[:, 0])
    pts_sorted = pts[sort_idx]
    gwp_min_run = np.minimum.accumulate(pts_sorted[:, 1])
    is_pareto = pts_sorted[:, 1] == gwp_min_run
    
    # Get Pareto front points
    pareto_idx = sort_idx[is_pareto]
    P = pts[pareto_idx]
    pareto_point_indices = point_indices[pareto_idx]
    
    # Sort Pareto points by cost
    pareto_sort = np.argsort(P[:, 0])
    P = P[pareto_sort]
    pareto_point_indices = pareto_point_indices[pareto_sort]
    
    # Create data for the plot
    data = []
    
    # All feasible points
    df_all = pd.DataFrame({
        'cost': pts[:, 0],
        'gwp': pts[:, 1],
    })
    
    # Add design variables for all points
    for i, idx in enumerate(point_indices):
        idx_tuple = tuple(idx)
        df_all.loc[i, 'power_split'] = splits[idx_tuple[0]]
        df_all.loc[i, 'toga_throttle'] = toga_throttle[idx_tuple[1]]
        df_all.loc[i, 'cruise_throttle'] = cruise_throttle[idx_tuple[2]]
        df_all.loc[i, 'n_pax'] = N_PAX[idx_tuple]
        
        # Add key design variables
        design_point = design_tensor[idx_tuple]
        df_all.loc[i, 'm_eps'] = design_point[0]  # Electric Propulsion System mass
        df_all.loc[i, 'm_fc'] = design_point[1]   # Fuel Cell mass
        df_all.loc[i, 'm_h2'] = design_point[2]   # Hydrogen mass
        df_all.loc[i, 'm_sto'] = design_point[3]  # Storage System mass
        df_all.loc[i, 'mtow'] = design_point[7]   # MTOW
    
    # All feasible points scatter
    scatter_all = go.Scatter(
        x=df_all['cost'],
        y=df_all['gwp'],
        mode='markers',
        marker=dict(
            size=5,
            color='blue',
            opacity=0.3,
        ),
        name='Feasible Designs',
        hovertemplate=(
            '<b>Feasible Design</b><br>'
            'Cost: %{x:.4f} M€/PAX<br>'
            'GWP: %{y:.4f} kg CO₂e/PAX<br>'
            'Power Split: %{customdata[0]:.2f}<br>'
            'TOGA Throttle: %{customdata[1]:.2f}<br>'
            'Cruise Throttle: %{customdata[2]:.2f}<br>'
            'PAX: %{customdata[3]:.1f}<br>'
            'MTOW: %{customdata[4]:.1f} kg<br>'
            'm_EPS: %{customdata[5]:.1f} kg<br>'
            'm_FC: %{customdata[6]:.1f} kg<br>'
            'm_H2: %{customdata[7]:.1f} kg<br>'
            'm_STO: %{customdata[8]:.1f} kg<br>'
        ),
        customdata=df_all[['power_split', 'toga_throttle', 'cruise_throttle', 'n_pax', 
                           'mtow', 'm_eps', 'm_fc', 'm_h2', 'm_sto']].values
    )
    data.append(scatter_all)
    
    # Create dataframe for Pareto points
    df_pareto = pd.DataFrame({
        'cost': P[:, 0],
        'gwp': P[:, 1],
    })
    
    # Add design variables for Pareto points
    for i, idx in enumerate(pareto_point_indices):
        idx_tuple = tuple(idx)
        df_pareto.loc[i, 'power_split'] = splits[idx_tuple[0]]
        df_pareto.loc[i, 'toga_throttle'] = toga_throttle[idx_tuple[1]]
        df_pareto.loc[i, 'cruise_throttle'] = cruise_throttle[idx_tuple[2]]
        df_pareto.loc[i, 'n_pax'] = N_PAX[idx_tuple]
        
        # Add key design variables
        design_point = design_tensor[idx_tuple]
        df_pareto.loc[i, 'm_eps'] = design_point[0]
        df_pareto.loc[i, 'm_fc'] = design_point[1]
        df_pareto.loc[i, 'm_h2'] = design_point[2]
        df_pareto.loc[i, 'm_sto'] = design_point[3]
        df_pareto.loc[i, 'mtow'] = design_point[7]
    
    # Pareto front line
    pareto_line = go.Scatter(
        x=df_pareto['cost'],
        y=df_pareto['gwp'],
        mode='lines',
        line=dict(color='red', width=2),
        name='Pareto Front',
        hoverinfo='skip'
    )
    data.append(pareto_line)
    
    # Pareto points
    pareto_points = go.Scatter(
        x=df_pareto['cost'],
        y=df_pareto['gwp'],
        mode='markers',
        marker=dict(
            size=8,
            color='red',
            line=dict(color='black', width=1)
        ),
        name='Pareto Optimal',
        hovertemplate=(
            '<b>Pareto Optimal Design</b><br>'
            'Cost: %{x:.4f} M€/PAX<br>'
            'GWP: %{y:.4f} kg CO₂e/PAX<br>'
            'Power Split: %{customdata[0]:.2f}<br>'
            'TOGA Throttle: %{customdata[1]:.2f}<br>'
            'Cruise Throttle: %{customdata[2]:.2f}<br>'
            'PAX: %{customdata[3]:.1f}<br>'
            'MTOW: %{customdata[4]:.1f} kg<br>'
            'm_EPS: %{customdata[5]:.1f} kg<br>'
            'm_FC: %{customdata[6]:.1f} kg<br>'
            'm_H2: %{customdata[7]:.1f} kg<br>'
            'm_STO: %{customdata[8]:.1f} kg<br>'
        ),
        customdata=df_pareto[['power_split', 'toga_throttle', 'cruise_throttle', 'n_pax',
                             'mtow', 'm_eps', 'm_fc', 'm_h2', 'm_sto']].values
    )
    data.append(pareto_points)
    
    # Add star_boy point if provided
    if star_boy is not None:
        star = go.Scatter(
            x=[star_boy[0]],
            y=[star_boy[1]],
            mode='markers',
            marker=dict(
                size=15,
                symbol='star',
                color='gold',
                line=dict(color='black', width=1)
            ),
            name='Sensitivity-Optimal Design',
            hoverinfo='name'
        )
        data.append(star)
    
    # Create the figure
    layout = go.Layout(
        title='Pareto Front: Cost vs GWP (Interactive)',
        xaxis=dict(title='Lifetime cost [M€/PAX]'),
        yaxis=dict(title='GWP [kg CO₂e/PAX]'),
        hovermode='closest',
        legend=dict(x=0.01, y=0.99, bordercolor='Black', borderwidth=1),
        template='plotly_white'
    )
    
    fig = go.Figure(data=data, layout=layout)
    
    # Save to HTML
    fig.write_html(fpath)
    
    return fig

def pareto_front(
        selection_tensor: np.ndarray,
        fpath: str = "data/plots/pareto_front.png",
        dpi: int = 1200,
        remove_outliers: bool = False,
        outlier_threshold: float = 1.5,
        bounds: np.ndarray = np.array([[0.1, 1.0],
                                       [0.1, 1.0],
                                       [0.1, 1.0]]),
        annotate_weights: bool = True,
        star_boy: np.ndarray | None = None,
        chosen_section: tuple = None,  # ((x1, y1), (x2, y2))
        focus_point: tuple = None,
    ) -> None:
    """
    Plot Pareto front (cost vs GWP) with optional highlight region.

    Parameters
    ----------
    selection_tensor : ndarray, shape (..., 2)
        Last axis = (cost, GWP) *per PAX* (any NaN => ignored).
    fpath : str
        Where to save the .png.
    remove_outliers : bool
        If True, remove points lying outside `outlier_threshold`·IQR
        in either objective.
    annotate_weights : bool
        Write the weight that makes each frontier vertex optimal.
    star_boy : ndarray | None
        If provided, highlights this point with a star marker.
    chosen_section : tuple
        If provided, highlights this section of the Pareto front with an 
        organic-looking shaded region. Format: ((x1, y1), (x2, y2))
    """

    cost   = selection_tensor[..., 0].ravel()
    gwp    = selection_tensor[..., 1].ravel()
    mask   = np.isfinite(cost) & np.isfinite(gwp)
    cost, gwp = cost[mask], gwp[mask]

    if remove_outliers:
        def _cut(x):
            q1, q3 = np.percentile(x, [25, 75])
            iqr = q3 - q1
            lo, hi = q1 - outlier_threshold*iqr, q3 + outlier_threshold*iqr
            return (x >= lo) & (x <= hi)
        keep = _cut(cost) & _cut(gwp)
        cost, gwp = cost[keep], gwp[keep]

    pts = np.column_stack((cost, gwp))
    sort_idx    = np.argsort(pts[:, 0])
    pts_sorted  = pts[sort_idx]
    gwp_min_run = np.minimum.accumulate(pts_sorted[:, 1])
    is_pareto   = pts_sorted[:, 1] == gwp_min_run
    pareto_idx  = sort_idx[is_pareto]
    P           = pts[pareto_idx]
    # Now strictly ascending cost, strictly descending GWP
    P = P[np.argsort(P[:, 0])]

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(pts[:, 0], pts[:, 1], s=8, alpha=.3, label="Feasible Design")
    ax.plot(P[:, 0], P[:, 1], "r-", lw=2, label="Pareto Front")
    ax.scatter(P[:, 0], P[:, 1], c="red", edgecolors="k", zorder=3)

    # Add a gold star for the head point
    if star_boy is not None:
        ax.scatter(star_boy[0], star_boy[1], s=350, marker="*", color="gold", edgecolors="k", zorder=5, label="Chosen Design Point")
    
    # Create an organic-looking shaded region for the chosen section
    if chosen_section is not None:
        (x1, y1), (x2, y2) = chosen_section
        
        # Find Pareto points within the x range
        mask = (P[:, 0] >= min(x1, x2)) & (P[:, 0] <= max(x1, x2))
        section_points = P[mask]
        
        if len(section_points) > 1:
            # Create organic shape around the section
            x_pareto = section_points[:, 0]
            y_pareto = section_points[:, 1]
            
            # Extend the x-range slightly for a more organic look
            x_range = max(x_pareto) - min(x_pareto)
            x_ext = [min(x_pareto) - x_range * 0.05, max(x_pareto) + x_range * 0.05]
            
            # Get smoothed versions of the curves with more points for organic look
            from scipy.interpolate import make_interp_spline
            
            if len(x_pareto) >= 4:  # Need at least 4 points for cubic spline
                x_smooth = np.linspace(min(x_pareto), max(x_pareto), 100)
                spl = make_interp_spline(x_pareto, y_pareto, k=3)  # cubic spline
                y_smooth = spl(x_smooth)
            else:
                x_smooth = x_pareto
                y_smooth = y_pareto
            
            # Create upper and lower boundaries for the shaded region
            y_offset = np.max(y_smooth) * 0.15  # 15% of max y for bubble height
            
            # Upper boundary (higher GWP) - create wavy organic top
            y_upper = y_smooth + y_offset * np.sin(np.linspace(0, 3*np.pi, len(x_smooth))) * 0.5 + y_offset
            
            # Lower boundary (follows Pareto curve closely with small buffer)
            y_buffer = y_offset * 0.3  # Smaller buffer near the Pareto front
            y_lower = y_smooth - y_buffer
            
            # Create polygon vertices by going around the region
            x_polygon = np.concatenate([x_smooth, x_smooth[::-1]])
            y_polygon = np.concatenate([y_upper, y_lower[::-1]])
            
            # Plot the shaded region with a nice gradient alpha
            ax.fill(x_polygon, y_polygon, color='#8A2BE2', alpha=0.15, 
                    zorder=2, label="Chosen Design Space")
    
            
            points = np.array([x_polygon, y_polygon]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # Create a custom colormap with varying transparency
            cmap = LinearSegmentedColormap.from_list(
                "custom_purple", [(0.54, 0.17, 0.89, 0.6), (0.54, 0.17, 0.89, 0.1)])
            
            lc = LineCollection(segments, cmap=cmap, alpha=0.5, linewidth=1.5)
            line_colors = np.linspace(0, 1, len(segments))
            lc.set_array(line_colors)
            ax.add_collection(lc)
    
    if focus_point is not None:
        center_x, center_y, radius = focus_point
        
        # Find Pareto points within the radius
        distances = np.sqrt((P[:, 0] - center_x)**2 + (P[:, 1] - center_y)**2)
        mask = distances <= radius
        section_points = P[mask]
        
        if len(section_points) > 1:
            # Sort points by x-coordinate for proper interpolation
            sort_idx = np.argsort(section_points[:, 0])
            x_pareto = section_points[sort_idx, 0]
            y_pareto = section_points[sort_idx, 1]
            
            # Handle duplicate x values by keeping only unique x values
            unique_indices = []
            seen_x = set()
            for i, x in enumerate(x_pareto):
                if x not in seen_x:
                    seen_x.add(x)
                    unique_indices.append(i)
            
            x_pareto = x_pareto[unique_indices]
            y_pareto = y_pareto[unique_indices]
            
            # Get smoothed versions of the curves
            from scipy.interpolate import make_interp_spline, interp1d
            
            if len(x_pareto) >= 4:  # Need at least 4 points for cubic spline
                x_smooth = np.linspace(min(x_pareto), max(x_pareto), 100)
                spl = make_interp_spline(x_pareto, y_pareto, k=3)  # cubic spline
                y_smooth = spl(x_smooth)
            elif len(x_pareto) >= 2:  # Fall back to linear interpolation
                x_smooth = np.linspace(min(x_pareto), max(x_pareto), 100)
                f = interp1d(x_pareto, y_pareto)
                y_smooth = f(x_smooth)
            else:
                x_smooth = x_pareto
                y_smooth = y_pareto
            
            # Create bubble-like shape around the points
            y_offset = radius * 0.5  # Scale the bubble height based on radius
            
            # Create an organic, bubble-like shape with varying offsets
            theta = np.linspace(0, 2*np.pi, len(x_smooth))
            y_upper = y_smooth + y_offset * (0.8 + 0.2 * np.sin(3*theta))
            y_lower = y_smooth - y_offset * 0.3  # Smaller buffer near the Pareto front
            
            # Create polygon vertices by going around the bubble
            x_polygon = np.concatenate([x_smooth, x_smooth[::-1]])
            y_polygon = np.concatenate([y_upper, y_lower[::-1]])
            
            # Plot the shaded region with gradient alpha
            ax.fill(x_polygon, y_polygon, color='#8A2BE2', alpha=0.15, 
                    zorder=2, label="Chosen Design Space")
            
            # Add a subtle outline with gradient transparency
            points = np.array([x_polygon, y_polygon]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # Create custom colormap for the outline
            cmap = LinearSegmentedColormap.from_list(
                "custom_purple", [(0.54, 0.17, 0.89, 0.6), (0.54, 0.17, 0.89, 0.1)])
            
            lc = LineCollection(segments, cmap=cmap, alpha=0.5, linewidth=1.5)
            line_colors = np.linspace(0, 1, len(segments))
            lc.set_array(line_colors)
            ax.add_collection(lc)
            
            # Draw a small semi-transparent circle at the center point
            ax.scatter([center_x], [center_y], color='#8A2BE2', alpha=0.3, 
                      s=80, zorder=4, edgecolors='none')

    ax.set_xlabel("Lifetime Cost [M€/PAX]")
    ax.set_ylabel(r"Global Warming Potential [kg CO$_\text{2,eq}$/PAX]")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(fpath, dpi=dpi)
    plt.show()
    
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
    
    print(tensor.shape)
        
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
    print(f"Percentage pruned by all constraints: {np.sum(~valid3d) / np.prod(design.shape[:-1]) * 100:.2f}%")
    
    # Calculate the scalar field and selection tensor
    scalar_field, selection_tensor = objective_function(weight, design, N_PAX)
    
    # Create a Pareto front plot
    # pareto_front(selection_tensor, fpath="data/plots/pareto_front_no_cons.png", dpi=600)
    
    # # Generate a score distribution plot
    score_distribution(selection_tensor, fpath="data/plots/cost_gwp_distrib.png", dpi=600, method='kde')
    constrained_selection_tensor = selection_tensor.copy()
    constrained_selection_tensor[~valid4d] = np.nan  # Apply the valid mask to the selection tensor
    
    # pareto_front(constrained_selection_tensor, fpath="data/plots/pareto_front.png", dpi=600)
    score_distribution(constrained_selection_tensor, fpath="data/plots/cost_gwp_distrib_constrained.png", dpi=600, method='kde')
    
    # Save the scalar field and selection tensor to a pickle file
    with open('data/logs/CostnSus.pkl', 'wb') as f:
        pickle.dump(selection_tensor, f)
    
    with open('data/logs/SWAG_score.pkl', 'wb') as f:
        pickle.dump(scalar_field, f)
    
    
    # — Simple scalar‑weighted optimum (for comparison) —
    # scalar_opt_flat = np.nanargmin(scalar_field)
    # scalar_opt_idx = np.unravel_index(scalar_opt_flat, selection_tensor.shape)
    
    # scalar_opt_flat = np.nanargmin(scalar_field)
    # scalar_opt_idx = np.unravel_index(scalar_opt_flat, selection_tensor.shape)
    
    pareto_boundary = _pareto_mask(constrained_selection_tensor)
    # Add dimension for broadcasting
    pareto_boundary = pareto_boundary[..., np.newaxis]
    # Now broadcast to match selection tensor shape
    pareto_boundary = np.broadcast_to(pareto_boundary, selection_tensor.shape)
    
    # Get the power split, TOGA throttle, and cruise throttle for 3 points:
    pareto_options = selection_tensor.copy()
    print(pareto_options.shape)
    pareto_options[~pareto_boundary] = np.nan  # Apply the Pareto mask
    
    head = np.argwhere(
        (pareto_options[..., 0] > 3.425) & (pareto_options[..., 0] < 3.43) &\
            (pareto_options[..., 1] > 10.12) & (pareto_options[..., 1] < 10.14)
    )
    # print(head)
        
    knee = np.argwhere(
        (pareto_options[..., 0] > 3.4450) & (pareto_options[..., 0] < 3.4475) &\
            (pareto_options[..., 1] > 23.598) & (pareto_options[..., 1] < 23.602)
    )
    
    # ankle = np.argwhere(
    #     (pareto_options[..., 0] > -np.inf) & (pareto_options[..., 0] < np.inf) &\
    #         (pareto_options[..., 1] > -np.inf) & (pareto_options[..., 1] < np.inf)
    # )
    
    toe = np.argwhere(
        (pareto_options[..., 0] > 3.505) & (pareto_options[..., 0] < 3.506) &\
            (pareto_options[..., 1] > 31.5) & (pareto_options[..., 1] < 31.55)
    )
    # print(toe)
    
    toenail = np.argwhere(
        (pareto_options[..., 0] > 3.4315) & (pareto_options[..., 0] < 3.4320) &\
            (pareto_options[..., 1] > 29.9) & (pareto_options[..., 1] < 30.0)
    )
    
    # print(knee)
    # print(constrained_selection_tensor[*knee[0], ...])
    
    # Plot the convergence of MTOW
    converged_mtow = convergence[0, 0, 0, :].ravel()
    plt.figure(figsize=(10, 6))
    plt.plot(converged_mtow, marker='o', linestyle='-', color='blue', markersize=4)
    # plt.title('Convergence of MTOW')
    plt.xlabel('Iteration')
    plt.ylabel('MTOW [kg]')
    plt.grid()
    plt.tight_layout()
    plt.savefig('data/plots/convergence_mtow.png', dpi=600)
    
    
    # print("\n=== PARETO OPTIMAL POINTS ===")
    
    # # For the head point
    # print_design_point("head", tensor[*head[0], ...], N_PAX[*head[0]])
    # print(f"Configuration: Split = {splits[head[0, 0]]}, TOGA = {toga_throttle[head[0, 1]]}, Cruise = {cruise_throttle[head[0, 2]]}")
    # print(f"Objectives: Cost = {constrained_selection_tensor[*head[0], 0]} M€/PAX, GWP = {constrained_selection_tensor[*head[0], 1]} kg CO₂e/PAX")
    
    # # For the knee point
    # print_design_point("knee", tensor[*knee[0], ...], N_PAX[*knee[0]])
    # print(f"Configuration: Split = {splits[knee[0, 0]]}, TOGA = {toga_throttle[knee[0, 1]]}, Cruise = {cruise_throttle[knee[0, 2]]}")
    # # print(f"Objectives: Cost = {constrained_selection_tensor[*knee[0], 0]:.2f} M€/PAX, GWP = {constrained_selection_tensor[*knee[0], 1]:.2f} kg CO₂e/PAX")

    # For the ankle point
    # print_design_point("ankle", tensor[*ankle[0], ...], N_PAX[*ankle[0]])
    # print(f"Configuration: Split = {splits[ankle[0, 0]]}, TOGA = {toga_throttle[ankle[0, 1]]}, Cruise = {cruise_throttle[ankle[0, 2]]}")
    # print(f"Objectives: Cost = {constrained_selection_tensor[*ankle[0], 0]} M€/PAX, GWP = {constrained_selection_tensor[*ankle[0], 1]} kg CO₂e/PAX")

    # # For the toe point
    # print_design_point("toe", tensor[*toe[0], ...], N_PAX[*toe[0]])
    # print(f"Configuration: Split = {splits[toe[0, 0]]}, TOGA = {toga_throttle[toe[0, 1]]}, Cruise = {cruise_throttle[toe[0, 2]]}")
    # print(f"Objectives: Cost = {constrained_selection_tensor[*toe[0], 0]} M€/PAX, GWP = {constrained_selection_tensor[*toe[0], 1]} kg CO₂e/PAX")
    
    # print_design_point("toenail", tensor[*toenail[0], ...], N_PAX[*toenail[0]])
    # print(f"Configuration: Split = {splits[toenail[0, 0]]:.2f}, TOGA = {toga_throttle[toenail[0, 1]]:.2f}, Cruise = {cruise_throttle[toenail[0, 2]]:.2f}")
    # print(f"Objectives: Cost = {constrained_selection_tensor[*toenail[0], 0]:.2f} M€/PAX, GWP = {constrained_selection_tensor[*toenail[0], 1]:.2f} kg CO₂e/PAX")
   
    # Get percentage change for pareto-optimal points
    x = pareto_options[..., 0].ravel()
    y = pareto_options[..., 1].ravel()

    # Filter out NaN values
    valid = np.isfinite(x) & np.isfinite(y)
    x_valid = x[valid]
    y_valid = y[valid]

    # Sort by x (cost) to ensure proper Pareto front ordering
    sort_idx = np.argsort(x_valid)
    x_sorted = x_valid[sort_idx]
    y_sorted = y_valid[sort_idx]

    # Calculate percentage changes and slopes
    x_percent_change = (x_sorted[3:] - x_sorted[:-3]) * 2 / (x_sorted[3:] + x_sorted[:-3])
    y_percent_change = (y_sorted[3:] - y_sorted[:-3]) * 2 / (y_sorted[3:] + y_sorted[:-3])

    # Calculate slopes, handling possible division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        slope = y_percent_change / x_percent_change
        
    # Filter out any remaining infinities or NaNs
    valid_slopes = np.isfinite(slope)
    clean_slopes = slope[valid_slopes]

    print(f"Number of valid Pareto points: {len(x_valid)}")
    print(f"Number of valid slopes: {len(clean_slopes)}")
    
    pareto_front(constrained_selection_tensor,
                 fpath="data/plots/pareto_front_no_star.png",
                 dpi=1200,
                #  star_boy=np.array([3.4413104293093415, 29.287180605687336]),
                # star_boy=np.array([3.425018682093658, 10.122170765960028])
                # star_boy=np.array([3.5051841306121516, 31.5128137833263])
                #  focus_point=(3.44, 29.29, 1)
                 )
    
    pareto_front(constrained_selection_tensor,
                 fpath="data/plots/pareto_front_star.png",
                 dpi=1200,
                 star_boy=np.array([3.4413104293093415, 29.287180605687336]),
                # star_boy=np.array([3.425018682093658, 10.122170765960028])
                # star_boy=np.array([3.5051841306121516, 31.5128137833263])
                #  focus_point=(3.44, 29.29, 1)
                 )
    
    # Use the corresponding x values for plotting slopes
    # Since slopes are calculated between adjacent points, 
    # we'll use the x values at the start of each segment
    # plt.plot(x_sorted[:-3][valid_slopes], clean_slopes, 'o', markersize=2, label='Pareto Points')
    # plt.xlabel('Cost [M€/PAX]')
    # plt.ylabel('Slope (GWP change / Cost change)')
    # plt.show()
        
    interactive_pareto_front(
        constrained_selection_tensor,
        tensor,
        splits,
        toga_throttle,
        cruise_throttle,
        N_PAX,
        fpath="data/plots/interactive_pareto.html",
        star_boy=np.array([3.44, 29.29])  # Optional: highlight a specific point
    )
    
    
if __name__ == "__main__":
    weight = 0.5  # Default weight for cost vs GWP
    main(
        weight=weight,
        # dim_bounds=np.array([
        #     [0.25, 0.35],
        #     [0.2, 0.35],
        #     [0.25, 0.4]
        #     ]),
        dim_bounds=np.array([
            [0.1, 0.9],
            [0.1, 1.0],
            [0.1, 1.0]
        ]),
        fpath_design='data/logs/result_tensor.pkl',
        fpath_loading='data/logs/loading_tensor.pkl'
    )