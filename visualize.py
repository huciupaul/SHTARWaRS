import numpy as np
import pickle
import os
import imageio.v2 as imageio
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib as mpl

# Local imports
from global_constants import NOx_pPax_TO, NOx_pPax_cruise
from performance.MTOW.mtow import power_and_wing_loading as pws
from performance.Integration.integration import main as cg_excursion
from cost import calc_cost as f_cost
from sus import get_sus as f_gwp

# Define a function to create the visualization with given mask
def create_visualization(valid_mask, output_filename, title):
    # Get indices
    x_idx, y_idx, z_idx = np.where(valid_mask)
    
    # Skip if no valid points
    if len(x_idx) == 0:
        print(f"No valid points for {title}. Skipping...")
        return
    
    # Convert indices to actual coordinate values
    x_coords = power_splits[x_idx]
    y_coords = toga_throttle[y_idx] 
    z_coords = cruise_throttle[z_idx]
    
    # Generate frames for the rotation
    print(f"Generating frames for {title}...")
    num_frames = 120  # For a smooth animation (120 frames = 4s at 30fps)
    
    # Set figure size (in inches) and DPI for higher quality
    fig_width, fig_height = 10, 10
    dpi = 150
    
    # Create color gradient based on position
    colormap = cm.viridis
    colors = colormap((x_coords - dim_bounds[0, 0]) / (dim_bounds[0, 1] - dim_bounds[0, 0]))
    
    # Define starting and ending angles for exactly one rotation
    start_angle = 0
    end_angle = 360  # Full 360 rotation
    
    # Create frames directory if it doesn't exist
    frames_dir = f"frames_{output_filename.replace('.gif', '')}"
    if not os.path.exists(frames_dir):
        os.makedirs(frames_dir)
    
    for i in range(num_frames):
        # Create a new figure for each frame with fixed size
        fig = plt.figure(figsize=(fig_width, fig_height), dpi=dpi, facecolor='black')
        ax = fig.add_subplot(111, projection='3d', facecolor='black')
        
        # Enhanced scatter plot
        scatter = ax.scatter(x_coords, y_coords, z_coords, 
                   c=colors,
                   s=50,  # Larger marker size for better visibility
                   alpha=0.8,
                   marker='o',  # Circular markers look more modern
                   edgecolors='white',  # Add white edges to markers
                   linewidth=0.5)
        
        # Add a subtle grid for better spatial perception
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        
        ax.xaxis.pane.set_edgecolor('white')
        ax.yaxis.pane.set_edgecolor('white')
        ax.zaxis.pane.set_edgecolor('white')
        
        ax.xaxis.pane.set_alpha(0.2)
        ax.yaxis.pane.set_alpha(0.2)
        ax.zaxis.pane.set_alpha(0.2)
        
        ax.grid(True, alpha=0.2)
        
        # Add title
        # plt.title(title, fontsize=30, color='white', pad=20)
        
        # Improved axis labels with custom styling
        ax.set_xlabel('Power Split', fontsize=30, labelpad=20, color='white')
        ax.set_ylabel('TOGA Throttle', fontsize=30, labelpad=20, color='white')
        ax.set_zlabel('Cruise Throttle', fontsize=30, labelpad=20, color='white')
        
        # Make tick labels white
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        ax.tick_params(axis='z', colors='white')
        
        # Set axis limits
        ax.set_xlim(dim_bounds[0, 0], dim_bounds[0, 1])
        ax.set_ylim(dim_bounds[1, 0], dim_bounds[1, 1])
        ax.set_zlim(dim_bounds[2, 0], dim_bounds[2, 1])
        
        # Set the view angle - simplified to ensure exactly one rotation
        # Linear interpolation from start_angle to end_angle
        azim = start_angle + i * (end_angle - start_angle) / (num_frames - 1)
        elev = 25
        ax.view_init(elev=elev, azim=azim)
        
        # Adjust axis aspect ratio for better cube shape
        ax.set_box_aspect([1, 1, 1])
        
        # Leave some padding around the plot
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
        
        # Save the frame WITHOUT bbox_inches='tight' to ensure consistent size
        frame_filename = f'{frames_dir}/frame_{i:03d}.png'
        plt.savefig(frame_filename, facecolor='black')
        plt.close(fig)  # Close the figure to free memory
        print(f"Saved frame {i+1}/{num_frames}")
    
    # Create GIF from frames
    print(f"Creating GIF for {title}...")
    images = []
    for i in range(num_frames):
        filename = f'{frames_dir}/frame_{i:03d}.png'
        images.append(imageio.imread(filename))
    
    # Save the GIF with smoother animation - set loop=0 for no looping
    imageio.mimsave(output_filename, images, duration=0.033, loop=1)
    print(f"GIF saved as {output_filename}")


# Set better visualization style
plt.style.use('dark_background')
mpl.rcParams['font.family'] = 'Avenir'
mpl.rcParams['font.size'] = 25

# Load data
print("Loading tensors...")
with open('data/logs/result_tensor.pkl', 'rb') as f:
    design_tensor = pickle.load(f)

with open('data/logs/loading_tensor.pkl', 'rb') as f:
    loading_tensor = pickle.load(f)

# Define dimension bounds
dim_bounds = np.array([[0.1, 0.9], [0.1, 1.0], [0.1, 1.0]])
N, M, P, Q = design_tensor.shape
power_splits = np.linspace(dim_bounds[0, 0], dim_bounds[0, 1], N)
toga_throttle = np.linspace(dim_bounds[1, 0], dim_bounds[1, 1], M)
cruise_throttle = np.linspace(dim_bounds[2, 0], dim_bounds[2, 1], P)

# Create unconstrained visualization
unconstrained_mask = np.ones(design_tensor.shape[:-1], dtype=bool)
create_visualization(unconstrained_mask, 'design_space_unconstrained.gif', 'Design Space (No Constraints)')

# Create constrained visualization
try:
    print("Applying constraints...")
    # The cg_excursion function returns a tuple where the first element is the mask
    valid_integration_tuple = cg_excursion(design_tensor)
    if isinstance(valid_integration_tuple, tuple):
        valid_integration = valid_integration_tuple[0]
    else:
        valid_integration = valid_integration_tuple
    
    # Ensure valid_integration is boolean
    valid_integration = np.array(valid_integration, dtype=bool)
    
    # The pws function might also return a tuple
    valid_mass_result = pws(loading_tensor)
    if isinstance(valid_mass_result, tuple):
        valid_mass = valid_mass_result[0]
    else:
        valid_mass = valid_mass_result
    
    # Ensure valid_mass is boolean
    valid_mass = np.array(valid_mass, dtype=bool)
    
    # Apply constraints
    constrained_mask = unconstrained_mask.copy()
    constrained_mask &= valid_integration
    constrained_mask &= valid_mass
    
    create_visualization(constrained_mask, 'design_space_constrained.gif', 'Design Space (With Constraints)')
    
except Exception as e:
    print(f"Error applying constraints: {e}")
    print("Only the unconstrained visualization was created.")