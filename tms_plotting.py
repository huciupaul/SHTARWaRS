import matplotlib.pyplot as plt
import numpy as np
import re
import pandas as pd
from matplotlib.colors import to_rgba


def fluids_table(fluids):
    # Filter fluids for hydrogen and coolant
    selected_fluids = fluids
    labels = []
    mass_flows = []
    pressures = []
    temperatures = []
    for key in fluids.keys():
        f = fluids[key]
        labels.append(key)
        mass_flows.append(f['massflow'])
        pressures.append(f['pressure'])
        temperatures.append(f['temperature'])
    # Convert pressure from Pascal to bar and round values to 4 significant figures
    mass_flows = [float(f"{v:.4g}") for v in mass_flows]
    pressures = [float(f"{v/1e5:.4g}") for v in pressures]  # 1 bar = 1e5 Pa
    temperatures = [float(f"{v:.4g}") for v in temperatures]
    df = pd.DataFrame({
        'Fluid': labels,
        'Mass Flow': mass_flows,
        'Pressure': pressures,
        'Temperature': temperatures
    })

    # Add a 'Type' column based on fluid name
    def fluid_type(name):
        if re.search(r'h2', name, re.IGNORECASE):
            return 'H2'
        elif re.search(r'cool', name, re.IGNORECASE):
            return 'Coolant'
        else:
            return 'Other'

    df['Type'] = df['Fluid'].apply(fluid_type)

    # Sort by Type then Fluid
    df = df.sort_values(['Type', 'Fluid']).reset_index(drop=True)

    # Prepare cell colors
    def get_colors(values):
        colors = []
        prev = None
        for v in values:
            if prev is None:
                colors.append('white')
            elif v > prev:
                colors.append('#b6fcb6')  # light green
            elif v < prev:
                colors.append('#fcb6b6')  # light red
            else:
                colors.append('#f0f0f0')  # gray
            prev = v
        return colors

    cell_colors = []
    for col in ['Mass Flow', 'Pressure', 'Temperature']:
        cell_colors.append(get_colors(df[col].values))

    # Add Fluid and Type columns as white
    cell_colors = [['white'] * len(df)] + cell_colors + [['white'] * len(df)]

    # Transpose to match table shape
    cell_colors = list(map(list, zip(*cell_colors)))

    # Plot table
    fig, ax = plt.subplots(figsize=(8, 0.5 + 0.5 * len(df)))
    ax.axis('off')
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellColours=cell_colors,
        loc='center'
    )
    # Make row height larger
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    for key, cell in table.get_celld().items():
        cell.set_height(0.04)  # Increase this value for larger rows
    table.auto_set_column_width(col=list(range(len(df.columns))))
    fig.suptitle('Fluid Properties Table', fontsize=14)
    plt.show()