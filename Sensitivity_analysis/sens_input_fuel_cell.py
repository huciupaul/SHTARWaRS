'''
Sensitivity analysis input template
DO NOT MODIFY THIS FILE, ADD YOUR TRADEOFFS TO A COPY OF THIS FILE
This file is a template for the sensitivity analysis input.

There are four arrays and two matrices which must be defined:
    Arrays:
    - design_option_names: A list of strings representing the names of the design options.
    - tradeoff_criteria_names: A list of strings representing the names of the tradeoff criteria.
    - tradeoff_criteria_weights: A list of floats representing the weights of the tradeoff criteria, given in % of total weight.
    - tradeoff_criteria_weight_margins: A list of floats representing the margins of the tradeoff criteria weights,
        in change in proportion of total weight. (if the weight is 0.1 (10%) and the margin is 0.05 (5%), the weight will be between 5% and 15%)
    Matrices:
    - tradeoff_criteria_values: A matrix of floats representing the values of the tradeoff criteria for each design option.
        The matrix should be a 2D list, where each row represents a design option and each column represents a tradeoff criterion.
    - tradeoff_criteria_value_margins: A matrix of floats representing the margins of the tradeoff criteria values,
        in change in % of the value. (if the value is 10 and the margin is 0.2 (20%), the value will be between 8 and 12)
        The matrix should be a 2D list, where each row represents a design option and each column represents a tradeoff criterion.

        
The size of the arrays and matrices must be consistent, but can be of size 3-5.
'''


import numpy as np



# Design option names
design_option_names = np.array([
    'LTPEM',  # Low Temperature Proton Exchange Membrane
    'HTPEM',  # High Temperature Proton Exchange Membrane
    'SOFC'    # Solid Oxide Fuel Cell
])

# Tradeoff criteria names
tradeoff_criteria_names = np.array([
    'Sustainability',           # Sustainability
    'Cost',                     # Cost
    'Efficiency',               # Tradeoff Criterion 3
    'Specific power',           # Tradeoff Criterion 4
    'Aircraft integration'      # Tradeoff Criterion 5
])

# Tradeoff criteria weights
tradeoff_criteria_weights = np.array([
    0.15, # Sustainability
    0.1,  # Cost
    0.2,  # Tradeoff Criterion 3
    0.45, # Tradeoff Criterion 4
    0.1   # Tradeoff Criterion 5
])

# Tradeoff criteria weight margins
tradeoff_criteria_weight_margins = np.array([
    0.1,  # Sustainability
    0.1,  # Cost
    0.15,  # Tradeoff Criterion 3
    0.2,  # Tradeoff Criterion 4
    0.1   # Tradeoff Criterion 5
])

# Tradeoff criteria values

tradeoff_criteria_values = np.array([
    [1, 0.886120996441281, .6, 0.677987421383648, .7],  # Design Option 1
    [1, 0.713467048710602, .6, 0.677987421383648, 1],  # Design Option 2
    [.5, 0.399038461538462, .65, 0.493081761006289,  .3]  # Design Option 3
])

# Tradeoff criteria value margins
tradeoff_criteria_value_margins = np.array([
    [0.2, 0.05, 0.1, 0.1, 0.2],  # Design Option 1
    [0.2, 0.15, 0.1, 0.1, 0.2],  # Design Option 2
    [0.2, 0.35, 0.1, 0.1, 0.2]  # Design Option 3
])