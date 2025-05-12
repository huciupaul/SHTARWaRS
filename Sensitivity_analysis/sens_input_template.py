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
design_option_names = [
    'Design Option 1',
    'Design Option 2',
    'Design Option 3',
    'Design Option 4',
    'Design Option 5'
]

# Tradeoff criteria names
tradeoff_criteria_names = [
    'Sustainability',           # Sustainability
    'Cost',                     # Cost
    'Tradeoff Criterion 3',     # Tradeoff Criterion 3
    'Tradeoff Criterion 4',     # Tradeoff Criterion 4
    'Tradeoff Criterion 5'      # Tradeoff Criterion 5
]

# Tradeoff criteria weights
tradeoff_criteria_weights = [
    0.15, # Sustainability
    0.1,  # Cost
    0.2,  # Tradeoff Criterion 3
    0.45, # Tradeoff Criterion 4
    0.1   # Tradeoff Criterion 5
]

# Tradeoff criteria weight margins
tradeoff_criteria_weight_margins = [
    0.05,  # Sustainability
    0.05,  # Cost
    0.05,  # Tradeoff Criterion 3
    0.05,  # Tradeoff Criterion 4
    0.05   # Tradeoff Criterion 5
]

# Tradeoff criteria values

tradeoff_criteria_values = [
    [0.8, 0.7, 0.6, 0.5, 0.4],  # Design Option 1
    [0.7, 0.6, 0.5, 0.4, 0.3],  # Design Option 2
    [0.6, 0.5, 0.4, 0.3, 0.2],  # Design Option 3
    [0.5, 0.4, 0.3, 0.2, 0.1],  # Design Option 4
    [0.4, 0.3, 0.2, 0.1, 0]     # Design Option 5
]

# Tradeoff criteria value margins
tradeoff_criteria_value_margins = [
    [0.1, 0.1, 0.1, 0.1, 0.1],  # Design Option 1
    [0.1, 0.1, 0.1, 0.1, 0.1],  # Design Option 2
    [0.1, 0.1, 0.1, 0.1, 0.1],  # Design Option 3
    [0.1, 0.1, 0.1, 0.1, 0.1],  # Design Option 4
    [0.1, 0.1, 0.1, 0.1, 0.1]   # Design Option 5
]