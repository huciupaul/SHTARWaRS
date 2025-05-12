'''
Sensitivity analysis for DSE Tradeoffs

Instructions:
1. Copy the sens_input_template.py file to a new file named sens_input<<your_tradeoff>>.py.
2. Modify the new file to include your tradeoff criteria and weights, following the instructions given in the template.
3. Ensure that the new file is in the same directory as this script.
4. Run this script to see the results.

'''


from sensitivity import SensitivityAnalyzer
import numpy as np

# Check if the input file is in the same directory
try:
    import sens_input_template as sens_input
except ImportError:
    raise ImportError("The input file 'sens_input_template.py' is not in the same directory. Please copy the template and modify it as instructed.")


# Check if the input file has the required variables andf if they are in the correct format
# Check if the input file has the required variables

required_variables = [
    'design_option_names',
    'tradeoff_criteria_names',
    'tradeoff_criteria_weights',
    'tradeoff_criteria_weight_margins',
    'tradeoff_criteria_values',
    'tradeoff_criteria_value_margins'
]

for var in required_variables:
    if not hasattr(sens_input, var):
        raise ValueError(f"The input file is missing the required variable: {var}. Please ensure all required variables are defined.")
    else:
        print(f"Variable '{var}' is present in the input file.")

# Check if the input arrays and matrices are the correct size

try:
    if len(sens_input.design_option_names) != len(sens_input.tradeoff_criteria_values):
        raise ValueError("The number of design options does not match the number of rows in the tradeoff criteria values matrix.")

    if len(sens_input.design_option_names) != len(sens_input.tradeoff_criteria_value_margins):
        raise ValueError("The number of design options does not match the number of rows in the tradeoff criteria value margins matrix.") 
    
    if len(sens_input.tradeoff_criteria_names) != len(sens_input.tradeoff_criteria_values[0]):
        raise ValueError("The number of tradeoff criteria does not match the number of columns in the tradeoff criteria values matrix.")

    if len(sens_input.tradeoff_criteria_names) != len(sens_input.tradeoff_criteria_value_margins[0]):
        raise ValueError("The number of tradeoff criteria does not match the number of columns in the tradeoff criteria value margins matrix.")
    
    if len(sens_input.tradeoff_criteria_weights) != len(sens_input.tradeoff_criteria_names):
        raise ValueError("The number of tradeoff criteria weights does not match the number of tradeoff criteria.")
    
    if len(sens_input.tradeoff_criteria_weight_margins) != len(sens_input.tradeoff_criteria_names):
        raise ValueError("The number of tradeoff criteria weight margins does not match the number of tradeoff criteria.")
except ValueError as e:
    raise ValueError(f"Input file validation error: {e}")

# Create output dictionary
output = {
    'design_option_names': sens_input.design_option_names,
    'tradeoff_criteria_names': sens_input.tradeoff_criteria_names,
    'tradeoff_criteria_weights': sens_input.tradeoff_criteria_weights,
    'tradeoff_criteria_weight_margins': sens_input.tradeoff_criteria_weight_margins,
    'tradeoff_criteria_values': sens_input.tradeoff_criteria_values,
    'tradeoff_criteria_value_margins': sens_input.tradeoff_criteria_value_margins
}

# Tradeoff Value calculation
# Calculate the tradeoff values for each design option
def calculate_tradeoff_values(tradeoff_values, tradeoff_weights):
    """
    Calculate the tradeoff values for each design option.
    
    Parameters:
    tradeoff_criteria_values (numpy.ndarray): The tradeoff criteria values matrix.
    tradeoff_criteria_weights (numpy.ndarray): The tradeoff criteria weights array.
    
    Returns:
    numpy.ndarray: The tradeoff values for each design option.
    """
    # Preprocess weights so they sum to 1
    tradeoff_weights = tradeoff_weights / np.sum(tradeoff_weights)


    return np.dot(tradeoff_values, tradeoff_weights)

output['Initial Tradeoff Values'] = calculate_tradeoff_values(sens_input.tradeoff_criteria_values, sens_input.tradeoff_criteria_weights)


# Analyze the sensitivity of the tradeoff values to changes in the tradeoff criteria weights
for i in range(len(sens_input.tradeoff_criteria_weights)):
    # 



print("Calculating tradeoff values...")
tradeoff_values = calculate_tradeoff_values(sens_input.tradeoff_criteria_values, sens_input.tradeoff_criteria_weights)
print(tradeoff_values)