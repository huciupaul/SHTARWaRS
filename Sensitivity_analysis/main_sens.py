'''
Sensitivity analysis for DSE Tradeoffs

Instructions:
1. Copy the sens_input_template.py file to a new file named sens_input<<your_tradeoff>>.py.
2. Modify the new file to include your tradeoff criteria and weights, following the instructions given in the template.
3. Ensure that the new file is in the same directory as this script.
4. Run this script to see the results.

'''


import numpy as np

# Check if the input file is in the same directory
try:
    import sens_input_fuel_cell as sens_input
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

# Initialize the tradeoff winner array to count how many times each design option wins
output['Tradeoff Wins'] = np.zeros((len(sens_input.design_option_names)))


# Function to calculate the winning design option
def find_winner(tradeoff_values):
    """
    Find the design option with the highest tradeoff value.
    
    Parameters:
    tradeoff_values (numpy.ndarray): The tradeoff values for each design option.
    
    Returns:
    int: The index of the winning design option.
    """
    return np.argmax(tradeoff_values)



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

    # Preprocess values so the maximum value of each column is 1
    tradeoff_values = tradeoff_values / np.max(tradeoff_values, axis=0)

    # print("tradeoff values: ", tradeoff_values)

    tradeoff_values = np.dot(tradeoff_values, tradeoff_weights)
    
    # Find the winning design option
    winning_design_option = find_winner(tradeoff_values)
    # Increment the tradeoff wins for the winning design option
    output['Tradeoff Wins'][winning_design_option] += 1

    return tradeoff_values




# Create initial tradeoff results
output['Initial Tradeoff Values'] = calculate_tradeoff_values(sens_input.tradeoff_criteria_values, sens_input.tradeoff_criteria_weights)
# Initialize array to store amount of tradeoff wins




# Initialize total tradeoff values
total_tradeoff_values = np.array([[0,0,0]])

# Analyze the sensitivity of the tradeoff values to changes in the tradeoff criteria weights
# Loop over each tradeoff criterion
for i in range(len(sens_input.tradeoff_criteria_weights)):

    value_array = np.linspace(sens_input.tradeoff_criteria_weights[i] - sens_input.tradeoff_criteria_weight_margins[i],
                               sens_input.tradeoff_criteria_weights[i] + sens_input.tradeoff_criteria_weight_margins[i], 21)
    
    # Initialize tradeoff values for the current criterion
    tradeoff_values_this_criterion = np.array([[0,0,0]])

    # Initialize wins array for the current criterion
    tradeoff_wins_this_criterion = np.zeros((len(sens_input.design_option_names)))

    # Loop over each weight margin
    for j in range(len(value_array)):
        # Create a copy of the weights
        weights = np.copy(sens_input.tradeoff_criteria_weights)
        # Set the current weight to the new value
        weights[i] = value_array[j]
        # Calculate the tradeoff values
        tradeoff_values = calculate_tradeoff_values(sens_input.tradeoff_criteria_values, weights)

        # Find the winning design option
        winning_design_option = find_winner(tradeoff_values)
        # Increment the tradeoff wins for the winning design option
        tradeoff_wins_this_criterion[winning_design_option] += 1
    
        # Append the tradeoff values to the tradeoff values for this criterion
        tradeoff_values_this_criterion = np.append(tradeoff_values_this_criterion, [tradeoff_values], axis = 0)
        
        # Append the tradeoff values to the total tradeoff values array
        total_tradeoff_values = np.append(total_tradeoff_values, [tradeoff_values], axis=0)

    # Store the tradeoff values in the output dictionary
    output[f'Tradeoff Values for {sens_input.tradeoff_criteria_names[i]}'] = tradeoff_values_this_criterion

    # Store the tradeoff wins in the output dictionary
    output[f'Tradeoff Wins for {sens_input.tradeoff_criteria_names[i]}'] = tradeoff_wins_this_criterion


# Analyze the sensitivity of the tradeoff values to changes in the tradeoff scores for each design option
# Loop over each design option
for i in range(len(sens_input.design_option_names)):
    # Create a copy of the tradeoff values
    tradeoff_values = np.copy(sens_input.tradeoff_criteria_values)
    
    # Loop over each tradeoff criterion
    for j in range(len(sens_input.tradeoff_criteria_names)):
        # Create a copy of the tradeoff values for this design option
        tradeoff_values_this_design_option = np.copy(tradeoff_values[i])
        
        # Create an array of values for the current design option
        value_array = np.linspace(tradeoff_values_this_design_option[j] * (1 - 2*sens_input.tradeoff_criteria_value_margins[i][j]),
                                   tradeoff_values_this_design_option[j] * (1 + 2*sens_input.tradeoff_criteria_value_margins[i][j]), 21)
        
        # Initialize tradeoff values for the current design option
        tradeoff_values_this_design_option = np.array([[0,0,0]])

        # Initialize wins array for the current design option
        tradeoff_wins_this_design_option = np.zeros((len(sens_input.design_option_names)))

        tradeoff_values[i][j] = 1
        
        # Loop over each value margin
        for k in range(len(value_array)):
            # Set the current value to the new value
            # print(f"Design Option: {i}, Tradeoff Criterion: {j}, Value: {value_array[k]}")
            tradeoff_values[i][j] = value_array[k]
            # Calculate the tradeoff values
            tradeoff_scores = calculate_tradeoff_values(tradeoff_values, sens_input.tradeoff_criteria_weights)

            # Find the winning design option
            winning_design_option = find_winner(tradeoff_scores)
            # Increment the tradeoff wins for the winning design option
            tradeoff_wins_this_design_option[winning_design_option] += 1
            
            # Append the tradeoff values to the tradeoff values for this design option
            tradeoff_values_this_design_option = np.append(tradeoff_values_this_design_option, [tradeoff_scores], axis=0)
            # Append the tradeoff values to the total tradeoff values array
            total_tradeoff_values = np.append(total_tradeoff_values, [tradeoff_scores], axis=0)

        # Store the tradeoff values in the output dictionary
        output[f'Tradeoff Values for {sens_input.design_option_names[i]}'] = tradeoff_values_this_design_option

        # Store the tradeoff wins in the output dictionary
        output[f'Tradeoff Wins for {sens_input.design_option_names[i]} and {sens_input.tradeoff_criteria_names[j]}'] = tradeoff_wins_this_design_option




# Store the total tradeoff values in the output dictionary
output['Total Tradeoff Values'] = total_tradeoff_values


print(output)

print("winning design options: ", output['Tradeoff Wins'])

