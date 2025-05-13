import numpy as np
import matplotlib.pyplot as plt
from docutils.parsers.rst.directives.tables import ListTable

# Check if the input file is in the same directory
try:
    import sens_input_fuel_cell as sens_input
except ImportError:
    raise ImportError("The input file 'sens_input_template.py' is not in the same directory. Please copy the template and modify it as instructed.")

try:
    import main_sens as main
except ImportError:
    raise ImportError("The input file 'main_sens' is not in the same directory. Please copy the template and modify it as instructed.")


# input is what design will be selected for the table to check how many wins it has
Which_to_check = 1


'''
dict_keys(['design_option_names', 'tradeoff_criteria_names', 'tradeoff_criteria_weights', 'tradeoff_criteria_weight_margins', 'tradeoff_criteria_values', 'tradeoff_criteria_value_margins', 'Tradeoff Wins',
           'Initial Tradeoff Values', 'Tradeoff Values for Sustainability', 'Tradeoff Wins for Sustainability', 'Tradeoff Values for Cost', 'Tradeoff Wins for Cost',
           'Tradeoff Values for Efficiency', 'Tradeoff Wins for Efficiency', 'Tradeoff Values for Specific power', 'Tradeoff Wins for Specific power',
           'Tradeoff Values for Aircraft integration', 'Tradeoff Wins for Aircraft integration', 'Tradeoff Values for LTPEM', 'Tradeoff Wins for LTPEM and Sustainability',
           'Tradeoff Wins for LTPEM and Cost', 'Tradeoff Wins for LTPEM and Efficiency', 'Tradeoff Wins for LTPEM and Specific power', 'Tradeoff Wins for LTPEM and Aircraft integration',
           'Tradeoff Values for HTPEM', 'Tradeoff Wins for HTPEM and Sustainability', 'Tradeoff Wins for HTPEM and Cost', 'Tradeoff Wins for HTPEM and Efficiency', 'Tradeoff Wins for HTPEM and Specific power',
           'Tradeoff Wins for HTPEM and Aircraft integration', 'Tradeoff Values for SOFC', 'Tradeoff Wins for SOFC and Sustainability', 'Tradeoff Wins for SOFC and Cost', 
           'Tradeoff Wins for SOFC and Efficiency', 'Tradeoff Wins for SOFC and Specific power', 'Tradeoff Wins for SOFC and Aircraft integration', 'Total Tradeoff Values'])


# total effect on winners
win_frac_sus = 100 * main.output['Tradeoff Wins for Sustainability'][1] /  sum(main.output['Tradeoff Wins for Sustainability'])
win_frac_eff = 100 * main.output['Tradeoff Wins for Efficiency'][1] /  sum(main.output['Tradeoff Wins for Efficiency'])
win_frac_spp = 100 * main.output['Tradeoff Wins for Specific power'][1] /  sum(main.output['Tradeoff Wins for Specific power'])
win_frac_cost = 100 * main.output['Tradeoff Wins for Cost'][1] /  sum(main.output['Tradeoff Wins for Cost'])
win_frac_int = 100 * main.output['Tradeoff Wins for Aircraft integration'][1] /  sum(main.output['Tradeoff Wins for Aircraft integration'])

'design_option_names': sens_input.design_option_names,
    'tradeoff_criteria_names': sens_input.tradeoff_criteria_names,
    'tradeoff_criteria_weights': sens_input.tradeoff_criteria_weights,
    'tradeoff_criteria_weight_margins': sens_input.tradeoff_criteria_weight_margins,
    'tradeoff_criteria_values': sens_input.tradeoff_criteria_values,
    'tradeoff_criteria_value_margins': sens_input.tradeoff_criteria_value_margins
}
'''


win_matrix = []
for j in range(len(sens_input.design_option_names)):
    d = []
    for i in range(len(sens_input.tradeoff_criteria_names)):
        x = main.output['Tradeoff Wins for ' + sens_input.design_option_names[j] + ' and ' + sens_input.tradeoff_criteria_names[i]][Which_to_check]
        d.append(x)
    win_matrix.append(d)
win_matrix = np.array(win_matrix)
print(win_matrix)



'''
wins = main.output['Tradeoff Wins']
namesw = ['LT-PEM', 'HT-PEM', 'SOFC']


#graph winner total and percentage per weight in bar chart
plt.figure()
plt.bar(namesw, wins, color='blue')
plt.title("W")
plt.ylabel("Values")
plt.show()

names = ['Sustainability', 'efficiency', 'Specific power', 'lifetime cost', 'integration']
percent = [win_frac_sus, win_frac_eff, win_frac_spp, win_frac_cost, win_frac_int]

plt.figure()
plt.bar(names, percent, color='green')
plt.title('Percentage of wins')
plt.ylabel("Values")
plt.show()
'''


def visualize_redness(data, row_labels, col_labels):
    data = np.array(data)

    if data.shape != (len(row_labels), len(col_labels)):
        raise ValueError("Shape of data does not match number of row and column labels")

    # Normalize values for color (lower values = more red)
    norm = (data - np.min(data)) / (np.max(data) - np.min(data) + 1e-5)
    cell_colors = [[(1, norm[i, j], norm[i, j]) for j in range(data.shape[1])] for i in range(data.shape[0])]

    fig, ax = plt.subplots()
    ax.set_axis_off()

    table = ax.table(
        cellText=data.astype(int),
        cellColours=cell_colors,
        rowLabels=row_labels,
        colLabels=col_labels,
        loc='center',
        cellLoc='center'
    )

    table.scale(1, 2)  # Adjust for better visibility
    plt.title("Redder = Lower Value", fontsize=14)
    plt.show()


data = win_matrix
names = []
criteria = []
for i in range(len(sens_input.design_option_names)):
    names.append(sens_input.design_option_names[i])

for i in range(len(sens_input.tradeoff_criteria_names)):
    criteria.append(sens_input.tradeoff_criteria_names[i])


sens_input.design_option_names[j]

row_labels = names
col_labels = criteria

visualize_redness(data, row_labels, col_labels)



#print(main.total_tradeoff_values)
#print(np.shape(main.total_tradeoff_values))
#print(main.total_tradeoff_values[1,0])






