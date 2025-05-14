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

names = []
criteria = []
for i in range(len(sens_input.design_option_names)):
    names.append(sens_input.design_option_names[i])

for i in range(len(sens_input.tradeoff_criteria_names)):
    criteria.append(sens_input.tradeoff_criteria_names[i])


'''
dict_keys(['design_option_names', 'tradeoff_criteria_names', 'tradeoff_criteria_weights', 'tradeoff_criteria_weight_margins', 'tradeoff_criteria_values', 'tradeoff_criteria_value_margins', 'Tradeoff Wins',
           'Initial Tradeoff Values', 'Tradeoff Values for Sustainability', 'Tradeoff Wins for Sustainability', 'Tradeoff Values for Cost', 'Tradeoff Wins for Cost',
           'Tradeoff Values for Efficiency', 'Tradeoff Wins for Efficiency', 'Tradeoff Values for Specific power', 'Tradeoff Wins for Specific power',
           'Tradeoff Values for Aircraft integration', 'Tradeoff Wins for Aircraft integration', 'Tradeoff Values for LTPEM', 'Tradeoff Wins for LTPEM and Sustainability',
           'Tradeoff Wins for LTPEM and Cost', 'Tradeoff Wins for LTPEM and Efficiency', 'Tradeoff Wins for LTPEM and Specific power', 'Tradeoff Wins for LTPEM and Aircraft integration',
           'Tradeoff Values for HTPEM', 'Tradeoff Wins for HTPEM and Sustainability', 'Tradeoff Wins for HTPEM and Cost', 'Tradeoff Wins for HTPEM and Efficiency', 'Tradeoff Wins for HTPEM and Specific power',
           'Tradeoff Wins for HTPEM and Aircraft integration', 'Tradeoff Values for SOFC', 'Tradeoff Wins for SOFC and Sustainability', 'Tradeoff Wins for SOFC and Cost', 
           'Tradeoff Wins for SOFC and Efficiency', 'Tradeoff Wins for SOFC and Specific power', 'Tradeoff Wins for SOFC and Aircraft integration', 'Total Tradeoff Values'])

'''

#compiles the amount of time a design has won for sensitivity analysis of criteria
win_matrix = []
for j in range(len(sens_input.design_option_names)):
    d = []
    for i in range(len(sens_input.tradeoff_criteria_names)):
        x = main.output['Tradeoff Wins for ' + sens_input.design_option_names[j] + ' and ' + sens_input.tradeoff_criteria_names[i]][Which_to_check]
        d.append(x)
    win_matrix.append(d)
win_matrix = np.array(win_matrix)

wins = main.output['Tradeoff Wins']

#graph winner total and percentage per weight in bar chart
plt.figure()
plt.bar(names, wins, color='blue')
plt.title("Number of wins per design option")
plt.show()


#compares the wins of one design to the total amount it could have won
win_frac_sus = 100 * main.output['Tradeoff Wins for Sustainability'][1] /  sum(main.output['Tradeoff Wins for Sustainability'])
percent = []
for i in range(len(sens_input.tradeoff_criteria_names)):
    percent.append(100 * main.output['Tradeoff Wins for ' + sens_input.tradeoff_criteria_names[i]][Which_to_check] /  sum(main.output['Tradeoff Wins for ' + sens_input.tradeoff_criteria_names[i]]))

plt.figure()
plt.bar(criteria, percent, color='green')
plt.title('Percentage of wins')
plt.ylabel("Values")
plt.show()


def visualize_red_yellow_green_above_half(data, row_labels, col_labels):
    data = np.array(data)

    if data.shape != (len(row_labels), len(col_labels)):
        raise ValueError("Shape of data does not match number of row and column labels")

    max_val = np.max(data)
    half_max = max_val / 2
    three_quarter_max = 0.75 * max_val

    if max_val == 0:
        raise ValueError("All data values are zero — can't compute color gradient.")

    def value_to_color(val):
        if val <= half_max:
            return (1, 0, 0)  # Solid red
        elif val <= three_quarter_max:
            # red → yellow (1,0,0) → (1,1,0)
            norm = (val - half_max) / (three_quarter_max - half_max + 1e-5)
            r = 1
            g = norm
            b = 0
            return (r, g, b)
        else:
            # yellow → green (1,1,0) → (0,1,0)
            norm = (val - three_quarter_max) / (max_val - three_quarter_max + 1e-5)
            r = 1 - norm
            g = 1
            b = 0
            return (r, g, b)

    cell_colors = [[value_to_color(data[i, j]) for j in range(data.shape[1])] for i in range(data.shape[0])]

    fig, ax = plt.subplots()
    ax.set_axis_off()

    table = ax.table(
        cellText=data.astype(int),
        cellColours=cell_colors,
        rowLabels=row_labels,
        colLabels=col_labels,
        loc='center',
        cellLoc='right',
    )

    table.scale(1, 2)
    plt.title("Color: ≤50% Max = Red, >50% = Red → Yellow → Green", fontsize=14)
    plt.show()

data = win_matrix

sens_input.design_option_names[j]

row_labels = names
col_labels = criteria

visualize_red_yellow_green_above_half(data, row_labels, col_labels)



#print(main.total_tradeoff_values)
#print(np.shape(main.total_tradeoff_values))
#print(main.total_tradeoff_values[1,0])






