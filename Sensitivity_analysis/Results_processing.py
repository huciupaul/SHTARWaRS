import numpy as np
import matplotlib.pyplot as plt

# Check if the input file is in the same directory
try:
    import main_sens as main
except ImportError:
    raise ImportError("The input file 'sens_input_template.py' is not in the same directory. Please copy the template and modify it as instructed.")
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
print(main.output['Tradeoff Wins for Specific power'])

# sustainability total effect on winners
win_frac_sus = 100 * main.output['Tradeoff Wins for Sustainability'][1] /  sum(main.output['Tradeoff Wins for Sustainability'])
win_frac_eff = 100 * main.output['Tradeoff Wins for Efficiency'][1] /  sum(main.output['Tradeoff Wins for Efficiency'])
win_frac_spp = 100 * main.output['Tradeoff Wins for Specific power'][1] /  sum(main.output['Tradeoff Wins for Specific power'])
win_frac_cost = 100 * main.output['Tradeoff Wins for Cost'][1] /  sum(main.output['Tradeoff Wins for Cost'])
win_frac_int = 100 * main.output['Tradeoff Wins for Aircraft integration'][1] /  sum(main.output['Tradeoff Wins for Aircraft integration'])

wins = main.output['Tradeoff Wins']
namesw = ['LT-PEM', 'HT-PEM', 'SOFC']

plt.figure()
plt.bar(namesw, wins, color='blue')
plt.title("Set 1")
plt.ylabel("Values")
plt.show()

names = ['Sustainability', 'efficiency', 'Specific power', 'lifetime cost', 'integration']
percent = [win_frac_sus, win_frac_eff, win_frac_spp, win_frac_cost, win_frac_int]

plt.figure()
plt.bar(names, percent, color='green')
plt.title('Percentage of wins')
plt.ylabel("Values")
plt.show()


#print(main.total_tradeoff_values)
#print(np.shape(main.total_tradeoff_values))
#print(main.total_tradeoff_values[1,0])






