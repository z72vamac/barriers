import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


datos = pd.read_excel('first_experiment.xlsx')

tabla_comparadora = datos.groupby(['n_N', 'Neighborhood', 'Preprocessing']).mean()[['n_B', 'Gap', 'Runtime']].round(2).reset_index()

print(tabla_comparadora)

tabla_comparadora.to_excel('summary_first_experiment.xlsx')