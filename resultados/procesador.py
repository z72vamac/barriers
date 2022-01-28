import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

with_prepro_segmentos = pd.read_csv('resultados_with_prepro_segmentos_visitar.csv')
with_prepro_segmentos['Preprocessing'] = True
with_prepro_segmentos['Neighborhood'] = 'Segment'

without_prepro_segmentos = pd.read_csv('resultados_without_prepro_segmentos_visitar.csv')
without_prepro_segmentos['Preprocessing'] = False
without_prepro_segmentos['Neighborhood'] = 'Segment'

with_prepro_bolas = pd.read_csv('resultados_with_prepro_bolas_visitar.csv')
with_prepro_bolas['Preprocessing'] = True
with_prepro_bolas['Neighborhood'] = 'Ball'

without_prepro_bolas = pd.read_csv('resultados_without_prepro_bolas_visitar.csv')
without_prepro_bolas['Preprocessing'] = False
without_prepro_bolas['Neighborhood'] = 'Ball'

datos = pd.concat([without_prepro_segmentos, with_prepro_segmentos, without_prepro_bolas, with_prepro_bolas])

datos.to_excel('first_experiment.xlsx')

sns.set(style="darkgrid")

g = sns.catplot(x='n_N', y='Runtime', hue='Preprocessing', col='Neighborhood', kind='box', data=datos, aspect=0.5, sharey=True)

g.set(yscale='log')

import tikzplotlib

matplotlib.rcParams['axes.unicode_minus'] = False

tikzplotlib.save('runtime_boxplot_first_experiment.tex', encoding='utf-8')

plt.savefig('runtime_boxplot_first_experiment.png')
