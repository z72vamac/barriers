import gurobipy as gp
import pdb
from HTSPS_with_prepro import HTSPS_with_prepro
from entorno import Circulo

import numpy as np
import pandas as pd

init = False

if init:
    dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])
else:
    dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'NodeCount', 'ObjVal'])

for nP in [5, 10, 20, 30, 50, 80]:
    for instance in range(10):
        
        print('\n\nResolviendo la instancia ' + str(instance) + ' con un numero ' + str(nP) + ' de entornos.\n\n')
        
        segments = np.genfromtxt('./instancias/segmentos'+ str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        barriers = []
        for lista in segments:
            barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])
            
        bolas = np.genfromtxt('./instancias/bolas' + str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        N = [Circulo(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

        resultados = HTSPS_with_prepro(barriers, N, timeLimit = 3600)
        
        serie = pd.Series([instance] + resultados, index = dataframe.columns)
        
        dataframe = dataframe.append(serie, ignore_index=True)
        dataframe.to_csv('./resultados/resultados_with_prepro_updated.csv')

        
