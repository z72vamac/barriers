import gurobipy as gp
import pdb
from tspn_b import HTSPS_new
from neighborhood import Circle

import numpy as np
import pandas as pd

init = False

if init:
    dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])
else:
    dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'NodeCount', 'ObjVal'])

for nP in range(5, 31, 5):
    for instance in range(5):
        
        print('Resolviendo la instancia ' + str(instance) + ' con un numero ' + str(nP) + ' de neighborhoods.')
        
        segments = np.genfromtxt('./instancias/segmentos'+ str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        barriers = []
        for lista in segments:
            barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])
            
        bolas = np.genfromtxt('./instancias/bolas' + str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        N = [Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

        resultados = HTSPS_new(barriers, N)
        
        serie = pd.Series([instance] + resultados, index = dataframe.columns)
        
        dataframe = dataframe.append(serie, ignore_index=True)
        dataframe.to_csv('./resultados/resultados_noinit.csv')

        
