import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product, permutations, chain
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
import neighborhood as e
import copy
import estimacion_M as eM
from auxiliar_functions import *
import networkx as nx
from HTSPS import HTSPS
# from HTSPS_CPLEX import HTSPS_CPLEX

# semilla 1 hay punto de corte
np.random.seed(5)

# Lectura de datos
segmentos = np.genfromtxt('segments.csv', delimiter = ',')

# Numero total de segmentos que hay en los datos
nT = len(segmentos)

# Numero de segmentos que queremos visitar
nS = 5
segmentos_visitar = []

# Numero de barreras que queremos tener
nB = 5
barreras = []


todos_los_segmentos = segmentos[np.random.choice(nT, size = nS + nB, replace = False)]


for i in range(nS):
    row = todos_los_segmentos[i, :]
    poligono = e.Poligono(V = np.array([[row[0], row[1]], [row[2], row[3]]]), col = 'blue')
    segmentos_visitar.append(poligono)
    
for i in range(nS, nB+nS):
    row = todos_los_segmentos[i, :]
    poligono = e.Poligono(V = np.array([[row[0], row[1]], [row[2], row[3]], [row[0], row[1]]]), col = 'red')
    barreras.append(poligono)

# barreras = []

# barrera1 = Poligono(V = np.array([[7, 11], [9, 8]]))
# barreras.append(barrera1)
#
# barrera2 = Poligono(V = np.array([[12, 7], [11, 2]]))
# barreras.append(barrera2)
#
# barrera3 = Poligono(V = np.array([[9, 6], [8, 2]]))
# barreras.append(barrera3)
#
# barrera4 = Poligono(V = np.array([[11, 10], [16, 4]]))
# barreras.append(barrera4)
# datos = Data([], m = 10, r = 4, capacity = 100, modo = 2, tmax = 20, init = False, show = False)
# datos.generar_muestra()

print([barreras[i].V for i in range(len(barreras))])

print([segmentos_visitar[i].V for i in range(len(segmentos_visitar))])

parametros = Data(data = [], m = 100, init = True, tmax = 30)

HTSPS(segmentos_visitar, barreras, parametros)
# HTSPS_CPLEX(segmentos_visitar, barreras, parametros)