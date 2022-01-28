import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product, permutations, chain
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
from entorno import Circulo
import copy
import estimacion_M as eM
from auxiliar_functions import *
import networkx as nx
# from HTSPS_with_prepro2 import HTSPS_with_prepro2
# from HTSPS_with_prepro import HTSPS_with_prepro
# from HTSPS_with_prepro3 import HTSPS_with_prepro3
# from HTSPS_with_prepro4 import HTSPS_with_prepro4
# from HTSPS_new_ven import HTSPS
from HTSPS_new import HTSPS
from HSPPS_new import HSPPS
from HTSPS_new_visible import HTSPS_ven

# from HTSPS_without_prepro import HTSPS_without_prepro


# segments = np.genfromtxt('./instancias/segmentos5-8.csv', delimiter = ',')
#
# barriers = []
# for lista in segments:
#     barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])
#
# bolas = np.genfromtxt('./instancias/bolas5-8.csv', delimiter = ',')
# N = [Circulo(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]
# # segmentos_visitar = np.genfromtxt('./instancias/segmentos_visitar50-8.csv', delimiter = ',')
# # N = [Poligonal(V = [np.array([lista[0], lista[1]]), np.array([lista[2], lista[3]])]) for lista in segmentos_visitar] # 105.164
# # resultados = HTSPS(barriers, N, prepro = False, log = True, timeLimit = 30)
# resultados = HTSPS_ven(barriers, N, prepro = False, log = True, timeLimit = 60)


# print(N)


barrier1 = [[20, 80], [40, 30]]
barrier2 = [[70, 95], [40, 70]]
barrier3 = [[95, 60], [60, 70]]
barrier4 = [[60, 50], [90, 10]]
barrier5 = [[10, 70], [20, 50]]
# barrier6 = [[30, 70], [70, 20]]

barriers = [barrier1, barrier2, barrier3, barrier4, barrier5]
# barriers = [barrier1, barrier3, barrier4, barrier5]
# barriers = [barrier1, barrier4]

N1 = Circulo(center=[20, 10], radii=10)
N2 = Circulo(center=[90, 90], radii=5)
N3 = Circulo(center=[35, 85], radii=9)
N4 = Circulo(center=[85, 40], radii=11)

N = [N1, N2, N3, N4]
# N = [N1, N4]

resultados = HTSPS(barriers, N, prepro=False, log=1, timeLimit=10)

# resultados = HTSPS_with_prepro(barriers, N, log = 1, timeLimit = 600)

print(resultados)

# HTSPS_with_prepro(barriers, N)
