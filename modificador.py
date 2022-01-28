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
# from HTSPS_with_prepro import HTSPS_with_prepro
# from HTSPS_without_prepro import HTSPS_without_prepro

np.random.seed(10)



for nP in [5, 10, 20, 30, 50, 80]:
    for instance in range(10):
        
        segmentos_visitar = []

        
        bolas = np.genfromtxt('./instancias/bolas' + str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        for centro1, centro2, radio in bolas:
            angle = np.random.uniform(0, 360) 
            
            P1 = [centro1 + radio*np.cos(angle), centro2 + radio*np.sin(angle)]
            P2 = [centro1 + radio*np.cos(angle+ 180), centro2 + radio*np.sin(angle+180)]
            
            segmentos_visitar.append([P1[0], P1[1], P2[0], P2[1]])
            
        np.savetxt('instancias/segmentos_visitar' + str(nP) + '-' + str(instance) + '.csv', segmentos_visitar, delimiter = ",")
