"""Tenemos un conjunto E de neighborhoods ugt un conjunto de poligonales P de las que queremos recorrer un porcentaje alfa p . Buscamos
   un tour de mínima distancia que alterne poligonal-neighborhood ugt que visite todas las poligonales"""


# Incluimos primero los paquetes
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product, permutations, chain
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
from neighborhood import *
import copy
import estimacion_M as eM
from auxiliar_functions import path2matrix

# Definicion de los datos
""" P: conjunto de poligonales a agrupar
    E: conjunto de neighborhoods
    T: sucesion de etapas
    C(e): centro del neighborhood e
    R(e): radio del neighborhood e
    p: indice de las poligonales
    e: indice de los neighborhoods
    t: indice de las etapas
    n: dimension del problema
"""

np.random.seed(129)

segmento1 = Poligonal(np.array([[5, 6], [4, 9]]), alpha = 1)
segmento2 = Poligonal(np.array([[18, 7], [17, 9]]), alpha = 1)

barrera = Poligonal(np.array([[10, 4], [13, 11]]), alpha = 1)
# barrera = Poligonal(np.array([[8, 11], [13, 11]]), alpha = 1)

# Creacion del modelo
MODEL = gp.Model('barriers')


# Variable binaria ugt = 1 si se entra por el grafo g en la etapa t
yL1 = MODEL.addVar(vtype = GRB.BINARY, name = 'yL1')
dL1 = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dL1')
difL1 = MODEL.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'difL1')

yL1p = MODEL.addVar(vtype = GRB.BINARY, name = 'yL1p')
dL1p = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dL1p')
difL1p = MODEL.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'difL1p')

y1R = MODEL.addVar(vtype = GRB.BINARY, name = 'y1R')
d1R = MODEL.addVar(vtype = GRB.CONTINUOUS, lb =0.0, name = 'd1R')
dif1R = MODEL.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dif1R')

y1pR = MODEL.addVar(vtype = GRB.BINARY, name = 'y1pR')
d1pR = MODEL.addVar(vtype = GRB.CONTINUOUS, lb =0.0, name = 'd1pR')
dif1pR = MODEL.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dif1Rp')

yLR = MODEL.addVar(vtype = GRB.BINARY, name = 'yLR')
dLR = MODEL.addVar(vtype = GRB.CONTINUOUS, name = 'dLR')
difLR = MODEL.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'difLR')

gammaL = MODEL.addVar(vtype = GRB.BINARY, name = 'gammaL')
deltaL = MODEL.addVar(vtype = GRB.BINARY, name = 'deltaL')


xR = MODEL.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 30.0, name = 'xR')
xL = MODEL.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 30.0, name = 'xL')

# alpha1 = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0,  ub = 1000, name = 'alpha1')
# alpha1p = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1000, name = 'alpha1p')

alpha1 = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'alpha1')
alpha1p = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'alpha1p')
beta1 = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'beta1')


mu1 = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mu1')
mu2 = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mu2')

MODEL.update()

MODEL.addConstr(yL1 + yL1p + yLR == 1)
MODEL.addConstr(yL1 - y1R == 0)
MODEL.addConstr(yL1p - y1pR == 0)
MODEL.addConstr(-y1R - y1pR - yLR == -1)
MODEL.addConstr(yLR <= gammaL)
MODEL.addConstr(deltaL <= gammaL)
# MODEL.addConstr(dLR >= (1 - gammaL) * 10000)

puntomedio = np.mean(barrera.V, axis = 0)

MODEL.addConstrs(xR[dim] <= xL[dim] + alpha1*(barrera.V[0][dim] - xL[dim]) - beta1*(puntomedio[dim] - xL[dim]) + 10000*(1 - gammaL + deltaL) for dim in range(2))
MODEL.addConstrs(xR[dim] >= xL[dim] + alpha1*(barrera.V[0][dim] - xL[dim]) + beta1*(puntomedio[dim] - xL[dim]) - 10000*(1 - gammaL + deltaL) for dim in range(2))

MODEL.addConstrs(xR[dim] <= xL[dim] + alpha1*(barrera.V[0][dim] - xL[dim]) + alpha1p*(barrera.V[1][dim] - xL[dim]) + 10000*gammaL for dim in range(2))
MODEL.addConstrs(xR[dim] >= xL[dim] + alpha1*(barrera.V[0][dim] - xL[dim]) + alpha1p*(barrera.V[1][dim] - xL[dim]) - 10000*gammaL for dim in range(2))

MODEL.addConstrs(xR[dim] <= xL[dim] + alpha1p*(barrera.V[1][dim] - xL[dim]) - beta1*(puntomedio[dim] - xL[dim]) + 10000*(2 - gammaL - deltaL) for dim in range(2))
MODEL.addConstrs(xR[dim] >= xL[dim] + alpha1p*(barrera.V[1][dim] - xL[dim]) - beta1*(puntomedio[dim] - xL[dim]) - 10000*(2 - gammaL - deltaL) for dim in range(2))

MODEL.addConstrs(xL[dim] == mu1*segmento1.V[0][dim] + (1-mu1)*segmento1.V[1][dim] for dim in range(2))
MODEL.addConstrs(xR[dim] == mu2*segmento2.V[0][dim] + (1-mu2)*segmento2.V[1][dim] for dim in range(2))

MODEL.addConstrs(difL1[dim] >=  xL[dim] - barrera.V[0][dim] for dim in range(2))
MODEL.addConstrs(difL1[dim] >= -xL[dim] + barrera.V[0][dim] for dim in range(2))
MODEL.addConstr(gp.quicksum(difL1[dim]*difL1[dim] for dim in range(2)) <= dL1*dL1)

MODEL.addConstrs(difL1p[dim] >=  xL[dim] - barrera.V[1][dim] for dim in range(2))
MODEL.addConstrs(difL1p[dim] >= -xL[dim] + barrera.V[1][dim] for dim in range(2))
MODEL.addConstr(gp.quicksum(difL1p[dim]*difL1p[dim] for dim in range(2)) <= dL1p*dL1p)

MODEL.addConstrs(dif1R[dim] >=  xR[dim] - barrera.V[0][dim] for dim in range(2))
MODEL.addConstrs(dif1R[dim] >= -xR[dim] + barrera.V[0][dim] for dim in range(2))
MODEL.addConstr(gp.quicksum(dif1R[dim]*dif1R[dim] for dim in range(2)) <= d1R*d1R)

MODEL.addConstrs(dif1pR[dim] >=  xR[dim] - barrera.V[1][dim] for dim in range(2))
MODEL.addConstrs(dif1pR[dim] >= -xR[dim] + barrera.V[1][dim] for dim in range(2))
MODEL.addConstr(gp.quicksum(dif1pR[dim]*dif1pR[dim] for dim in range(2)) <= d1pR*d1pR)

MODEL.addConstrs(difLR[dim] >=  xL[dim] - xR[dim] for dim in range(2))
MODEL.addConstrs(difLR[dim] >= -xL[dim] + xR[dim] for dim in range(2))
MODEL.addConstr(gp.quicksum(difLR[dim]*difLR[dim] for dim in range(2)) <= dLR*dLR)


MODEL.setObjective(dL1*yL1 + dL1p*yL1p + d1R*y1R + d1pR*y1pR + dLR*yLR, GRB.MINIMIZE)

MODEL.update()

MODEL.Params.NonConvex = 2

# Optimizamos
MODEL.optimize()

if MODEL.Status == 3:

    MODEL.computeIIS()
    MODEL.write('infactible.ilp')

MODEL.write('solucion_nivel1.sol')
