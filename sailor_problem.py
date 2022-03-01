import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product, permutations, chain
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
from neighborhood import Circle, Poligonal
import copy
import estimacion_M as eM
from auxiliar_functions import *
import networkx as nx
import time

beta = 0
num_angles = 4
num_puntos = 20

L = 5
s = [0, 0]
t = [30, 0]

num_angles += 1
angle = np.linspace(0, 2*np.pi, num_angles)

u_j = np.zeros((num_angles, 2))

u_j[:, 0] = np.cos(angle)
u_j[:, 1] = np.sin(angle)

print(u_j)

# c_j = (1- (1/np.pi)*np.abs(np.abs(angle-beta)-np.pi))
c_j = 1 - np.cos(abs(angle-beta))


print(c_j)




MODEL = gp.Model('Sailor_Problem')

P = MODEL.addVars(num_puntos, 2, vtype = GRB.CONTINUOUS, name = 'P')
alpha = MODEL.addVars(num_puntos, num_angles, vtype = GRB.BINARY, name = 'alpha')
r = MODEL.addVars(num_puntos, vtype = GRB.CONTINUOUS, lb = 0, ub = L, name = 'r')
p = MODEL.addVars(num_puntos, vtype = GRB.CONTINUOUS, lb = 0, name = 'p')
aux1 = MODEL.addVars(num_puntos, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'aux1')
aux2 = MODEL.addVars(num_puntos, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'aux2')



MODEL.update()

# MODEL.addConstr(alpha[0, 0] == 1)

MODEL.addConstrs(P[0, dim] == s[dim] for dim in range(2))
MODEL.addConstrs(P[i+1, dim] == P[i, dim] + gp.quicksum(r[i]*alpha[i, j]*u_j[j, dim] for j in range(num_angles-1)) for i in range(num_puntos-1) for dim in range(2))
MODEL.addConstrs(P[num_puntos-1, dim] == t[dim] for dim in range(2))

MODEL.addConstrs(alpha.sum(i, '*') == 1 for i in range(num_puntos))

MODEL.addConstrs(aux1[i] >= gp.quicksum(j*(alpha[i+1, j] - alpha[i, j]) for j in range(num_angles-1)) for i in range(num_puntos-1))
MODEL.addConstrs(aux1[i] >= gp.quicksum(j*(alpha[i, j] - alpha[i+1, j]) for j in range(num_angles-1)) for i in range(num_puntos-1))

MODEL.addConstrs(aux2[i] >= beta - gp.quicksum(j*alpha[i, j] for j in range(num_angles-1)) for i in range(num_puntos-1))
MODEL.addConstrs(aux2[i] >= gp.quicksum(j*alpha[i, j] for j in range(num_angles-1)) - beta for i in range(num_puntos-1))


objective =  gp.quicksum(gp.quicksum(r[i]*c_j[j]*alpha[i, j] for j in range(num_angles-1)) for i in range(num_puntos-1)) + gp.quicksum(aux1[i] for i in range(num_puntos-1))

# objective = gp.quicksum(aux1[i] for i in range(num_puntos-1))
# 
#  +
            
MODEL.setObjective(objective, GRB.MINIMIZE)

MODEL.update()

MODEL.Params.NonConvex = 2

MODEL.optimize()

MODEL.write('solucion_sailor.sol')

fig, ax = plt.subplots()

for i in range(num_puntos):
    ax.scatter(P[i, 0].X, P[i, 1].X, s = 10, c = 'black')
    
plt.show()
# for i in range(num_puntos):
#     ax.plot([P[i, 0], b[1][0]], [b[0][1], b[1][1]], c = 'red')
# print(MODEL.getAttr('X', alpha))


