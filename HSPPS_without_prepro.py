# Con la formulacion nueva
# Shortest path por neighborhoods con barriers

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product, permutations, chain
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
from neighborhood import Circle
import copy
import estimacion_M as eM
from auxiliar_functions import *
import networkx as nx

barrier1 = [[20, 80], [40, 30]]
barrier2 = [[70, 100], [40, 70]]
barrier3 = [[100, 60], [60, 70]]
barrier4 = [[60, 50], [90, 10]]
barrier5 = [[10, 70], [20, 50]]
barrier6 = [[30, 70], [70, 20]]

barriers = [barrier1, barrier2, barrier3, barrier4, barrier5]
barriers.append(barrier6)

NS = Circle(center = [20, 10], radii = 10)
NT = Circle(center = [90, 90], radii = 5)

N = [NS, NT]

VS = [-1]
ES = list(product(VS, range(len(barriers)), range(2))) #+ list(product(range(len(barriers)), range(2), VS))

def nocortan(barrier1, barrier2):
    det1 = (barrier1[0][0] - barrier2[0][0])*(barrier1[1][1] - barrier2[0][1]) - (barrier1[1][0] - barrier2[0][0])*(barrier1[0][1] - barrier2[0][1])
    det2 = (barrier1[0][0] - barrier2[1][0])*(barrier1[1][1] - barrier2[1][1]) - (barrier1[1][0] - barrier2[1][0])*(barrier1[0][1] - barrier2[1][1])

    det3 = (barrier2[0][0] - barrier1[0][0])*(barrier2[1][1] - barrier1[0][1]) - (barrier2[1][0] - barrier1[0][0])*(barrier2[0][1] - barrier1[0][1])
    det4 = (barrier2[0][0] - barrier1[1][0])*(barrier2[1][1] - barrier1[1][1]) - (barrier2[1][0] - barrier1[1][0])*(barrier2[0][1] - barrier1[1][1])

    return ((det1*det2) >= 0) or ((det3*det4) >= 0)

VB = list(product(range(len(barriers)), range(2)))

EB = []
for v, i in VB:
    for w, j in VB:
        if v != w:
            barrier = [barriers[v][i], barriers[w][j]]
            
            if all([nocortan(barrieri, barrier) for barrieri in barriers]):
                EB.append((v, i, w, j))

# EB = list(product(range(len(barriers)), range(2), range(len(barriers)), range(2)))

VT = [-2]
ET = list(product(VT, range(len(barriers)), range(2))) #+ list(product(range(len(barriers)), range(2), VT))

V = VS + VB + VT
E = ES + EB + ET

print("V = " + str(V))
print("E = " + str(E))

# print(barriers)



def determinant(P, Q, R):
    a11 = Q[0]-P[0]
    a12 = R[0]-P[0]
    a21 = Q[1]-P[1]
    a22 = R[1]-P[1]
    
    return a11 * a22 - a12 * a21
    
    
# epsilon(S / T, B, i) = 1 si (P_{S/T}, P_B^i)\in E_{S/T}
epsilon_index = ES + ET
p_index = ES + ET
print("epsilon = " + str(epsilon_index))

y_index = E
print("y = " + str(y_index))

dist_index = E
print("dist = " + str(dist_index))

dif_index = []

for a in product(E, range(2)):
    tupla = tuple([element for element in a[0]] + [a[1]])
    dif_index.append(tupla)
    
print("dif = " + str(dif_index))


# delta(S / T, B, i, B') = 1 si (P_{S/T}, P_B^i) y B' do not intersect.
delta_index = list(product(VS + VT, range(len(barriers)), range(2), range(len(barriers))))
print("delta = " + str(delta_index))

# gamma(S / T, B, i, B') = alpha(S / T, B, i, B', 0, order2)*alpha(S / T, B, i, B', 1, order2)
gamma_index = list(product(VS + VT, range(len(barriers)), range(2), range(len(barriers)), range(2)))
print("gamma = " + str(gamma_index))

# beta(S/T, B, i, B', order) = 1 si sign(P_{S/T} P_B^i | B') is the same (order = 0); sign(B' | P_{S/T} P_B^i) is the same (order = 1);
beta_index = list(product(VS + VT, range(len(barriers)), range(2), range(len(barriers)), range(2)))
print("beta = " + str(beta_index))

# alpha(S/T, B, i, B', order1, order2) = 1 si sign(P_{S/T} | B') > 0 (order1 = 0); sign(P_B^i | B') > 0 (order1) = 1;
alpha_index = list(product(VS + VT, range(len(barriers)), range(2), range(len(barriers)), range(2), range(2)))
print("alpha = " + str(alpha_index))

# P_S and P_T: indices of the points in the neighborhoods
P_index = list(product(VS + VT, range(2)))
print("P_index = " + str(P_index))

# socp variables:
d_inside_index = VS + VT
dif_inside_index = list(product(VS + VT, range(2)))

print("d_inside_index = " + str(d_inside_index))
print("dif_inside_index = " + str(dif_inside_index))

MODEL = gp.Model('HSPPS_Model')

epsilon = MODEL.addVars(epsilon_index, vtype = GRB.BINARY, name = 'epsilon')
p = MODEL.addVars(p_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'p')
y = MODEL.addVars(y_index, vtype = GRB.BINARY, name = 'y')
dist = MODEL.addVars(dist_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dist')
dif = MODEL.addVars(dif_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dif')
delta = MODEL.addVars(delta_index, vtype = GRB.BINARY, name = 'delta')
gamma = MODEL.addVars(gamma_index, vtype = GRB.BINARY, name = 'gamma')
beta = MODEL.addVars(beta_index, vtype = GRB.BINARY, name = 'beta')
alpha = MODEL.addVars(alpha_index, vtype = GRB.BINARY, name = 'alpha')
P = MODEL.addVars(P_index, vtype = GRB.CONTINUOUS, name = 'P')
d_inside = MODEL.addVars(d_inside_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'd_inside')
dif_inside = MODEL.addVars(dif_inside_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dif_inside')

MODEL.update()

# Alpha constraint
L = -10000
U = 10000

# alpha-C
for a, b, c, d, e, f in alpha.keys():
    if e == 0 and f == 0:
        MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*L <= determinant([P[a, 0], P[a, 1]], barriers[d][0], barriers[d][1]))
        MODEL.addConstr(determinant([P[a, 0], P[a, 1]], barriers[d][0], barriers[d][1]) <= U*alpha[a, b, c, d, e, f])
    
    # elif e == 1 and f == 0:
    #     MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*L <= determinant(barriers[b][c], barriers[d][0], barriers[d][1]))
    #     MODEL.addConstr(determinant(barriers[b][c], barriers[d][0], barriers[d][1]) <= U*alpha[a, b, c, d, e, f])

    elif e == 1 and f == 0:
        MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*abs(determinant(barriers[b][c], barriers[d][0], barriers[d][1]))*(-1) <= determinant(barriers[b][c], barriers[d][0], barriers[d][1]))
        MODEL.addConstr(determinant(barriers[b][c], barriers[d][0], barriers[d][1]) <= abs(determinant(barriers[b][c], barriers[d][0], barriers[d][1]))*alpha[a, b, c, d, e, f])
       
    else:
        MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*L <= determinant(barriers[d][e], [P[a, 0], P[a, 1]], barriers[b][c]))
        MODEL.addConstr(determinant(barriers[d][e], [P[a, 0], P[a, 1]], barriers[b][c]) <= U*alpha[a, b, c, d, e, f])
    
    # if e == 1 and f == 1:
    #     MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*L <= determinant(barriers[d][0], [P[a, 0], P[a, 1]], barriers[b][c]))
    #     MODEL.addConstr(determinant(barriers[b][c], barriers[d][0], barriers[d][1]) <= M*alpha[a, b, c, d, e, f])


# beta-C
for a, b, c, d, f in beta.keys():
    MODEL.addConstr(beta[a, b, c, d, f] == 2*gamma[a, b, c, d, f] - alpha[a, b, c, d, 0, f]-alpha[a, b, c, d, 1, f]+1)


# gamma-C
for a, b, c, d, f in gamma.keys():
    MODEL.addConstr(gamma[a, b, c, d, f] <= alpha[a, b, c, d, 0, f])
    MODEL.addConstr(gamma[a, b, c, d, f] <= alpha[a, b, c, d, 1, f])
    MODEL.addConstr(gamma[a, b, c, d, f] >= alpha[a, b, c, d, 0, f] + alpha[a, b, c, d, 1, f] - 1)
    
    

# delta-C
for a, b, c, d in delta.keys():
    MODEL.addConstr((beta[a, b, c, d, 0] + beta[a, b, c, d, 1] <= 2*delta[a, b, c, d]))
    MODEL.addConstr(delta[a, b, c, d] <= 2*(beta[a, b, c, d, 0] + beta[a, b, c, d, 1]))
    
# epsilon-C and y-C
for a, b, c in epsilon.keys():
    MODEL.addConstr((delta.sum(a, b, c, '*')-len(barriers)) + 1 <= epsilon[a, b, c])
    MODEL.addConstr(len(barriers)*epsilon[a, b, c] <= delta.sum(a, b, c, '*'))
    
    MODEL.addConstr(y[a, b, c] <= epsilon[a, b, c])



# NS and NT constraints
for a, dim in P.keys():
    MODEL.addConstr(dif_inside[a, dim] >= P[a, dim] - N[abs(a)-1].center[dim])
    MODEL.addConstr(dif_inside[a, dim] >= N[abs(a)-1].center[dim] - P[a, dim])
    
for a in VS+VT:
    MODEL.addConstr(gp.quicksum(dif_inside[a, dim]*dif_inside[a, dim] for dim in range(2)) <= d_inside[a]*d_inside[a])
    MODEL.addConstr(d_inside[a] <= N[abs(a)-1].radii)

# MODEL.addConstr(y[-1, 2, 1] >= 0.5)
# MODEL.addConstr(y[-2, 2, 1] >= 0.5)
# MODEL.addConstr(P[-1, 0] == 30)
# MODEL.addConstr(P[-1, 1] == 10)


# dist constraints
for index in dif_index:
    if len(index) == 4:
        a, b, c, dim = index
        print((a, b, c))
        # if a < 0:
        MODEL.addConstr(dif[index] >= P[a, dim] - barriers[b][c][dim])
        MODEL.addConstr(dif[index] >= barriers[b][c][dim] - P[a, dim])
        #
        # if a >= 0:
        #     MODEL.addConstrs(dif[index, dim] >= P[c, dim] - barriers[a][b][dim] for dim in range(2))
        #     MODEL.addConstrs(dif[index, dim] >= barriers[a][b][dim] - P[c, dim] for dim in range(2))
                 
    if len(index) == 5:
        a, b, c, d, e = index
        dist[a, b, c, d] = np.linalg.norm(np.array(barriers[a][b])-np.array(barriers[c][d]))

MODEL.addConstrs(gp.quicksum(dif[a, b, c, dim]*dif[a, b, c, dim] for dim in range(2)) <= dist[a, b, c]*dist[a, b, c] for a, b, c in ES + ET)


L_out = 0
U_out = 10000

# p constraints
for a, b, c in p.keys():
    MODEL.addConstr(p[a, b, c] >= L_out*y[a, b, c])
    MODEL.addConstr(p[a, b, c] >= dist[a, b, c] - U_out*(1 - y[a, b, c]))
    


# flow conservation constraints
# for index in y_index:
#     if len(index) == 3:
MODEL.addConstr(gp.quicksum(y[tupla] for tupla in E if len(tupla) == 3 and tupla[0] == -1) == 1)

for v, i in VB:
    tuplas_salen = gp.quicksum([y[tupla] for tupla in E if ((len(tupla) == 4 and tupla[0] == v and tupla[1] == i) or (len(tupla) == 3 and tupla[0] == -2 and tupla[1] == v and tupla[2] == i))])
    tuplas_entran = gp.quicksum([y[tupla] for tupla in E if ((len(tupla) == 4 and tupla[2] == v and tupla[3] == i) or (len(tupla) == 3 and tupla[0] == -1 and tupla[1] == v and tupla[2] == i))])

    MODEL.addConstr(tuplas_salen - tuplas_entran == 0)
    
MODEL.addConstr(- gp.quicksum(y[tupla] for tupla in E if len(tupla) == 3 and tupla[0] == -2) == -1)

MODEL.update()

objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(dist[index]*y[index] for index in EB)
MODEL.setObjective(objective, GRB.MINIMIZE)

MODEL.update()

MODEL.write('prueba.lp')

MODEL.optimize()

try:
    MODEL.write('solucion.sol')
except:
    MODEL.computeIIS()
    MODEL.write('prueba.ilp')


y_indices = []

for index in E:
    if y[index].X > 0.5:
        y_indices.append(index)
        
print(y_indices)


# MODEL.addConstr(gp.quicksum( y[v, w] for v, w in edges_tuple.select(0, '*')) - gp.quicksum( y[w, v] for w, v in edges_tuple.select('*', 0)) ==  1)
# MODEL.write('prueba.lp')
