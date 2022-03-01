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
import time


def HTSPS_without_prepro(barriers, N, log = False, timeLimit = 7200, init = False):
    # N = [N1, N2]
    
    first_time = time.time()
    
    VN = [-i for i in range(1, 1+len(N))]
    EN = list(product(VN, range(len(barriers)), range(2))) + list(product(range(len(barriers)), range(2), VN))
    
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
    
    # VT = [-2]
    # ET = list(product(VT, range(len(barriers)), range(2))) #+ list(product(range(len(barriers)), range(2), VT))
    
    V = VN + VB
    E = EN + EB
    
    if log:
        print("VN = " + str(VN))
        print("VB = " + str(VB))
        
        print("EN = " + str(EN))
        print("EB = " + str(EB))
    
    
    # print(barriers)
    
    
    
    def determinant(P, Q, R):
        a11 = Q[0]-P[0]
        a12 = R[0]-P[0]
        a21 = Q[1]-P[1]
        a22 = R[1]-P[1]
        
        return a11 * a22 - a12 * a21
        
        
    # epsilon(S / T, B, i) = 1 si (P_{S/T}, P_B^i)\in E_{S/T}
    epsilon_index = list(product(VN, range(len(barriers)), range(2)))
    p_index = EN
    if log:
        print("epsilon = " + str(epsilon_index))
    
    y_index = E
    
    if log:
        print("y = " + str(y_index))
    
    dist_index = E
    
    if log:
        print("dist = " + str(dist_index))
    
    dif_index = []
    
    for a in product(E, range(2)):
        tupla = tuple([element for element in a[0]] + [a[1]])
        dif_index.append(tupla)
    
    if log:    
        print("dif = " + str(dif_index))
    
    
    # delta(S / T, B, i, B') = 1 si (P_{S/T}, P_B^i) y B' do not intersect.
    delta_index = list(product(VN, range(len(barriers)), range(2), range(len(barriers)))) 
    
    if log:
        print("delta = " + str(delta_index))
    
    # gamma(S / T, B, i, B') = alpha(S / T, B, i, B', 0, order2)*alpha(S / T, B, i, B', 1, order2)
    gamma_index = list(product(VN, range(len(barriers)), range(2), range(len(barriers)), range(2)))
    
    if log:
        print("gamma = " + str(gamma_index))
    
    # beta(S/T, B, i, B', order) = 1 si sign(P_{S/T} P_B^i | B') is the same (order = 0); sign(B' | P_{S/T} P_B^i) is the same (order = 1);
    beta_index = list(product(VN, range(len(barriers)), range(2), range(len(barriers)), range(2)))
    
    if log:
        print("beta = " + str(beta_index))
    
    # alpha(S/T, B, i, B', order1, order2) = 1 si sign(P_{S/T} | B') > 0 (order1 = 0); sign(P_B^i | B') > 0 (order1) = 1;
    alpha_index = list(product(VN, range(len(barriers)), range(2), range(len(barriers)), range(2), range(2)))
    
    if log:
        print("alpha = " + str(alpha_index))
    
    # P_S and P_T: indices of the points in the neighborhoods
    P_index = list(product(VN, range(2)))
    
    if log:
        print("P_index = " + str(P_index))
    
    # socp variables:
    d_inside_index = VN
    dif_inside_index = list(product(VN, range(2)))
    
    if log:
        print("d_inside_index = " + str(d_inside_index))
        print("dif_inside_index = " + str(dif_inside_index))
    
    # z variables:
    # z_index = list(product(VN, VN))
    # print("z_index = " + str(z_index))
    
    
    # f variables:
    g_index = y_index
    
    if log:
        print("g_index = " + str(g_index))
    
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
    landa = MODEL.addVars(d_inside_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'landa')

    # z = MODEL.addVars(z_index, vtype = GRB.BINARY, name = 'z')
    g = MODEL.addVars(g_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'g')
    
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
            if determinant(barriers[b][c], barriers[d][0], barriers[d][1]) != 0:
                MODEL.addConstr(alpha[a, b, c, d, e, f] >= determinant(barriers[b][c], barriers[d][0], barriers[d][1])/abs(determinant(barriers[b][c], barriers[d][0], barriers[d][1])))
           
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
        
        MODEL.addConstr(y[a, b, c] + y[b, c, a] <= 2*epsilon[a, b, c])
    
    
    
    # NS and NT constraints
    for a, dim in P.keys():
        neighborhood = N[abs(a)-1]
        
        if type(neighborhood) is Circle:
            MODEL.addConstr(dif_inside[a, dim] >= P[a, dim] - N[abs(a)-1].center[dim])
            MODEL.addConstr(dif_inside[a, dim] >= N[abs(a)-1].center[dim] - P[a, dim])
        
            MODEL.addConstr(gp.quicksum(dif_inside[a, dim]*dif_inside[a, dim] for dim in range(2)) <= d_inside[a]*d_inside[a])
            MODEL.addConstr(d_inside[a] <= N[abs(a)-1].radii)
        
        if type(neighborhood) is Poligonal:
            MODEL.addConstrs(P[a, dim] == landa[a]*neighborhood.V[0][dim] + (1-landa[a])*neighborhood.V[1][dim] for dim in range(2))
    
    # MODEL.addConstr(y[-1, 2, 1] >= 0.5)
    # MODEL.addConstr(y[-2, 2, 1] >= 0.5)
    # MODEL.addConstr(P[-1, 0] == 30)
    # MODEL.addConstr(P[-1, 1] == 10)
    
    # MODEL.addConstr(epsilon[-4, 3, 1] >= 0.5)
    # MODEL.addConstr(epsilon[3, 1, -4] >= 0.5)
    
    # MODEL.addConstr(y[3, 1, -1] >= 0.5)
    
    # dist constraints
    for index in dif_index:
        if len(index) == 4:
            a, b, c, dim = index
            # print((a, b, c))
            if a < 0:
                MODEL.addConstr(dif[index] >= P[a, dim] - barriers[b][c][dim])
                MODEL.addConstr(dif[index] >= barriers[b][c][dim] - P[a, dim])
            
            elif c < 0:
                MODEL.addConstr(dif[index] >= P[c, dim] - barriers[a][b][dim])
                MODEL.addConstr(dif[index] >= barriers[a][b][dim] - P[c, dim])
            #
            # if a >= 0:
            #     MODEL.addConstrs(dif[index, dim] >= P[c, dim] - barriers[a][b][dim] for dim in range(2))
            #     MODEL.addConstrs(dif[index, dim] >= barriers[a][b][dim] - P[c, dim] for dim in range(2))
                     
        if len(index) == 5:
            a, b, c, d, e = index
            dist[a, b, c, d] = np.linalg.norm(np.array(barriers[a][b])-np.array(barriers[c][d]))
    
    MODEL.addConstrs(gp.quicksum(dif[a, b, c, dim]*dif[a, b, c, dim] for dim in range(2)) <= dist[a, b, c]*dist[a, b, c] for a, b, c in EN)
    
    
    L_out = 0
    U_out = 10000        
        
    
    def estima_L(neighborhood, punto):
        if type(neighborhood) is Circle:
            centro = neighborhood.center
            dist = np.linalg.norm(np.array(centro) - np.array(punto)) - neighborhood.radii
            
            # dist = 0
        
        if type(neighborhood) is Poligonal:
            dr = np.array(neighborhood.V[1]) - np.array(neighborhood.V[0])
            
            A, B = -dr[1], dr[0]
            
            C = -A*neighborhood.V[0][0] - B*neighborhood.V[0][1]
            
            dist = abs(A*punto[0] + B*punto[1] + C)/np.sqrt(A**2 + B**2)
            
        return dist
    
    def estima_U(neighborhood, punto):
        if type(neighborhood) is Circle:
            
            dist = np.linalg.norm(np.array(neighborhood.center) - np.array(punto))+ neighborhood.radii
            
            # dist = 10000
            
        if type(neighborhood) is Poligonal:
            
            dist = max([np.linalg.norm(np.array(neighborhood.V[i]) - np.array(punto)) for i in range(2)])
        
        return dist
    
    # p constraints
    for a, b, c in p.keys():
        if a < 0:
            neighborhood = N[abs(a) - 1]
            punto = barriers[b][c]
            L_out = estima_L(neighborhood, punto)
            U_out = estima_U(neighborhood, punto)
        
        if c < 0:
            neighborhood = N[abs(c) - 1]
            punto = barriers[a][b]
            L_out = estima_L(neighborhood, punto)
            U_out = estima_U(neighborhood, punto)
                        
        MODEL.addConstr(p[a, b, c] >= L_out*y[a, b, c])
        MODEL.addConstr(p[a, b, c] >= dist[a, b, c] - U_out*(1 - y[a, b, c]))
        
    
    # MODEL.addConstrs(z[v, v] == 0 for v in VN)
    
    for index in g_index:
        if len(index) == 3:
            vn, v, i = index
            MODEL.addConstr(g[vn, v, i] <= (len(N)-1)*y[vn, v, i])
        
        if len(index) == 4:
            v, i, w, j = index
            MODEL.addConstr(g[v, i, w, j] <= (len(N)-1)*y[v, i, w, j])
        
    
    
    MODEL.addConstrs(gp.quicksum(y[v, i, vn] for v, i in VB) >= 1 for vn in VN)
    
    for v in V:
        if v in VN:
            vn = v
            MODEL.addConstr(gp.quicksum(y[v, i, vn] for v, i in VB) == gp.quicksum(y[vn, v, i] for v, i in VB))
            if vn <= -2:
                MODEL.addConstr(gp.quicksum(g[vn, v, i] for v, i in VB) - gp.quicksum(g[v, i, vn] for v, i in VB) == 1)
        else:
            vb, i = v
            MODEL.addConstr(gp.quicksum(y[(w, j, vb, i)] for w, j in VB if (w, j, vb, i) in E) + gp.quicksum(y[(vn, vb, i)] for vn in VN if (vn, vb, i) in E) == gp.quicksum(y[(vb, i, w, j)] for w, j in VB if (vb, i, w, j) in E) + gp.quicksum(y[(vb, i, vn)] for vn in VN if (vb, i, vn) in E))
            MODEL.addConstr(gp.quicksum(g[(vb, i, w, j)] for w, j in VB if (vb, i, w, j) in E) + gp.quicksum(g[(vb, i, vn)] for vn in VN if (vb, i, vn) in E) - gp.quicksum(g[w, j, vb, i] for w, j in VB if (w, j, vb, i) in E) - gp.quicksum(g[(vn, vb, i)] for vn in VN if (vn, vb, i) in E) == 0)
    
    # MODEL.addConstrs(gp.quicksum(y[v, i, vn] for v, i in VB) == 1 for vn in VN)
    # MODEL.addConstrs(gp.quicksum(y[vn, v, i] for v, i in VB) == 1 for vn in VN)
    #
    # MODEL.addConstrs(gp.quicksum(g[v, i, vn] for v, i in VB) - gp.quicksum(g[vn, v, i] for v, i in VB) == 1 for vn in VN)
    # MODEL.addConstrs(gp.quicksum(g[index] for index in EB if index[2] == v and index[3] == i) + gp.quicksum(g[index] for index in EN if index[1] == v and index[2] == i) - gp.quicksum(g[index] for index in EB if index[0] == v and index[1] == i) - gp.quicksum(g[index] for index in EN if index[0] == v and index[1] == i) == 0 for v, i in VB)
    
    
    # MODEL.addConstrs(gp.quicksum(f[-1, w, k] for w in VN if w <= -2) == 1 for k in VN if k <= -2)
    # MODEL.addConstrs(gp.quicksum(f[v, w, w] ))
    
    # MODEL.addConstrs(gp.quicksum(z[w, v] for v in VN if w != v) == 1 for w in VN)
    
    
    # flow conservation constraints
    # for index in y_index:
    #     if len(index) == 3:
    # MODEL.addConstrs(gp.quicksum(y[tupla] for tupla in EN if tupla[0] == v) == 1 for v in VN)
    #
    # for v, i in VB:
    #     tuplas_salen = gp.quicksum([y[tupla] for tupla in EB if tupla[0] == v and tupla[1] == i]) + gp.quicksum([y[tupla] for tupla in EN if tupla[0] == v and tupla[1] == i])
    #     tuplas_entran = gp.quicksum([y[tupla] for tupla in EB if tupla[2] == v and tupla[3] == i]) + gp.quicksum([y[tupla] for tupla in EN if tupla[1] == v and tupla[2] == i])
    #
    #     MODEL.addConstr(tuplas_salen - tuplas_entran == 0)
    #
    # for v in VN:
    #     tuplas_salen = gp.quicksum([y[tupla] for tupla in EN if tupla[0] == v])
    #     tuplas_entran = gp.quicksum([y[tupla] for tupla in EN if tupla[2] == v])
    #
    #     MODEL.addConstr(tuplas_salen - tuplas_entran == 0)
    #
    # MODEL.addConstrs(gp.quicksum(y[tupla] for tupla in EN if tupla[2] == w) == 1 for w in VN)
    
    MODEL.update()
    
    objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(dist[index]*y[index] for index in EB)
    MODEL.setObjective(objective, GRB.MINIMIZE)
    
    second_time = time.time()
    
    time_elapsed = second_time - first_time
    MODEL.update()
    
    MODEL.Params.Threads = 6
    MODEL.Params.timeLimit = timeLimit - time_elapsed
    
    MODEL.write('prueba.lp')
    
    MODEL.optimize()
    
    # try:
    #     MODEL.write('solucion.sol')
    # except:
    #     MODEL.computeIIS()
    #     MODEL.write('prueba.ilp')
    
    
    # y_indices = []
    
    # for index in E:
    #     if y[index].X > 0.5:
    #         y_indices.append(index)
            
    # print(y_indices)
    
    # g_indices = []
    
    # for index in g_index:
    #     if g[index].X > 0.5:
    #         g_indices.append(g[index])
    
    # print(g_indices)
    
    resultados = []
    
    resultados.append(len(N))
    resultados.append(len(barriers))
    
    if MODEL.Status == 3:
        MODEL.computeIIS()
        MODEL.write('casa.ilp')
        if init:
            try:
                resultados.append(resultado_h[0])
                resultados.append(resultado_h[1])
            except:
                resultados.append(np.nan)
                resultados.append(np.nan)
        
        resultados = resultados + [np.nan, np.nan, np.nan, np.nan]
        
        return resultados
    
    if MODEL.SolCount == 0:
        if init:
            try:
                resultados.append(resultado_h[0])
                resultados.append(resultado_h[1])
            except:
                resultados.append(np.nan)
                resultados.append(np.nan)
    
        resultados = resultados + [np.nan, np.nan, np.nan, np.nan]

        return resultados
    
    
    if init:
        try:
            resultados.append(resultado_h[0])
            resultados.append(resultado_h[1])
        except:
            resultados.append(np.nan)
            resultados.append(np.nan)
            print('El heuristico no ha encontrado initial sol')
  
    resultados.append(MODEL.getAttr('MIPGap'))
    resultados.append(MODEL.Runtime + time_elapsed)
    resultados.append(MODEL.getAttr('NodeCount'))
    resultados.append(MODEL.ObjVal)
    
    if log:
        fig, ax = plt.subplots()
        
        for b in barriers:
            ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c = 'red')
        
        for n in N:
            ax.add_artist(n.figura)
        
        P_vals = MODEL.getAttr('x', P)
        # print(P_vals)
        
        points = []
        for keys, vals in P_vals.items():
            points.append(vals)
            
        points = np.array(points).reshape((len(N), 2))
        
        for i in points:
            ax.scatter(i[0], i[1], s = 10, c = 'black')
        
        # print(points)
        
        segmentos = []
        
        for indice in y_indices:
            if len(indice) == 3:
                a, b, c = indice
                if a < 0:
                    segmentos.append([points[abs(a)-1][0], barriers[b][c][0], points[abs(a)-1][1], barriers[b][c][1]])
                
                if c < 0:
                    segmentos.append([barriers[a][b][0], points[abs(c)-1][0], barriers[a][b][1], points[abs(c)-1][1]])
                
            else:
                a, b, c, d = indice
                segmentos.append([barriers[a][b][0], barriers[c][d][0], barriers[a][b][1], barriers[c][d][1]])
                
            
        # print(segmentos)
        for segmento in segmentos:
            ax.arrow(segmento[0], segmento[2], segmento[1]-segmento[0], segmento[3]-segmento[2], width = 0.1, head_width = 1, length_includes_head = True, color = 'black')
        
        
        
        plt.axis([0, 100, 0, 100])
        
        ax.set_aspect('equal')
        plt.show()
        
    return resultados

# barrier1 = [[20, 80], [40, 30]]
# barrier2 = [[70, 100], [40, 70]]
# barrier3 = [[100, 60], [60, 70]]
# barrier4 = [[60, 50], [90, 10]]
# barrier5 = [[10, 70], [20, 50]]
# # barrier6 = [[30, 70], [70, 20]]
#
# barriers = [barrier1, barrier2, barrier3, barrier4, barrier5]
# # barriers.append(barrier6)
#
# N1 = Circle(center = [20, 10], radii = 10)
# N2 = Circle(center = [90, 90], radii = 5)
# N3 = Circle(center = [40, 90], radii = 10)
# N4 = Circle(center = [85, 40], radii = 13)
#
# N = [N1, N2, N3, N4]
#
# HTSPS_new(barriers, N)
# MODEL.addConstr(gp.quicksum( y[v, w] for v, w in edges_tuple.select(0, '*')) - gp.quicksum( y[w, v] for w, v in edges_tuple.select('*', 0)) ==  1)
# MODEL.write('prueba.lp')
