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
from neighborhood import Circle, Poligonal
import copy
import estimacion_M as eM
from auxiliar_functions import *
import networkx as nx
import time


def HTSPS_with_prepro3(barriers, N, log = False, timeLimit = 7200, init = False):
    # N = [N1, N2]
    
    first_time = time.time()
    def nocortan(barrier1, barrier2):
        det1 = (barrier1[0][0] - barrier2[0][0])*(barrier1[1][1] - barrier2[0][1]) - (barrier1[1][0] - barrier2[0][0])*(barrier1[0][1] - barrier2[0][1])
        det2 = (barrier1[0][0] - barrier2[1][0])*(barrier1[1][1] - barrier2[1][1]) - (barrier1[1][0] - barrier2[1][0])*(barrier1[0][1] - barrier2[1][1])
    
        det3 = (barrier2[0][0] - barrier1[0][0])*(barrier2[1][1] - barrier1[0][1]) - (barrier2[1][0] - barrier1[0][0])*(barrier2[0][1] - barrier1[0][1])
        det4 = (barrier2[0][0] - barrier1[1][0])*(barrier2[1][1] - barrier1[1][1]) - (barrier2[1][0] - barrier1[1][0])*(barrier2[0][1] - barrier1[1][1])
    
        return ((det1*det2) >= 0) or ((det3*det4) >= 0)
    
    # def ve(punto, segmento):
    #     barrera1 = [punto, segmento[0]]
    #     flag1 = all([nocortan(barrera1, barrera) for barrera in barriers])
    #     barrera2 = [punto, segmento[1]]
    #     flag2 = all([nocortan(barrera2, barrera) for barrera in barriers])
    #
    #     return flag1 | flag2
    
    def no_ve(neighborhood, punto):
        if type(neighborhood) is Circle:
            centro = neighborhood.center
            
            dr = np.array(centro) - np.array(punto)
            dr_u = dr / np.linalg.norm(dr)
            
            nr_u = np.array([-dr_u[1], dr_u[0]])
            
            lambdas = np.linspace(-neighborhood.radii, neighborhood.radii, 20)
            
            # flag1 representa si corta alguna o no
            # flag2 representa si todos los puntos tienen algun corte
            
            flag2 = True

            for landa in lambdas:
                
                for barrera in barriers:
                    flag1 = False
                    if not(nocortan([punto, np.array(centro)+landa*nr_u], [barrera[0], barrera[1]])):
                        flag1 = True
                        break
                
                flag2 = flag2 & flag1
                
                if not(flag2):
                    return flag2
                
            return flag2
                           
            # punto1 = np.array(centro) + neighborhood.radii*nr_u
            # segmento1 = [punto, punto1]
            # flag1 = any([not(nocortan(segmento1, barrera)) for barrera in barriers])
            #
            # punto2 = np.array(centro) - neighborhood.radii*nr_u
            # segmento2 = [punto, punto2]
            # flag2 = any([not(nocortan(segmento2, barrera)) for barrera in barriers])
            
        if type(neighborhood) is Poligonal:
            
            nr = np.array(neighborhood.V[1]) - np.array(neighborhood.V[0])
            
            nr_u = nr / np.linalg.norm(nr)
            
            lambdas = np.linspace(0, neighborhood.longitud, 20)

            flag2 = True

            for landa in lambdas:
                
                for barrera in barriers:
                    flag1 = False
                    if not(nocortan([punto, np.array(neighborhood.V[0]) + landa*nr_u], [barrera[0], barrera[1]])):
                        flag1 = True
                        break
                
                flag2 = flag2 & flag1
                
                if not(flag2):
                    return flag2
                
            return flag2        
            
        
    def ve(neighborhood, punto):
        return not(no_ve(neighborhood, punto))
        
        
    
    VN = [-i for i in range(1, 1+len(N))]
    EN = [] #list(product(VN, range(len(barriers)), range(2))) + list(product(range(len(barriers)), range(2), VN))
    
    
    for a in VN:
        for b in range(len(barriers)):
            for c in range(2):
                # if ve(N[abs(a)-1], barriers[b][c]): #or ve(N[abs(a)-1], barriers[b][1]):
                    EN.append((a, b, c))
                    EN.append((b, c, a))     

    # for a, b, c in EN:
    
    # print(EN)
        
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

    # print("VN = " + str(VN))
    # print("VB = " + str(VB))
    #
    # print("EN = " + str(EN))
    # print("EB = " + str(EB))    
    
    # print(barriers)
    
    
    
    def determinant(P, Q, R):
        a11 = Q[0]-P[0]
        a12 = R[0]-P[0]
        a21 = Q[1]-P[1]
        a22 = R[1]-P[1]
        
        return a11 * a22 - a12 * a21
        
        
    # epsilon(S / T, B, i) = 1 si (P_{S/T}, P_B^i)\in E_{S/T}
    epsilon_index = []
    
    for a, b, c in EN:
        if a < 0:
            epsilon_index.append((a, b, c))
    
    # print(epsilon_index)
    
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
    delta_index = [] #list(product(VN, range(len(barriers)), range(2), range(len(barriers)))) 
    
    
    for a, b, c in epsilon_index:
        for d in range(len(barriers)):
            delta_index.append((a, b, c, d))

    
    if log:
        print("delta = " + str(delta_index))
    
    # gamma(S / T, B, i, B') = alpha(S / T, B, i, B', 0, order2)*alpha(S / T, B, i, B', 1, order2)
    alpha_index = []
    
    for a, b, c, d in delta_index:
        for e in range(2):
            alpha_index.append((a, b, c, d, e))
    
    # list(product(VN, range(len(barriers)), range(2), range(len(barriers)), range(2)))
    
    beta_index = delta_index
    
    gamma_index = delta_index
    # for a, b, c, d in delta_index:
    #     for e in range(2):
    #         beta_index.append((a, b, c, d, e))

    
    # beta(S/T, B, i, B', order) = 1 si sign(P_{S/T} P_B^i | B') is the same (order = 0); sign(B' | P_{S/T} P_B^i) is the same (order = 1);
    
    if log:
        print("beta = " + str(beta_index))
    
    # alpha(S/T, B, i, B', order1, order2) = 1 si sign(P_{S/T} | B') > 0 (order1 = 0); sign(P_B^i | B') > 0 (order1) = 1;
    # alpha_index = [] # list(product(VN, range(len(barriers)), range(2), range(len(barriers)), range(2), range(2)))
    #
    # for a, b, c, d, e in beta_index:
    #     for f in range(2):
    #         alpha_index.append((a, b, c, d, e, f))
                    
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
    dif = MODEL.addVars(dif_index, vtype = GRB.CONTINUOUS, name = 'dif')
    delta = MODEL.addVars(delta_index, vtype = GRB.BINARY, name = 'delta')
    gamma = MODEL.addVars(gamma_index, vtype = GRB.BINARY, name = 'gamma')
    aux = MODEL.addVars(gamma_index, vtype = GRB.BINARY, name = 'aux')
    beta = MODEL.addVars(beta_index, vtype = GRB.BINARY, name = 'beta')
    alpha = MODEL.addVars(alpha_index, vtype = GRB.BINARY, name = 'alpha')
    P = MODEL.addVars(P_index, vtype = GRB.CONTINUOUS, name = 'P')
    d_inside = MODEL.addVars(d_inside_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'd_inside')
    dif_inside = MODEL.addVars(dif_inside_index, vtype = GRB.CONTINUOUS, name = 'dif_inside')
    landa = MODEL.addVars(d_inside_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'landa')
    # z = MODEL.addVars(z_index, vtype = GRB.BINARY, name = 'z')
    g = MODEL.addVars(g_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'g')
    
    MODEL.update()
    
    # Alpha constraint
    L = -1e6
    U = 1e6
    
    # param2[-1,0,0,3] = -4 / 3
    
    # alpha-C
    for a, b, c, d in beta.keys():
        x1 = P[a, 0]
        y1 = P[a, 1]
        
        x2 = barriers[b][c][0]
        y2 = barriers[b][c][1]
        
        x3 = barriers[d][0][0]
        y3 = barriers[d][0][1]
        
        x4 = barriers[d][1][0]
        y4 = barriers[d][1][1]
        
        # MODEL.addConstr(((y3-y1)*(x2-x1)-(x3-x1)*(y2-y1)) >= L*(1-beta[a, b, c, d])*((y3-y1)*(x2-x1)-(x3-x1)*(y2-y1)))
        # MODEL.addConstr(((y3-y1)*(x2-x1)-(x3-x1)*(y2-y1)) <= U*(2-gamma[a, b, c, d])*((y3-y1)*(x2-x1)-(x3-x1)*(y2-y1)))
        

        # MODEL.addConstr((y3-y1)*(x2-x1)-(x3-x1)*(y2-y1) >= L*(1 - alpha[a, b, c, d, 0]))
        MODEL.addConstr((y3-y1)*(x2-x1)-(x3-x1)*(y2-y1) <= U*alpha[a, b, c, d, 0])
        
        # MODEL.addConstr((x4-x3)*(y2-y1)-(y4-y3)*(x2-x1) >= L*(1 - alpha[a, b, c, d, 1]))
        MODEL.addConstr((x4-x3)*(y2-y1)-(y4-y3)*(x2-x1) <= U*alpha[a, b, c, d, 1])
        
        MODEL.addConstr(beta[a, b, c, d] == alpha[a, b, c, d, 0]*alpha[a, b, c, d, 1] + (1-alpha[a, b, c, d, 0])*(1-alpha[a, b, c, d, 1]))
        
        
        # MODEL.addConstr((y3-y1)*(x2-x1)-(x3-x1)*(y2-y1) - (x4-x3)*(y2-y1)-(y4-y3)*(x2-x1) <= U*(1-gamma[a, b, c, d]))
        MODEL.addConstr((y3-y1)*(x2-x1)-(x3-x1)*(y2-y1) - (x4-x3)*(y2-y1)-(y4-y3)*(x2-x1) >= L*(1 - gamma[a, b, c, d]))
        #
        MODEL.addConstr(2*delta[a, b, c, d] >= beta[a, b, c, d] + gamma[a, b, c, d] - 1)
        MODEL.addConstr(delta[a, b, c, d] <= beta[a, b, c, d])
        MODEL.addConstr(delta[a, b, c, d] <= gamma[a, b, c, d])
        
        
        

        
        
        
        
            
    
    # beta-C
    # for a, b, c, d, f in beta.keys():
    #     MODEL.addConstr(beta[a, b, c, d, f] == 2*gamma[a, b, c, d, f] - alpha[a, b, c, d, 0, f]-alpha[a, b, c, d, 1, f]+1)
    
    
    # gamma-C
    # for a, b, c, d, e in beta.keys():
    #     MODEL.addConstr(beta[a, b, c, d, e] <=  alpha[a, b, c, d, 0, e] + alpha[a, b, c, d, 1, e])
    #     MODEL.addConstr(beta[a, b, c, d, e] >=  alpha[a, b, c, d, 0, e] - alpha[a, b, c, d, 1, e])
    #     MODEL.addConstr(beta[a, b, c, d, e] >= -alpha[a, b, c, d, 0, e] + alpha[a, b, c, d, 1, e])
    #     MODEL.addConstr(beta[a, b, c, d, e] <= 2 - alpha[a, b, c, d, 0, e] - alpha[a, b, c, d, 1, e])
                
        
    
    # delta-C
    # for a, b, c, d in delta.keys():
    #     MODEL.addConstr(delta[a, b, c, d] >= beta[a, b, c, d, 0] + beta[a, b, c, d, 1] - 1)
        
    # epsilon-C and y-C
    for a, b, c in epsilon.keys():
        # MODEL.addConstr((delta.sum(a, b, c, '*')-len(barriers)) + 1 <= epsilon[a, b, c])
        MODEL.addConstr(len(barriers)*(1 - epsilon[a, b, c]) >= delta.sum(a, b, c, '*'))
        
        # MODEL.addConstr(delta.sum(a, b, c, '*') <= len(barriers)*(1- epsilon[a, b, c]))
        
        # Si algun delta es 1, entonces epsilon tiene que ser 0
        # MODEL.addConstr(len(barriers)*epsilon[a, b, c] <= len(barriers) - delta.sum(a, b, c, '*'))
        # MODEL.addConstr(1 - delta.sum(a, b, c, '*') <= len(barriers)*epsilon[a, b, c])
        
        # Si todos los deltas son 0, entonces epsilon tiene que ser 1
        
        # 1 - epsilon <= suma
        # epsilon <= suma
        # MODEL.addConstr(len(barriers)*epsilon[a, b, c] <= len(barriers) - delta.sum(a, b, c, '*'))
        
        # MODEL.addConstr(len(barriers)*epsilon[a, b, c] <= len(barriers) - delta.sum(a, b, c, '*'))
        
        # MODEL.addConstr(1 - epsilon[a, b, c] <= delta.sum(a, b, c, '*'))
        MODEL.addConstr(y[a, b, c] + y[b, c, a] <= 2*epsilon[a, b, c])
        
        # if no_ve(N[abs(a)-1], barriers[b][0]) and no_ve(N[abs(a)-1], barriers[b][1]):
        #     MODEL.addConstr(epsilon[a, b, c] <= 0)
    
    
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
        
    
    
    MODEL.addConstrs(gp.quicksum(y[v, i, vn] for v, i in VB if (v, i, vn) in EN) >= 1 for vn in VN)
    
    for v in V:
        if v in VN:
            vn = v
            MODEL.addConstr(gp.quicksum(y[v, i, vn] for v, i in VB  if (v, i, vn) in EN) == gp.quicksum(y[vn, v, i] for v, i in VB if (vn, v, i) in EN))
            if vn <= -2:
                MODEL.addConstr(gp.quicksum(g[vn, v, i] for v, i in VB if (vn, v, i) in EN) - gp.quicksum(g[v, i, vn] for v, i in VB if (v, i, vn) in EN) == 1)
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
    
    objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(dist[index]*y[index] for index in EB) + gp.quicksum(10000*(delta[index]) for index in delta.keys())
    MODEL.setObjective(objective, GRB.MINIMIZE)
    
    second_time = time.time()
    
    time_elapsed = second_time - first_time
    
    MODEL.update()
    
    MODEL.Params.Threads = 6
    MODEL.Params.timeLimit = timeLimit - time_elapsed
    # MODEL.Params.NonConvex = 2
    # MODEL.presolve()
    
    MODEL.write('prueba.lp')
    


    MODEL.optimize()

    
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

        return resultados
    
    if MODEL.SolCount == 0:
        if init:
            try:
                resultados.append(resultado_h[0])
                resultados.append(resultado_h[1])
            except:
                resultados.append(np.nan)
                resultados.append(np.nan)
    
        return resultados
    
    y_indices = []
    
    for index in E:
        if y[index].X > 0.5:
            y_indices.append(index)
            
    print(y_indices)
    
    g_indices = []
    
    for index in g_index:
        if g[index].X > 0.5:
            g_indices.append(g[index])
    
    print(g_indices)
    
    
    MODEL.write('solucion.sol')
    
    resultados = []
    
    resultados.append(len(N))
    resultados.append(len(barriers))
    
    resultados.append(MODEL.getAttr('MIPGap'))
    resultados.append(MODEL.Runtime + time_elapsed)
    resultados.append(MODEL.getAttr('NodeCount'))
    resultados.append(MODEL.ObjVal)
    
    if init:
        try:
            resultados.append(resultado_h[0])
            resultados.append(resultado_h[1])
        except:
            resultados.append(np.nan)
            resultados.append(np.nan)
            print('El heuristico no ha encontrado initial sol')
    
    if log:
        fig, ax = plt.subplots()
        
        for b in barriers:
            ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c = 'red')
        
        for n in N:
            ax.add_artist(n.artist)
        
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
