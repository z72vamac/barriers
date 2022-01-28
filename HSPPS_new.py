# Con la formulacion nueva
# Shortest path por entornos con barriers

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product, permutations, chain
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
from entorno import Circulo, Poligonal
import copy
import estimacion_M as eM
from auxiliar_functions import *
import networkx as nx
import time


#
# VS = [-1]
# ES = list(product(VS, range(len(barriers)), range(2))) #+ list(product(range(len(barriers)), range(2), VS))

def HSPPS(barriers, N, prepro = True, log = False, timeLimit = 7200, init = False):
    
    first_time = time.time()
    
    def nocortan(barrier1, barrier2):
        det1 = (barrier1[0][0] - barrier2[0][0])*(barrier1[1][1] - barrier2[0][1]) - (barrier1[1][0] - barrier2[0][0])*(barrier1[0][1] - barrier2[0][1])
        det2 = (barrier1[0][0] - barrier2[1][0])*(barrier1[1][1] - barrier2[1][1]) - (barrier1[1][0] - barrier2[1][0])*(barrier1[0][1] - barrier2[1][1])
    
        det3 = (barrier2[0][0] - barrier1[0][0])*(barrier2[1][1] - barrier1[0][1]) - (barrier2[1][0] - barrier1[0][0])*(barrier2[0][1] - barrier1[0][1])
        det4 = (barrier2[0][0] - barrier1[1][0])*(barrier2[1][1] - barrier1[1][1]) - (barrier2[1][0] - barrier1[1][0])*(barrier2[0][1] - barrier1[1][1])
    
        return ((det1*det2) >= 0) or ((det3*det4) >= 0)

    def no_ve(entorno, punto):
        if type(entorno) is Circulo:
            centro = entorno.center
            
            dr = np.array(centro) - np.array(punto)
            dr_u = dr / np.linalg.norm(dr)
            
            nr_u = np.array([-dr_u[1], dr_u[0]])
            
            lambdas = np.linspace(-entorno.radii, entorno.radii, 20)
            
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
                           
            # punto1 = np.array(centro) + entorno.radii*nr_u
            # segmento1 = [punto, punto1]
            # flag1 = any([not(nocortan(segmento1, barrera)) for barrera in barriers])
            #
            # punto2 = np.array(centro) - entorno.radii*nr_u
            # segmento2 = [punto, punto2]
            # flag2 = any([not(nocortan(segmento2, barrera)) for barrera in barriers])
            
        if type(entorno) is Poligonal:
            
            dr = np.array(entorno.V[1]) - np.array(entorno.V[0])
            
            dr_u = dr / np.linalg.norm(dr)
            
            lambdas = np.linspace(0, entorno.longitud, 20)

            flag2 = True

            for landa in lambdas:
                
                for barrera in barriers:
                    flag1 = False
                    if not(nocortan([punto, np.array(entorno.V[0]) + landa*dr_u], [barrera[0], barrera[1]])):
                        flag1 = True
                        break
                
                flag2 = flag2 & flag1
                
                if not(flag2):
                    return flag2
                
            return flag2        
            
        
    def ve(entorno, punto):
        return not(no_ve(entorno, punto))
        
    VN = [(-1, 0), (-2, 0)]
    EN = []
    
    for (a, b) in VN:
        for c in range(len(barriers)):
            for d in range(2):
                if prepro:
                    if ve(N[abs(a)-1], barriers[c][0]) or ve(N[abs(a)-1], barriers[c][1]):
                        EN.append((a, b, c, d))
                        EN.append((c, d, a, b))
                else:
                    EN.append((a, b, c, d))
                    EN.append((c, d, a, b))       
    
    
    VB = list(product(range(len(barriers)), range(2)))
    
    EB = []
    for v, i in VB:
        for w, j in VB:
            if v != w:
                barrier = [barriers[v][i], barriers[w][j]]
                
                if all([nocortan(barrieri, barrier) for barrieri in barriers]):
                    EB.append((v, i, w, j))
    
    # EB = list(product(range(len(barriers)), range(2), range(len(barriers)), range(2)))
        
    V = VN + VB
    E = EN + EB
    
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
    epsilon_index = []
    
    for a, b, c, d in EN:
        if a < 0:
            epsilon_index.append((a, b, c, d))
            
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
    
    for a, b, c, d in E:
        for dim in range(2):
            dif_index.append((a, b, c, d, dim))
    
    if log:    
        print("dif = " + str(dif_index))
    
    
    # delta(S / T, B, i, B') = 1 si (P_{S/T}, P_B^i) y B' do not intersect.
    delta_index = [] #list(product(VN, range(len(barriers)), range(2), range(len(barriers)))) 
    
    
    for a, b, c, d in epsilon_index:
        for e in range(len(barriers)):
            delta_index.append((a, b, c, d, e))

    
    if log:
        print("delta = " + str(delta_index))
    
    # gamma(S / T, B, i, B') = alpha(S / T, B, i, B', 0, order2)*alpha(S / T, B, i, B', 1, order2)
    gamma_index = [] # list(product(VN, range(len(barriers)), range(2), range(len(barriers)), range(2)))

    for a, b, c, d, e in delta_index:
        for f in range(2):
            gamma_index.append((a, b, c, d, e, f))
    
    
    if log:
        print("gamma = " + str(gamma_index))
    
    # beta(S/T, B, i, B', order) = 1 si sign(P_{S/T} P_B^i | B') is the same (order = 0); sign(B' | P_{S/T} P_B^i) is the same (order = 1);
    beta_index = gamma_index
    
    beta0_index = []
    beta1_index = []
    
    for a, b, c, d, e, h in beta_index:
        if h == 0:
            beta0_index.append((a, c, d, e))
        else:
            beta1_index.append((a, c, d, e))
    
    beta0_index = list(set(beta0_index))
    beta1_index = list(set(beta1_index))
    
    if log:
        print("beta = " + str(beta_index))
    
    # alpha(S/T, B, i, B', order1, order2) = 1 si sign(P_{S/T} | B') > 0 (order1 = 0); sign(P_B^i | B') > 0 (order1) = 1;
    alpha_index = [] # list(product(VN, range(len(barriers)), range(2), range(len(barriers)), range(2), range(2)))

    for a, b, c, d, e, f in gamma_index:
        for g in range(2):
            alpha_index.append((a, b, c, d, e, f, g))
                    
    if log:
        print("alpha = " + str(alpha_index))
        
    alpha00_index = []
    alpha10_index = []
    alpha_1_index = []

    for a, b, c, d, e, f, h in alpha_index:
        if f == 0 and h == 0:
            alpha00_index.append((a, e))
        
        if f == 1 and h == 0:
            alpha10_index.append((c, d, e))
            
        if h == 1:
            alpha_1_index.append((a, c, d, e, f))
            
    alpha00_index = list(set(alpha00_index))
    alpha10_index = list(set(alpha10_index))
    alpha_1_index = list(set(alpha_1_index))
    
    
    # P_S and P_T: indices of the points in the neighborhoods
    P_index = []
    for index in VN:
        for dim in range(2):
            P_index.append((index[0], index[1], dim))
    
    if log:
        print("P_index = " + str(P_index))
    
    # socp variables:
    d_inside_index = VN
    
    dif_inside_index = []
    
    for (a, b) in VN:
        for dim in range(2):
            dif_inside_index.append((a, b, dim))
    
    if log:
        print("d_inside_index = " + str(d_inside_index))
        print("dif_inside_index = " + str(dif_inside_index))
    
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
    # beta = MODEL.addVars(beta_index, vtype = GRB.BINARY, name = 'beta')
    # alpha = MODEL.addVars(alpha_index, vtype = GRB.BINARY, name = 'alpha')
    
    alpha00 = MODEL.addVars(alpha00_index, vtype = GRB.BINARY, name = 'alpha00')
    alpha10 = MODEL.addVars(alpha10_index, vtype = GRB.BINARY, name = 'alpha10')
    alpha_1 = MODEL.addVars(alpha_1_index, vtype = GRB.BINARY, name = 'alpha_1')
    
    beta0 = MODEL.addVars(beta0_index, vtype = GRB.BINARY, name = 'beta0')
    beta1 = MODEL.addVars(beta1_index, vtype = GRB.BINARY, name = 'beta1')
    
    
    # P = MODEL.addVars(P_index, vtype = GRB.CONTINUOUS, lb = 0.1, ub = 99.9, name = 'P')
    P = MODEL.addVars(P_index, vtype = GRB.CONTINUOUS, name = 'P')
    
    d_inside = MODEL.addVars(d_inside_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'd_inside')
    dif_inside = MODEL.addVars(dif_inside_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dif_inside')
    landa = MODEL.addVars(d_inside_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'landa')
    # z = MODEL.addVars(z_index, vtype = GRB.BINARY, name = 'z')
    
    MODEL.update()
    
    # Alpha constraint
    L = -10000
    U = 10000
    
    # alpha-C
    for a, e in alpha00.keys():
        L, U = estima_det(N[abs(a)-1], barriers[e])
        
        if U <= 0:
            alpha00[a, e] = 0
        
        elif L >= 0:
            alpha00[a, e] = 1
        
        else:
        
            MODEL.addConstr((1 - alpha00[a, e])*L <= determinant([P[a, 0, 0], P[a, 0, 1]], barriers[e][0], barriers[e][1]))
            MODEL.addConstr(determinant([P[a, 0, 0], P[a, 0, 1]], barriers[e][0], barriers[e][1]) <= U*alpha00[a, e])
     
    for c, d, e in alpha10.keys():
        if determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0:
            alpha10[c, d, e] = (determinant(barriers[c][d], barriers[e][0], barriers[e][1]) / abs(determinant(barriers[c][d], barriers[e][0], barriers[e][1])) >= 0)*1
        
        # MODEL.addConstrs((1 - alpha10[c, d, e])*abs(determinant(barriers[c][d], barriers[e][0], barriers[e][1]))*(-1) <= determinant(barriers[c][d], barriers[e][0], barriers[e][1]) for c, d, e in alpha10.keys() if determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0)
        # MODEL.addConstrs(abs(determinant(barriers[c][d], barriers[e][0], barriers[e][1]))*alpha10[c, d, e] >= determinant(barriers[c][d], barriers[e][0], barriers[e][1]) for c, d, e in alpha10.keys() if determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0) 
    
    # MODEL.addConstrs((1 - alpha10[c, d, e])*abs(determinant(barriers[c][d], barriers[e][0], barriers[e][1]))*(-1) <= determinant(barriers[c][d], barriers[e][0], barriers[e][1]) for c, d, e in alpha10.keys() if determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0)
    # MODEL.addConstrs(abs(determinant(barriers[c][d], barriers[e][0], barriers[e][1]))*alpha10[c, d, e] >= determinant(barriers[c][d], barriers[e][0], barriers[e][1]) for c, d, e in alpha10.keys() if determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0) 

    # MODEL.addConstrs((1 - alpha_1[a, c, d, e, f])*estima_det(N[abs(a)-1], [barriers[e][f], barriers[c][d]])[0] <= determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]) for a, c, d, e, f in alpha_1.keys())
    # MODEL.addConstrs(determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]) <= estima_det(N[abs(a)-1], [barriers[e][f], barriers[c][d]])[1]*alpha_1[a, c, d, e, f] for a, c, d, e, f in alpha_1.keys())

    L = -20000
    U = 20000
    
    MODEL.addConstrs((1 - alpha_1[a, c, d, e, f])*L <= determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]) for a, c, d, e, f in alpha_1.keys())
    MODEL.addConstrs(determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]) <= U*alpha_1[a, c, d, e, f] for a, c, d, e, f in alpha_1.keys())

    # MODEL.addConstrs((1 - alpha_1[a, c, d, e, f])*estima_det(N[abs(a)-1], [barriers[e][f], barriers[c][d]])[0] <= determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]) for a, c, d, e, f in alpha_1.keys())
    # MODEL.addConstrs(determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]) <= estima_det(N[abs(a)-1], [barriers[e][f], barriers[c][d]])[1]*alpha_1[a, c, d, e, f] for a, c, d, e, f in alpha_1.keys())

    
    # MODEL.addConstr(alpha00[-2,0] == 0)
    # MODEL.addConstr(alpha00[-2,1] == 1)
    
    
    # beta-C
    # for a, b, c, d, e, h in beta.keys():
    #     MODEL.addConstr(beta[a, b, c, d, e, h] == 2*gamma[a, b, c, d, e, h] - alpha[a, b, c, d, e, 0, h]-alpha[a, b, c, d, e, 1, h]+1)
    # for a, b, c, d, e, h in beta.keys():
    #     if h == 0:
    #         MODEL.addConstr(beta[a, b, c, d, e, h] == 2*gamma[a, b, c, d, e, h] - alpha00[a, e] - alpha10[c, d, e] + 1)
    #     else:
    #         MODEL.addConstr(beta[a, b, c, d, e, h] == 2*gamma[a, b, c, d, e, h] - alpha_1[a, c, d, e, 0] - alpha_1[a, c, d, e, 1] + 1)

    # MODEL.addConstrs(beta0[a, c, d, e] == 2*alpha00[a, e]*alpha10[c, d, e] - alpha00[a, e] - alpha10[c, d, e] + 1 for a, c, d, e in beta0.keys())
    # MODEL.addConstrs(beta1[a, c, d, e] == 2*alpha_1[a, c, d, e, 0]*alpha_1[a, c, d, e, 1] - alpha_1[a, c, d, e, 0] - alpha_1[a, c, d, e, 1] + 1 for a, c, d, e in beta1.keys())

    MODEL.addConstrs(beta0[a, c, d, e] == 2*gamma[a, 0, c, d, e, 0] - alpha00[a, e] - alpha10[c, d, e] + 1 for a, c, d, e in beta0.keys())
    MODEL.addConstrs(beta1[a, c, d, e] == 2*gamma[a, 0, c, d, e, 1] - alpha_1[a, c, d, e, 0] - alpha_1[a, c, d, e, 1] + 1 for a, c, d, e in beta1.keys())
    
    # gamma-C
    for a, b, c, d, e, h in gamma.keys():
        if h == 0:
            MODEL.addConstr(gamma[a, b, c, d, e, h] <= alpha00[a, e])
            MODEL.addConstr(gamma[a, b, c, d, e, h] <= alpha10[c, d, e])
            MODEL.addConstr(gamma[a, b, c, d, e, h] >= alpha00[a, e] + alpha10[c, d, e] - 1)
    
        if h == 1:
            MODEL.addConstr(gamma[a, b, c, d, e, h] <= alpha_1[a, c, d, e, 0])
            MODEL.addConstr(gamma[a, b, c, d, e, h] <= alpha_1[a, c, d, e, 1])
            MODEL.addConstr(gamma[a, b, c, d, e, h] >= alpha_1[a, c, d, e, 0] + alpha_1[a, c, d, e, 1] - 1)            
        
        
    
    # delta-C
    for a, c, d, e in beta0.keys():
        MODEL.addConstr((beta0[a, c, d, e] + beta1[a, c, d, e] <= 2*delta[a, 0, c, d, e]))
        MODEL.addConstr(delta[a, 0, c, d, e] <= 2*(beta0[a, c, d, e] + beta1[a, c, d, e]))
        
    # epsilon-C and y-C
    for a, b, c, d in epsilon.keys():
        # print((a, b, c, d))
        MODEL.addConstr((delta.sum(a, b, c, d, '*')-len(barriers)) + 1 <= epsilon[a, b, c, d])
        MODEL.addConstr(len(barriers)*epsilon[a, b, c, d] <= delta.sum(a, b, c, d, '*'))
        
        # if a < 0 and b >= 0 and c >= 0:
        MODEL.addConstr(y[a, b, c, d] + y[c, d, a, b] <= 2*epsilon[a, b, c, d])
        
        # if a < 0 and b < 0 and c == 0:
        #     MODEL.addConstr(y[a, b, c] + y[b, a, c] <= 2*epsilon[a, b, c])
            
        
            
        
        # if no_ve(N[abs(a)-1], barriers[b][0]) and no_ve(N[abs(a)-1], barriers[b][1]):
        #     MODEL.addConstr(epsilon[a, b, c] <= 0)
    
    
    # NS and NT constraints
    for a, b, dim in P.keys():
        entorno = N[abs(a)-1]
        
        if type(entorno) is Circulo:
            MODEL.addConstr(dif_inside[a, b, dim] >= P[a, b, dim] - N[abs(a)-1].center[dim])
            MODEL.addConstr(dif_inside[a, b, dim] >= N[abs(a)-1].center[dim] - P[a, b, dim])
        
            MODEL.addConstr(gp.quicksum(dif_inside[a, b, dim]*dif_inside[a, b, dim] for dim in range(2)) <= d_inside[a, b]*d_inside[a, b])
            MODEL.addConstr(d_inside[a, b] <= N[abs(a)-1].radii)
        
        if type(entorno) is Poligonal:
            MODEL.addConstrs(P[a, b, dim] == landa[a, b]*entorno.V[0][dim] + (1-landa[a, b])*entorno.V[1][dim] for dim in range(2))
        
    
    # MODEL.addConstr(y[-1, 2, 1] >= 0.5)
    # MODEL.addConstr(y[-2, 2, 1] >= 0.5)
    # MODEL.addConstr(P[-1, 0] == 30)
    # MODEL.addConstr(P[-1, 1] == 10)
    
    # MODEL.addConstr(epsilon[-4, 3, 1] >= 0.5)
    # MODEL.addConstr(epsilon[3, 1, -4] >= 0.5)
    
    # MODEL.addConstr(y[3, 1, -1] >= 0.5)
    
    # dist constraints
    for a, b, c, d, dim in dif_index:
        if (a, b, c, d) in EB:
            MODEL.addConstr(dist[a, b, c, d] == np.linalg.norm(np.array(barriers[a][b])-np.array(barriers[c][d])))
            
        if (a, b, c, d) in EN:
            if a < 0:
                # MODEL.addConstr(dist[a, b, c, d] <= estima_L(N[abs(a)-1], barriers[c][d]))
                MODEL.addConstr(dif[a, b, c, d, dim] >=   P[a, b, dim] - barriers[c][d][dim])
                MODEL.addConstr(dif[a, b, c, d, dim] >= - P[a, b, dim] + barriers[c][d][dim])
                MODEL.addConstr(gp.quicksum(dif[a, b, c, d, dim]*dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d]*dist[a, b, c, d])
            if c < 0:
                # MODEL.addConstr(dist[a, b, c, d] <= estima_L(N[abs(c)-1], barriers[a][b]))

                MODEL.addConstr(dif[c, d, a, b, dim] >=   P[c, d, dim] - barriers[a][b][dim])
                MODEL.addConstr(dif[c, d, a, b, dim] >= - P[c, d, dim] + barriers[a][b][dim])
                MODEL.addConstr(gp.quicksum(dif[c, d, a, b, dim]*dif[c, d, a, b, dim] for dim in range(2)) <= dist[c, d, a, b]*dist[c, d, a, b])

        # if (a, b, c, d) in ENN:
            # MODEL.addConstr(dif[a, b, c, d, dim] >=   P[a, b, dim] - P[c, d, dim])
            # MODEL.addConstr(dif[a, b, c, d, dim] >= - P[a, b, dim] + P[c, d, dim])
            # MODEL.addConstr(gp.quicksum(dif[a, b, c, d, dim]*dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d]*dist[a, b, c, d])
    
    L_out = 0
    U_out = 1000000
            
    # p constraints
    for a, b, c, d in p.keys():
            
        if a < 0:
            entorno = N[abs(a) - 1]
            punto = barriers[c][d]
            L_out = estima_L(entorno, punto)
            U_out = estima_U(entorno, punto)
    
        if c < 0:
            entorno = N[abs(c) - 1]
            punto = barriers[a][b]
            L_out = estima_L(entorno, punto)
            U_out = estima_U(entorno, punto)

        
        # print((L_out, U_out))
        MODEL.addConstr(p[a, b, c, d] >= L_out*y[a, b, c, d])
        MODEL.addConstr(p[a, b, c, d] >= dist[a, b, c, d] - U_out*(1 - y[a, b, c, d]))

        
    
    
    # flow conservation constraints
    # for index in y_index:
    #     if len(index) == 3:
    # MODEL.addConstr(gp.quicksum(y[tupla] for tupla in E if tupla[0] == -1) == 1)
    #
    # for v, i in VB:
    #     tuplas_salen = gp.quicksum([y[a, b, c, d] for a, b, c, d in E if (a == v and b == i)])
    #     tuplas_entran = gp.quicksum([y[a, b, c, d] for a, b, c, d in E if (c == v and d == i)])
    #
    #     MODEL.addConstr(tuplas_salen - tuplas_entran == 0)
    #
    # MODEL.addConstr(- gp.quicksum(y[a, b, c, d] for a, b, c, d in E if c == -2) == -1)
    
    for v, i in V:
        if v == -1:
            MODEL.addConstr(gp.quicksum(y[v, i, w, j] for w, j in V if (v, i, w, j) in E) - gp.quicksum(y[w, j, v, i] for w, j in V if (w, j, v, i) in E) == 1) 
        elif v == -2:
            MODEL.addConstr(gp.quicksum(y[v, i, w, j] for w, j in V if (v, i, w, j) in E) - gp.quicksum(y[w, j, v, i] for w, j in V if (w, j, v, i) in E) == -1) 
        
        else:
            MODEL.addConstr(gp.quicksum(y[v, i, w, j] for w, j in V if (v, i, w, j) in E) - gp.quicksum(y[w, j, v, i] for w, j in V if (w, j, v, i) in E) == 0) 
           
          
        
    
    MODEL.update()
    
    objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(dist[index]*y[index] for index in EB)
    MODEL.setObjective(objective, GRB.MINIMIZE)
    
    second_time = time.time()
    
    time_elapsed = second_time - first_time
    
    MODEL.update()
    
    MODEL.Params.Threads = 6
    MODEL.Params.timeLimit = timeLimit - time_elapsed
    # MODEL.Params.LazyConstraints = 1
    MODEL.Params.NumericFocus = 1
    # MODEL.Params.NonConvex = 2
        
    MODEL.write('prueba.lp')
    MODEL.write('prueba.mps')
    
    


    MODEL.optimize()
    
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
    
    try:
        MODEL.write('solucion.sol')
    except:
        MODEL.computeIIS()
        MODEL.write('prueba.ilp')


    y_indices = []
    
    for index in E:
        if y[index].X > 0.5:
            y_indices.append(index)
    
    if log:
        print(y_indices)
    
    # g_indices = []
    #
    # for index in g_index:
    #     if g[index].X > 0.5:
    #         g_indices.append(g[index])
    #
    # if log: 
    #     print(g_indices)

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
            ax.add_artist(n.artist)
        
        P_vals = MODEL.getAttr('x', P)
        print(P_vals)
        
        points = []
        for keys, vals in P_vals.items():
            points.append(vals)
            
        points = np.array(points).reshape((len(N), 2))
        print(points)
        
        for i in points:
            ax.scatter(i[0], i[1], s = 10, c = 'black')
        
        # print(points)
        
        segmentos = []
        
        for a, b, c, d in y_indices:
            if (a, b, c, d) in EN:
                if a < 0:
                    segmentos.append([points[abs(a)-1][0], barriers[c][d][0], points[abs(a)-1][1], barriers[c][d][1]])
                if c < 0:
                    segmentos.append([barriers[a][b][0], points[abs(c)-1][0], barriers[a][b][1], points[abs(c)-1][1]])
            if (a, b, c, d) in EB:
                segmentos.append([barriers[a][b][0], barriers[c][d][0], barriers[a][b][1], barriers[c][d][1]])
                
            # if (a, b, c, d) in ENN:
            #         segmentos.append([points[abs(a)-1][0], points[abs(c)-1][0], points[abs(a)-1][1], points[abs(c)-1][1]])
                
                
                
            
        # print(segmentos)
        for segmento in segmentos:
            ax.arrow(segmento[0], segmento[2], segmento[1]-segmento[0], segmento[3]-segmento[2], width = 0.1, head_width = 1, length_includes_head = True, color = 'black')
        
        
        
        
        # plt.axis([-5, 105, -5, 105])
        plt.axis([0, 100, 0, 100])
        
        ax.set_aspect('equal')
        plt.show()

    return resultados
# barrier1 = [[20, 80], [40, 30]]
# barrier2 = [[70, 95], [40, 70]]
# barrier3 = [[95, 60], [60, 70]]
# barrier4 = [[60, 50], [90, 10]]
# barrier5 = [[10, 70], [20, 50]]
# # barrier6 = [[30, 70], [70, 20]]
#
# barriers = [barrier1, barrier2, barrier3, barrier4, barrier5]
# # barriers.append(barrier6)
#
# NS = Circulo(center = [20, 10], radii = 10)
# NT = Circulo(center = [90, 90], radii = 5)
#
# N = [NS, NT]
#
# HSPPS(barriers, N, log = True)
    
    # MODEL.addConstr(gp.quicksum( y[v, w] for v, w in edges_tuple.select(0, '*')) - gp.quicksum( y[w, v] for w, v in edges_tuple.select('*', 0)) ==  1)
    # MODEL.write('prueba.lp')
