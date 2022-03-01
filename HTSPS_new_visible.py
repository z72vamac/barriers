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


def HTSPS_ven(barriers, N, prepro=True, log=False, timeLimit=7200, init=False):
    # N = [N1, N2]
    
    first_time = time.time()

    def nocortan(barrier1, barrier2):
        det1 = (barrier1[0][0] - barrier2[0][0]) * (barrier1[1][1] - barrier2[0][1]) - (barrier1[1][0] - barrier2[0][0]) * (barrier1[0][1] - barrier2[0][1])
        det2 = (barrier1[0][0] - barrier2[1][0]) * (barrier1[1][1] - barrier2[1][1]) - (barrier1[1][0] - barrier2[1][0]) * (barrier1[0][1] - barrier2[1][1])
    
        det3 = (barrier2[0][0] - barrier1[0][0]) * (barrier2[1][1] - barrier1[0][1]) - (barrier2[1][0] - barrier1[0][0]) * (barrier2[0][1] - barrier1[0][1])
        det4 = (barrier2[0][0] - barrier1[1][0]) * (barrier2[1][1] - barrier1[1][1]) - (barrier2[1][0] - barrier1[1][0]) * (barrier2[0][1] - barrier1[1][1])
    
        return ((det1 * det2) >= 0) or ((det3 * det4) >= 0)
    
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
                    if not(nocortan([punto, np.array(centro) + landa * nr_u], [barrera[0], barrera[1]])):
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
                    if not(nocortan([punto, np.array(neighborhood.V[0]) + landa * nr_u], [barrera[0], barrera[1]])):
                        flag1 = True
                        break
                
                flag2 = flag2 & flag1
                
                if not(flag2):
                    return flag2
                
            return flag2        
        
    def ve(neighborhood, punto):
        return not(no_ve(neighborhood, punto))
    
    VN = list(product(range(-len(N), 0), range(1)))
    VN = VN[::-1]
    EN = []  # list(product(VN, range(len(barriers)), range(2))) + list(product(range(len(barriers)), range(2), VN))
    
    for (a, b) in VN:
        for c in range(len(barriers)):
            for d in range(2):
                if prepro:
                    if ve(N[abs(a) - 1], barriers[c][0]) or ve(N[abs(a) - 1], barriers[c][1]):
                        EN.append((a, b, c, d))
                        EN.append((c, d, a, b))
                else:
                    EN.append((a, b, c, d))
                    EN.append((c, d, a, b))                     

    ENN = []
    
    for (a, b) in VN:
        for (c, d) in VN:
            if a != c:
                ENN.append((a, b, c, d))
                # ENN.append((b, 0, a))
            
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
    E = EN + EB + ENN
    
    if log:
        print("VN = " + str(VN))
        print("VB = " + str(VB))
        
        print("EN = " + str(EN))
        print("EB = " + str(EB))
        
        print("ENN = " + str(ENN))

    # print("VN = " + str(VN))
    # print("VB = " + str(VB))
    #
    # print("EN = " + str(EN))
    # print("EB = " + str(EB))    
    
    # print(barriers)
    
    def determinant(P, Q, R):
        a11 = Q[0] - P[0]
        a12 = R[0] - P[0]
        a21 = Q[1] - P[1]
        a22 = R[1] - P[1]
        
        return a11 * a22 - a12 * a21
        
    # epsilon(S / T, B, i) = 1 si (P_{S/T}, P_B^i)\in E_{S/T}
    epsilon_index = []
    
    for a, b, c, d in EN + ENN:
        if a < 0:
            epsilon_index.append((a, b, c, d))
            
    # print(epsilon_index)
    
    p_index = EN + ENN
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
    delta_index = []  # list(product(VN, range(len(barriers)), range(2), range(len(barriers)))) 
    
    for a, b, c, d in epsilon_index:
        for e in range(len(barriers)):
            delta_index.append((a, b, c, d, e))
    
    if log:
        print("delta = " + str(delta_index))
    
    # gamma(S / T, B, i, B') = alpha(S / T, B, i, B', 0, order2)*alpha(S / T, B, i, B', 1, order2)
    gamma_index = []  # list(product(VN, range(len(barriers)), range(2), range(len(barriers)), range(2)))

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
    alpha_index = []  # list(product(VN, range(len(barriers)), range(2), range(len(barriers)), range(2), range(2)))

    for a, b, c, d, e, f in gamma_index:
        for g in range(2):
            alpha_index.append((a, b, c, d, e, f, g))

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

    if log:
        print("alpha00_index = " + str(alpha00_index))
        print("alpha10_index = " + str(alpha10_index))
        print("alpha_1_index = " + str(alpha_1_index))
        
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
    
    # z variables:
    # z_index = list(product(VN, VN))
    # print("z_index = " + str(z_index))
    
    # f variables:
    g_index = y_index
    
    if log:
        print("g_index = " + str(g_index))
    
    MODEL = gp.Model('HSPPS_Model')
    
    epsilon = MODEL.addVars(epsilon_index, vtype=GRB.BINARY, name='epsilon')
    p = MODEL.addVars(p_index, vtype=GRB.CONTINUOUS, lb=0.0, name='p')
    y = MODEL.addVars(y_index, vtype=GRB.BINARY, name='y')
    dist = MODEL.addVars(dist_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dist')
    dif = MODEL.addVars(dif_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')
    delta = MODEL.addVars(delta_index, vtype=GRB.BINARY, name='delta')
    gamma = MODEL.addVars(gamma_index, vtype=GRB.BINARY, name='gamma')
    # beta = MODEL.addVars(beta_index, vtype = GRB.BINARY, name = 'beta')
    # alpha = MODEL.addVars(alpha_index, vtype = GRB.BINARY, name = 'alpha')
    
    alpha00 = MODEL.addVars(alpha00_index, vtype=GRB.BINARY, name='alpha00')
    alpha10 = MODEL.addVars(alpha10_index, vtype=GRB.BINARY, name='alpha10')
    alpha_1 = MODEL.addVars(alpha_1_index, vtype=GRB.BINARY, name='alpha_1')
    
    beta0 = MODEL.addVars(beta0_index, vtype=GRB.BINARY, name='beta0')
    beta1 = MODEL.addVars(beta1_index, vtype=GRB.BINARY, name='beta1')
        
    P = MODEL.addVars(P_index, vtype=GRB.CONTINUOUS, name='P')
    
    d_inside = MODEL.addVars(d_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='d_inside')
    dif_inside = MODEL.addVars(dif_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif_inside')
    landa = MODEL.addVars(d_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landa')
    # z = MODEL.addVars(z_index, vtype = GRB.BINARY, name = 'z')
    g = MODEL.addVars(g_index, vtype=GRB.CONTINUOUS, lb=0.0, name='g')
    
    MODEL.update()
    
    # Alpha constraint
    # print(alpha00.keys())
    # alpha-C
    for a, e in alpha00.keys():
        
        L, U = estima_det(N[abs(a) - 1], barriers[e])
        
        if U <= 0:
            alpha00[a, e] = 0
        
        elif L >= 0:
            alpha00[a, e] = 1
        
        else:
            MODEL.addConstr((1 - alpha00[a, e]) * L <= determinant([P[a, 0, 0], P[a, 0, 1]], barriers[e][0], barriers[e][1]))
            MODEL.addConstr(determinant([P[a, 0, 0], P[a, 0, 1]], barriers[e][0], barriers[e][1]) <= U * alpha00[a, e])
    
    L = -20000
    U = 20000
    
    for c, d, e in alpha10.keys():
        if (a, 0, c, d) in EN:
            if determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0:
                alpha10[c, d, e] = (determinant(barriers[c][d], barriers[e][0], barriers[e][1]) / abs(determinant(barriers[c][d], barriers[e][0], barriers[e][1])) >= 0) * 1

        elif (a, 0, c, d) in ENN:
            MODEL.addConstr((1 - alpha10[c, d, e]) * L <= determinant([P[c, 0, 0], P[c, 0, 1]], barriers[e][0], barriers[e][1]))
            MODEL.addConstr(determinant([P[c, 0, 0], P[c, 0, 1]], barriers[e][0], barriers[e][1]) <= U * alpha10[c, d, e])   
    
    for a, c, d, e, f in alpha_1.keys():
        if (a, 0, c, d) in EN:
            MODEL.addConstr((1 - alpha_1[a, c, d, e, f]) * L <= determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]))
            MODEL.addConstr(determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]) <= U * alpha_1[a, c, d, e, f])
        
        if (a, 0, c, d) in ENN:
            MODEL.addConstr((1 - alpha_1[a, c, d, e, f]) * L <= determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], [P[c, 0, 0], P[c, 0, 1]]))
            MODEL.addConstr(determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], [P[c, 0, 0], P[c, 0, 1]]) <= U * alpha_1[a, c, d, e, f])    
            
    # for a, b, c, d, e, f, h in alpha.keys():
    #     # print((a, b, c, d, e, f, h))
    #     if (a, b, c, d) in EN:
    #         if f == 0 and h == 0:
    #             MODEL.addConstr((1 - alpha[a, b, c, d, e, f, h])*L <= determinant([P[a, 0, 0], P[a, 0, 1]], barriers[e][0], barriers[e][1]))
    #             MODEL.addConstr(determinant([P[a, 0, 0], P[a, 0, 1]], barriers[e][0], barriers[e][1]) <= U*alpha[a, b, c, d, e, f, h])
    #
    #         # elif e == 1 and f == 0:
    #         #     MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*L <= determinant(barriers[b][c], barriers[d][0], barriers[d][1]))
    #         #     MODEL.addConstr(determinant(barriers[b][c], barriers[d][0], barriers[d][1]) <= U*alpha[a, b, c, d, e, f])
    #
    #         elif f == 1 and h == 0:
    #             MODEL.addConstr((1 - alpha[a, b, c, d, e, f, h])*abs(determinant(barriers[c][d], barriers[e][0], barriers[e][1]))*(-1) <= determinant(barriers[c][d], barriers[e][0], barriers[e][1]))
    #             if determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0:
    #                 MODEL.addConstr(alpha[a, b, c, d, e, f, h] >= determinant(barriers[c][d], barriers[e][0], barriers[e][1])/abs(determinant(barriers[c][d], barriers[e][0], barriers[e][1])))
    #
    #         else:
    #             MODEL.addConstr((1 - alpha[a, b, c, d, e, f, h])*L <= determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]))
    #             MODEL.addConstr(determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], barriers[c][d]) <= U*alpha[a, b, c, d, e, f, h])
    #
    #     if (a, b, c, d) in ENN:
    #         # if b < 0:
    #             if f == 0 and h == 0:
    #                 MODEL.addConstr((1 - alpha[a, b, c, d, e, f, h])*L <= determinant([P[a, 0, 0], P[a, 0, 1]], barriers[e][0], barriers[e][1]))
    #                 MODEL.addConstr(determinant([P[a, 0, 0], P[a, 0, 1]], barriers[e][0], barriers[e][1]) <= U*alpha[a, b, c, d, e, f, h])            
    #
    #             elif f == 1 and h == 0:
    #                 MODEL.addConstr((1 - alpha[a, b, c, d, e, f, h])*L <= determinant([P[c, 0, 0], P[c, 0, 1]], barriers[e][0], barriers[e][1]))
    #                 MODEL.addConstr(determinant([P[c, 0, 0], P[c, 0, 1]], barriers[e][0], barriers[e][1]) <= U*alpha[a, b, c, d, e, f, h])    
    #
    #             else:
    #                 MODEL.addConstr((1 - alpha[a, b, c, d, e, f, h])*L <= determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], [P[c, 0, 0], P[c, 0, 1]]))
    #                 MODEL.addConstr(determinant(barriers[e][f], [P[a, 0, 0], P[a, 0, 1]], [P[c, 0, 0], P[c, 0, 1]]) <= U*alpha[a, b, c, d, e, f, h])                        
            
            # if b == 0:
            #     if e == 0 and f == 0:
            #         MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*L <= determinant([P[a, 0], P[a, 1]], barriers[d][0], barriers[d][1]))
            #         MODEL.addConstr(determinant([P[a, 0], P[a, 1]], barriers[d][0], barriers[d][1]) <= U*alpha[a, b, c, d, e, f])            
            #
            #     elif e == 1 and f == 0:
            #         MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*L <= determinant([P[c, 0], P[c, 1]], barriers[d][0], barriers[d][1]))
            #         MODEL.addConstr(determinant([P[c, 0], P[c, 1]], barriers[d][0], barriers[d][1]) <= U*alpha[a, b, c, d, e, f])    
            #
            #     else:
            #         MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*L <= determinant(barriers[d][e], [P[a, 0], P[a, 1]], [P[c, 0], P[c, 1]]))
            #         MODEL.addConstr(determinant(barriers[d][e], [P[a, 0], P[a, 1]], [P[c, 0], P[c, 1]]) <= U*alpha[a, b, c, d, e, f])                        
                
                # if determinant(barriers[b][c], barriers[d][0], barriers[d][1]) != 0:
                #     MODEL.addConstr(alpha[a, b, c, d, e, f] >= determinant(barriers[b][c], barriers[d][0], barriers[d][1])/abs(determinant(barriers[b][c], barriers[d][0], barriers[d][1])))
       
        # if e == 1 and f == 1:
        #     MODEL.addConstr((1 - alpha[a, b, c, d, e, f])*L <= determinant(barriers[d][0], [P[a, 0], P[a, 1]], barriers[b][c]))
        #     MODEL.addConstr(determinant(barriers[b][c], barriers[d][0], barriers[d][1]) <= M*alpha[a, b, c, d, e, f])
    
    MODEL.addConstrs(beta0[a, c, d, e] == 2 * gamma[a, 0, c, d, e, 0] - alpha00[a, e] - alpha10[c, d, e] + 1 for a, c, d, e in beta0.keys())
    MODEL.addConstrs(beta1[a, c, d, e] == 2 * gamma[a, 0, c, d, e, 1] - alpha_1[a, c, d, e, 0] - alpha_1[a, c, d, e, 1] + 1 for a, c, d, e in beta1.keys())
    
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
        MODEL.addConstr((beta0[a, c, d, e] + beta1[a, c, d, e] <= 2 * delta[a, 0, c, d, e]))
        MODEL.addConstr(delta[a, 0, c, d, e] <= 2 * (beta0[a, c, d, e] + beta1[a, c, d, e]))
        
    # epsilon-C and y-C
    for a, b, c, d in epsilon.keys():
        # print((a, b, c, d))
        MODEL.addConstr((delta.sum(a, b, c, d, '*') - len(barriers)) + 1 <= epsilon[a, b, c, d])
        MODEL.addConstr(len(barriers) * epsilon[a, b, c, d] <= delta.sum(a, b, c, d, '*'))
        
        # if a < 0 and b >= 0 and c >= 0:
        MODEL.addConstr(y[a, b, c, d] + y[c, d, a, b] <= 2 * epsilon[a, b, c, d])
        
        # if a < 0 and b < 0 and c == 0:
        #     MODEL.addConstr(y[a, b, c] + y[b, a, c] <= 2*epsilon[a, b, c])
        
        # if no_ve(N[abs(a)-1], barriers[b][0]) and no_ve(N[abs(a)-1], barriers[b][1]):
        #     MODEL.addConstr(epsilon[a, b, c] <= 0)
    
    # NS and NT constraints
    for a, b, dim in P.keys():
        neighborhood = N[abs(a) - 1]
        
        if type(neighborhood) is Circle:
            MODEL.addConstr(dif_inside[a, b, dim] >= P[a, b, dim] - N[abs(a) - 1].center[dim])
            MODEL.addConstr(dif_inside[a, b, dim] >= N[abs(a) - 1].center[dim] - P[a, b, dim])
        
            MODEL.addConstr(gp.quicksum(dif_inside[a, b, dim] * dif_inside[a, b, dim] for dim in range(2)) <= d_inside[a, b] * d_inside[a, b])
            MODEL.addConstr(d_inside[a, b] <= N[abs(a) - 1].radii)
        
        if type(neighborhood) is Poligonal:
            MODEL.addConstrs(P[a, b, dim] == landa[a] * neighborhood.V[0][dim] + (1 - landa[a]) * neighborhood.V[1][dim] for dim in range(2))
    
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
            MODEL.addConstr(dist[a, b, c, d] == np.linalg.norm(np.array(barriers[a][b]) - np.array(barriers[c][d])))
            
        if (a, b, c, d) in EN:
            if a < 0:
                MODEL.addConstr(dif[a, b, c, d, dim] >= P[a, b, dim] - barriers[c][d][dim])
                MODEL.addConstr(dif[a, b, c, d, dim] >= -P[a, b, dim] + barriers[c][d][dim])
                MODEL.addConstr(gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] * dist[a, b, c, d])
            if c < 0:
                MODEL.addConstr(dif[c, d, a, b, dim] >= P[c, d, dim] - barriers[a][b][dim])
                MODEL.addConstr(dif[c, d, a, b, dim] >= -P[c, d, dim] + barriers[a][b][dim])
                MODEL.addConstr(gp.quicksum(dif[c, d, a, b, dim] * dif[c, d, a, b, dim] for dim in range(2)) <= dist[c, d, a, b] * dist[c, d, a, b])

        if (a, b, c, d) in ENN:
            MODEL.addConstr(dif[a, b, c, d, dim] >= P[a, b, dim] - P[c, d, dim])
            MODEL.addConstr(dif[a, b, c, d, dim] >= -P[a, b, dim] + P[c, d, dim])
            MODEL.addConstr(gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] * dist[a, b, c, d])
    
    L_out = 0
    U_out = 10000
    
    def estima_L(neighborhood, punto):

        if type(punto) is Circle:
            centro1 = neighborhood.center
            centro2 = punto.center
            
            dist = np.linalg.norm(np.array(centro1) - np.array(centro2)) - neighborhood.radii - punto.radii
        
        else:
            if type(neighborhood) is Circle:
                centro = neighborhood.center
                dist = np.linalg.norm(np.array(centro) - np.array(punto)) - neighborhood.radii
                
                # dist = 0
            
            if type(neighborhood) is Poligonal:
                dr = np.array(neighborhood.V[1]) - np.array(neighborhood.V[0])
                
                A, B = -dr[1], dr[0]
                
                C = -A * neighborhood.V[0][0] - B * neighborhood.V[0][1]
                
                dist = abs(A * punto[0] + B * punto[1] + C) / np.sqrt(A ** 2 + B ** 2)
        
        return dist
    
    def estima_U(neighborhood, punto):
        
        if type(punto) is Circle:
            centro1 = neighborhood.center
            centro2 = punto.center
            
            dist = np.linalg.norm(np.array(centro1) - np.array(centro2)) + neighborhood.radii + punto.radii
        
        else:
            if type(neighborhood) is Circle:
                
                dist = np.linalg.norm(np.array(neighborhood.center) - np.array(punto)) + neighborhood.radii
                
                # dist = 10000
                
            if type(neighborhood) is Poligonal:
                
                dist = max([np.linalg.norm(np.array(neighborhood.V[i]) - np.array(punto)) for i in range(2)])
            
        return dist
    
    # p constraints
    for a, b, c, d in p.keys():
        if (a, b, c, d) in EN:
            
            if a < 0:
                neighborhood = N[abs(a) - 1]
                punto = barriers[c][d]
                L_out = estima_L(neighborhood, punto)
                U_out = estima_U(neighborhood, punto)
        
            if c < 0:
                neighborhood = N[abs(c) - 1]
                punto = barriers[a][b]
                L_out = estima_L(neighborhood, punto)
                U_out = estima_U(neighborhood, punto)
        
        else:
            neighborhood1 = N[abs(a) - 1]
            neighborhood2 = N[abs(c) - 1]
            L_out = estima_L(neighborhood1, neighborhood2)
            U_out = estima_U(neighborhood1, neighborhood2)
            
        MODEL.addConstr(p[a, b, c, d] >= L_out * y[a, b, c, d])
        MODEL.addConstr(p[a, b, c, d] >= dist[a, b, c, d] - U_out * (1 - y[a, b, c, d]))
    
    # MODEL.addConstrs(z[v, v] == 0 for v in V
    
    # Restriccion 1
    MODEL.addConstrs(gp.quicksum(y[v, i, vn, j] for v, i in V if (v, i, vn, j) in EN + ENN) >= 1 for vn, j in VN)
    
    # Restriccion 2
    for v, i in V:
        MODEL.addConstr(gp.quicksum(y[v, i, vn, j] for vn, j in V if (v, i, vn, j) in E) == gp.quicksum(y[vn, j, v, i] for vn, j in V if (vn, j, v, i) in E))
    
    # Restriccion 3
    for vn, i in VN:
        if vn <= -2:
            MODEL.addConstr(gp.quicksum(g[vn, i, v, j] for v, j in V if (vn, i, v, j) in EN + ENN) - gp.quicksum(g[v, j, vn, i] for v, j in V if (v, j, vn, i) in EN + ENN) == 1)
    
    # Restriccion 4
    for vb, i in VB:
        MODEL.addConstr(gp.quicksum(g[(w, j, vb, i)] for w, j in V if (w, j, vb, i) in E) - gp.quicksum(g[(vb, i, w, j)] for w, j in V if (vb, i, w, j) in E) == 0)
  
    # Restriccion 5
    MODEL.addConstrs(g[a, b, c, d] <= (len(N) - 1) * y[a, b, c, d] for a, b, c, d in g.keys())
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
    
    objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(dist[index] * y[index] for index in EB)
    MODEL.setObjective(objective, GRB.MINIMIZE)
    
    second_time = time.time()
    
    time_elapsed = second_time - first_time
    
    MODEL.update()
    
    MODEL.Params.Threads = 6
    MODEL.Params.timeLimit = timeLimit - time_elapsed
    MODEL.Params.NonConvex = 2
    
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
    
    g_indices = []
    
    for index in g_index:
        if g[index].X > 0.5:
            g_indices.append(g[index])
    
    if log: 
        print(g_indices)

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
            ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c='red')
        
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
            ax.scatter(i[0], i[1], s=10, c='black')
        
        # print(points)
        
        segmentos = []
        
        for a, b, c, d in y_indices:
            if (a, b, c, d) in EN:
                if a < 0:
                    segmentos.append([points[abs(a) - 1][0], barriers[c][d][0], points[abs(a) - 1][1], barriers[c][d][1]])
                if c < 0:
                    segmentos.append([barriers[a][b][0], points[abs(c) - 1][0], barriers[a][b][1], points[abs(c) - 1][1]])
            if (a, b, c, d) in EB:
                segmentos.append([barriers[a][b][0], barriers[c][d][0], barriers[a][b][1], barriers[c][d][1]])
                
            if (a, b, c, d) in ENN:
                    segmentos.append([points[abs(a) - 1][0], points[abs(c) - 1][0], points[abs(a) - 1][1], points[abs(c) - 1][1]])
            
        # print(segmentos)
        for segmento in segmentos:
            ax.arrow(segmento[0], segmento[2], segmento[1] - segmento[0], segmento[3] - segmento[2], width=0.1, head_width=1, length_includes_head=True, color='black')
        
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
