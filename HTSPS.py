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
from entorno import *
import copy
import estimacion_M as eM
from auxiliar_functions import *
import networkx as nx
from HSPPS import HSPPS



def HTSPS(segmentos, barreras, datos):
    nS = len(segmentos)
    S_index = range(1, nS+1)
    S_indexprima = range(1, nS+2)
    
    nB = len(barreras)
    
    no_vistos_segm = []
    vistos_c_segm = []
    
    def cortan(A, B, C, D):
        det1 = (A[0] - C[0])*(B[1] - C[1]) - (B[0] - C[0])*(A[1] - C[1])
        det2 = (A[0] - D[0])*(B[1] - D[1]) - (B[0] - D[0])*(A[1] - D[1])
    
        det3 = (C[0] - A[0])*(D[1] - A[1]) - (D[0] - A[0])*(C[1] - A[1])
        det4 = (C[0] - B[0])*(D[1] - B[1]) - (D[0] - B[0])*(C[1] - B[1])
    
        return ((det1*det2)/abs(det1*det2) <= 0) and ((det3*det4)/abs(det3*det4) <= 0)
    
    
    for s in S_index:
        for i in range(1, nB+1):
            barrera = barreras[i-1]
    
            Ci = barrera.V[0]
            Cip = barrera.V[1]
    
            segmento = segmentos[s-1]
    
            A = segmento.V[0]
            B = segmento.V[1]
    
            comb1 = []
            comb2 = []
    
            alls = []
    
            for j in range(1, nB+1):
                if j != i:
                    barreraj = barreras[j-1]
    
                    Cj = barreraj.V[0]
                    Cjp = barreraj.V[1]
    
                    comb1.append(cortan(A, Ci, Cj, Cjp))
                    alls.append(not(cortan(A, Ci, Cj, Cjp)))
    
                    comb2.append(cortan(B, Ci, Cj, Cjp))
                    alls.append(not(cortan(B, Ci, Cj, Cjp)))
    
            if any(comb1) and any(comb2):
                no_vistos_segm.append((s, i*100))
    
            comb1 = []
            comb2 = []
    
            for j in range(1, nB+1):
                if j != i:
                    barreraj = barreras[j-1]
    
                    Cj = barreraj.V[0]
                    Cjp = barreraj.V[1]
    
                    comb1.append(cortan(A, Cip, Cj, Cjp))
                    alls.append(not(cortan(A, Cip, Cj, Cjp)))
    
                    comb2.append(cortan(B, Cip, Cj, Cjp))
                    alls.append(not(cortan(B, Cip, Cj, Cjp)))
    
            if any(comb1) and any(comb2):
                no_vistos_segm.append((s, i*100+1))
    
            if all(alls):
                vistos_c_segm.append((s, i))
    
    vertices = []
    
    # vertices.append(0)
    # vertices.append(5)
    
    for c in range(1, nB+1):
        for v in range(barreras[c-1].num_segmentos):
            vertices.append(c*100+v)
    
    vistos_segm = [(s, i) for i in vertices for s in S_index if (s, i) not in no_vistos_segm]
    
    barreras_vistas = list(set([(s, i // 100) for s, i in vistos_segm]))
    
    print(barreras_vistas)
    
    # vistos_c_segm1 = [i for i in range(1, nB+1) if i * 100 in vistos_segm1 and i * 100 + 1 in vistos_segm1]
    # vistos_c_segm2 = [i for i in range(1, nB+1) if i * 100 in vistos_segm2 and i * 100 + 1 in vistos_segm2]
    
    no_vistos_c_segm = [(s, i) for (s, i) in vistos_segm if (s, i // 100) not in vistos_c_segm]
    # no_vistos_c_segm2 = [i for i in vistos_segm2 if i // 100 not in vistos_c_segm2]
    
    # print(vertices)
    
    vertices_pos = []
    
    # vertices_pos.append(segmento1.V[0])
    # vertices_pos.append(segmento1.V[1])
    # vertices_pos.append(segmento2.V[0])
    # vertices_pos.append(segmento2.V[1])
    
    for barrera in barreras:
        vertices_pos.append(barrera.V[0])
        vertices_pos.append(barrera.V[1])
    
    # print(vertices_pos)
    
    edges = []
    
    for v in vertices:
        for w in vertices:
            if v <= w:
                barrera1 = v // 100
                v1 = v % 100
    
                barrera2 = w // 100
                v2 = w % 100
    
                if barrera1 >= 1 and barrera2 >= 1:
                    if barrera1 != barrera2:
    
                        # No cruce ninguna otra superficie
                        Ci = barreras[barrera1-1].V[v1]
                        Cip = barreras[barrera2-1].V[v2]
    
                        comb = []
                        for j in range(1, nB + 1):
                            if j != barrera1 and j != barrera2:
    
                                Cj = barreras[j-1].V[0]
                                Cjp = barreras[j-1].V[1]
    
                                comb.append(not(cortan(Cj, Cjp, Ci, Cip)))
    
                        if all(comb):
                            edges.append((v, w))
    
    
    Ar = np.zeros((len(vertices), len(vertices)))
    
    for i, v in zip(range(len(vertices)), vertices):
        for j, w in zip(range(len(vertices)), vertices):
            if (v, w) in edges:
                # print(i, j)
                Ar[i, j] = 1
    
    grafo = Grafo(vertices_pos, Ar, 1)
    
    for s, i in vistos_segm:
            edges.append((s, i))
    
    
    # print(edges)
    # fig, ax = plt.subplots()
    # plt.axis([0, 100, 0, 100])
    #
    #
    # for b in range(1, nB+1):
    #     ax.add_artist(barreras[b-1].artist)
    #     ax.annotate(str(b), xy = (barreras[b-1].V[0][0], barreras[b-1].V[0][1]))
    #
    # ax.add_artist(segmento1.artist)
    # ax.add_artist(segmento2.artist)
    #
    # plt.show()
    
    # orig = -2
    # dest = -1
    
    
    for i in S_index:
        for j in S_index:
            if i > j:
                edges.append((i, j))
    
    edges_total = []
    for i, j in edges:
        edges_total.append((i, j))
        edges_total.append((j, i))
    
    
    edges_tuple = gp.tuplelist(edges_total)
    print(edges_tuple)
    
    print('Vertices no vistos por s: ' + str(no_vistos_segm))
    print('Vertices vistos por s: ' + str(vistos_segm))
    print('Barreras vistas por s: ' + str(barreras_vistas))
    print('Barreras vistas completamente por s: ' + str(vistos_c_segm))
    print('Vertices no vistos completamente por s: ' + str(no_vistos_c_segm))
    
    print()
    
    # x_index = edges
    
    zSij_index = []
    zSijp_index = []
    zSjjp_index = []
    zSijjp_index = []
    
    
    for s1 in S_index:
        for u, j in barreras_vistas:
            if u == s1:
                for s2, i in vistos_segm:#no_vistos_c_segm1:
                    if i // 100 != j and s2 == s1:
                        zSij_index.append((s1, i, j))
                        zSijp_index.append((s1, i, j))
                        zSjjp_index.append((s1, i, j))
                        zSijjp_index.append((s1, i, j))
    
    # print(zSij_index)
    
    y_index = []
    
    for s in S_index:
        for i, j in edges_tuple:
            y_index.append((i, j, s))
    
    dS_index = []
    
    for s, i in vistos_segm:#vertices:
        # if i >= 100:
            dS_index.append((s, i))
    
    # print(dS_index)
    dSS_index = []
    
    for s1 in S_index:
        for s2 in S_index:
            if s1 != s2:
                dSS_index.append((s1, s2))
    
    
    for s in S_index:
        vertices.append(s)
    
    # fig, ax = plt.subplots()
    # plt.axis([0, 100, 0, 100])
    #
    # for b in range(1, nB+1):
    #     ax.add_artist(barreras[b-1].artist)
    #     ax.annotate(str(b), xy = (barreras[b-1].V[0][0], barreras[b-1].V[0][1]))
    #
    # for s in S_index:
    #     ax.add_artist(segmentos[s-1].artist)
    #     ax.annotate(str(s), xy = (segmentos[s-1].V[0][0], segmentos[s-1].V[0][1]))
    #
    # nx.draw(grafo.G, grafo.pos, node_size=40, edge_color = 'grey')
    #
    # plt.savefig('barreras_nivel3.png')
    #
    # plt.show()
    
    MODEL = gp.Model('TSP with Barriers')
    
    xS = MODEL.addVars(S_index, 2, vtype = GRB.CONTINUOUS, name = 'xS')
    betaS = MODEL.addVars(S_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'betaS')
    
    zSij = MODEL.addVars(zSij_index, vtype = GRB.BINARY, name = 'zSij')
    zSijp = MODEL.addVars(zSijp_index, vtype = GRB.BINARY, name = 'zSijp')
    zSjjp = MODEL.addVars(zSjjp_index, vtype = GRB.BINARY, name = 'zSjjp')
    zSijjp = MODEL.addVars(zSijjp_index, vtype = GRB.BINARY, name = 'zSijjp')
    deltaS = MODEL.addVars(zSij_index, 2, vtype = GRB.BINARY, name = 'deltaS')
    gammaS = MODEL.addVars(zSij_index, vtype = GRB.BINARY, name ='gammaS')
    alphaS = MODEL.addVars(vistos_segm, vtype = GRB.BINARY, name = 'alphaS')
    
    dS = MODEL.addVars(dS_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dS')
    difS = MODEL.addVars(dS_index, 2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'difS')
    pS = MODEL.addVars(dS_index, S_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'pS')
    
    y = MODEL.addVars(y_index, vtype = GRB.BINARY, name = 'y')
    
    u = MODEL.addVars(vertices, S_index, vtype = GRB.BINARY, name = 'u')
    v = MODEL.addVars(vertices, S_index, vtype = GRB.BINARY, name = 'v')
    mu = MODEL.addVars(vertices, S_index, vtype = GRB.BINARY, name = 'mu')
    sit = MODEL.addVars(vertices, S_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = len(vertices)-1, name = 'sit')
    
    
    dSS = MODEL.addVars(dSS_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dSS')
    difSS = MODEL.addVars(dSS_index, 2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'difSS')
    pSS = MODEL.addVars(dSS_index, S_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'pSS')
    zSSij = MODEL.addVars(dSS_index, range(1, nB+1), vtype = GRB.BINARY, name = 'zSSij')
    zSSijp = MODEL.addVars(dSS_index, range(1, nB+1), vtype = GRB.BINARY, name = 'zSSijp')
    zSSjjp = MODEL.addVars(dSS_index, range(1, nB+1), vtype = GRB.BINARY, name = 'zSSjjp')
    zSSijjp = MODEL.addVars(dSS_index, range(1, nB+1), vtype = GRB.BINARY, name = 'zSSijjp')
    deltaSS = MODEL.addVars(dSS_index, range(1, nB+1), 2, vtype = GRB.BINARY, name = 'deltaSS')
    gammaSS = MODEL.addVars(dSS_index, range(1, nB+1), vtype = GRB.BINARY, name = 'gammaSS')
    alphaSS = MODEL.addVars(dSS_index, vtype = GRB.BINARY, name = 'alphaSS')
    
    
    MODEL.update()
    
    # init = True
    if datos.init:
        tour = HSPPS(segmentos, barreras)
        
        for s, t in zip(S_index, S_index):
            u[tour[t-1]+1, t].start = 1
            if t < nS:
                v[tour[t]+1, t].start = 1
        
        v[tour[0]+1, nS].start = 1
        
        # for s, t in zip(S_index, S_index):
        #     MODEL.addConstr(u[tour[t-1]+1, t] >= .5)
        #     if t < nS:
        #         MODEL.addConstr(v[tour[t]+1, t] >= .5)
        #
        # MODEL.addConstr(v[tour[0]+1, nS] >= 0.5)    
        # [1, 5, 10, 7, 3, 6, 9, 13, 15, 2, 14, 12, 8, 4, 11]
        # [12, 4, 8, 11, 1, 5, 10, 7, 3, 6, 13, 9, 15, 2, 14]
        
    # MODEL.addConstr(muR == 1)
    for t in S_index:
        for i in vertices:
            if i >= 100:
                MODEL.addConstr(gp.quicksum(y[i, j, t] for i, j in edges_tuple.select(i, '*')) - gp.quicksum(y[j, i, t] for j, i in edges_tuple.select('*', i)) == 0, 'node%s' % i)
    
    MODEL.addConstrs(gp.quicksum(u[i, t] for i in S_index) == 1 for t in S_index)
    MODEL.addConstrs(gp.quicksum(u[i, t] for i in vertices if i  >= 100) == 0 for t in S_index)
    
    MODEL.addConstrs(u.sum(s, '*') == 1 for s in S_index)
    
    MODEL.addConstrs(gp.quicksum(v[i, t] for i in S_index) == 1 for t in S_index)
    MODEL.addConstrs(gp.quicksum(v[i, t] for i in vertices if i >= 100) == 0 for t in S_index)
    
    MODEL.addConstrs(v.sum(s, '*') == 1 for s in S_index)
    
    MODEL.addConstrs(v[s, t] == u[s, t+1] for s in S_index for t in S_index[:-1])
    MODEL.addConstrs(v[s, nS] == u[s, 1] for s in S_index)
    
    
    # MODEL.addConstrs(mu.sum(i, '*') == 1 for i in S_index)
    # MODEL.addConstrs(gp.quicksum(mu[i, t] for i in S_index) == 1 for t in S_index)
    
    
    MODEL.addConstrs(y.sum('*', i, t) + u[i, t] == mu[i, t] for i, t in mu.keys())
    MODEL.addConstrs(y.sum(i, '*', t) + v[i, t] == mu[i, t] for i, t in mu.keys())
    
    # MODEL.addConstr(y[3, 1, 3] == 1)
    # MODEL.addConstr(y[1, 201, 4] == 1)
    
    # MODEL.addConstr(y[4, 3, 3] == 1)
    # MODEL.addConstr(y[1, 3, 1] == 1)
    # MODEL.addConstr(y[3, 4, 2] == 1)
    # MODEL.addConstr(y[4, 201, 3] == 1)
    # MODEL.addConstr(y[4, 1, 2] == 1)
    # MODEL.addConstr(y[1, 3, 3] == 1)
    # MODEL.addConstr(y[3, 201, 4] == 1)
    # MODEL.addConstr(y[201, 2, 4] == 1)
    
    
    
    
    
    # MODEL.addConstrs(u.sum(s, '*') == 1 for s in S_index)
    # MODEL.addConstrs(u.sum('*', s) == 1 for s in S_index)
    
    def powerset(iterable):
        "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
        s = list(iterable)
        return chain.from_iterable(combinations(s, r) for r in range(2, len(s)+1))
    
    
    for t in S_index:
        for i, j in edges_tuple:
            MODEL.addConstr(sit[i, t] - sit[j, t] + len(vertices)*y[i, j, t] <= len(vertices)-1)
    # for set in powerset(list(edges_tuple)):
    #     print(set)
    
    # print(edges_tuple)
    # for t in S_index:
    #     restricciones = MODEL.addConstrs(gp.quicksum(y[i, j, t] for i, j in set) <= len(set) - 1 for set in powerset(list(edges_tuple)) if len(set) < nS)
    #     restricciones.Lazy = 3
    print(dS.keys())
    print(difS.keys())
    #
    for s, i in dS.keys():
        MODEL.addConstrs(difS[s, i, dim] >=  xS[s, dim] - barreras[i // 100 - 1].V[i % 100][dim] for dim in range(2))
        MODEL.addConstrs(difS[s, i, dim] >= -xS[s, dim] + barreras[i // 100 - 1].V[i % 100][dim] for dim in range(2))
        MODEL.addConstr(difS[s, i, 0]*difS[s, i, 0] + difS[s, i, 1]*difS[s, i, 1] <= dS[s, i]*dS[s, i])
    
    
    MODEL.addConstrs(difSS[s1, s2, dim] >=  xS[s1, dim] - xS[s2, dim] for dim in range(2) for s1, s2 in dSS.keys())
    MODEL.addConstrs(difSS[s1, s2, dim] >= -xS[s1, dim] + xS[s2, dim] for dim in range(2) for s1, s2 in dSS.keys())
    MODEL.addConstrs(difSS[s1, s2, 0]*difSS[s1, s2, 0] + difSS[s1, s2, 1]*difSS[s1, s2, 1] <= dSS[s1, s2]*dSS[s1, s2] for s1, s2 in dSS.keys())
    
    
    BigM = 1e5
    SmallM = -1e5
    
    for s, i, j in zSij.keys():
        Mj = (barreras[j - 1].V[0] + barreras[j - 1].V[1])/2
    
        Ci = barreras[i // 100 - 1].V[i % 100]
    
        Cj = barreras[j - 1].V[0]
        Cjp = barreras[j - 1].V[1]
    
        detSij = (xS[s, 0] - Cj[0])*(Ci[1] - Cj[1]) - (Ci[0] - Cj[0])*(xS[s, 1] - Cj[1])
    
        MODEL.addConstr(detSij >= SmallM*zSij[s, i, j])
        MODEL.addConstr(detSij <= BigM*(1 - zSij[s, i, j]))
        # MODEL.addConstr(zSij[s, i, j] == 1 - bLij[s, i, j])
    
        detSijp = (xS[s, 0] - Cjp[0])*(Ci[1] - Cjp[1]) - (Ci[0] - Cjp[0])*(xS[s, 1] - Cjp[1])
    
        MODEL.addConstr(detSijp >= SmallM*zSijp[s, i, j])
        MODEL.addConstr(detSijp <= BigM*(1 - zSijp[s, i, j]))
    
        MODEL.addConstr(deltaS[s, i, j, 0] <=   zSijp[s, i, j] + zSij[s, i, j])
        MODEL.addConstr(deltaS[s, i, j, 0] >=   zSijp[s, i, j] - zSij[s, i, j])
        MODEL.addConstr(deltaS[s, i, j, 0] >= - zSijp[s, i, j] + zSij[s, i, j])
        MODEL.addConstr(deltaS[s, i, j, 0] <= 2 - zSijp[s, i, j] - zSij[s, i, j])
    
    
        # MODEL.addConstr( 1 - abs(zSijp[s, i, j] + zSij[s, i, j]) + 1 <= gammaS[s, i, j])
    
        # MODEL.addConstr(zSijp[s, i, j] == 1 - bLijp[s, i, j])
    
        detSjjp = (Cj[0] - xS[s, 0])*(Cjp[1] - xS[s, 1]) - (Cjp[0] - xS[s, 0])*(Cj[1] - xS[s, 1])
    
        MODEL.addConstr(detSjjp >= SmallM*zSjjp[s, i, j])
        MODEL.addConstr(detSjjp <= BigM*(1 - zSjjp[s, i, j]))
        # MODEL.addConstr(zSij[s, i, j] == 1 - bLij[s, i, j])
    
        detSijjp = (Cj[0] - Ci[0])*(Cjp[1] - Ci[1]) - (Cjp[0] - Ci[0])*(Cj[1] - Ci[1])
    
        MODEL.addConstr(- BigM*zSijjp[s, i, j] <= detSijjp)
        MODEL.addConstr( BigM*(1 - zSijjp[s, i, j]) >= detSijjp)
    
        MODEL.addConstr(deltaS[s, i, j, 1] <=   zSjjp[s, i, j] + zSijjp[s, i, j])
        MODEL.addConstr(deltaS[s, i, j, 1] >=   zSjjp[s, i, j] - zSijjp[s, i, j])
        MODEL.addConstr(deltaS[s, i, j, 1] >= - zSjjp[s, i, j] + zSijjp[s, i, j])
        MODEL.addConstr(deltaS[s, i, j, 1] <= 2 - zSjjp[s, i, j] - zSijjp[s, i, j])
    
        MODEL.addConstr(gammaS[s, i, j] <= deltaS[s, i, j, 0])
        MODEL.addConstr(gammaS[s, i, j] <= deltaS[s, i, j, 1])
        MODEL.addConstr(gammaS[s, i, j] >= deltaS[s, i, j, 0] + deltaS[s, i, j, 1] - 1)
    
    
    
        # MODEL.addConstr(zSij[s, i, j] == 1 - bLij[s, i, j])
    
        # MODEL.addConstr(deltaS[s, i, j] >= zetaL[s, i, j])
        # MODEL.addConstr(zetaL[s, i, j] >= etaL[s, i, j])
    
        # MODEL.addConstr(zetaL[s, i, j] <= deltaS[s, i, j])
    
    # Si sum(gammaS) >= 1 ==> yLi == 0
    # sum(gammaS) <= datos.m*(1 - zSi) \/ yLi <= zSi
    
    for s, i in vistos_segm:
        MODEL.addConstr(gp.quicksum(gammaS.select(s, i, '*')) <= nB*(1- alphaS[s, i]))
        MODEL.addConstr(y.sum(s, i, '*') <= alphaS[s, i])
        MODEL.addConstr(y.sum(i, s, '*') <= alphaS[s, i])
    
    
    for s1, s2, j in gammaSS.keys():
        Cj = barreras[j - 1].V[0]
        Cjp = barreras[j - 1].V[1]
    
        detSSij = (xS[s1, 0] - Cj[0])*(xS[s2, 1] - Cj[1]) - (xS[s2, 0] - Cj[0])*(xS[s1, 1] - Cj[1])
    
        MODEL.addConstr(detSSij >= SmallM*zSSij[s1, s2, j])
        MODEL.addConstr(detSSij <= BigM*(1 - zSSij[s1, s2, j]))
        # MODEL.addConstr(zRij[s, i, j] == 1 - bRij[s, i, j])
    
        detSSijp = (xS[s1, 0] - Cjp[0])*(xS[s2, 1] - Cjp[1]) - (xS[s2, 0] - Cjp[0])*(xS[s1, 1] - Cjp[1])
    
        MODEL.addConstr(detSSijp >= SmallM*zSSijp[s1, s2, j])
        MODEL.addConstr(detSSijp <= BigM*(1 - zSSijp[s1, s2, j]))
    
        MODEL.addConstr(deltaSS[s1, s2, j, 0] <=   zSSijp[s1, s2, j] + zSSij[s1, s2, j])
        MODEL.addConstr(deltaSS[s1, s2, j, 0] >=   zSSijp[s1, s2, j] - zSSij[s1, s2, j])
        MODEL.addConstr(deltaSS[s1, s2, j, 0] >= - zSSijp[s1, s2, j] + zSSij[s1, s2, j])
        MODEL.addConstr(deltaSS[s1, s2, j, 0] <= 2 - zSSijp[s1, s2, j] - zSSij[s1, s2, j])
    
        detSSjjp = (Cj[0] - xS[s1, 0])*(Cjp[1] - xS[s1, 1]) - (Cjp[0] - xS[s1, 0])*(Cj[1] - xS[s1, 1])
    
        MODEL.addConstr(detSSjjp >= SmallM*zSSjjp[s1, s2, j])
        MODEL.addConstr(detSSjjp <= BigM*(1 - zSSjjp[s1, s2, j]))
    
        detSSijjp = (Cj[0] - xS[s2, 0])*(Cjp[1] - xS[s2, 1]) - (Cjp[0] - xS[s2, 0])*(Cj[1] - xS[s2, 1])
    
        MODEL.addConstr(- BigM*zSSijjp[s1, s2, j] <= detSSijjp)
        MODEL.addConstr( BigM*(1 - zSSijjp[s1, s2, j]) >= detSSijjp)
    
        MODEL.addConstr(deltaSS[s1, s2, j, 1] <=   zSSjjp[s1, s2, j] + zSSijjp[s1, s2, j])
        MODEL.addConstr(deltaSS[s1, s2, j, 1] >=   zSSjjp[s1, s2, j] - zSSijjp[s1, s2, j])
        MODEL.addConstr(deltaSS[s1, s2, j, 1] >= - zSSjjp[s1, s2, j] + zSSijjp[s1, s2, j])
        MODEL.addConstr(deltaSS[s1, s2, j, 1] <= 2 - zSSjjp[s1, s2, j] - zSSijjp[s1, s2, j])
    
        MODEL.addConstr(gammaSS[s1, s2, j] <= deltaSS[s1, s2, j, 0])
        MODEL.addConstr(gammaSS[s1, s2, j] <= deltaSS[s1, s2, j, 1])
        MODEL.addConstr(gammaSS[s1, s2, j] >= deltaSS[s1, s2, j, 0] + deltaSS[s1, s2, j, 1] - 1)
    
    
    MODEL.addConstrs(gammaSS.sum(s1, s2, '*') <= nB*(1- alphaSS[s1, s2]) for s1, s2 in alphaSS.keys())
    #
    MODEL.addConstrs(y.sum(s1, s2, '*') <= alphaSS[s1, s2] for s1, s2 in alphaSS.keys())# for i in vistos_segm2:
    #     # if i >= 100:
    #         MODEL.addConstr(y[i, 5] <= 1 - gp.quicksum(gammaR[u, v] for u, v in gammaR.keys() if u == i))
    
    
    # MODEL.addConstr(y[1, 3, 1] == 1)
    # MODEL.addConstr(y[3, 4, 2] == 1)
    # MODEL.addConstr(y[4, 400, 3] == 1)
    # MODEL.addConstr(y[400, 2, 3] == 1)
    # MODEL.addConstr(y[2, 400, 4] == 1)
    # MODEL.addConstr(y[400, 1, 4] == 1)
    
    # for i in vistos_segm2:#no_vistos_c_segm1:
    #     # if i != j*100:
    #         if not(i // 100 in vistos_c_segm2):
    #             MODEL.addConstr(gammaR.sum(i, '*') <= 1)
    
    # for i in no_vistos_segm2:
    #     # if i in no_vistos_segm1:
    #         MODEL.addConstr(y[i, 5] == 0)
    #
    objective = 0
    
    BigM = 10000
    
    MODEL.addConstrs(xS[s, dim] == betaS[s]*segmentos[s-1].V[0][dim] + (1-betaS[s])*segmentos[s-1].V[1][dim] for s, dim in xS.keys())
    
    # MODEL.addConstrs(pS[s, i, t] >= dS[s, i] - BigM*(1 - y[s, i, t]) for s, i, t in pS.keys())
    # MODEL.addConstrs(pS[s, i, t] >= dS[s, i] - BigM*(1 - y[i, s, t]) for s, i, t in pS.keys())
    
    # print(pS.keys())
    
    for s, i, t in pS.keys():
        
        SmallM = eM.estima_SmallM_local(segmentos[s-1], barreras[i // 100 -1])
        # BigM = eM.estima_BigM_local(segmentos[s-1], barreras[i // 100 -1])
        MODEL.addConstr(pS[s, i, t] >= SmallM*y[s, i, t])
        MODEL.addConstr(pS[s, i, t] >= SmallM*y[i, s, t])
        
        MODEL.addConstr(pS[s, i, t] >= dS[s, i] - BigM*(1 - y[s, i, t]))
        MODEL.addConstr(pS[s, i, t] >= dS[s, i] - BigM*(1 - y[i, s, t]))

    
    for s1, s2, t in y.keys(): 
        if s1 < 100 and s2 < 100:
            SmallM = eM.estima_SmallM_local(segmentos[s1-1], segmentos[s2-1])
            MODEL.addConstr(pSS[s1, s2, t] >= SmallM*y[s1, s2, t])
            MODEL.addConstr(pSS[s1, s2, t] >= dSS[s1, s2] - BigM*(1 - y[s1, s2, t]))
    
    
    objective = gp.quicksum(np.linalg.norm(barreras[i//100 - 1].V[i % 100] - barreras[j // 100 - 1].V[j % 100])*y[i, j, t] for i, j, t in y.keys() if i >= 100 and j >= 100) + gp.quicksum(pS[keys] for keys in pS.keys()) + gp.quicksum(pSS[keys] for keys in pSS.keys())
    
    # print(y.keys())
    # for i, j, t in y.keys():
    #     if i >= 100 and j >= 100:
    #         objective += np.linalg.norm(barreras[i//100 - 1].V[i % 100] - barreras[w // 100 - 1].V[j % 100])*y[i, j, t]
    #
    #     if i < 100 and j >= 100:
    #             MODEL.addConstr(pS[i, j] >= dS[i, j] - BigM*(1 - y[i, j, t]))
    #             MODEL.addConstr(pS[j, i] >= dS[j, i] - BigM*(1 - y[j, i, t]))
    #             objective += pS[i, j]/(2*nS)
    #             objective += pS[j, i]/(2*nS)
    #
    #
    #     # if i >= 100 and j < 100:
    #     #         MODEL.addConstr(pS[i, j] >= dS[i, j] - BigM*(1 - y[i, j, t]))
    #     #         objective += pS[i, j]
    #     #
    #     # if v >= 100 and w < 100:
    #     #     if (v, w) in y.ke
    #     #     MODEL.addConstr(pS[v, w] >= dS[v, w] - BigM*(1 - y[v, w]))
    #     #     objective += pL[v, w]
    #
    #     if i < 100 and j < 100:
    #         MODEL.addConstr(pSS[i, j] >= dSS[i, j] - BigM*(1 - y[i, j, t]))
    #         objective += pSS[i, j]/nS
    
    MODEL.setObjective(objective, GRB.MINIMIZE)
    
    MODEL.Params.NonConvex = 2
    MODEL.Params.TimeLimit = datos.tmax
    # MODEL.Params.LazyConstraints = 1
    # MODEL.Params.Heuristics = 0
    
    MODEL.update()
    
    MODEL.write('modelo_nivel3.lp')
    
    # MODEL.read('solucion_nivel3.sol')
    
    MODEL.optimize()
    
    if MODEL.Status == 3:
        MODEL.computeIIS()
        MODEL.write('infactible_nivel3.ilp')
    
    MODEL.write('solucion_nivel3.sol')
    
    
    # print(entro)
    # print(salgo)
    print('y')
    selected_y = gp.tuplelist((i, j, t) for i,j,t in y.keys() if y[i, j, t].X > 0.5)
    print(selected_y)
    
    print('u')
    selected_u = gp.tuplelist((i, t) for i,t in u.keys() if u[i, t].X > 0.5)
    print(selected_u)
    
    print('v')
    selected_v = gp.tuplelist((i, t) for i,t in v.keys() if v[i, t].X > 0.5)
    print(selected_v)
    
    
    # if MODEL.Status == GRB.Status.OPTIMAL:
    #     print('y:')
    #     for v, w in y.keys():
    #         if(y[v, w].X > 0):
    #             print(v, w)
    #
    #     print('u:')
    #     for v, w in u.keys():
    #         if(u[v, w].X > 0):
    #             print(v, w)
    
    
    fig, ax = plt.subplots()
    plt.axis([0, 100, 0, 100])
    
    path = []
    # path.append(0)
    
    for t in S_index:
        tripleta = selected_u.select('*', t)
        if tripleta:
            path.append(tripleta[0][0])
    
    # path.append(nG+1)
    print(path)
    
    ind = 0
    path_C = []
    path_D = []
    
    #path_C.append(orig)
    # path_C.append([xS[s, 0].X, xLt[0, 1].X])
    index_g = path[0]
    for t in S_index:
        s = path[t-1]
        #    if ind < datos.m:
        # path_C.append([xLt[t, 0].X, xLt[t, 1].X])
        # path_D = []
        path_D.append([xS[s, 0].X, xS[s, 1].X])
        # index_g = s
        # index_g = g
        # index_t = t
        # for g, ti in selected_ugt:
        #     if ti == t:
        #         index_g = g
        #         index_t = ti
        #
        count = 0
        # path_D.append([P[index_g-1, 0], P[index_g-1, 1]])
        # path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
        limite = len(selected_y.select('*', '*', t))
        print(selected_y.select('*', '*', t))
        while count < limite:
            for g1, g2, ti in selected_y.select('*', '*', t):
                # print(g1, g2, ti)
                # print(index_g)
                if index_g == g1:
                    count += 1
                    index_g = g2
                    # path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
                    if index_g >= 100:
                        path_D.append(barreras[index_g // 100 - 1].V[index_g % 100])
                    else:
                        path_D.append([xS[index_g, 0].X, xS[index_g, 1].X])
    
    
    
        # path_D.append([xRt[t, 0].X, xRt[t, 1].X])
        # paths_D.append(path_D)
        # path_C.append([xRt[t, 0].X, xRt[t, 1].X])
    
    # path_C.append([xLt[datos.m+1, 0].X, xLt[datos.m+1, 1].X])
    
    # for b in barreras:
    #     ax.add_artist(b.artist)
    
    # for b in barreras:
    #     for v, i in zip(b.V, range(b.num_puntos)):
    #         ax.annotate(str(i), xy = (v[0], v[1]))
    
    # red_edges = gp.tuplelist([(i,j) for i,j in y.keys() if y[i,j].X > 0.5])  #+ [(j,i) for i,j in edges_tuple if x[i,j].X >
    # print(red_edges)
    #
    # selected_u = gp.tuplelist([(s1, s2) for s1, s2 in u.keys() if u[s1, s2].X > 0.5])
    # print(selected_u)
    
    # tour = subtour(selected_u)
    # print(tour)
    #
    # index_v = orig
    #
    # ind = 1
    # while ind <= len(red_edges):
    #     for (v, w) in red_edges:
    #         if v == index_v:
    #             index_v = w
    #             ind += 1
    #             path.append(index_v)
    
    # print(path)
    
    # path_P = [[xS[0].X, xS[1].X]]
    #
    # for v in path[1:-1]:
    #     path_P.append(barreras[v//100 - 1].V[v%100])
    
    # path_P.append([xS[0].X, xS[1].X])
    
    for b in range(1, nB+1):
        ax.add_artist(barreras[b-1].artist)
        ax.annotate(str(b), xy = (barreras[b-1].V[0][0], barreras[b-1].V[0][1]))
    
    for s in S_index:
        ax.add_artist(segmentos[s-1].artist)
        ax.annotate(str(s), xy = (segmentos[s-1].V[0][0], segmentos[s-1].V[0][1]))
    
    for s in S_index:
        plt.plot(xS[s, 0].X, xS[s, 1].X, 'ko', markersize=5)
    
    ax.add_artist(Polygon(path_D, fill=False, closed = True, animated=False, linestyle='-', lw = 2, alpha=1, color='black'))
    
    nx.draw(grafo.G, grafo.pos, node_size=40, edge_color = 'grey')
    
    plt.savefig('barreras.png')
    
    plt.show()
    
    return None
