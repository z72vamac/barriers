# Incluimos primero los paquetes

""" Aqui estamos resolviendo el caso para solo dos segments, shortest-path entre dos segments con barreras"""

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
from auxiliar_functions import *
import networkx as nx


# # semilla 1 hay punto de corte
# np.random.seed(5)
#
# # Lectura de datos
# segments = np.genfromtxt('segments.csv', delimiter = ',')
#
# # Numero total de segments que hay en los datos
# nT = len(segments)
#
# # Numero de segments que queremos visitar
# nS = 8
# segments_visitar = []
#
# # Numero de barreras que queremos tener
# nB = 5
# barreras = []
#
#
# todos_los_segments = segments[np.random.choice(nT, size = nS + nB, replace = False)]
#
#
# for i in range(nS):
#     row = todos_los_segments[i, :]
#     poligono = e.Poligono(V = np.array([[row[0], row[1]], [row[2], row[3]]]), col = 'blue')
#     segments_visitar.append(poligono)
#
# for i in range(nS, nB+nS):
#     row = todos_los_segments[i, :]
#     poligono = e.Poligono(V = np.array([[row[0], row[1]], [row[2], row[3]], [row[0], row[1]]]), col = 'red')
#     barreras.append(poligono)

# barreras = []

# barrera1 = Poligono(V = np.array([[7, 11], [9, 8]]))
# barreras.append(barrera1)
#
# barrera2 = Poligono(V = np.array([[12, 7], [11, 2]]))
# barreras.append(barrera2)
#
# barrera3 = Poligono(V = np.array([[9, 6], [8, 2]]))
# barreras.append(barrera3)
#
# barrera4 = Poligono(V = np.array([[11, 10], [16, 4]]))
# barreras.append(barrera4)
# datos = Data([], m = 30, r = 2, capacity = 100, modo = 2, tmax = 3000, init = False, show = False)
# datos.generar_muestra()
#
# barreras = nBostrar_datos()
# print(barreras)

def HSPPS(segments, barreras):
    """
    Input: source: neighborhood from the drone departs.
           target: neighborhood where the drone finishes the tour.
           
    Output: tour: list of segments that the drone visits to make the complete tour.
    """
    
    nB = len(barreras)
    
    nS = len(segments)
    
    distancias = np.zeros((nS, nS))
    
    for segmento1, count1 in zip(segments, range(nS)):
        for segmento2, count2 in zip(segments, range(nS)):
            if count1 > count2:
                no_vistos_segm1 = []
                no_vistos_segm2 = []
                
                vistos_c_segm1 = []
                vistos_c_segm2 = []
                
                
                # for i in range(1, nB+1):
                #     barrera = barreras[i-1]
                #
                #     Ci = barrera.V[0]
                #     Cip = barrera.V[1]
                #
                #     A1 = segmento1.V[0]
                #     B1 = segmento1.V[1]
                #
                #     segm1 = Poligono(np.array([Ci, A1]))
                #     segm2 = Poligono(np.array([Ci, B1]))
                #
                #     comb1 = []
                #     comb2 = []
                #
                #     alls = []
                #
                #     for j in range(1, nB+1):
                #         if j != i:
                #             barreraj = barreras[j-1]
                #
                #             result = min_dist(segm1, barreraj)
                #             comb1.append(result[0] < 0.005)
                #             alls.append(result[0] > 0.005)
                #
                #             result = min_dist(segm2, barreraj)
                #             comb2.append(result[0] < 0.005)
                #             alls.append(result[0] > 0.005)
                #
                #     if any(comb1) and any(comb2):
                #         no_vistos_segm1.append(i*100)
                #
                #     segm1 = Poligono(np.array([Cip, A1]))
                #     segm2 = Poligono(np.array([Cip, B1]))
                #
                #     comb1 = []
                #     comb2 = []
                #
                #     for j in range(1, nB+1):
                #         if j != i:
                #             barreraj = barreras[j-1]
                #
                #             result = min_dist(segm1, barreraj)
                #             comb1.append(result[0] < 0.005)
                #             alls.append(result[0] > 0.005)
                #
                #             result = min_dist(segm2, barreraj)
                #             comb2.append(result[0] < 0.005)
                #             alls.append(result[0] > 0.005)
                #
                #     if any(comb1) and any(comb2):
                #         no_vistos_segm1.append(i*100+1)
                #
                #     if all(alls):
                #         vistos_c_segm1.append(i)
                #
                #     A2 = segmento2.V[0]
                #     B2 = segmento2.V[1]
                #
                #     segm1 = Poligono(np.array([Ci, A2]))
                #     segm2 = Poligono(np.array([Ci, B2]))
                #
                #     comb1 = []
                #     comb2 = []
                #
                #     alls = []
                #
                #     for j in range(1, nB+1):
                #         if j != i:
                #             barreraj = barreras[j-1]
                #
                #             result = min_dist(segm1, barreraj)
                #             comb1.append(result[0] < 0.005)
                #             alls.append(result[0] > 0.005)
                #
                #             result = min_dist(segm2, barreraj)
                #             comb2.append(result[0] < 0.005)
                #             alls.append(result[0] > 0.005)
                #
                #     if any(comb1) and any(comb2):
                #         no_vistos_segm2.append(i*100)
                #
                #     segm1 = Poligono(np.array([Cip, A2]))
                #     segm2 = Poligono(np.array([Cip, B2]))
                #
                #     comb1 = []
                #     comb2 = []
                #
                #     for j in range(1, nB+1):
                #         if j != i:
                #             barreraj = barreras[j-1]
                #
                #             result = min_dist(segm1, barreraj)
                #             comb1.append(result[0] < 0.005)
                #             alls.append(result[0] > 0.005)
                #
                #             result = min_dist(segm2, barreraj)
                #             comb2.append(result[0] < 0.005)
                #             alls.append(result[0] > 0.005)
                #
                #     if any(comb1) and any(comb2):
                #         no_vistos_segm2.append(i*100+1)
                #
                #     if all(alls):
                #         vistos_c_segm2.append(i)
                
                def cortan(A, B, C, D):
                    det1 = (A[0] - C[0])*(B[1] - C[1]) - (B[0] - C[0])*(A[1] - C[1])
                    det2 = (A[0] - D[0])*(B[1] - D[1]) - (B[0] - D[0])*(A[1] - D[1])
                
                    det3 = (C[0] - A[0])*(D[1] - A[1]) - (D[0] - A[0])*(C[1] - A[1])
                    det4 = (C[0] - B[0])*(D[1] - B[1]) - (D[0] - B[0])*(C[1] - B[1])
                
                    return ((det1*det2)/abs(det1*det2) <= 0) and ((det3*det4)/abs(det3*det4) <= 0)
                
                for i in range(1, nB+1):
                    barrera = barreras[i-1]
                
                    Ci = barrera.V[0]
                    Cip = barrera.V[1]
                
                    A1 = segmento1.V[0]
                    B1 = segmento1.V[1]
                
                    comb1 = []
                    comb2 = []
                
                    alls = []
                
                    for j in range(1, nB+1):
                        if j != i:
                            barreraj = barreras[j-1]
                
                            Cj = barreraj.V[0]
                            Cjp = barreraj.V[1]
                
                            comb1.append(cortan(A1, Ci, Cj, Cjp))
                            alls.append(not(cortan(A1, Ci, Cj, Cjp)))
                
                            comb2.append(cortan(B1, Ci, Cj, Cjp))
                            alls.append(not(cortan(B1, Ci, Cj, Cjp)))
                
                    if any(comb1) and any(comb2):
                        no_vistos_segm1.append(i*100)
                
                    comb1 = []
                    comb2 = []
                
                    for j in range(1, nB+1):
                        if j != i:
                            barreraj = barreras[j-1]
                
                            Cj = barreraj.V[0]
                            Cjp = barreraj.V[1]
                
                            comb1.append(cortan(A1, Cip, Cj, Cjp))
                            alls.append(not(cortan(A1, Cip, Cj, Cjp)))
                
                            comb2.append(cortan(B1, Cip, Cj, Cjp))
                            alls.append(not(cortan(B1, Cip, Cj, Cjp)))
                
                    if any(comb1) and any(comb2):
                        no_vistos_segm1.append(i*100+1)
                
                    if all(alls):
                        vistos_c_segm1.append(i)
                
                    A1 = segmento2.V[0]
                    B1 = segmento2.V[1]
                
                    comb1 = []
                    comb2 = []
                
                    alls = []
                
                    for j in range(1, nB+1):
                        if j != i:
                            barreraj = barreras[j-1]
                
                            Cj = barreraj.V[0]
                            Cjp = barreraj.V[1]
                
                            comb1.append(cortan(A1, Ci, Cj, Cjp))
                            alls.append(not(cortan(A1, Ci, Cj, Cjp)))
                
                            comb2.append(cortan(B1, Ci, Cj, Cjp))
                            alls.append(not(cortan(B1, Ci, Cj, Cjp)))
                
                    if any(comb1) and any(comb2):
                        no_vistos_segm2.append(i*100)
                
                    comb1 = []
                    comb2 = []
                
                    for j in range(1, nB+1):
                        if j != i:
                            barreraj = barreras[j-1]
                
                            Cj = barreraj.V[0]
                            Cjp = barreraj.V[1]
                
                            comb1.append(cortan(A1, Cip, Cj, Cjp))
                            alls.append(not(cortan(A1, Cip, Cj, Cjp)))
                
                            comb2.append(cortan(B1, Cip, Cj, Cjp))
                            alls.append(not(cortan(B1, Cip, Cj, Cjp)))
                
                    if any(comb1) and any(comb2):
                        no_vistos_segm2.append(i*100+1)
                
                    if all(alls):
                        vistos_c_segm2.append(i)
                
                vertices = []
                
                # vertices.append(0)
                # vertices.append(5)
                
                for c in range(1, nB+1):
                    for v in range(barreras[c-1].num_segmentos):
                        vertices.append(c*100+v)
                
                vistos_segm1 = [i for i in vertices if i not in no_vistos_segm1]
                vistos_segm2 = [i for i in vertices if i not in no_vistos_segm2]
                
                barreras_vistas1 = list(set([i // 100 for i in vistos_segm1]))
                barreras_vistas2 = list(set([i // 100 for i in vistos_segm2]))
                
                
                # vistos_c_segm1 = [i for i in range(1, nB+1) if i * 100 in vistos_segm1 and i * 100 + 1 in vistos_segm1]
                # vistos_c_segm2 = [i for i in range(1, nB+1) if i * 100 in vistos_segm2 and i * 100 + 1 in vistos_segm2]
                
                no_vistos_c_segm1 = [i for i in vistos_segm1 if i // 100 not in vistos_c_segm1]
                no_vistos_c_segm2 = [i for i in vistos_segm2 if i // 100 not in vistos_c_segm2]
                
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
                
                for v in vistos_segm1:
                        edges.append((0, v))
                
                for v in vistos_segm2:
                        edges.append((v, 5))
                
                
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
                
                edges_total = []
                MODEL = gp.Model('Shortest-Path with Barriers')
                
                edges.append((0, 5))
                
                
                for i, j in edges:
                    edges_total.append((i, j))
                    edges_total.append((j, i))
                
                
                edges_tuple = gp.tuplelist(edges_total)
                edges_tuple
                
                print('Vertices no vistos por s1: ' + str(no_vistos_segm1))
                print('Vertices vistos por s1: ' + str(vistos_segm1))
                print('Barreras vistas por s1: ' + str(barreras_vistas1))
                print('Barreras vistas completamente por s1: ' + str(vistos_c_segm1))
                print('Vertices no vistos completamente por s1: ' + str(no_vistos_c_segm1))
                
                print()
                
                print('Vertices no vistos por s2: ' + str(no_vistos_segm2))
                print('Vertices vistos por s2: ' + str(vistos_segm2))
                print('Barreras vistas por s2: ' + str(barreras_vistas2))
                print('Barreras vistas completamente por s2: ' + str(vistos_c_segm2))
                print('Vertices no vistos completamente por s2: ' + str(no_vistos_c_segm2))
                
                print()
                
                x_index = edges
                
                zLij_index = []
                zLijp_index = []
                zLjjp_index = []
                zLijjp_index = []
                
                
                for j in barreras_vistas1:
                    for i in vistos_segm1:#no_vistos_c_segm1:
                        if i // 100 != j:
                            zLij_index.append((i, j))
                            zLijp_index.append((i, j))
                            zLjjp_index.append((i, j))
                            zLijjp_index.append((i, j))
                
                zRij_index = []
                zRijp_index = []
                zRjjp_index = []
                zRijjp_index = []
                
                for j in barreras_vistas2:
                    for i in vistos_segm2:#no_vistos_c_segm1:
                        if i // 100 != j:
                            zRij_index.append((i, j))
                            zRijp_index.append((i, j))
                            zRjjp_index.append((i, j))
                            zRijjp_index.append((i, j))
                
                dL_index = []
                
                for i in vistos_segm1:#vertices:
                    # if i >= 100:
                        dL_index.append(i)
                
                dR_index = []
                
                for i in vistos_segm2:#vertices:
                    # if i >= 100:
                        dR_index.append(i)
                
                
                xL = MODEL.addVars(2, vtype = GRB.CONTINUOUS, name = 'xL')
                muL = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muL')
                
                zLij = MODEL.addVars(zLij_index, vtype = GRB.BINARY, name = 'zLij')
                zLijp = MODEL.addVars(zLijp_index, vtype = GRB.BINARY, name = 'zLijp')
                zLjjp = MODEL.addVars(zLjjp_index, vtype = GRB.BINARY, name = 'zLjjp')
                zLijjp = MODEL.addVars(zLijjp_index, vtype = GRB.BINARY, name = 'zLijjp')
                deltaL = MODEL.addVars(zLij_index, 2, vtype = GRB.BINARY, name = 'deltaL')
                gammaL = MODEL.addVars(zLij_index, vtype = GRB.BINARY, name ='gammaL')
                alphaL = MODEL.addVars(vistos_segm1, vtype = GRB.BINARY, name = 'alphaL')
                
                dL = MODEL.addVars(dL_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dL')
                difL = MODEL.addVars(dL_index, 2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'difL')
                
                xR = MODEL.addVars(2, vtype = GRB.CONTINUOUS, name = 'xR')
                muR = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muR')
                
                zRij = MODEL.addVars(zRij_index, vtype = GRB.BINARY, name = 'zRij')
                zRijp = MODEL.addVars(zRijp_index, vtype = GRB.BINARY, name = 'zRijp')
                zRjjp = MODEL.addVars(zRjjp_index, vtype = GRB.BINARY, name = 'zRjjp')
                zRijjp = MODEL.addVars(zRijjp_index, vtype = GRB.BINARY, name = 'zRijjp')
                deltaR = MODEL.addVars(zRij_index, 2, vtype = GRB.BINARY, name = 'deltaR')
                gammaR = MODEL.addVars(zRij_index, vtype = GRB.BINARY, name ='gammaR')
                alphaR = MODEL.addVars(vistos_segm2, vtype = GRB.BINARY, name = 'alphaR')
                
                dR = MODEL.addVars(dR_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dR')
                difR = MODEL.addVars(dR_index, 2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'difR')
                
                y = MODEL.addVars(edges_tuple, vtype = GRB.BINARY, name = 'y')
                pL = MODEL.addVars(edges_tuple, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'pL')
                pR = MODEL.addVars(edges_tuple, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'pR')
                
                dLR = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dLR')
                difLR = MODEL.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'difLR')
                pLR = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'pLR')
                zLRij = MODEL.addVars(range(1, nB+1), vtype = GRB.BINARY, name = 'zLRij')
                zLRijp = MODEL.addVars(range(1, nB+1), vtype = GRB.BINARY, name = 'zLRijp')
                zLRjjp = MODEL.addVars(range(1, nB+1), vtype = GRB.BINARY, name = 'zLRjjp')
                zLRijjp = MODEL.addVars(range(1, nB+1), vtype = GRB.BINARY, name = 'zLRijjp')
                deltaLR = MODEL.addVars(range(1, nB+1), 2, vtype = GRB.BINARY, name = 'deltaLR')
                gammaLR = MODEL.addVars(range(1, nB+1), vtype = GRB.BINARY, name = 'gammaLR')
                alphaLR = MODEL.addVar(vtype = GRB.BINARY, name = 'alphaLR')
                
                
                orig = 0
                dest = 5
                
                MODEL.update()
                
                # MODEL.addConstr(muR == 1)
                # for v in vertices:
                #     MODEL.addConstr(gp.quicksum(y[v, w] for v, w in edges_tuple.select(v, '*')) - gp.quicksum(y[w, v] for w, v in edges_tuple.select('*', v)) == 0, 'node%s' % v)
                
                for v in vertices:
                    MODEL.addConstr(gp.quicksum(y[v, w] for v, w in edges_tuple.select(v, '*')) - gp.quicksum(y[w, v] for w, v in edges_tuple.select('*', v)) == (1 if v == orig else -1 if v == dest else 0), 'node%s' % v)
                
                MODEL.addConstr(gp.quicksum( y[v, w] for v, w in edges_tuple.select(0, '*')) - gp.quicksum( y[w, v] for w, v in edges_tuple.select('*', 0)) ==  1)
                
                MODEL.addConstr(gp.quicksum( y[v, w] for v, w in edges_tuple.select(5, '*')) - gp.quicksum( y[w, v] for w, v in edges_tuple.select('*', 5)) ==  -1)
                
                
                
                for i in dL.keys():
                    MODEL.addConstrs(difL[i, dim] >=  xL[dim] - barreras[i // 100 - 1].V[i % 100][dim] for dim in range(2))
                    MODEL.addConstrs(difL[i, dim] >= -xL[dim] + barreras[i // 100 - 1].V[i % 100][dim] for dim in range(2))
                    MODEL.addConstr(difL[i, 0]*difL[i, 0] + difL[i, 1]*difL[i, 1] <= dL[i]*dL[i])
                
                for i in dR.keys():
                    MODEL.addConstrs(difR[i, dim] >=  xR[dim] - barreras[i // 100 - 1].V[i % 100][dim] for dim in range(2))
                    MODEL.addConstrs(difR[i, dim] >= -xR[dim] + barreras[i // 100 - 1].V[i % 100][dim] for dim in range(2))
                    MODEL.addConstr(difR[i, 0]*difR[i, 0] + difR[i, 1]*difR[i, 1] <= dR[i]*dR[i])
                
                MODEL.addConstrs(difLR[dim] >=  xL[dim] - xR[dim] for dim in range(2))
                MODEL.addConstrs(difLR[dim] >= -xL[dim] + xR[dim] for dim in range(2))
                MODEL.addConstr(difLR[0]*difLR[0] + difLR[1]*difLR[1] <= dLR*dLR)
                
                
                BigM = 1e4
                SmallM = -1e4
                
                for i, j in zLij.keys():
                    Mj = (barreras[j - 1].V[0] + barreras[j - 1].V[1])/2
                
                    Ci = barreras[i // 100 - 1].V[i % 100]
                
                    Cj = barreras[j - 1].V[0]
                    Cjp = barreras[j - 1].V[1]
                
                    detLij = (xL[0] - Cj[0])*(Ci[1] - Cj[1]) - (Ci[0] - Cj[0])*(xL[1] - Cj[1])
                
                    MODEL.addConstr(detLij >= SmallM*zLij[i, j])
                    MODEL.addConstr(detLij <= BigM*(1 - zLij[i, j]))
                    # MODEL.addConstr(zLij[i, j] == 1 - bLij[i, j])
                
                    detLijp = (xL[0] - Cjp[0])*(Ci[1] - Cjp[1]) - (Ci[0] - Cjp[0])*(xL[1] - Cjp[1])
                
                    MODEL.addConstr(detLijp >= SmallM*zLijp[i, j])
                    MODEL.addConstr(detLijp <= BigM*(1 - zLijp[i, j]))
                
                    MODEL.addConstr(deltaL[i, j, 0] <=   zLijp[i, j] + zLij[i, j])
                    MODEL.addConstr(deltaL[i, j, 0] >=   zLijp[i, j] - zLij[i, j])
                    MODEL.addConstr(deltaL[i, j, 0] >= - zLijp[i, j] + zLij[i, j])
                    MODEL.addConstr(deltaL[i, j, 0] <= 2 - zLijp[i, j] - zLij[i, j])
                
                
                    # MODEL.addConstr( 1 - abs(zLijp[i, j] + zLij[i, j]) + 1 <= gammaL[i, j])
                
                    # MODEL.addConstr(zLijp[i, j] == 1 - bLijp[i, j])
                
                    detLjjp = (Cj[0] - xL[0])*(Cjp[1] - xL[1]) - (Cjp[0] - xL[0])*(Cj[1] - xL[1])
                
                    MODEL.addConstr(detLjjp >= SmallM*zLjjp[i, j])
                    MODEL.addConstr(detLjjp <= BigM*(1 - zLjjp[i, j]))
                    # MODEL.addConstr(zLij[i, j] == 1 - bLij[i, j])
                
                    detLijjp = (Cj[0] - Ci[0])*(Cjp[1] - Ci[1]) - (Cjp[0] - Ci[0])*(Cj[1] - Ci[1])
                
                    MODEL.addConstr(- BigM*zLijjp[i, j] <= detLijjp)
                    MODEL.addConstr( BigM*(1 - zLijjp[i, j]) >= detLijjp)
                
                    MODEL.addConstr(deltaL[i, j, 1] <=   zLjjp[i, j] + zLijjp[i, j])
                    MODEL.addConstr(deltaL[i, j, 1] >=   zLjjp[i, j] - zLijjp[i, j])
                    MODEL.addConstr(deltaL[i, j, 1] >= - zLjjp[i, j] + zLijjp[i, j])
                    MODEL.addConstr(deltaL[i, j, 1] <= 2 - zLjjp[i, j] - zLijjp[i, j])
                
                    MODEL.addConstr(gammaL[i, j] <= deltaL[i, j, 0])
                    MODEL.addConstr(gammaL[i, j] <= deltaL[i, j, 1])
                    MODEL.addConstr(gammaL[i, j] >= deltaL[i, j, 0] + deltaL[i, j, 1] - 1)
                
                
                
                    # MODEL.addConstr(zLij[i, j] == 1 - bLij[i, j])
                
                    # MODEL.addConstr(deltaL[i, j] >= zetaL[i, j])
                    # MODEL.addConstr(zetaL[i, j] >= etaL[i, j])
                
                    # MODEL.addConstr(zetaL[i, j] <= deltaL[i, j])
                
                # Si sum(gammaL) >= 1 ==> yLi == 0
                # sum(gammaL) <= nB*(1 - zLi) \/ yLi <= zLi
                
                for i in vistos_segm1:
                    MODEL.addConstr(gp.quicksum(gammaL[u, v] for u, v in gammaL.keys() if u == i) <= nB*(1- alphaL[i]))
                    MODEL.addConstr(y[0, i] <= alphaL[i])
                    # if i >= 100:
                        # MODEL.addConstr(1 - y[0, i] >= gp.quicksum(gammaL[u, v] for u, v in gammaL.keys() if u == i))
                
                for i, j in zRij.keys():
                    Mj = (barreras[j - 1].V[0] + barreras[j - 1].V[1])/2
                
                    Ci = barreras[i // 100 - 1].V[i % 100]
                
                    Cj = barreras[j - 1].V[0]
                    Cjp = barreras[j - 1].V[1]
                
                    detRij = (xR[0] - Cj[0])*(Ci[1] - Cj[1]) - (Ci[0] - Cj[0])*(xR[1] - Cj[1])
                
                    MODEL.addConstr(detRij >= SmallM*zRij[i, j])
                    MODEL.addConstr(detRij <= BigM*(1 - zRij[i, j]))
                    # MODEL.addConstr(zRij[i, j] == 1 - bRij[i, j])
                
                    detRijp = (xR[0] - Cjp[0])*(Ci[1] - Cjp[1]) - (Ci[0] - Cjp[0])*(xR[1] - Cjp[1])
                
                    MODEL.addConstr(detRijp >= SmallM*zRijp[i, j])
                    MODEL.addConstr(detRijp <= BigM*(1 - zRijp[i, j]))
                
                    MODEL.addConstr(deltaR[i, j, 0] <=   zRijp[i, j] + zRij[i, j])
                    MODEL.addConstr(deltaR[i, j, 0] >=   zRijp[i, j] - zRij[i, j])
                    MODEL.addConstr(deltaR[i, j, 0] >= - zRijp[i, j] + zRij[i, j])
                    MODEL.addConstr(deltaR[i, j, 0] <= 2 - zRijp[i, j] - zRij[i, j])
                
                    detRjjp = (Cj[0] - xR[0])*(Cjp[1] - xR[1]) - (Cjp[0] - xR[0])*(Cj[1] - xR[1])
                
                    MODEL.addConstr(detRjjp >= SmallM*zRjjp[i, j])
                    MODEL.addConstr(detRjjp <= BigM*(1 - zRjjp[i, j]))
                
                    detRijjp = (Cj[0] - Ci[0])*(Cjp[1] - Ci[1]) - (Cjp[0] - Ci[0])*(Cj[1] - Ci[1])
                
                    MODEL.addConstr(- BigM*zRijjp[i, j] <= detRijjp)
                    MODEL.addConstr( BigM*(1 - zRijjp[i, j]) >= detRijjp)
                
                    MODEL.addConstr(deltaR[i, j, 1] <=   zRjjp[i, j] + zRijjp[i, j])
                    MODEL.addConstr(deltaR[i, j, 1] >=   zRjjp[i, j] - zRijjp[i, j])
                    MODEL.addConstr(deltaR[i, j, 1] >= - zRjjp[i, j] + zRijjp[i, j])
                    MODEL.addConstr(deltaR[i, j, 1] <= 2 - zRjjp[i, j] - zRijjp[i, j])
                
                    MODEL.addConstr(gammaR[i, j] <= deltaR[i, j, 0])
                    MODEL.addConstr(gammaR[i, j] <= deltaR[i, j, 1])
                    MODEL.addConstr(gammaR[i, j] >= deltaR[i, j, 0] + deltaR[i, j, 1] - 1)
                
                for i in vistos_segm2:
                    MODEL.addConstr(gp.quicksum(gammaR[u, v] for u, v in gammaR.keys() if u == i) <= nB*(1- alphaR[i]))
                    MODEL.addConstr(y[i, 5] <= alphaR[i])
                
                
                for j in gammaLR.keys():
                    Cj = barreras[j - 1].V[0]
                    Cjp = barreras[j - 1].V[1]
                
                    detLRij = (xL[0] - Cj[0])*(xR[1] - Cj[1]) - (xR[0] - Cj[0])*(xL[1] - Cj[1])
                
                    MODEL.addConstr(detLRij >= SmallM*zLRij[j])
                    MODEL.addConstr(detLRij <= BigM*(1 - zLRij[j]))
                    # MODEL.addConstr(zRij[i, j] == 1 - bRij[i, j])
                
                    detLRijp = (xL[0] - Cjp[0])*(xR[1] - Cjp[1]) - (xR[0] - Cjp[0])*(xL[1] - Cjp[1])
                
                    MODEL.addConstr(detLRijp >= SmallM*zLRijp[j])
                    MODEL.addConstr(detLRijp <= BigM*(1 - zLRijp[j]))
                
                    MODEL.addConstr(deltaLR[j, 0] <=   zLRijp[j] + zLRij[j])
                    MODEL.addConstr(deltaLR[j, 0] >=   zLRijp[j] - zLRij[j])
                    MODEL.addConstr(deltaLR[j, 0] >= - zLRijp[j] + zLRij[j])
                    MODEL.addConstr(deltaLR[j, 0] <= 2 - zLRijp[j] - zLRij[j])
                
                    detLRjjp = (Cj[0] - xL[0])*(Cjp[1] - xL[1]) - (Cjp[0] - xL[0])*(Cj[1] - xL[1])
                
                    MODEL.addConstr(detLRjjp >= SmallM*zLRjjp[j])
                    MODEL.addConstr(detLRjjp <= BigM*(1 - zLRjjp[j]))
                
                    detLRijjp = (Cj[0] - xR[0])*(Cjp[1] - xR[1]) - (Cjp[0] - xR[0])*(Cj[1] - xR[1])
                
                    MODEL.addConstr(- BigM*zLRijjp[j] <= detLRijjp)
                    MODEL.addConstr( BigM*(1 - zLRijjp[j]) >= detLRijjp)
                
                    MODEL.addConstr(deltaLR[j, 1] <=   zLRjjp[j] + zLRijjp[j])
                    MODEL.addConstr(deltaLR[j, 1] >=   zLRjjp[j] - zLRijjp[j])
                    MODEL.addConstr(deltaLR[j, 1] >= - zLRjjp[j] + zLRijjp[j])
                    MODEL.addConstr(deltaLR[j, 1] <= 2 - zLRjjp[j] - zLRijjp[j])
                
                    MODEL.addConstr(gammaLR[j] <= deltaLR[j, 0])
                    MODEL.addConstr(gammaLR[j] <= deltaLR[j, 1])
                    MODEL.addConstr(gammaLR[j] >= deltaLR[j, 0] + deltaLR[j, 1] - 1)
                
                
                MODEL.addConstr(gp.quicksum(gammaLR[j] for j in gammaLR.keys()) <= nB*(1- alphaLR))
                MODEL.addConstr(y[0, 5] <= alphaLR)# for i in vistos_segm2:
                #     # if i >= 100:
                #         MODEL.addConstr(y[i, 5] <= 1 - gp.quicksum(gammaR[u, v] for u, v in gammaR.keys() if u == i))
                
                
                # for i in vistos_segm2:#no_vistos_c_segm1:
                #     # if i != j*100:
                #         if not(i // 100 in vistos_c_segm2):
                #             MODEL.addConstr(gammaR.sum(i, '*') <= 1)
                
                # for i in no_vistos_segm2:
                #     # if i in no_vistos_segm1:
                #         MODEL.addConstr(y[i, 5] == 0)
                
                objective = 0
                
                BigM = 10000
                
                for v, w in edges_tuple:
                    if v >= 100 and w >= 100:
                        objective += np.linalg.norm(barreras[v//100 - 1].V[v % 100] - barreras[w // 100 - 1].V[w % 100])*y[v, w]
                
                    if v == 0 and w in vistos_segm1:
                        MODEL.addConstr(pL[v, w] >= dL[w] - BigM*(1 - y[v, w]))
                        objective += pL[v, w]
                
                    if v in vistos_segm1 and w == 0:
                        MODEL.addConstr(pL[v, w] >= dL[v] - BigM*(1 - y[v, w]))
                        objective += pL[v, w]
                
                    if v == 5 and w in vistos_segm2:
                        MODEL.addConstr(pR[v, w] >= dR[w] - BigM*(1 - y[v, w]))
                        objective += pR[v, w]
                
                    if v in vistos_segm2 and w == 5:
                        MODEL.addConstr(pR[v, w] >= dR[v] - BigM*(1 - y[v, w]))
                        objective += pR[v, w]
                
                    if v == 0 and w == 5:
                        MODEL.addConstr(pLR >= dLR - BigM*(1 - y[0, 5]))
                        objective += pLR
                
                MODEL.addConstrs(xL[dim] == muL*segmento1.V[0][dim] + (1-muL)*segmento1.V[1][dim] for dim in range(2))
                MODEL.addConstrs(xR[dim] == muR*segmento2.V[0][dim] + (1-muR)*segmento2.V[1][dim] for dim in range(2))
                
                
                MODEL.setObjective(objective, GRB.MINIMIZE)
                
                MODEL.Params.NonConvex = 2
                # MODEL.Params.Heuristics = 0
                
                MODEL.update()
                
                MODEL.write('modelo_nivel2.lp')
                
                MODEL.optimize()
                
                if MODEL.Status == 3:
                    MODEL.computeIIS()
                    MODEL.write('infactible.ilp')
                
                MODEL.write('solucion_nivel2.sol')
                
                distancias[count1, count2] = MODEL.ObjVal
                
                # if MODEL.Status == GRB.Status.OPTIMAL:
                #    print('The final solution is:')
                #    for v, w in edges_total:
                #        if(y[v, w].X > 0):
                #            print(v, w, y[v, w].X)
                #
                # fig, ax = plt.subplots()
                # plt.axis([0, 100, 0, 100])
                #
                #
                # for b in barreras:
                #     ax.add_artist(b.artist)
                
                # for b in barreras:
                #     for v, i in zip(b.V, range(b.num_puntos)):
                #         ax.annotate(str(i), xy = (v[0], v[1]))
                
                # red_edges = gp.tuplelist([(i,j) for i,j in edges_tuple if y[i,j].X > 0.5])  #+ [(j,i) for i,j in edges_tuple if x[i,j].X > 0]
                #
                #
                # path = [orig]
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
                #
                # print(path)
                #
                # path_P = [[xL[0].X, xL[1].X]]
                #
                # for v in path[1:-1]:
                #     path_P.append(barreras[v//100 - 1].V[v%100])
                #
                # path_P.append([xR[0].X, xR[1].X])
                #
                # for b in range(1, nB+1):
                #     ax.add_artist(barreras[b-1].artist)
                #     ax.annotate(str(b), xy = (barreras[b-1].V[0][0], barreras[b-1].V[0][1]))
                #
                # ax.add_artist(segmento1.artist)
                # ax.add_artist(segmento2.artist)
                #
                # ax.add_artist(Polygon(path_P, fill=False, closed = False, animated=False, linestyle='-', lw = 2, alpha=1, color='black'))
                #
                # nx.draw(grafo.G, grafo.pos, node_size=40, edge_color = 'grey')
                #
                # plt.savefig('barreras.png')
                
                # plt.show()
    

        # Callback - use lazy constraints to eliminate sub-tours
    def subtourelim(model, where):
        if where == GRB.Callback.MIPSOL:
            # make a list of edges selected in the solution
            vals = model.cbGetSolution(model._vars)
            selected = gp.tuplelist((i, j) for i, j in model._vars.keys()
                                    if vals[i, j] > 0.5)
            # find the shortest cycle in the selected edge list
            tour = subtour(selected)
            if len(tour) < nS:
                # add subtour elimination constr. for every pair of cities in tour
                model.cbLazy(gp.quicksum(model._vars[i, j]
                                         for i, j in combinations(tour, 2))
                             <= len(tour)-1)
    
    
    # Given a tuplelist of edges, find the shortest subtour
    
    def subtour(edges):
        unvisited = list(range(nS))
        cycle = range(nS+1)  # initial length has 1 more city
        while unvisited:  # true if list is non-empty
            thiscycle = []
            neighbors = unvisited
            while neighbors:
                current = neighbors[0]
                thiscycle.append(current)
                unvisited.remove(current)
                neighbors = [j for i, j in edges.select(current, '*')
                             if j in unvisited]
            if len(cycle) > len(thiscycle):
                cycle = thiscycle
        return cycle
    
    distancias = distancias + distancias.T
    m = gp.Model()
    
    # Create variables
    
    index = []
    dist = {}
    
    for i in range(nS):
        for j in range(nS):
            if i != j:
                index.append((i, j))
                dist[(i, j)] = distancias[i, j]
                
                
    vars = m.addVars(index, obj=dist, vtype=GRB.BINARY, name='e')
    for i, j in vars.keys():
        vars[j, i] = vars[i, j]  # edge in opposite direction
    
    # You could use Python looping constructs and m.addVar() to create
    # these decision variables instead.  The following would be equivalent
    # to the preceding m.addVars() call...
    #
    # vars = tupledict()
    # for i,j in dist.keys():
    #   vars[i,j] = m.addVar(obj=dist[i,j], vtype=GRB.BINARY,
    #                        name='e[%d,%d]'%(i,j))
    
    
    # Add degree-2 constraint
    
    m.addConstrs(vars.sum(i, '*') == 2 for i in range(nS))
    
    # Using Python looping constructs, the preceding would be...
    #
    # for i in range(n):
    #   m.addConstr(sum(vars[i,j] for j in range(n)) == 2)
    
    
    # Optimize model
    
    m._vars = vars
    m.Params.lazyConstraints = 1
    m.optimize(subtourelim)
    
    vals = m.getAttr('x', vars)
    selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)
    
    tour = subtour(selected)
    assert len(tour) == nS
    
    print('')
    print('Optimal tour: %s' % str(tour))
    print('Optimal cost: %g' % m.objVal)
    print('')
    
    return tour

# distancias = HSPPS(segments_visitar, barreras)


# [0, 3, 7, 1, 5, 2, 6, 4]
# [1, 5, 7, 3, 6, 2, 8, 4]
