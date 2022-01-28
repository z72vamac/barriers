# Este documento es para definir funciones que se utilizan en varios ficheros.

import numpy as np
import gurobipy as gp
from gurobipy import GRB
from entorno import Elipse, Poligono, Poligonal, Circulo
import estimacion_M as eM
import copy
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines
import networkx as nx

def path2matrix(path):
    "Toma un camino y lo devuelve como matriz de adyacencia"
    m = len(path)
    zcc = np.zeros([m, m])
    for i in range(m-1):
        zcc[path[i]][path[i+1]]=1
    zcc[path[m-1]][path[0]] = 1
    return zcc

def matrix2path(matrix):
    " Toma una matriz y lo devuelve como camino "
    matrix = np.array(matrix, int)
    ind = 0
    path = []
    while ind not in path:
        path.append(ind)
        lista = matrix[ind]
        counter = 0
        for i in lista:
            if i == 1:
                ind = counter
                break
            counter += 1
    return path


def subtour(edges):
    "Genera un subtour de una lista de aristas"
    m = len(edges)
    unvisited = list(range(m))
    cycle = range(m+1)  # initial length has 1 more city
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

def subtour_cplex(edges):
    "Genera un subtour de una lista de aristas"
    m = len(edges)
    unvisited = list(range(m))
    cycle = range(m+1)  # initial length has 1 more city
    while unvisited:  # true if list is non-empty
        thiscycle = []
        neighbors = unvisited
        while neighbors:
            current = neighbors[0]
            thiscycle.append(current)
            unvisited.remove(current)
            edges_selected = [(i, j) for i, j in edges if i == current]
            neighbors = [j for i, j in edges_selected
                         if j in unvisited]
        if len(cycle) > len(thiscycle):
            cycle = thiscycle
    return cycle

def subtours(edges):
    "Genera un subtour de una lista de aristas"
    m = len(edges)
    unvisited = edges
    cycle = range(m+1)  # initial length has 1 more city
    while unvisited:  # true if list is non-empty
        thiscycle = []
        neighbors = unvisited
        while neighbors:
            current = neighbors[0]
            thiscycle.append(current)
            unvisited.remove(current)
            neighbors = [j for i, j in edges.select(current, '*')
                         if j in unvisited] + [i for i, j in edges.select('*', current) if i in unvisited]
        if len(cycle) > len(thiscycle):
            cycle = thiscycle
    return cycle

# def subtours(edges):
#     "Genera un subtour de una lista de aristas"
#     m = len(edges)
#     unvisited = list(range(m))
#     cycles = []  # initial length has 1 more city
#     while unvisited:  # true if list is non-empty
#         thiscycle = []
#         neighbors = unvisited
#         while neighbors:
#             current = neighbors[0]
#             thiscycle.append(current)
#             unvisited.remove(current)
#             neighbors = [j for i, j in edges.select(current, '*')
#                          if j in unvisited]
#         # if len(cycle) > len(thiscycle):
#         cycles.append(thiscycle)
#     return cycles
#
#
# def subtour_s(edges):
#     m = int(len(edges)/2)
#     unvisited = list(range(m))
#     cycle = list(range(m)) # Dummy - guaranteed to be replaced
#     while unvisited:  # true if list is non-empty
#         thiscycle = []
#         neighbors = unvisited
#         while neighbors:
#             current = neighbors[0]
#             thiscycle.append(current)
#             unvisited.remove(current)
#             neighbors = [j for i, j in edges.select(current, '*')
#                          if j in unvisited]
#         if len(thiscycle) <= len(cycle):
#             cycle = thiscycle # New shortest subtour
#     return cycle

def min_dist(comp0, comp1):

        MODEL = gp.Model('minima_distancia')

        x0 = MODEL.addVars(2, vtype = GRB.CONTINUOUS, name = 'x0')
        x1 = MODEL.addVars(2, vtype = GRB.CONTINUOUS, name = 'x1')

        if type(comp0) is Poligono:
            mu0 = MODEL.addVars(comp0.num_puntos, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mu0')

        if type(comp0) is Poligonal:
            landa_index = []
            sublanda_index = []
            s_index = []

            # landa de la variable de entrada en la poligonal c
            landa_index.append(0)
            # landa de la variable de salida en la poligonal c
            for segm in range(comp0.num_segmentos):
                s_index.append(segm)
            for punto in range(comp0.num_puntos):
                sublanda_index.append(punto)

            landa0 = MODEL.addVar(vtype=GRB.CONTINUOUS, name='landa')
            sublanda0 = MODEL.addVars(sublanda_index, vtype=GRB.CONTINUOUS,
                                 lb=0.0, ub=1.0, name='sublanda')
            s0 = MODEL.addVars(s_index, vtype=GRB.BINARY, name='s')

        if type(comp1) is Poligono:
            mu1 = MODEL.addVars(comp1.num_puntos, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mu1')

        if type(comp1) is Poligonal:
            landa_index = []
            sublanda_index = []
            s_index = []

            # landa de la variable de salida en la poligonal c
            for segm in range(comp1.num_segmentos):
                s_index.append(segm)
            for punto in range(comp1.num_puntos):
                sublanda_index.append(punto)

            landa1 = MODEL.addVar(vtype=GRB.CONTINUOUS, name='landa')
            sublanda1 = MODEL.addVars(sublanda_index, vtype=GRB.CONTINUOUS,
                                 lb=0.0, ub=1.0, name='sublanda')
            s1 = MODEL.addVars(s_index, vtype=GRB.BINARY, name='s')

        dif01 = MODEL.addVars(2, vtype = GRB.CONTINUOUS, name = 'dif01')
        d01 = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'd01')

        MODEL.update()

        MODEL.addConstr(dif01[0] >=  x0[0] - x1[0])
        MODEL.addConstr(dif01[0] >= -x0[0] + x1[0])
        MODEL.addConstr(dif01[1] >=  x0[1] - x1[1])
        MODEL.addConstr(dif01[1] >= -x0[1] + x1[1])
        MODEL.addConstr(dif01[0] * dif01[0] + dif01[1] * dif01[1] <= d01 * d01)
        # MODEL.addConstr(d01 <= estima_BigM_local(comp0, comp1))

        if type(comp0) is Poligono:
            MODEL.addConstr(gp.quicksum(mu0[j] for j in range(comp0.num_puntos)) == 1, name = 'envConv')
            for j in range(2):
                MODEL.addConstr(x0[j] == gp.quicksum(mu0[v]*comp0.V[v][j] for v in range(comp0.num_puntos)), name = 'inP1')
        if type(comp0) is Elipse:
            MODEL.addConstr(comp0.P[0, 0] * x0[0] * x0[0] + comp0.P[1, 0] * x0[0] * x0[1] +
                            comp0.P[0, 1] * x0[0] * x0[1] + comp0.P[1, 1] * x0[1] * x0[1] +
                            comp0.q[0] * x0[0] + comp0.q[1] * x0[1] + comp0.r <= 0, name='inC1')
        if type(comp0) is Poligonal:
            for i in range(2):
                for punto in range(1, comp0.num_puntos):
                    MODEL.addConstr(landa0 - punto >= sublanda0[punto] - comp0.num_puntos * (1 - s0[punto - 1]))
                    MODEL.addConstr(landa0 - punto <= sublanda0[punto] + comp0.num_puntos * (1 - s0[punto - 1]))
                MODEL.addConstr(sublanda0[0] <= s0[0])
                MODEL.addConstr(sublanda0[comp0.num_puntos - 1] <= s0[comp0.num_puntos - 2])
                for punto in range(1, comp0.num_puntos - 1):
                    MODEL.addConstr(sublanda0[punto] <= s0[punto - 1] + s0[punto])
                MODEL.addConstr(s0.sum('*') == 1)
                MODEL.addConstr(sublanda0.sum('*') == 1)
                for j in range(2):
                    MODEL.addConstr(x0[j] == gp.quicksum(sublanda0[punto] * comp0.V[punto][j] for punto in range(comp0.num_puntos)), name='seg1')
        if type(comp1) is Poligono:
            MODEL.addConstr(gp.quicksum(mu1[j] for j in range(comp1.num_puntos)) == 1, name = 'envConv')
            for j in range(2):
                MODEL.addConstr(x1[j] == gp.quicksum(mu1[v]*comp1.V[v][j] for v in range(comp1.num_puntos)), name = 'inP1')
        if type(comp1) is Elipse:
            MODEL.addConstr(comp1.P[0, 0] * x1[0] * x1[0] + comp1.P[1, 0] * x1[0] * x1[1] +
                            comp1.P[0, 1] * x1[0] * x1[1] + comp1.P[1, 1] * x1[1] * x1[1] +
                            comp1.q[0] * x1[0] + comp1.q[1] * x1[1] + comp1.r <= 0, name='inC1')
        if type(comp1) is Poligonal:
            for i in range(2):
                for punto in range(1, comp1.num_puntos):
                    MODEL.addConstr(landa1 - punto >= sublanda1[punto] - comp1.num_puntos * (1 - s1[punto - 1]))
                    MODEL.addConstr(landa1 - punto <= sublanda1[punto] + comp1.num_puntos * (1 - s1[punto - 1]))
                MODEL.addConstr(sublanda1[0] <= s1[0])
                MODEL.addConstr(sublanda1[comp1.num_puntos - 1] <= s1[comp1.num_puntos - 2])
                for punto in range(1, comp1.num_puntos - 1):
                    MODEL.addConstr(sublanda1[punto] <= s1[punto - 1] + s1[punto])
                MODEL.addConstr(s1.sum('*') == 1)
                MODEL.addConstr(sublanda1.sum('*') == 1)
                for j in range(2):
                    MODEL.addConstr(x1[j] == gp.quicksum(sublanda1[punto] * comp1.V[punto][j] for punto in range(comp1.num_puntos)), name='seg1')

        # MODEL.setParam('OutputFlag', 1)

        MODEL.setObjective(d01, GRB.MINIMIZE)
        # MODEL.Params.FeasibilityTol = 1e-2
        MODEL.Params.OutputFlag = 0

        MODEL.update()

        MODEL.optimize()

        x_0 = [x0[0].X, x0[1].X]
        x_1 = [x1[0].X, x1[1].X]

        return d01.X, x_0, x_1

def dist_grafo(point, graph):
    """Short summary.

    Parameters
    ----------
    point : list of two coordinates
        Point to calculate the distance.
    graph : entorno graph
        Graph from which we want to search the closest point to the origin.

    Returns
    -------
    list
        Closest point to punto

    """
    minimum = np.inf
    # x = []

    for segm in graph.aristas:
        first_V = graph.V[segm // 100 - 1]
        second_V = graph.V[segm % 100]

        M = gp.Model('min_dist')

        landa = M.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'landa')
        x_segm = M.addVars(2, vtype = GRB.CONTINUOUS, name = 'x_segm')
        dif = M.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dif')
        d = M.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'd')

        M.update()

        M.addConstrs(x_segm[dim] == landa*first_V[dim] + (1-landa)*second_V[dim] for dim in range(2))
        M.addConstrs(dif[dim] >=  x_segm[dim] - point[dim] for dim in range(2))
        M.addConstrs(dif[dim] >= -x_segm[dim] + point[dim] for dim in range(2))
        M.addConstr(dif[0]*dif[0] + dif[1]*dif[1] <= d*d)

        M.update()

        M.setObjective(d, GRB.MINIMIZE)

        M.setParam('OutputFlag', 0)
        M.update()

        M.optimize()

        if d.X < minimum:
            x = [x_segm[0].X, x_segm[1].X]
            minimum = d.X
    # else:
    #     for v in graph.V:
    #         if np.linalg.norm(point - v) <= minimum:
    #             x = v
    #             minimum = np.linalg.norm(point - v)

    return minimum, x

def dist_point_point_polygon(point1, polygon1, point2, polygon2):
    """Short summary.

    Parameters
    ----------
    point : list of two coordinates
        Point to calculate the distance.
    graph : entorno graph
        Graph from which we want to search the closest point to the origin.

    Returns
    -------
    list
        Closest point to punto

    """

    M = gp.Model('min_dist')

    landa1 = M.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'landa1')
    landa2 = M.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'landa2')
    x1 = M.addVars(2, vtype = GRB.CONTINUOUS, name = 'x1')
    x2 = M.addVars(2, vtype = GRB.CONTINUOUS, name = 'x2')
    dif = M.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dif')
    d = M.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'd')
    mu1 = M.addVars(polygon1.num_puntos, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mu1')
    mu2 = M.addVars(polygon2.num_puntos, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mu2')

    M.update()

    M.addConstrs(x1[dim] == landa1*point1[dim] + (1-landa1)*point2[dim] for dim in range(2))
    M.addConstrs(x2[dim] == landa2*point1[dim] + (1-landa2)*point2[dim] for dim in range(2))

    M.addConstr(gp.quicksum(mu1[v] for v in range(polygon1.num_puntos)) == 1, name='envConv1')
    M.addConstr(gp.quicksum(mu2[v] for v in range(polygon2.num_puntos)) == 1, name='envConv2')

    for dim in range(2):
        M.addConstr(x1[dim] == gp.quicksum(mu1[v] * polygon1.V[v][dim] for v in range(polygon1.num_puntos)))
        M.addConstr(x2[dim] == gp.quicksum(mu2[v] * polygon2.V[v][dim] for v in range(polygon2.num_puntos)))

    M.addConstrs(dif[dim] >=  x1[dim] - x2[dim] for dim in range(2))
    M.addConstrs(dif[dim] >= -x1[dim] + x2[dim] for dim in range(2))
    M.addConstr(dif[0]*dif[0] + dif[1]*dif[1] <= d*d)

    M.update()

    M.setObjective(d, GRB.MINIMIZE)

    M.setParam('OutputFlag', 0)
    M.update()

    M.optimize()
    # else:
    #     for v in graph.V:
    #         if np.linalg.norm(point - v) <= minimum:
    #             x = v
    #             minimum = np.linalg.norm(point - v)

    return d.X

def dist_point_polygon(point1, point2, polygon2):
    """Short summary.

    Parameters
    ----------
    point : list of two coordinates
        Point to calculate the distance.
    graph : entorno graph
        Graph from which we want to search the closest point to the origin.

    Returns
    -------
    list
        Closest point to punto

    """

    M = gp.Model('min_dist')

    landa1 = M.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'landa1')
    landa2 = M.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'landa2')
    x1 = M.addVars(2, vtype = GRB.CONTINUOUS, name = 'x1')
    x2 = M.addVars(2, vtype = GRB.CONTINUOUS, name = 'x2')
    dif = M.addVars(2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dif')
    d = M.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'd')
    mu1 = M.addVars(polygon1.num_puntos, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mu1')
    mu2 = M.addVars(polygon2.num_puntos, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mu2')

    M.update()

    M.addConstrs(x1[dim] == landa1*point1[dim] + (1-landa1)*point2[dim] for dim in range(2))
    M.addConstrs(x2[dim] == landa2*point1[dim] + (1-landa2)*point2[dim] for dim in range(2))

    M.addConstr(gp.quicksum(mu1[v] for v in range(polygon1.num_puntos)) == 1, name='envConv1')
    M.addConstr(gp.quicksum(mu2[v] for v in range(polygon2.num_puntos)) == 1, name='envConv2')

    for dim in range(2):
        M.addConstr(x1[dim] == gp.quicksum(mu1[v] * polygon1.V[v][dim] for v in range(polygon1.num_puntos)))
        M.addConstr(x2[dim] == gp.quicksum(mu2[v] * polygon2.V[v][dim] for v in range(polygon2.num_puntos)))

    M.addConstrs(dif[dim] >=  x1[dim] - x2[dim] for dim in range(2))
    M.addConstrs(dif[dim] >= -x1[dim] + x2[dim] for dim in range(2))
    M.addConstr(dif[0]*dif[0] + dif[1]*dif[1] <= d*d)

    M.update()

    M.setObjective(d, GRB.MINIMIZE)

    M.setParam('OutputFlag', 0)
    M.update()

    M.optimize()
    # else:
    #     for v in graph.V:
    #         if np.linalg.norm(point - v) <= minimum:
    #             x = v
    #             minimum = np.linalg.norm(point - v)

    return d.X

def dist_point_segment(point, segment):
    
    MODEL = gp.Model('Dist_Point_Segment')
    
    landa = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0)
    P_segm = MODEL.addVars(2, vtype = GRB.CONTINUOUS)
    dif = MODEL.addVars(2, vtype = GRB.CONTINUOUS)
    d = MODEL.addVar(vtype = GRB.CONTINUOUS, lb = 0.0)
    
    MODEL.update()
    
    MODEL.addConstrs(P_segm[dim] == landa*segment.V[0][dim] + (1-landa)*segment.V[1][dim] for dim in range(2))
    
    MODEL.addConstrs(dif[dim] >=  P_segm[dim] - point[dim] for dim in range(2))
    MODEL.addConstrs(dif[dim] >= -P_segm[dim] + point[dim] for dim in range(2))
    
    MODEL.addConstr(dif[0]*dif[0] + dif[1]*dif[1] <= d*d)
    
    MODEL.setObjective(d, GRB.MINIMIZE)
    
    MODEL.update()
    
    MODEL.Params.OutputFlag = 0
    
    MODEL.optimize()
    
    return d.X


def determinant(P, Q, R):
    a11 = Q[0]-P[0]
    a12 = R[0]-P[0]
    a21 = Q[1]-P[1]
    a22 = R[1]-P[1]
    
    return Q[0]*R[1] - Q[0]*P[1] - P[0]*R[1] - R[0]*Q[1] + R[0]*P[1] + P[0]*Q[1]  

def estima_det(entorno, barrera):
        
    PB1x = barrera[0][0]
    PB2x = barrera[1][0]
    PB1y = barrera[0][1]
    PB2y = barrera[1][1]
    
    if type(entorno) is Circulo:
        
        # print((entorno.center, entorno.radii, barrera))
        centro = entorno.center
        radio = entorno.radii
        
        a = 1
        b = 1
        c = 0
        d = -2*centro[0]
        e = -2*centro[1]
        f = centro[0]**2 + centro[1]**2 - radio**2
        
        # print((a, b, c, d, e, f))
        # Ecuacion seria entonces (x - centro[0])^2 + (y- centro[1])^2 <= radio^2
        # x^2 + centro[0]^2 - 2x·centro[0] + y^2 + centro[1]^2 - 2y·centro[1] - radio^2 <= 0
        # a = 1; b = 1; c = 0; d = -2·centro[0]; e = -2·centro[1]; f = centro[0]^2 + centro[1]^2 - radio^2
        
        if abs(PB1y - PB2y) >= 1e-1:
            
            m = (2*a*(PB1x - PB2x) - c*(PB2y - PB1y))/(2*b*(PB2y - PB1y) - c*(PB1x - PB2x))
            n = (d*(PB1x - PB2x) - e*(PB2y - PB1y))/(2*b*(PB2y - PB1y) - c*(PB1x - PB2x))
            
            
            # print((m, n))
            
            x_mas = (-(2*b*m*n+c*n+d+e*m)+np.sqrt((2*b*m*n + c*n + d + e*m)**2 - 4*(a + b*m*m + c*m)*(n*n*b + e*n + f)))/(2*(a + b*m*m + c*m))
            x_menos = (-(2*b*m*n+c*n+d+e*m)-np.sqrt((2*b*m*n + c*n + d + e*m)**2 - 4*(a + b*m*m + c*m)*(n*n*b + e*n + f)))/(2*(a + b*m*m + c*m))
            
            # print((x_mas, x_menos))


            y_mas = m*x_mas + n
            y_menos = m*x_menos + n

            # print((x_mas, y_mas))
            # print((x_menos, y_menos))
                            
            eval_mas = x_mas*PB1y - x_mas*PB2y + y_mas*PB2x - y_mas*PB1x + PB1x*PB2y - PB1y*PB2x
            
            eval_mas = determinant([x_mas, y_mas], [PB1x, PB1y], [PB2x, PB2y])
            eval_menos = determinant([x_menos, y_menos], [PB1x, PB1y], [PB2x, PB2y])
            
            # eval_menos = x_menos*PB1y - x_menos*PB2y + y_menos*PB2x - y_menos*PB1x + PB1x*PB2y - PB1y*PB2x
            
            L = min(eval_mas, eval_menos)
            # U = max(eval_mas, eval_menos)
            U = max(eval_mas, eval_menos)
        
        else:
            x = centro[0]
            y_mas = centro[1] + radio
            y_menos = centro[1] - radio
            
            eval_mas = x*PB1y - x*PB2y + y_mas*PB2x - y_mas*PB1x + PB1x*PB2y - PB1y*PB2x
            eval_menos = x*PB1y - x*PB2y + y_menos*PB2x - y_menos*PB1x + PB1x*PB2y - PB1y*PB2x
            
            L = min(eval_mas, eval_menos)
            U = max(eval_mas, eval_menos)                
        
        # print(L)
        # print(U)
    if type(entorno) is Poligonal:
        x_mas = entorno.V[0][0]
        y_mas = entorno.V[0][1]
        
        x_menos = entorno.V[1][0]
        y_menos = entorno.V[1][1]

        eval_mas = determinant([x_mas, y_mas], [PB1x, PB1y], [PB2x, PB2y])
        eval_menos = determinant([x_menos, y_menos], [PB1x, PB1y], [PB2x, PB2y])
        
        L = min(eval_mas, eval_menos)
        U = max(eval_mas, eval_menos)              
        
    # L = -20000
    # U = 20000
    
    return L, U
    
def estima_L(entorno, punto):

    if type(punto) is Circulo:
        centro1 = entorno.center
        centro2 = punto.center
        
        dis = np.linalg.norm(np.array(centro1) - np.array(centro2)) - entorno.radii - punto.radii
    
    else:
        if type(entorno) is Circulo:
            centro = entorno.center
            dis = np.linalg.norm(np.array(centro) - np.array(punto)) - entorno.radii
            
            # dist = 0
        
        if type(entorno) is Poligonal:
            
            dis = dist_point_segment(punto, entorno)
            # dr = np.array(entorno.V[1]) - np.array(entorno.V[0])
            #
            # A, B = -dr[1], dr[0]
            #
            # C = -A*entorno.V[0][0] - B*entorno.V[0][1]
            #
            # dis = abs(A*punto[0] + B*punto[1] + C)/np.sqrt(A**2 + B**2)
            
    
    return dis

def estima_U(entorno, punto):
    
    if type(punto) is Circulo:
        centro1 = entorno.center
        centro2 = punto.center
        
        dis = np.linalg.norm(np.array(centro1) - np.array(centro2)) + entorno.radii + punto.radii
    
    else:
        if type(entorno) is Circulo:
            
            dis = np.linalg.norm(np.array(entorno.center) - np.array(punto)) + entorno.radii
            
            # dist = 10000
            
        if type(entorno) is Poligonal:
            
            dis = max([np.linalg.norm(np.array(entorno.V[i]) - np.array(punto)) for i in range(2)])
        
    return dis
