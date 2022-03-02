# Este documento es para definir funciones que se utilizan en varios ficheros.

import numpy as np
import gurobipy as gp
from gurobipy import GRB
import neighborhood as neigh
import estimacion_M as eM
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines
import networkx as nx
import itertools
import sympy
# from sympy.abc import x, y

def path2matrix(path):
    "Toma un camino y lo devuelve como matriz de adyacencia"
    m = len(path)
    zcc = np.zeros([m, m])
    for i in range(m - 1):
        zcc[path[i]][path[i + 1]] = 1
    zcc[path[m - 1]][path[0]] = 1
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
    cycle = range(m + 1)  # initial length has 1 more city
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
    cycle = range(m + 1)  # initial length has 1 more city
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
    cycle = range(m + 1)  # initial length has 1 more city
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

    x0 = MODEL.addVars(2, vtype=GRB.CONTINUOUS, name='x0')
    x1 = MODEL.addVars(2, vtype=GRB.CONTINUOUS, name='x1')

    if type(comp0) is neigh.Poligono:
        mu0 = MODEL.addVars(comp0.num_puntos, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='mu0')

    if type(comp0) is neigh.Poligonal:
        landa_index = []
        sublanda_index = []
        s_index = []

        # landa de la variable de entrada en la neigh.Poligonal c
        landa_index.append(0)
        # landa de la variable de salida en la neigh.Poligonal c
        for segm in range(comp0.num_segmentos):
            s_index.append(segm)
        for punto in range(comp0.num_puntos):
            sublanda_index.append(punto)

        landa0 = MODEL.addVar(vtype=GRB.CONTINUOUS, name='landa')
        sublanda0 = MODEL.addVars(sublanda_index, vtype=GRB.CONTINUOUS,
                                  lb=0.0, ub=1.0, name='sublanda')
        s0 = MODEL.addVars(s_index, vtype=GRB.BINARY, name='s')

    if type(comp1) is neigh.Poligono:
        mu1 = MODEL.addVars(comp1.num_puntos, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='mu1')

    if type(comp1) is neigh.Poligonal:
        landa_index = []
        sublanda_index = []
        s_index = []

        # landa de la variable de salida en la neigh.Poligonal c
        for segm in range(comp1.num_segmentos):
            s_index.append(segm)
        for punto in range(comp1.num_puntos):
            sublanda_index.append(punto)

        landa1 = MODEL.addVar(vtype=GRB.CONTINUOUS, name='landa')
        sublanda1 = MODEL.addVars(sublanda_index, vtype=GRB.CONTINUOUS,
                                  lb=0.0, ub=1.0, name='sublanda')
        s1 = MODEL.addVars(s_index, vtype=GRB.BINARY, name='s')

    dif01 = MODEL.addVars(2, vtype=GRB.CONTINUOUS, name='dif01')
    d01 = MODEL.addVar(vtype=GRB.CONTINUOUS, lb=0.0, name='d01')

    MODEL.update()

    MODEL.addConstr(dif01[0] >= x0[0] - x1[0])
    MODEL.addConstr(dif01[0] >= -x0[0] + x1[0])
    MODEL.addConstr(dif01[1] >= x0[1] - x1[1])
    MODEL.addConstr(dif01[1] >= -x0[1] + x1[1])
    MODEL.addConstr(dif01[0] * dif01[0] + dif01[1] * dif01[1] <= d01 * d01)
    # MODEL.addConstr(d01 <= estima_BigM_local(comp0, comp1))

    if type(comp0) is neigh.Poligono:
        MODEL.addConstr(gp.quicksum(mu0[j] for j in range(comp0.num_puntos)) == 1, name='envConv')
        for j in range(2):
            MODEL.addConstr(x0[j] == gp.quicksum(mu0[v] * comp0.V[v][j] for v in range(comp0.num_puntos)), name='inP1')
    if type(comp0) is Elipse:
        MODEL.addConstr(comp0.P[0, 0] * x0[0] * x0[0] + comp0.P[1, 0] * x0[0] * x0[1] +
                        comp0.P[0, 1] * x0[0] * x0[1] + comp0.P[1, 1] * x0[1] * x0[1] +
                        comp0.q[0] * x0[0] + comp0.q[1] * x0[1] + comp0.r <= 0, name='inC1')
    if type(comp0) is neigh.Poligonal:
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
                MODEL.addConstr(
                    x0[j] == gp.quicksum(sublanda0[punto] * comp0.V[punto][j] for punto in range(comp0.num_puntos)),
                    name='seg1')
    if type(comp1) is neigh.Poligono:
        MODEL.addConstr(gp.quicksum(mu1[j] for j in range(comp1.num_puntos)) == 1, name='envConv')
        for j in range(2):
            MODEL.addConstr(x1[j] == gp.quicksum(mu1[v] * comp1.V[v][j] for v in range(comp1.num_puntos)), name='inP1')
    if type(comp1) is Elipse:
        MODEL.addConstr(comp1.P[0, 0] * x1[0] * x1[0] + comp1.P[1, 0] * x1[0] * x1[1] +
                        comp1.P[0, 1] * x1[0] * x1[1] + comp1.P[1, 1] * x1[1] * x1[1] +
                        comp1.q[0] * x1[0] + comp1.q[1] * x1[1] + comp1.r <= 0, name='inC1')
    if type(comp1) is neigh.Poligonal:
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
                MODEL.addConstr(
                    x1[j] == gp.quicksum(sublanda1[punto] * comp1.V[punto][j] for punto in range(comp1.num_puntos)),
                    name='seg1')

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
    graph : neighborhood graph
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

        landa = M.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landa')
        x_segm = M.addVars(2, vtype=GRB.CONTINUOUS, name='x_segm')
        dif = M.addVars(2, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')
        d = M.addVar(vtype=GRB.CONTINUOUS, lb=0.0, name='d')

        M.update()

        M.addConstrs(x_segm[dim] == landa * first_V[dim] + (1 - landa) * second_V[dim] for dim in range(2))
        M.addConstrs(dif[dim] >= x_segm[dim] - point[dim] for dim in range(2))
        M.addConstrs(dif[dim] >= -x_segm[dim] + point[dim] for dim in range(2))
        M.addConstr(dif[0] * dif[0] + dif[1] * dif[1] <= d * d)

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
    graph : neighborhood graph
        Graph from which we want to search the closest point to the origin.

    Returns
    -------
    list
        Closest point to punto

    """

    M = gp.Model('min_dist')

    landa1 = M.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landa1')
    landa2 = M.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landa2')
    x1 = M.addVars(2, vtype=GRB.CONTINUOUS, name='x1')
    x2 = M.addVars(2, vtype=GRB.CONTINUOUS, name='x2')
    dif = M.addVars(2, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')
    d = M.addVar(vtype=GRB.CONTINUOUS, lb=0.0, name='d')
    mu1 = M.addVars(polygon1.num_puntos, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='mu1')
    mu2 = M.addVars(polygon2.num_puntos, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='mu2')

    M.update()

    M.addConstrs(x1[dim] == landa1 * point1[dim] + (1 - landa1) * point2[dim] for dim in range(2))
    M.addConstrs(x2[dim] == landa2 * point1[dim] + (1 - landa2) * point2[dim] for dim in range(2))

    M.addConstr(gp.quicksum(mu1[v] for v in range(polygon1.num_puntos)) == 1, name='envConv1')
    M.addConstr(gp.quicksum(mu2[v] for v in range(polygon2.num_puntos)) == 1, name='envConv2')

    for dim in range(2):
        M.addConstr(x1[dim] == gp.quicksum(mu1[v] * polygon1.V[v][dim] for v in range(polygon1.num_puntos)))
        M.addConstr(x2[dim] == gp.quicksum(mu2[v] * polygon2.V[v][dim] for v in range(polygon2.num_puntos)))

    M.addConstrs(dif[dim] >= x1[dim] - x2[dim] for dim in range(2))
    M.addConstrs(dif[dim] >= -x1[dim] + x2[dim] for dim in range(2))
    M.addConstr(dif[0] * dif[0] + dif[1] * dif[1] <= d * d)

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
    graph : neighborhood graph
        Graph from which we want to search the closest point to the origin.

    Returns
    -------
    list
        Closest point to punto

    """

    M = gp.Model('min_dist')

    landa1 = M.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landa1')
    landa2 = M.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landa2')
    x1 = M.addVars(2, vtype=GRB.CONTINUOUS, name='x1')
    x2 = M.addVars(2, vtype=GRB.CONTINUOUS, name='x2')
    dif = M.addVars(2, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')
    d = M.addVar(vtype=GRB.CONTINUOUS, lb=0.0, name='d')
    mu1 = M.addVars(polygon1.num_puntos, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='mu1')
    mu2 = M.addVars(polygon2.num_puntos, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='mu2')

    M.update()

    M.addConstrs(x1[dim] == landa1 * point1[dim] + (1 - landa1) * point2[dim] for dim in range(2))
    M.addConstrs(x2[dim] == landa2 * point1[dim] + (1 - landa2) * point2[dim] for dim in range(2))

    M.addConstr(gp.quicksum(mu1[v] for v in range(polygon1.num_puntos)) == 1, name='envConv1')
    M.addConstr(gp.quicksum(mu2[v] for v in range(polygon2.num_puntos)) == 1, name='envConv2')

    for dim in range(2):
        M.addConstr(x1[dim] == gp.quicksum(mu1[v] * polygon1.V[v][dim] for v in range(polygon1.num_puntos)))
        M.addConstr(x2[dim] == gp.quicksum(mu2[v] * polygon2.V[v][dim] for v in range(polygon2.num_puntos)))

    M.addConstrs(dif[dim] >= x1[dim] - x2[dim] for dim in range(2))
    M.addConstrs(dif[dim] >= -x1[dim] + x2[dim] for dim in range(2))
    M.addConstr(dif[0] * dif[0] + dif[1] * dif[1] <= d * d)

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

    landa = MODEL.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0)
    P_segm = MODEL.addVars(2, vtype=GRB.CONTINUOUS)
    dif = MODEL.addVars(2, vtype=GRB.CONTINUOUS)
    d = MODEL.addVar(vtype=GRB.CONTINUOUS, lb=0.0)

    MODEL.update()

    MODEL.addConstrs(P_segm[dim] == landa * segment.V[0][dim] + (1 - landa) * segment.V[1][dim] for dim in range(2))

    MODEL.addConstrs(dif[dim] >= P_segm[dim] - point[dim] for dim in range(2))
    MODEL.addConstrs(dif[dim] >= -P_segm[dim] + point[dim] for dim in range(2))

    MODEL.addConstr(dif[0] * dif[0] + dif[1] * dif[1] <= d * d)

    MODEL.setObjective(d, GRB.MINIMIZE)

    MODEL.update()

    MODEL.Params.OutputFlag = 0

    MODEL.optimize()

    return d.X


def determinant(P, Q, R):
    a11 = Q[0] - P[0]
    a12 = R[0] - P[0]
    a21 = Q[1] - P[1]
    a22 = R[1] - P[1]

    return a11 * a22 - a12 * a21


def intersect(barrier1, barrier2):
    """
        Check if two segments intersect or not
        :param barrier1:
        :param barrier2:
        :return: True, if they intersect; False, otherwise.
        """
    det1 = determinant(barrier1[0], barrier2[0], barrier2[1])
    det2 = determinant(barrier1[1], barrier2[0], barrier2[1])

    det3 = determinant(barrier2[0], barrier1[0], barrier1[1])
    det4 = determinant(barrier2[1], barrier1[0], barrier1[1])

    expr1 = 0 <= det1*det2
    expr2 = 0 <= det3*det4

    expr = not((expr1 or expr2))

    return expr


def cansee(point, neighborhood, barriers):
    """
    Function that determines if a point can see a neighborhood
    :param point: Point to test if it is visible.
    :param neighborhood: Neighborhood to test if it is visible.
    :param barriers: Set of barriers that can avoid that the point is visible by the neighborhood.
    :return: True, if there exists a point in neighborhood that can be joined with point; False, otherwise.
    """
    if type(neighborhood) is neigh.Circle:
        # Define the center of the neigh.Circle
        center = np.array(neighborhood.center)

        # dr = vector joining the center and the point of the barrier
        dr = center - np.array(point)

        # dr_u = unitary vector
        dr_u = dr / np.linalg.norm(dr)

        # nr_u = normal vector to dr
        nr_u = np.array([-dr_u[1], dr_u[0]])

        # Generate a discretization of the diameter that is normal to the segment joining the center and the point
        mus = np.linspace(-neighborhood.radii, neighborhood.radii, 20)

        return any(
            [not (any([(intersect([point, center + mu * nr_u], barrier)) for barrier in barriers])) for mu in mus])

    if type(neighborhood) is neigh.Poligonal:
        # Define the extreme point of the segment
        extreme_point = np.array(neighborhood.V[0])

        # drs = vector joining the extreme points of the segment
        drs = np.array(neighborhood.V[1]) - np.array(neighborhood.V[0])

        # drs_u = unitary vector
        drs_u = drs / np.linalg.norm(drs)

        # Generate a discretization of the segment
        mus = np.linspace(0, neighborhood.longitud, 20)

        return any([not (
            any([intersect([point, extreme_point + mu * drs_u], [barrier[0], barrier[1]]) for barrier in barriers])) for
                    mu in mus])


def canseeN(neighborhood1, neighborhood2, barriers):
    """
    Extends the definition of cansee function.
    :param neighborhood1:
    :param neighborhood2:
    :param barriers:
    :return: True, if the exists a point joining neighborhood1 and neighborhood2; False, otherwhise
    """
    if type(neighborhood1) is neigh.Circle and type(neighborhood2) is neigh.Circle:
        # Define the centers of both neighborhoods
        center1 = np.array(neighborhood1.center)
        center2 = np.array(neighborhood2.center)

        # dr = vector joining the centers of the neighborhoods
        dr = center2 - center1

        # dr_u = unitary vector
        dr_u = dr / np.linalg.norm(dr)

        # nr_u = normal vector to dr
        nr_u = np.array([-dr_u[1], dr_u[0]])

        # Generate a discretization of the neighborhoods that is normal to the segment joining the center and the point
        mus1 = np.linspace(-neighborhood1.radii, neighborhood1.radii, 20)
        mus2 = np.linspace(-neighborhood2.radii, neighborhood2.radii, 20)

        return any(
            [not (any([intersect([center1 + mu1 * nr_u, center2 + mu2 * nr_u], [barrier[0], barrier[1]]) for barrier in
                       barriers]))
             for mu1 in mus1 for mu2 in mus2])

    if type(neighborhood1) is neigh.Circle and type(neighborhood2) is neigh.Poligonal:
        # Define the centers of both neighborhoods
        center = np.array(neighborhood1.center)
        extreme_point = np.array(neighborhood2.V[0])

        # dr = vector joining the centers of the neighborhoods
        dr = extreme_point - center

        # dr_u = unitary vector
        dr_u = dr / np.linalg.norm(dr)

        # nr_u = normal vector to dr
        nr_u = np.array([-dr_u[1], dr_u[0]])

        # drs = vector joining the extreme points of the segment
        drs = np.array(neighborhood2.V[1]) - np.array(neighborhood2.V[0])

        # drs_u = unitary vector
        drs_u = drs / np.linalg.norm(drs)

        # Generate a discretization of the neighborhoods that is normal to the segment joining the center and the point
        mus1 = np.linspace(-neighborhood1.radii, neighborhood1.radii, 20)
        mus2 = np.linspace(0, neighborhood2.longitud, 20)

        return any(
            [not (
                any([intersect([center + mu1 * nr_u, extreme_point + mu2 * drs_u], [barrier[0], barrier[1]]) for barrier
                     in barriers]))
             for mu1 in mus1 for mu2 in mus2])

    if type(neighborhood1) is neigh.Poligonal and type(neighborhood2) is neigh.Circle:
        neighborhood1, neighborhood2 = neighborhood2, neighborhood1

        # Define the centers of both neighborhoods
        center = np.array(neighborhood1.center)
        extreme_point = np.array(neighborhood2.V[0])

        # dr = vector joining the centers of the neighborhoods
        dr = extreme_point - center

        # dr_u = unitary vector
        dr_u = dr / np.linalg.norm(dr)

        # nr_u = normal vector to dr
        nr_u = np.array([-dr_u[1], dr_u[0]])

        # drs = vector joining the extreme points of the segment
        drs = np.array(neighborhood2.V[1]) - np.array(neighborhood2.V[0])

        # drs_u = unitary vector
        drs_u = drs / np.linalg.norm(drs)

        # Generate a discretization of the neighborhoods that is normal to the segment joining the center and the point
        mus1 = np.linspace(-neighborhood1.radii, neighborhood1.radii, 20)
        mus2 = np.linspace(0, neighborhood2.longitud, 20)

        return any(
            [not (
                any([intersect([center + mu1 * nr_u, extreme_point + mu2 * drs_u], [barrier[0], barrier[1]]) for barrier
                     in barriers]))
             for mu1 in mus1 for mu2 in mus2])

    if type(neighborhood1) is neigh.Poligonal and type(neighborhood2) is neigh.Poligonal:
        # Define the extreme points of both neighborhoods
        extreme_point1 = np.array(neighborhood1.V[0])
        extreme_point2 = np.array(neighborhood2.V[0])

        # drs = vector joining the extreme points of the segments
        drs1 = np.array(neighborhood1.V[1]) - np.array(neighborhood1.V[0])
        drs2 = np.array(neighborhood2.V[1]) - np.array(neighborhood2.V[0])

        # drs_u = unitary vectors
        drs_u1 = drs1 / np.linalg.norm(drs1)
        drs_u2 = drs2 / np.linalg.norm(drs2)

        # Generate a discretization of the neighborhoods that is normal to the segment joining the center and the point
        mus1 = np.linspace(0, neighborhood1.longitud, 20)
        mus2 = np.linspace(0, neighborhood2.longitud, 20)

        return any(
            [not (
                any([intersect([extreme_point1 + mu1 * drs_u1, extreme_point2 + mu2 * drs_u2], [barrier[0], barrier[1]])
                     for barrier in barriers]))
             for mu1 in mus1 for mu2 in mus2])


def estima_det(neighborhood, barrera):
    PB1x = barrera[0][0]
    PB2x = barrera[1][0]
    PB1y = barrera[0][1]
    PB2y = barrera[1][1]

    if type(neighborhood) is neigh.Circle:

        # print((neighborhood.center, neighborhood.radii, barrera))
        center = neighborhood.center
        radio = neighborhood.radii

        a = 1
        b = 1
        c = 0
        d = -2 * center[0]
        e = -2 * center[1]
        f = center[0] ** 2 + center[1] ** 2 - radio ** 2

        x, y = sympy.symbols("x y", real=True)

        eq1 = sympy.Eq((2*a*(PB2x-PB1x)+c*(PB2y-PB1y))*x+(c*(PB2x-PB1x)+2*b*(PB2y-PB1y))*y+(d*(PB2x-PB1x)+e*(PB2y-PB1y)), 0)
        eq2 = sympy.Eq(a*x**2+b*y**2+c*x*y+d*x+e*y+f, 0)

        sols = sympy.solve([eq1, eq2])

        x_mas = float(sols[0][x])
        y_mas = float(sols[0][y])

        x_menos = float(sols[1][x])
        y_menos = float(sols[1][y])

        eval_mas = determinant([x_mas, y_mas], [PB1x, PB1y], [PB2x, PB2y])
        eval_menos = determinant([x_menos, y_menos], [PB1x, PB1y], [PB2x, PB2y])

        L = min(eval_mas, eval_menos)
        U = max(eval_mas, eval_menos)
        # print((a, b, c, d, e, f))
        # Ecuacion seria entonces (x - center[0])^2 + (y- center[1])^2 <= radio^2
        # x^2 + center[0]^2 - 2x·center[0] + y^2 + center[1]^2 - 2y·center[1] - radio^2 <= 0
        # a = 1; b = 1; c = 0; d = -2·center[0]; e = -2·center[1]; f = center[0]^2 + center[1]^2 - radio^2

        # if abs(PB1y - PB2y) >= 1e-1:
        #
        #     m = (2 * a * (PB1x - PB2x) - c * (PB2y - PB1y)) / (2 * b * (PB2y - PB1y) - c * (PB1x - PB2x))
        #     n = (d * (PB1x - PB2x) - e * (PB2y - PB1y)) / (2 * b * (PB2y - PB1y) - c * (PB1x - PB2x))
        #
        #     # print((m, n))
        #
        #     x_mas = (-(2 * b * m * n + c * n + d + e * m) + np.sqrt(
        #         (2 * b * m * n + c * n + d + e * m) ** 2 - 4 * (a + b * m * m + c * m) * (n * n * b + e * n + f))) / (
        #                         2 * (a + b * m * m + c * m))
        #     x_menos = (-(2 * b * m * n + c * n + d + e * m) - np.sqrt(
        #         (2 * b * m * n + c * n + d + e * m) ** 2 - 4 * (a + b * m * m + c * m) * (n * n * b + e * n + f))) / (
        #                           2 * (a + b * m * m + c * m))
        #
        #     # print((x_mas, x_menos))
        #
        #     y_mas = m * x_mas + n
        #     y_menos = m * x_menos + n
        #
        #     # print((x_mas, y_mas))
        #     # print((x_menos, y_menos))
        #
        #     eval_mas = x_mas * PB1y - x_mas * PB2y + y_mas * PB2x - y_mas * PB1x + PB1x * PB2y - PB1y * PB2x
        #
        #     eval_mas = determinant([x_mas, y_mas], [PB1x, PB1y], [PB2x, PB2y])
        #     eval_menos = determinant([x_menos, y_menos], [PB1x, PB1y], [PB2x, PB2y])
        #
        #     # eval_menos = x_menos*PB1y - x_menos*PB2y + y_menos*PB2x - y_menos*PB1x + PB1x*PB2y - PB1y*PB2x
        #
        #     L = min(eval_mas, eval_menos)
        #     # U = max(eval_mas, eval_menos)
        #     U = max(eval_mas, eval_menos)
        #
        # else:
        #     x = center[0]
        #     y_mas = center[1] + radio
        #     y_menos = center[1] - radio
        #
        #     eval_mas = x * PB1y - x * PB2y + y_mas * PB2x - y_mas * PB1x + PB1x * PB2y - PB1y * PB2x
        #     eval_menos = x * PB1y - x * PB2y + y_menos * PB2x - y_menos * PB1x + PB1x * PB2y - PB1y * PB2x
        #
        #     L = min(eval_mas, eval_menos)
        #     U = max(eval_mas, eval_menos)

            # print(L)
        # print(U)
    if type(neighborhood) is neigh.Poligonal:
        x_mas = neighborhood.V[0][0]
        y_mas = neighborhood.V[0][1]

        x_menos = neighborhood.V[1][0]
        y_menos = neighborhood.V[1][1]

        eval_mas = determinant([x_mas, y_mas], [PB1x, PB1y], [PB2x, PB2y])
        eval_menos = determinant([x_menos, y_menos], [PB1x, PB1y], [PB2x, PB2y])

        L = min(eval_mas, eval_menos)
        U = max(eval_mas, eval_menos)

        # L = -20000
    # U = 20000

    return L, U


def estima_L(neighborhood, punto):
    if type(punto) is neigh.Circle:
        centro1 = neighborhood.center
        centro2 = punto.center

        dis = np.linalg.norm(np.array(centro1) - np.array(centro2)) - neighborhood.radii - punto.radii

    else:
        if type(neighborhood) is neigh.Circle:
            centro = neighborhood.center
            dis = np.linalg.norm(np.array(centro) - np.array(punto)) - neighborhood.radii

            # dist = 0

        if type(neighborhood) is neigh.Poligonal:
            dis = dist_point_segment(punto, neighborhood)
            # dr = np.array(neighborhood.V[1]) - np.array(neighborhood.V[0])
            #
            # A, B = -dr[1], dr[0]
            #
            # C = -A*neighborhood.V[0][0] - B*neighborhood.V[0][1]
            #
            # dis = abs(A*punto[0] + B*punto[1] + C)/np.sqrt(A**2 + B**2)

    return dis


def estima_U(neighborhood, punto):
    if type(punto) is neigh.Circle:
        centro1 = neighborhood.center
        centro2 = punto.center

        dis = np.linalg.norm(np.array(centro1) - np.array(centro2)) + neighborhood.radii + punto.radii

    else:
        if type(neighborhood) is neigh.Circle:
            dis = np.linalg.norm(np.array(neighborhood.center) - np.array(punto)) + neighborhood.radii

            # dist = 10000

        if type(neighborhood) is neigh.Poligonal:
            dis = max([np.linalg.norm(np.array(neighborhood.V[i]) - np.array(punto)) for i in range(2)])

    return dis

def dominant_set(neighborhoods, barriers):
    """
    Function that provides the dominant set of points candidates to be selected in the neighborhoods in the optimal
    solution.
    :param neighborhoods: neighborhoods for which we want to compute the dominant set
    :param barriers: barriers of the model
    :return: list of points that conforms the dominant set
    """

    vertices_neighborhood = list(itertools.product(range(-len(neighborhoods), 0), range(1)))
    vertices_neighborhood = vertices_neighborhood[::-1]
    # list(product(vertices_neighborhood, range(len(barriers)), range(2))) + list(product(range(len(barriers)),
    # range(2), vertices_neighborhood))

    dom_set = {}

    for (a, b) in vertices_neighborhood:
        for c in range(len(barriers)):
            for d in range(2):
                for e in range(len(barriers)):
                    for f in range(2):
                        if c <= e:
                            # Point of the barrier to check if is visible by the neighborhood
                            point1 = barriers[c][d]
                            point2 = barriers[e][f]

                            # Neighborhood to check if it is visible by the point
                            neighborhood = neighborhoods[abs(a) - 1]

                            if cansee(point1, neighborhood, barriers) and cansee(point2, neighborhood, barriers):
                                # Appending the feasible edges to edges_neighborhood
                                point = dominant_point(neighborhood, point1, point2, barriers)
                                if not(np.isnan(point[0])):
                                    dom_set[(a, c, d, e, f)] = point
                                    dom_set[(a, e, f, c, d)] = point

    print(dom_set)
    # # Second, we filter which points can be seen by the neighborhood
    # points = [point for point in points if cansee(point, neighborhood, barriers)]
    #
    # # Third, we compute for each pair of points, the dominant point in neighborhood
    # dom_set = []
    # for a in range(len(points)):
    #     for b in range(len(points)):
    #         if a <= b:
    #             point1 = points[a]
    #             point2 = points[b]
    #             dom_set.append(dominant_point(neighborhood, point1, point2))

    return dom_set

def dominant_point(neighborhood, point1, point2, barriers):
    """
    Model that computes the dominant point in its neighborhood with respect to the points point1 and point2
    :param neighborhood: neighborhood where the dominant point is located
    :param point1: list of coordinates
    :param point2: list of coordinates
    :param barriers: list of barriers that can not be crossed
    :return: dominant point
    """

    print(point1)
    print(point2)

    if type(neighborhood) is neigh.Circle:

        # We create the model
        model = gp.Model('Dominant point of a neighborhood')

        # Defining the variables
        # Dominant point
        point = model.addVars(2, vtype=GRB.CONTINUOUS, name="point")

        # Distance to point1
        dist1 = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="dist1")
        dif1 = model.addVars(2, lb=0.0, vtype=GRB.CONTINUOUS, name="dif1")

        # Distance to point2
        dist2 = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="dist2")
        dif2 = model.addVars(2, lb=0.0, vtype=GRB.CONTINUOUS, name="dif2")

        # Distance to center
        dist_inside = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="dist_inside")
        dif_inside = model.addVars(2, lb=0.0, vtype=GRB.CONTINUOUS, name="dif_inside")

        # alpha_index:
        alpha_index = [(i, j, b) for i in range(4) for j in range(2) for b in range(len(barriers))]

        # mu_index
        mu_index = [(i, j, k, b) for i in range(2) for j in range(2) for k in range(2) for b in range(len(barriers))]

        # z_index
        z_index = [(i, j, b) for i in range(2) for j in range(2) for b in range(len(barriers))]

        # gamma_index
        gamma_index = [(k, b) for k in range(2) for b in range(len(barriers))]

        # delta_index
        delta_index = [(k, b) for k in range(2) for b in range(len(barriers))]

        # Discrete variables
        alpha = model.addVars(alpha_index, vtype=GRB.BINARY, name="alpha") # first index: determinant, second index: barrier 1/2.
        mu = model.addVars(mu_index, vtype=GRB.BINARY, name="mu") # first index: max/min, second index: 1/2, third_index:barrier 1/2.
        z = model.addVars(z_index, vtype=GRB.BINARY, name="z") # first index: 12/34, second index: barrier 1/2.
        gamma = model.addVars(gamma_index, vtype=GRB.BINARY, name="gamma")
        delta = model.addVars(delta_index, vtype=GRB.BINARY, name="delta")

        # Updating the model
        model.update()

        # Constraints
        # Distance to point1
        model.addConstrs(dif1[dim] >=  point[dim] - point1[dim] for dim in range(2))
        model.addConstrs(dif1[dim] >= -point[dim] + point1[dim] for dim in range(2))
        model.addConstr(gp.quicksum(dif1[dim]*dif1[dim] for dim in range(2)) <= dist1*dist1)

        # Distance to point2
        model.addConstrs(dif2[dim] >=  point[dim] - point2[dim] for dim in range(2))
        model.addConstrs(dif2[dim] >= -point[dim] + point2[dim] for dim in range(2))
        model.addConstr(gp.quicksum(dif2[dim]*dif2[dim] for dim in range(2)) <= dist2*dist2)

        # The point is in neighborhood
        model.addConstrs(dif_inside[dim] >=  point[dim] - neighborhood.center[dim] for dim in range(2))
        model.addConstrs(dif_inside[dim] >= -point[dim] + neighborhood.center[dim] for dim in range(2))
        model.addConstr(gp.quicksum(dif_inside[dim]*dif_inside[dim] for dim in range(2)) <= dist_inside*dist_inside)
        model.addConstr(dist_inside <= neighborhood.radii)


        # BigM and SmallM
        bigM = 1e6
        smallM = -1e6

        # Points can not cross any barrier
        for u in range(len(barriers)):
            barrier = barriers[u]
            two_barriers = [[point, point1], [point, point2]]

            for b in range(2):
                barrier1 = two_barriers[b]

                det1 = determinant(barrier1[0], barrier[0], barrier[1])
                det2 = determinant(barrier1[1], barrier[0], barrier[1])

                det3 = determinant(barrier[0], barrier1[0], barrier1[1])
                det4 = determinant(barrier[1], barrier1[0], barrier1[1])

                model.addConstr(det1 <= bigM*alpha[0, b, u])
                model.addConstr(det1 >= smallM*(1 - alpha[0, b, u]))

                model.addConstr(bigM * alpha[1, b, u] - det2 >= 0)
                model.addConstr(smallM * (1 - alpha[1, b, u]) - det2 <= 0)

                model.addConstr(det3 <= bigM * alpha[2, b, u])
                model.addConstr(det3 >= smallM * (1 - alpha[2, b, u]))

                model.addConstr(det4 <= bigM * alpha[3, b, u])
                model.addConstr(det4 >= smallM * (1 - alpha[3, b, u]))

                model.addConstr(alpha[0, b, u] - alpha[1, b, u] == mu[0, 0, b, u] - mu[1, 0, b, u])
                model.addConstr(mu[0, 0, b, u] <= 1 - z[0, b, u])
                model.addConstr(mu[1, 0, b, u] <= z[0, b, u])

                model.addConstr(alpha[2, b, u] - alpha[3, b, u] == mu[0, 1, b, u] - mu[1, 1, b, u])
                model.addConstr(mu[0, 1, b, u] <= 1 - z[1, b, u])
                model.addConstr(mu[1, 1, b, u] <= z[1, b, u])

                # model.addConstr(mu[0, 0, b, u] + mu[1, 0, b, u] <= gamma[b, u])
                # model.addConstr(mu[0, 1, b, u] + mu[1, 1, b, u] <= delta[b, u])

                model.addConstr(mu[0, 0, b, u] + mu[1, 0, b, u] + mu[0, 1, b, u] + mu[1, 1, b, u] <= 1)
                # model.addConstr(gamma[b, u] + delta[b, u] <= 1)

            # model.addConstr(det1 * det2 <= 1e5 * u)
            # model.addConstr(det3 * det4 <= 1e5 * v)

            # constraint = model.addConstr(not(intersect([[point[0], point[1]], [point1[0], point1[1]]], barrier)))


        # Updating the model
        model.update()

        # Objective function
        model.setObjective(dist1+dist2, sense=GRB.MINIMIZE)

        # Updating the model
        model.update()

        # Parameters
        model.Params.OutputFlag = 1
        # model.Params.NonConvex = 2

        # Optimizing the model
        model.optimize()

        if model.Status == 3:
            model.computeIIS()
            model.write('infeasible_constraints.ilp')

            return [np.nan, np.nan]

        # Updating the model
        model.update()

    if type(neighborhood) is neigh.Poligonal:

        # We create the model
        model = gp.Model('Dominant point of a neighborhood')

        # Defining the variables
        # Dominant point
        point = model.addVars(2, vtype=GRB.CONTINUOUS, name="point")

        # Distance to point1
        dist1 = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="dist1")
        dif1 = model.addVars(2, lb=0.0, vtype=GRB.CONTINUOUS, name="dif1")

        # Distance to point2
        dist2 = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="dist2")
        dif2 = model.addVars(2, lb=0.0, vtype=GRB.CONTINUOUS, name="dif2")

        # Segment constraint
        mu_inside = model.addVar(lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="mu_inside")

        # Updating the model
        model.update()

        # Constraints
        # Distance to point1
        model.addConstrs(dif1[dim] >=  point[dim] - point1[dim] for dim in range(2))
        model.addConstrs(dif1[dim] >= -point[dim] + point1[dim] for dim in range(2))
        model.addConstr(gp.quicksum(dif1[dim]*dif1[dim] for dim in range(2)) <= dist1*dist1)

        # Distance to point2
        model.addConstrs(dif2[dim] >=  point[dim] - point2[dim] for dim in range(2))
        model.addConstrs(dif2[dim] >= -point[dim] + point2[dim] for dim in range(2))
        model.addConstr(gp.quicksum(dif2[dim]*dif2[dim] for dim in range(2)) <= dist2*dist2)

        # The point is in neighborhood
        model.addConstrs(point[dim] == mu_inside*neighborhood.V[0][dim] + (1-mu_inside)*neighborhood.V[1][dim] for dim in range(2))

        # Updating the model
        model.update()

        # Objective function
        model.setObjective(dist1+dist2, sense=GRB.MINIMIZE)

        # Updating the model
        model.update()

        # Parameters
        model.Params.OutputFlag = 0

        # Optimizing the model
        model.optimize()

        # Updating the model
        model.update()

    return [point[0].X, point[1].X]
