
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
from neighborhood import *
import copy
import estimacion_M as eM
from auxiliar_functions import *
import networkx as nx


np.random.seed(33)

point1 = Punto([5, 50])
point2 = Punto([95, 50])

datos = Data([], 2, grid = True, tmax = 600, alpha = False, init = True, show = True, mode = 2)
datos.generar_muestra()

# Creo la lista de vertices
vertices = []

for c in range(len(datos.data)):
    for v in range(datos.data[c].num_segmentos):
        vertices.append(c*10+v)

# print(vertices)

vertices_pos = []

for comp in datos.data:
    for v in range(len(comp.V)-1):
        vertices_pos.append(comp.V[v])

# vertices_pos = np.array(vertices_pos)
# print(vertices_pos)

edges = []

for v in vertices:
    for w in vertices:
        if v < w:
            pol1 = v // 10
            v1 = v % 10

            pol2 = w // 10
            v2 = w % 10

            if pol1 == pol2:
                if v2 - v1 == 1:
                    edges.append((v, w))

                if v2 - v1 == datos.data[pol1].num_segmentos-1:
                    edges.append((v, w))

            if pol1 != pol2:
                if dist_point_point_polygon(datos.data[pol1].V[v1], datos.data[pol1], datos.data[pol2].V[v2], datos.data[pol2]) >= np.linalg.norm(datos.data[pol1].V[v1] - datos.data[pol2].V[v2]) - 0.5:
                    edges.append((v, w))


Ar = np.zeros((len(vertices), len(vertices)))

for i, v in zip(range(len(vertices)), vertices):
    for j, w in zip(range(len(vertices)), vertices):
        if (v, w) in edges:
            # print(i, j)
            Ar[i, j] = 1

grafo = Grafo(vertices_pos, Ar, 1)

fig, ax = plt.subplots()
plt.axis([0, 100, 0, 100])


for c in datos.data:
    ax.add_artist(c.artist)

for comp in datos.data:
    for v, i in zip(comp.V, range(comp.num_segmentos)):
        ax.annotate(str(i), xy = (v[0], v[1]))

nx.draw(grafo.G, grafo.pos)

plt.show()
