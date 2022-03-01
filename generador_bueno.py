import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
import networkx as nx
from neighborhood import Circle
from matplotlib.patches import Circle
from scipy.spatial import Voronoi, voronoi_plot_2d

np.random.seed(1)

nN = 20
r = 10

Vx = np.random.uniform(0+r, 100-r, nN)
Vy = np.random.uniform(0+r, 100-r, nN)

V = np.array([[xi, yi] for xi, yi in zip(Vx, Vy)])

vor = Voronoi(V)

adyacencias = []
barreras = []

for key, vals in vor.ridge_dict.items():
    if -1 not in vals:
        adyacencias.append(key)
        barreras.append(vals)
        
print(adyacencias)
print(barreras)

fig = voronoi_plot_2d(vor)
plt.show()

def genera_separacion(punto1, punto2):
    dr = np.array(punto2)-np.array(punto1)

    dr_u = dr/np.sqrt(np.sum(dr**2))

    nr_u = np.array([-dr_u[1], dr_u[0]])

    distancia_centros = np.linalg.norm(dr)
    
    M = np.array(punto1) + 0.5*dr
    
    # alpha = np.random.rand()
    # beta = np.random.rand()
    
    alpha = 4
    beta = 4
    
    P1 = M + alpha*nr_u
    P2 = M - beta*nr_u
    
    return P1, P2


fig, ax = plt.subplots()

plt.axis([0, 100, 0, 100])

bolas = []

ax.set_aspect('equal')

for c in V:
    Circle = Circle(center = c, radii = np.random.uniform(0, 0.5))
    ax.add_artist(Circle.figura)
    bolas.append([Circle.center[0], Circle.center[1], Circle.radii])

# plt.scatter(Vx, Vy, s = 20)

np.savetxt('bolas.csv', bolas, delimiter = ",")

import operator

barreras = []

for i in range(len(V)):
    distancias = {}

    for j in range(len(V)):
        if (i, j) in adyacencias:
    #         distancias[j] = np.linalg.norm(np.array(subconjuntos[i].center) - np.array(subconjuntos[j].center))
    #
    # distancias = sorted(distancias.items(), key = operator.itemgetter(1))[0:5]
    # for key, value in distancias:
            P1, P2 = genera_separacion(V[i], V[j])

# data = [(lista[0][0], lista[0][1]), (lista[1][0], lista[1][1])]

            ax.plot([P1[0], P2[0]], [P1[1], P2[1]], '-', color = 'red')
    
            barreras.append([P1[0], P1[1], P2[0], P2[1]])
# ax.add_artist(Circle(center = [50, 50], radii = 50).figura)

np.savetxt('segments.csv', barreras, delimiter = ",")

plt.show()

# contador = 0
#
# Circles = []
# while contador <= nN:
#     vx = np.random.uniform(0+r, 100-r)
#     vy = np.random.uniform(0+r, 100-r)
#
#     radio = np.random.uniform(0, r)
#
#     if all([np.linalg.norm(np.array([vx, vy])-np.array(Circle.center)) >= Circle.radii+radio+5 for Circle in Circles]):
#         Circle_new = Circle(center = [vx, vy], radii = radio)
#         Circles.append(Circle_new)
#         contador += 1
#
#
#
# def genera_separacion(Circle1, Circle2):
#     C1 = Circle1.center
#     r1 = Circle1.radii
#
#     C2 = Circle2.center
#     r2 = Circle2.radii
#
#     # print((r1, r2))
#
#     r_min = min([r1, r2])
#
#     if r_min == r2:
#         copia_c = C1
#         copia_r = r1
#
#         C1 = C2
#         r1 = r2
#
#         C2 = copia_c
#         r2 = copia_r
#
#     r_max = max([r1, r2])
#
#     dr = np.array(C2)-np.array(C1)
#
#     dr_u = dr/np.sqrt(np.sum(dr**2))
#
#     nr_u = np.array([-dr_u[1], dr_u[0]])
#
#     distancia_centros = np.linalg.norm(dr)
#
#     # print((r1, r2, distancia_centros))
#
#     if r1 > r2:
#         print('ERROR')
#
#     # np.random.uniform(r_max, 1-r_max)
#     M = np.array(C1) + dr_u*(r1 + (distancia_centros - r1 - r2)/2)
#
#     # ax.plot(M[0], M[1])
#
#
#     r_dist1 = np.random.uniform(r_max, r_max)
#     r_dist2 = np.random.uniform(r_max, r_max)
#
#     P1 = M + r_dist1*nr_u
#     P2 = M - r_dist2*nr_u
#
#     return P1, P2
#
#
# subconjuntos = Circles[0:10]
#
# fig, ax = plt.subplots()
#
# plt.axis([0, 100, 0, 100])
#
# bolas = []
#
# ax.set_aspect('equal')
# for c in subconjuntos:
#     ax.add_artist(c.figura)
#     bolas.append([c.center[0], c.center[1], c.radii])
#
# np.savetxt('bolas.csv', bolas, delimiter = ",")
#
# import operator
#
# barreras = []
#
# for i in range(len(subconjuntos)):
#     distancias = {}
#
#     for j in range(len(subconjuntos)):
#         if j != i:
#             distancias[j] = np.linalg.norm(np.array(subconjuntos[i].center) - np.array(subconjuntos[j].center))
#
#     distancias = sorted(distancias.items(), key = operator.itemgetter(1))[0:5]
#     for key, value in distancias:
#         P1, P2 = genera_separacion(subconjuntos[i], subconjuntos[key])
#
# # data = [(lista[0][0], lista[0][1]), (lista[1][0], lista[1][1])]
#
#         ax.plot([P1[0], P2[0]], [P1[1], P2[1]], '-', color = 'red')
#
#         barreras.append([P1[0], P1[1], P2[0], P2[1]])
# # ax.add_artist(Circle(center = [50, 50], radii = 50).figura)
#
# np.savetxt('segments.csv', barreras, delimiter = ",")
#
# plt.show()

    
# genera_separacion(Circle1, Circle2) 
