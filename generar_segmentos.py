import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
import networkx as nx

np.random.seed(2)

nV = 100

Vx = np.random.uniform(0, 100, nV)
Vy = np.random.uniform(0, 100, nV)

V = np.array([[xi, yi] for xi, yi in zip(Vx, Vy)])

tri = Delaunay(V)

G = nx.Graph()

for v in range(nV):
    G.add_node(v)

for path in tri.simplices:
    nx.add_path(G, path)

first = list(set([x for x, y in G.edges]))

edges_G = list(G.edges).copy()

edges = []

H = nx.Graph()
acum = []

for i in first:
    possibles = [y for x, y in edges_G if x == i and x not in acum and y not in acum]
    if len(possibles) > 0:
        selected = np.random.choice(possibles)
        H.add_node(i, pos=V[i])
        H.add_node(selected, pos = V[selected])
        H.add_edge(i, selected)
        
        acum.append(i)
        acum.append(selected)
    
    # print(acum)
    # edges.append((i, selected))
    

# H.pos = nx.get_node_attributes(H, 'pos')
# print(H.pos)

# print(len(H.edges))
nEdges = len(H.edges)
# G.pos = nx.get_node_attributes(G, 'pos')

H.pos = nx.get_node_attributes(H, 'pos')

positions = np.zeros((nEdges, 4))

ind = 0
lista = []
i = 0
for arrays in H.pos.values():
    lista.append(arrays)
    ind += 1
    if ind >= 2:
        positions[i, :] = np.reshape(lista, (1, 4))
        i += 1
        ind = 0
        lista = []

np.savetxt('segments.csv', positions, delimiter = ",")

# print(H.pos)
# nx.draw(H, H.pos, node_size=10)
# plt.triplot(Vx, Vy, tri.simplices)

# plt.show()