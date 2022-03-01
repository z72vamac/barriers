import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
import networkx as nx
from neighborhood import Circle
from matplotlib.patches import Circle


np.random.seed(2)

nV = 300

Vx = np.random.uniform(0, 100, nV)
Vy = np.random.uniform(0, 100, nV)

V = np.array([[xi, yi] for xi, yi in zip(Vx, Vy)])

tri = Delaunay(V)

# print(tri.simplices)

G = nx.Graph()

for v in range(nV):
    G.add_node(v, pos = V[v])

for path in tri.simplices:
    nx.add_path(G, path)


first = list(set([x for x, y in G.edges]))

edges_G = list(G.edges).copy()
G.pos = nx.get_node_attributes(G, 'pos')


edges = []

from shapely.geometry import Polygon

# pgon = Polygon([Vx[68], Vx[88], Vx[58]], [Vy[68], Vy[88], Vy[58]])
def genera_Circle(P1, P2, P3):    
    pgon = Polygon([P1, P2, P3])
    area = pgon.area
    semi_perimetro = pgon.exterior.length / 2
    
    a = np.linalg.norm(np.array(P2) - np.array(P3))
    b = np.linalg.norm(np.array(P3) - np.array(P1))
    c = np.linalg.norm(np.array(P1) - np.array(P2))
    
    incentro = (a*np.array(P1) + b*np.array(P2) + c*np.array(P3))/(a+b+c)
    radio_circunferencia = area / semi_perimetro 
    
    radio_aleatorio = np.random.uniform(0, radio_circunferencia)
    
    Circle = Circle(center = incentro, radii = radio_aleatorio)
    
    return Circle

Circles = []

for lista in tri.simplices:
    a, b, c = tuple(lista)
    Circle = genera_Circle([Vx[a], Vy[a]], [Vx[b], Vy[b]], [Vx[c], Vy[c]])
    
    Circles.append(Circle)


def no_cortan(barrier1, barrier2):
    det1 = (barrier1[0][0] - barrier2[0][0])*(barrier1[1][1] - barrier2[0][1]) - (barrier1[1][0] - barrier2[0][0])*(barrier1[0][1] - barrier2[0][1])
    det2 = (barrier1[0][0] - barrier2[1][0])*(barrier1[1][1] - barrier2[1][1]) - (barrier1[1][0] - barrier2[1][0])*(barrier1[0][1] - barrier2[1][1])
    
    det3 = (barrier2[0][0] - barrier1[0][0])*(barrier2[1][1] - barrier1[0][1]) - (barrier2[1][0] - barrier1[0][0])*(barrier2[0][1] - barrier1[0][1])
    det4 = (barrier2[0][0] - barrier1[1][0])*(barrier2[1][1] - barrier1[1][1]) - (barrier2[1][0] - barrier1[1][0])*(barrier2[0][1] - barrier1[1][1])

    return (det1*det2 >= 0) or (det3*det4 >= 0)

def cortan(barrier1, barrier2):
    return not(no_cortan(barrier1, barrier2))

def se_ven(Circle1, Circle2):
    
    # true si se ven, false si no se ven
    
    def genera_circunferencia(Circle):
        centro = Circle.center
        radio = Circle.radii
        angle = np.linspace(0, 2*np.pi, 10)
        
        x = centro[0] + radio*np.cos(angle)
        y = centro[1] + radio*np.sin(angle)
        
        return zip(x, y)
    
    circunferencia1 = genera_circunferencia(Circle1)
    circunferencia2 = genera_circunferencia(Circle2)

    # flag = True indica que siguen cortando
    flag = True
    for punto1 in circunferencia1:    

        for punto2 in circunferencia2:
            # print(punto1[0])
            barrera1 = [punto1, punto2]
            
            # subflag = True si cortan
            subflag = False
            
            for i, j in B.edges():
                barrera2 = [V[i], V[j]]
                
                if cortan(barrera1, barrera2):
                    subflag = True
                    break
        
            flag = flag*subflag
     
    return not(flag)       
                
def sistema(lista):
    # filtro que dice true si todas las bolas que lo forman no se ven entre ellas
    return all([not(se_ven(Circles[i], Circles[j])) for i in lista for j in lista if i != j])

# print(se_ven(Circles[i], Circles[2]))
# for i in range(4000):
#     lista = np.random.choice(len(Circles), 10, replace = False) 
#     if sistema(lista):
#         print(lista)
        # break

lista = np.random.choice(len(Circles), 10, replace = False)

N = nx.Graph()

for i in lista:
    Circle = Circles[i]
    N.add_node(i, pos = Circle.center, rad = Circle.radii)
    
# for Circle, i in zip(Circles, range(len(Circles))):
#     N.add_node(i, pos = Circle.center, rad = Circle.radii)

N.pos = nx.get_node_attributes(N, 'pos')
N.rad = nx.get_node_attributes(N, 'rad')

fig, ax = plt.subplots()

for node in N.nodes():
    ax.add_artist(Circle(xy=N.pos[node], radius = N.rad[node]))
# nx.draw(N, N.pos, node_size = N.rad)

# Generacion del grafo B
B = nx.Graph()
acum = []

for i in first:
    possibles = [y for x, y in edges_G if x == i and x not in acum and y not in acum]
    if len(possibles) > 0:
        selected = np.random.choice(possibles)
        B.add_node(i, pos=V[i])
        B.add_node(selected, pos = V[selected])
        B.add_edge(i, selected)
        
        acum.append(i)
        acum.append(selected)
    
    # print(acum)
    # edges.append((i, selected))

# H.pos = nx.get_node_attributes(H, 'pos')
# print(H.pos)

# print(len(H.edges))
nEdges = len(B.edges)

B.pos = nx.get_node_attributes(B, 'pos')

positions = np.zeros((nEdges, 4))

ind = 0
lista = []
i = 0
for arrays in B.pos.values():
    lista.append(arrays)
    ind += 1
    if ind >= 2:
        positions[i, :] = np.reshape(lista, (1, 4))
        i += 1
        ind = 0
        lista = []
        
np.savetxt('segments.csv', positions, delimiter = ",")

plt.axis([0, 100, 0, 100])
ax.set_aspect('equal')
# print(H.pos)
# fig, ax = plt.subplots()
nx.draw(G, G.pos, node_size=10, node_color = 'red', edge_color = 'red')

# nx.draw(N, N.pos, node_size = 1, with_labels = True)
# ax.triplot(Vx, Vy, tri.simplices)

# ax.add_artist(Circles[0].figura)



# plt.triplot()
plt.show()