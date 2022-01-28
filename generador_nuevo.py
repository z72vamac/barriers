import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
import networkx as nx
from entorno import Circulo
from matplotlib.patches import Circle


nP = 10

V = np.random.uniform(0, 100, (nP, 2))

r = 1

radii = np.random.uniform(r*5, (r+1)*5, nP)

circulos = [Circulo(center = [V[i][0], V[i][1]], radii = radii[i]) for i in range(nP)]

print(circulos)