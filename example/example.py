import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product, permutations, chain
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
import neighborhood as e
import copy
import estimacion_M as eM
from auxiliar_functions import *
import networkx as nx
from HSPPS import HSPPS

help(HSPPS)

source = e.Elipse(P = np.identity(2), q = -8*np.ones(2), r = 4**2 + 4**2 - 2**2)

destination = e.Elipse(P = np.identity(2), q = np.array([-38, -24]), r = 19**2 + 12**2 - 3**2)

# example = e.Elipse(P = np.identity(2), q = np.array([-100, -100]), r = 2*(50**2)-2500)
fig, ax = plt.subplots()

ax.add_artist(source.artist)
ax.add_artist(destination.artist)
# ax.add_artist(example.artist)

plt.xlim(0, 100)
plt.ylim(0, 100)

ax.set_aspect('equal')




plt.show()