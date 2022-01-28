from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(5)

points = np.random.uniform(0, 1, (10, 2))

vor = Voronoi(points)

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

def recorta_barrera(barrera):
    """ Input: Lista de dos vertices """
    
    contador = 0
    
    tupla = ()
    for key, values in vor.ridge_dict.items():
        contador += 1
        if values[0] == barrera[0] and values[1] == barrera[1]:
            tupla = key
    
    if contador > len(vor.ridge_dict)+1:
        P1 = [0, 0]
        P2 = [0, 0]
        
        return P1, P2
    
    PB1 = vor.vertices[barrera[0]]
    PB2 = vor.vertices[barrera[1]]
    
    dB12 = PB2 - PB1
    mod_dB12 = np.linalg.norm(dB12)
    
    dB12_u = dB12/mod_dB12
    
    PN1 = vor.points[tupla[0]]
    PN2 = vor.points[tupla[1]]
    
    M = (PN1 + PN2) / 2
    
    dM1 = M - PB1
    mod_dM1 = np.linalg.norm(dM1)
    
    landa = mod_dM1/mod_dB12
    
    deviation = min(landa, 1-landa)
    
    P1 = M + dB12*deviation
    P2 = M - dB12*deviation
    
    return P1, P2

def genera_separacion(adyacencia):
    # dr = np.array(punto2)-np.array(punto1)
    
    # dr_u = dr/np.sqrt(np.sum(dr**2))
    
    barrera = vor.ridge_dict[adyacencia]
    # contador = 0
    #
    # tupla = ()
    # for key, values in vor.ridge_dict.items():
    #     contador += 1
    #     if values[0] == barrera[0] and values[1] == barrera[1]:
    #         tupla = key
    #
    # if contador > len(vor.ridge_dict)+1:
    #     P1 = [0, 0]
    #     P2 = [0, 0]
    #
    #     return P1, P2
    
    PB1 = vor.vertices[barrera[0]]
    PB2 = vor.vertices[barrera[1]]
    
    dB12 = PB2 - PB1
    mod_dB12 = np.linalg.norm(dB12)
    
    dB12_u = dB12/mod_dB12
    
    PN1 = vor.points[adyacencia[0]]
    PN2 = vor.points[adyacencia[1]]
    
    M = (PN1 + PN2) / 2
    
    dM1 = M - PB1
    mod_dM1 = np.linalg.norm(dM1)
    
    landa = mod_dM1/mod_dB12
    
    # alpha = np.random.rand()
    # beta = np.random.rand()
    
    alpha = 0.1
    beta = 0.1
    
    deviation = min(landa, 1-landa)
    
    P1 = M + dB12*deviation
    P2 = M - dB12*deviation
    
    return P1, P2


fig, ax = plt.subplots()

plt.axis([0, 1, 0, 1])

for i in points:
    plt.scatter(i[0], i[1], s = 5, c = 'black')

# plt.scatter(points[:, 0], points[:, 1], s = 1)


for a in adyacencias:
    P1, P2 = genera_separacion(a)
    ax.plot([P1[0], P2[0]], [P1[1], P2[1]], '-', color = 'red')
    # ax.annotate(str(b), xy = (P1[0], P1[1]))
    
    # PB1 = vor.vertices[b[0]]
    # PB2 = vor.vertices[b[1]]
    #
    # ax.plot([PB1[0], PB2[0]], [PB1[1], PB2[1]], '-', color = 'red')
    
    

plt.show()

    
    
    
    
# import matplotlib.pyplot as plt
#
# fig, ax = plt.subplots()
#
# # plt.axis([0, 10, 0, 10])
# fig = voronoi_plot_2d(vor, ax = ax)
#
# ax.set_xlim(left = -2, right = 2)
# ax.set_ylim(bottom = -2, top = 2)
#
# plt.show()

