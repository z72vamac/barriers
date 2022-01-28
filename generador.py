import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
import networkx as nx
from entorno import Circulo


def no_cortan(barrier1, barrier2):
    det1 = (barrier1[0][0] - barrier2[0][0])*(barrier1[1][1] - barrier2[0][1]) - (barrier1[1][0] - barrier2[0][0])*(barrier1[0][1] - barrier2[0][1])
    det2 = (barrier1[0][0] - barrier2[1][0])*(barrier1[1][1] - barrier2[1][1]) - (barrier1[1][0] - barrier2[1][0])*(barrier1[0][1] - barrier2[1][1])
    
    det3 = (barrier2[0][0] - barrier1[0][0])*(barrier2[1][1] - barrier1[0][1]) - (barrier2[1][0] - barrier1[0][0])*(barrier2[0][1] - barrier1[0][1])
    det4 = (barrier2[0][0] - barrier1[1][0])*(barrier2[1][1] - barrier1[1][1]) - (barrier2[1][0] - barrier1[1][0])*(barrier2[0][1] - barrier1[1][1])

    return (det1*det2 >= 0) or (det3*det4 >= 0)

def cortan(barrier1, barrier2):
    return not(no_cortan(barrier1, barrier2))

def genera_separacion(punto1, punto2, r = 8):
    dr = np.array(punto2)-np.array(punto1)

    dr_u = dr/np.sqrt(np.sum(dr**2))

    nr_u = np.array([-dr_u[1], dr_u[0]])

    # distancia_centros = np.linalg.norm(dr)
    
    M = np.array(punto1) + 0.5*dr
    
    # alpha = np.random.rand()
    # beta = np.random.rand()
    
    # alpha = np.random.uniform(0, r)
    # beta = np.random.uniform(0, r)
    
    P1 = M + r*nr_u
    P2 = M - r*nr_u
    
    return P1, P2

def d_PB(punto, segmento):
    
    B1 = np.array(segmento[0])
    B2 = np.array(segmento[1])
    dr = B2 - B1
    A = -dr[1]
    B = dr[0]
    C = -A*B1[0] - B*B1[1]
    
    return np.abs(A*punto[0] + B*punto[1] + C)/np.sqrt(A**2 + B**2)
    
        
def interseccion_long(punto, segmento):
    
    B1 = np.array(segmento[0])
    B2 = np.array(segmento[1])
    
    dr = B2 - B1
    nr = np.array([-dr[1], dr[0]])
    
    nr_u = nr / np.linalg.norm(nr)
    
    A = nr_u[0]
    B = nr_u[1]
    C = -A*B1[0] - B*B1[1]
    
    # r : Ax + By + C = 0

    # Sustituyendo: A*(punto[0] + landa*nr_u[0]) + B*(punto[1] + landa*nr_u[1]) + C = 0
    
    # landa(A*nr_u[0] + B*nr_u[1]) = -C - A*punto[0] - B*punto[1]
    
    # landa = (-C - A*punto[0] - B*punto[1])/(A*nr_u[0] + B*nr_u[1])
    P = punto + (-C-A*punto[0] - B*punto[1])/(A*nr_u[0] + B*nr_u[1])*nr_u
    
    landa1 = (np.linalg.norm(P-B1))/(np.linalg.norm(B1-B2))
    landa2 = (np.linalg.norm(P-B2))/(np.linalg.norm(B1-B2))
    
    long = min(landa1, landa2)*np.linalg.norm(dr)
    return long

# nP = 30

# r_dentro= 5
    

def generador(nP, index):

    print((nP, index))
    V = np.random.uniform(0, 100, (nP, 2))
        
    barreras = [[[0, 0], [100, 0]], [[100, 0], [100, 100]], [[100, 100], [0, 100]], [[0, 100], [0, 0]]]
    
    r_init = 10
    
    for i in range(10):
        # print(i)
        for punto1 in V:
            for punto2 in V:
                if any(punto1 != punto2):
                    if all([no_cortan([punto1, punto2], barrera) for barrera in barreras]):
                        r = r_init
                        P1, P2 = genera_separacion(punto1, punto2, r)
                        while any([not(no_cortan([P1, P2], barrera)) for barrera in barreras]):
                            r = r / 2
                            # print(r)
                            P1, P2 = genera_separacion(punto1, punto2, r)
                        barreras.append([P1, P2])
    
    # print(barreras)
    # print(r)
    
    # fig, ax = plt.subplots()
    
    # plt.axis([0, 100, 0, 100])
    
    circulos = []
    bolas = []
    
    def ve(punto, segmento):
        barrera1 = [punto, segmento[0]]
        flag1 = all([no_cortan(barrera1, barrera) for barrera in barreras])
        barrera2 = [punto, segmento[1]]
        flag2 = all([no_cortan(barrera2, barrera) for barrera in barreras])
        
        return flag1 | flag2
        
    for v in V:
                
        r_min = min([interseccion_long(v, barrera) for barrera in barreras if ve(v, barrera)])
        upper_bound = min(r_min, min([d_PB(v, barrera) for barrera in barreras]))
        radii = np.random.uniform(upper_bound, upper_bound)
        circulos.append(Circulo(center = v, radii = radii))
        bolas.append([v[0], v[1], radii])
    
    np.savetxt('instancias/bolas' + str(nP) + '-' + str(index) + '.csv', bolas, delimiter = ",")
    
        
    # for c in circulos:
    #     ax.add_artist(c.figura)
    
    
    segmentos = []
    for b in barreras:
        # ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c = 'red')
        segmentos.append([b[0][0], b[0][1], b[1][0], b[1][1]])
    
    np.savetxt('instancias/segmentos' + str(nP) + '-' + str(index) + '.csv', segmentos, delimiter = ",")
    
    # ax.set_aspect('equal')

np.random.seed(10)

for nP in range(5, 101, 5):
    for index in range(10):
        generador(nP, index)


# plt.show()