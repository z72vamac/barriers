import numpy as np
from itertools import chain, combinations

np.random.seed(5)

m = 25

puntos = np.random.uniform(0, 100, (m, 2))

print(puntos)

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(2, (len(s)+1)// 2))



# print(np.mean(puntos[(0, 3), :], axis = 0))

max = 0
conjunto1 = []

for tupla in powerset(range(m)):
    set1 = tupla
    set2 = tuple([i for i in range(m) if i not in tupla])
    
    c_masas1 = np.mean(puntos[tuple(set1), :])
    c_masas2 = np.mean(puntos[tuple(set2), :])
    
    distancia = np.linalg.norm(c_masas1-c_masas2) 
    print('Tupla: ' + str(tupla) + ", distancia: " + str(distancia))
    if distancia > max:
        max = np.linalg.norm(c_masas1-c_masas2)
        conjunto1 = set1
        

print(max)
print(conjunto1)