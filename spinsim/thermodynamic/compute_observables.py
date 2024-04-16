import numpy as np
import numpy.linalg as la
from mpmath import *
np.seterr(all='raise') 

dtype = 'float64' # max value709
boltz = 8.617333262e-5 #eV/K

""" 
Funcion para calcular la probabilidad de cada uno de los estados
input: 
    - ee (numpy array): Arreglo con los valores de energia ordenados de menor a mayor
    - t (float): Valor de temperatura
output:
    - Arreglo con las probabilidades asociadas a cada valor de energia
"""
def partition_function(ee: np.array, t: float):
    # Prueba de que el exponente supere el maximo permitido de la exponencial
    # antes de que ocurra un overflow
    if np.divide(-ee[0], t*boltz, dtype=dtype)>709:
        # Aca se ocupa la libreria mpmath
        mp.dps=100
        partition = [ exp( (-e)/(t*boltz) ) for e in ee ]
        Z = sum(partition)
        partition = [ float(p/Z) if p/Z >= 1e-8 else 0.0 for p in partition ]
    else:
        # En caso de que el entero sea menor, se utiliza las funcionalidades de numpy
        partition = np.exp( np.divide(-ee, t*boltz, dtype=dtype), dtype=dtype )
        Z = np.sum(partition, dtype=dtype)
        partition = [val if val>=1e-8  else 0.0 for val in np.divide( partition, Z, dtype=dtype )]
    return np.array(partition)


""" 
Funcion que calcula el calor especifico para un conjunto de temperaturas
input: 
    - O (numpy array): Operador hermitiano al que se le quiere calcular el
    calor especifico
    - temp (numpy array): Arreglo con las temperaturas
output:
    - Valor del calor especifico en cada temperatura
"""
def hermitian_specific_heat(O: np.array, temp: np.array) -> np.array:
    # Calcular valores propios del operador O hermitiano
    ee = np.linalg.eigvalsh(O)

    # Funcion para calcular el calor espeficico
    def specific_heat(t):
        partition = partition_function(ee, t)
        # Calculo del calor especifico usando la funcion de particion
        aux1 = np.sum(partition*ee, dtype=dtype)
        aux2 = np.sum(partition*(ee**2), dtype=dtype)
        return np.divide(aux2- (aux1**2), (t*t*boltz), dtype=dtype)

    return np.array( [ specific_heat(t) for t in temp ] )


op = -0.89*1e-3*np.array([[0, 1, 0], [1,0,1],[0,1,0]])
print( hermitian_specific_heat(op, [0.0001, 0.001, 0.01, 0.1, 1]) )