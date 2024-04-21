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
def prob_states(ee: np.array, t: float):
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
Funcion para calcular el logaritmo de la funcion de particion Z
input: 
    - ee (numpy array): Arreglo con los valores de energia ordenados de menor a mayor
    - t (float): Valor de temperatura
output:
    - Logaritmo natural de Z
"""
def Z_function(ee: np.array, t: float) -> np.array:
    # Prueba de que el exponente supere el maximo permitido de la exponencial
    # antes de que ocurra un overflow
    if np.divide(-ee[0], t*boltz, dtype=dtype)>709:
        # Aca se ocupa la libreria mpmath
        mp.dps=100
        partition = [ exp( (-e)/(t*boltz) ) for e in ee ]
        Z = sum(partition)
        Z = log(Z)
    else:
        # En caso de que el entero sea menor, se utiliza las funcionalidades de numpy
        partition = np.exp( np.divide(-ee, t*boltz, dtype=dtype), dtype=dtype )
        Z = np.sum(partition, dtype=dtype)
        Z =np.log(Z)
    return float(Z)



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
    def specific_heat(t: float) -> float:
        partition = prob_states(ee, t)
        aux1 = np.sum(partition*ee, dtype=dtype)
        aux2 = np.sum(partition*(ee**2), dtype=dtype)
        return np.divide(aux2- (aux1**2), (t*t*boltz), dtype=dtype)

    return np.array( [ specific_heat(t) for t in temp ] )


""" 
Funcion que calcula la entropia para un conjunto de temperaturas
input: 
    - O (numpy array): Operador hermitiano al que se le quiere calcular la entropia
    - temp (numpy array): Arreglo con las temperaturas
output:
    - Valor de la entropia en cada temperatura
"""
def hermitian_entropy(O: np.array, temp: np.array) -> np.array:
    # Calcular valores propios del operador O hermitiano
    ee = np.linalg.eigvalsh(O)

    # Funcion para calcular la entropia
    def entropy(t: float) -> float:
        partition = prob_states(ee, t)
        termal = np.sum( ee*partition, dtype=dtype)
        free_energy = -boltz*Z_function(ee, t)
        return np.divide(termal, t, dtype=dtype) - free_energy 
    
    return np.array( [ entropy(t) for t in temp ] )


""" 
Funcion que calcula el valor esperado de un operador
input: 
    - O (numpy array): Operador hermitiano al que se le toman los valores y vectores propios
    - Op (numpy array): Operador al que se le calcula el valor esperado
    - temp (numpy array): Arreglo con las temperaturas
output:
    - Arreglo de los valores esperados a diferentes temperaturas
"""
def hermitian_expected_value(O: np.array, Op: np.array, temp: np.array) -> np.array:
    # Calcular valores propios del operador O hermitiano
    ee, vv = np.linalg.eigh(O)
    proyeccion= np.array( [ ((vv[:,k]).T.conj()).dot(Op).dot(vv[:,k]) for k in range(len(ee))] )

    # Funcion para calcular el valor esperado generico
    def valor_esperado(t: float) -> float: 
        partition = prob_states(ee, t)
        exp_val = np.sum( proyeccion*partition, dtype=dtype)
        return exp_val   
    return np.array( [ valor_esperado(t) for t in temp ] )
