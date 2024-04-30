import numpy as np
import mpmath as mp
np.seterr(all='raise') 
dtype = 'float64'
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
    beta = 1.0/(t*boltz)
    ee_var = -ee*beta
    try:
        partition = np.exp( ee_var, dtype=dtype )
        Z = np.sum(partition, dtype=dtype)
        partition = np.divide( partition, Z, dtype=dtype )
    except FloatingPointError:
        mp.dps=70
        partition = [ mp.exp( e ) for e in ee_var ]
        Z = mp.fdiv( 1.0, mp.fsum(partition) )
        partition = [ float( mp.fmul(p, Z) ) for p in partition ]
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
    beta = 1.0/(t*boltz)
    ee_var = -ee*beta
    try:
        partition = np.exp( ee_var, dtype=dtype )
        Z = np.sum(partition, dtype=dtype)
        Z = np.log(Z)
    except FloatingPointError:
        mp.dps=70
        partition = [ mp.exp( e ) for e in ee_var ]
        Z = mp.fsum(partition)
        Z = mp.log(Z)
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
        aux1 = np.sum( partition*ee, dtype=dtype )
        aux2 = np.sum( partition*(ee**2), dtype=dtype )
        aux3 = aux2- (aux1**2)
        return np.divide(aux3, (t*t*boltz), dtype=dtype)

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


