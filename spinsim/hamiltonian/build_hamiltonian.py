import numpy as np
import scipy as sc

"""
Kronecker delta funcion
"""
def delta_function(a: int, b: int) -> int:
    if a == b:
        return 1
    return 0


""" 
Construir las matrices de pauli para spin S
input: 
    - size (float): tamaño de la matriz 2*S + 1
output:
    - I: pauli X de spin S en formato sparse
    - X: pauli I de spin S en formato sparse
    - Y: pauli Y de spin S en formato sparse
    - Z: pauli Z de spin S en formato sparse
"""
def x_matrix( size: int ) -> sc.sparse:
    spin = (size-1)/2.0
    base = np.zeros( (size, size) )
    for a in range(size):
        for b in range(size):
            base[a][b] = ( delta_function(a+2,b+1) + delta_function(a+1,b+2) )*np.sqrt( (spin + 1)*(a + b + 1) - (a+1)*(b+1) )
    return 0.5*sc.sparse.csr_matrix(base)

def y_matrix( size: int ) -> sc.sparse:
    spin = (size-1)/2.0
    base = np.zeros( (size, size) )
    for a in range(size):
        for b in range(size):
            base[a][b] = ( delta_function(a+1,b+2) - delta_function(a+2,b+1) )*np.sqrt( (spin + 1)*(a + b + 1) - (a+1)*(b+1) )
    return 0.5j*sc.sparse.csr_matrix(base)

def z_matrix( size: int ) -> sc.sparse:
    spin = (size-1)/2.0
    base = np.zeros( (size, size) )
    for a in range(size):
        for b in range(size):
            base[a][b] = delta_function(a+1,b+1)*(spin +1 - (a+1))
    return sc.sparse.csr_matrix(base)

def i_matrix( size: int ) -> sc.sparse:
    base = np.eye( size )
    return sc.sparse.csr_matrix(base)


""" 
Construir las matrices de pauli para spin S
input: 
    - spin (float): indica el valor del espin de cada sitio
output:
    - I: pauli X de spin S en formato sparse
    - X: pauli I de spin S en formato sparse
    - Y: pauli Y de spin S en formato sparse
    - Z: pauli Z de spin S en formato sparse
"""
def pauli_matrices(spin: float) -> dict:
    size = int( 2*spin+1 )
    return { "I": i_matrix(size), "X": x_matrix(size), "Y": y_matrix(size), "Z": z_matrix(size) }



"""
Construir matriz asociada al operador ingresado
input:
    - operator (string): String que representa el operador, esto debe estar ordenado segun la definicion del problema
    - spines ([float]): Lista del valor del spin en cada uno de los sitios
output:
    - Arreglo de numpy que representa el operador ingresado
"""
def construct_term( operator: str, spins: list) -> np.array:
    # Diccionario con las matrices de los espines no repetidos
    spin_system = list( set(spins) )
    matrices = { s : pauli_matrices(s) for s in spin_system }

    # Construccion del operador mediante una aplicacion consecutiva del
    # producto kronecker
    result = matrices[ spins[0] ][operator[0]]
    for i, op in enumerate(operator[1:]):
        result = sc.sparse.kron(result, matrices[spins[i+1]][op] )
    
    # Filtro para limpiar la matriz de terminos imaginarios cuando la cantidad de
    # operadores 'Y' es par.
    if operator.count('Y')%2 == 0:
        result = np.real(result)
    return result.toarray()


"""
Construir matriz asociada al hamiltoniano ingresado
input:
    - operator ([string]): Lista de string, cada uno representa un operador del hamiltoniano
    - spines ([float]): Lista del valor del spin en cada uno de los sitios
output:
    - Arreglo de numpy que representa el hamiltoniano
"""
def construct_hamiltonian(list_operators: list, spins: list) -> np.array:
    size = np.prod( 2*np.array(spins)+1 )
    base = np.zeros( (int(size), int(size)) )

    # Loop para construir cada termino del hamiltoniano
    for (exchange, op) in list_operators:
        base = base + exchange*construct_term(op, spins)
    return base
