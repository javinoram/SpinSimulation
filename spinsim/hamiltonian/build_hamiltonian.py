import numpy as np
import scipy as sc
import jax

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
def X_matrix( size: int ) -> sc.sparse:
    spin = (size-1)/2.0
    base = np.zeros( (size, size) )
    for a in range(size):
        for b in range(size):
            base[a][b] = ( delta_function(a+2,b+1) + delta_function(a+1,b+2) )*np.sqrt( (spin + 1)*(a + b + 1) - (a+1)*(b+1) )
    return 0.5*sc.sparse.csr_matrix(base)

def Y_matrix( size: int ) -> sc.sparse:
    spin = (size-1)/2.0
    base = np.zeros( (size, size) )
    for a in range(size):
        for b in range(size):
            base[a][b] = ( delta_function(a+1,b+2) - delta_function(a+2,b+1) )*np.sqrt( (spin + 1)*(a + b + 1) - (a+1)*(b+1) )
    return 0.5j*sc.sparse.csr_matrix(base)

def Z_matrix( size: int ) -> sc.sparse:
    spin = (size-1)/2.0
    base = np.zeros( (size, size) )
    for a in range(size):
        for b in range(size):
            base[a][b] = delta_function(a+1,b+1)*(spin +1 - (a+1))
    return sc.sparse.csr_matrix(base)

def I_matrix( size: int ) -> sc.sparse:
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
    return { "I": I_matrix(size), "X": X_matrix(size), "Y": Y_matrix(size), "Z": Z_matrix(size) }



"""
Construir matriz asociada al operador ingresado
input:
    - operadores (string): String que representa el operador, esto debe estar ordenado segun la definicion del problema
    - spines ([float]): Lista del valor del spin en cada uno de los sitios
output:
    - Matriz en formato sparse del operador ingresado
"""
def construct_term( operator: str, spins: list) -> sc.sparse:
    # Diccionario con las matrices de los espines no repetidos
    spin_system = list( set(spins) )
    matrices = { s : pauli_matrices(s) for s in spin_system }

    # Construccion del operador mediante una aplicacion consecutiva del
    # producto kronecker
    result = matrices[ spins[0] ][operator[0]]
    for i, op in enumerate(operator[1:]):
        result = sc.sparse.kron(result, matrices[spins[i]][op] )
    return result


def construct_hamiltonian() -> sc.sparse:
    return


#for i in pauli_matrices(1.0):
#    print( i.toarray() )

#print( construct_term("II", [0.5, 0.5]).toarray() )