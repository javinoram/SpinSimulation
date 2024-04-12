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
def pauli_matrices(spin: float):
    size = int( 2*spin+1 )
    return [ I_matrix(size), X_matrix(size), Y_matrix(size), Z_matrix(size) ]



def construct_term() -> sc.sparse:
    return


def construct_hamiltonian() -> sc.sparse:
    return


for i in pauli_matrices(0.5):
    print( i.toarray() )