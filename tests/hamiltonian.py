import spinsim as ss
import numpy as np

def test_hamiltonian_matrix():
    ##Matrix using the library
    indices = [ (0,1) ]
    exchanges = [ 1 ]
    espines = [ 0.5, 0.5 ]
    SiSjvectores = ss.operators.set_sij_vector(indices, 2) 
    terminos = []
    for J, op in zip(exchanges, SiSjvectores):
        terminos += [ [J, op[0]], [J, op[1]], [J, op[2]] ] 
    H = ss.hamiltonian.construct_hamiltonian(terminos, espines) 

    ##By-Hand Matrix
    sx = 0.5*np.array([[0,1], [1,0]])
    sy = 0.5*np.array([[0,-1j], [1j,0]])
    sz = 0.5*np.array([[1,0], [0,-1]])
    H2 = np.real( np.kron( sx, sx ) + np.kron( sy, sy ) + np.kron( sz, sz ) )
    assert np.array_equal( H, H2 )