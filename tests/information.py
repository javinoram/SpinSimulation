import spinsim as ss
import numpy as np

def test_partial_trace():
    original=np.array([[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]])
    new_matrix=ss.information.partial_trace_lr( original, [0.5, 0.5], 1 )
    assert np.array_equal( new_matrix, np.array([[10, 12],[18, 20]]) )


def test_vn_entropy():
    rho = np.array([1, 0])
    assert np.abs( ss.information.von_neumann_entropy(rho) ) <= 1e-7


def test_concurrence():
    rho = (0.5)*np.array([[1,1], [1,1]])
    assert np.abs( ss.information.concurrence(rho) ) <= 1e-7


def test_fidelity1():
    state1 = (1/np.sqrt(2))*np.array([1,1])
    state2 = (1/np.sqrt(2))*np.array([1,1])
    assert np.abs( 1 - ss.information.fidelity_states(state1, state2) ) <= 1e-7

def test_fidelity2(): 
    state1 = np.array([0,1])
    state2 = np.array([1,0])
    assert np.abs( ss.information.fidelity_states(state1, state2) ) <= 1e-7
