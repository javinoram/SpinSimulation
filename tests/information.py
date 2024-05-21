import spinsim as ss
import numpy as np

def test_partial_trace():
    original=np.array([[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]])
    new_matrix=ss.information.partial_trace_lr( original, [0.5, 0.5], 1 )
    assert np.array_equal( new_matrix, np.array([[10, 12],[18, 20]]) )
