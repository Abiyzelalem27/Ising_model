



import numpy as np
from itertools import combinations

def dicke_state(N, k):
    """
    Create the symmetric Dicke state |D_k^N>.
    
    N = total number of spins
    k = number of spins in state |1>
    """
    dimension = 2**N
    state = np.zeros(dimension, dtype=complex)

    # Find every possible arrangement containing k ones
    for positions in combinations(range(N), k):
        bits = ["0"] * N

        for position in positions:
            bits[position] = "1"

        basis_index = int("".join(bits), 2)
        state[basis_index] = 1

    # Normalize the state
    state = state / np.linalg.norm(state)

    return state