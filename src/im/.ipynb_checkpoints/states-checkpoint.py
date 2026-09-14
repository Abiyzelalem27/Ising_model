



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

def LMG_hamiltonian(N, J, Omega):
    """
    Construct the Lipkin-Meshkov-Glick (LMG) Hamiltonian
    in the symmetric Dicke basis.

    The Hamiltonian is

        H = -(J/N) Sz^2 - Omega Sx

    where the Dicke-state index k counts the number of spins
    in state |1>.

    Parameters
    ----------
    N :  Total number of spins.

    J : Collective spin-spin interaction strength along the
        z-direction.

    Omega :  Strength of the transverse field acting along the
        x-direction.

    Returns
    -------
    H : LMG Hamiltonian with shape (N + 1, N + 1).

    Sx :Collective spin operator in the x-direction.

    Sy : Collective spin operator in the y-direction.

    Sz : Collective spin operator in the z-direction.

    Notes
    -----
    The Hamiltonian is constructed in the permutation-symmetric
    Dicke subspace. Its dimension is N + 1 instead of 2**N.

    For this normalization, the ground-state quantum phase
    transition occurs at Omega/J = 1 in the large-N limit.
    """

    Sx, Sy, Sz, S_plus, S_minus = (
        im.operators.collective_spin_dicke(N)
    )

    H = -(J / N) * (Sz @ Sz) - Omega * Sx

    return H, Sx, Sy, Sz, S_plus, S_minus 