

import numpy as np

# Pauli matrices
X = np.array([[0, 1],
                    [1, 0]], dtype=complex)

Y = np.array([[0, -1j],
                    [1j, 0]], dtype=complex)

Z = np.array([[1, 0],
                    [0, -1]], dtype=complex)

I = np.eye(2, dtype=complex)


def collective_spin(N, alpha):
    """
    Construct the collective spin operator S_alpha(Sa).

    N: number of spins
    alpha : 'x', 'y', or 'z'
    """
    
    pauli = { "x": X, "y": Y, "z": Z}[alpha] # α(alpha) is only a label for a direction.
    S_alpha = np.zeros((2**N, 2**N), dtype=complex) # Sα(S_alpha) the collective spin operator in that direction.

    for i in range(N):
        operator = np.array([[1]], dtype=complex)
        for j in range(N):
            if j == i:
                operator = np.kron(operator, pauli)
            else:
                operator = np.kron(operator, I)
        S_alpha += 0.5 * operator
    return S_alpha

def collective_spin_dicke(N):
    """
    Construct collective spin operators in the Dicke basis.

    The basis is:
        |D_0^N>, |D_1^N>, ..., |D_N^N>

    Here, k counts the number of spins in state |1>.
    """

    dimension = N + 1

    Sz = np.zeros((dimension, dimension), dtype=complex)
    S_plus = np.zeros((dimension, dimension), dtype=complex)
    S_minus = np.zeros((dimension, dimension), dtype=complex)

    for k in range(dimension):

        # Sz |D_k^N> = (N/2 - k) |D_k^N>
        Sz[k, k] = N / 2 - k

        # S+ decreases k by one
        if k > 0:
            S_plus[k - 1, k] = np.sqrt(k * (N - k + 1))

        # S- increases k by one
        if k < N:
            S_minus[k + 1, k] = np.sqrt((N - k) * (k + 1))

    Sx = (S_plus + S_minus) / 2
    Sy = (S_plus - S_minus) / (2j)

    return Sx, Sy, Sz, S_plus, S_minus 

    

def commutator(A, B):
    """
    Calculate the commutator [A, B] = AB - BA.

    Parameters
    ----------
    A, B : numpy.ndarray
        Operators represented as matrices.

    Returns
    -------
    numpy.ndarray
        The commutator of A and B.
    """

    return A @ B - B @ A 

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

    Sx, Sy, Sz, S_plus, S_minus = (collective_spin_dicke(N))
    H = -(J / N) * (Sz @ Sz) - Omega * Sx 
    return H, Sx, Sy, Sz, S_plus, S_minus 

    