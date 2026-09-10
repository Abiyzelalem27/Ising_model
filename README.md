

# Ising_model

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A high-performance Python package for constructing **Collective Spin Operators** and the **All-to-All Transverse-Field Ising Model**. 

Designed for quantum dynamics simulations, this package supports exact operator representations in both full Hilbert space ($2^N \times 2^N$) and the permutation-symmetric **Dicke basis** ($(N+1) \times (N+1)$), enabling simulation of large spin systems ($N \gg 100$).

---

## 📌 Mathematical Background

### 1. Collective Spin Operators
The collective spin operator in direction $\alpha \in \{x,y,z\}$ for an $N$-qubit system is:
$$S_\alpha = \sum_{i=1}^{N}s_\alpha^{(i)} = \frac{1}{2}\sum_{i=1}^{N}\sigma_\alpha^{(i)}$$

In the symmetric Dicke basis $\{\vert{}D_0^N\rangle, \vert{}D_1^N\rangle, \dots, \vert{}D_N^N\rangle\}$, matrix dimensions drop exponentially from $2^N \times 2^N$ to $(N+1) \times (N+1)$:
* **$S_z$ diagonal action:** $S_z\vert{}D_k^N\rangle = \left(\frac{N}{2}-k\right)\vert{}D_k^N\rangle$
* **Ladder operators ($S_\pm$):**
  $$S_+\vert{}D_k^N\rangle = \sqrt{k(N-k+1)}\,\vert{}D_{k-1}^N\rangle, \qquad S_-\vert{}D_k^N\rangle = \sqrt{(N-k)(k+1)}\,\vert{}D_{k+1}^N\rangle$$

### 2. All-to-All Transverse-Field Ising Hamiltonian
The Hamiltonian with pairwise interactions and an external transverse field $\Omega$ is:
$$H = J\sum_{i<j}\sigma_z^{(i)}\sigma_z^{(j)} + \Omega\sum_i\sigma_x^{(i)}$$

Using collective operator identities:
* **Distinct pair convention ($\sum_{i<j}$):** $H = 2J S_z^2 + 2\Omega S_x - \frac{JN}{2}I \;\equiv\; 2J S_z^2 + 2\Omega S_x$
* **Double-counted convention ($\sum_{i\neq j}$):** $H = 4J S_z^2 + 2\Omega S_x - JNI \;\equiv\; 4J S_z^2 + 2\Omega S_x$

---

## 🚀 Installation 

Clone the repository and install locally in editable mode:

```bash
git clone [https://github.com/Abiyzelalem27/Ising_model.git](https://github.com/Abiyzelalem27/Ising_model.git)
cd Ising_model
pip install -e .