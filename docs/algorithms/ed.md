# Exact Diagonalization

In this section, we discuss how we can leverage the sparsity of the Hamiltonian to optimize the construction of the Hamiltonian matrix. To this end, consider the $1d$ long-range XY model:

$$
H = -J \sum_{i < j} J_{i,j} \left( S_i^x S_j^x + S_i^y S_j^y \right) = -\frac{J}{2} \sum_{i < j} J_{i,j} \left( S_i^+ S_j^- + S_i^- S_j^+ \right)
$$

For our ED calculator, to be able to compute the superfluid density, we also consider a boundary phase twist. This yields the following Hamiltonian: 

$$
H = -\frac{J}{2} \sum_{i=1}^{N} \sum_{r=1}^{(N-1)/2} J_{i,j} \left( e^{r \theta} S_i^+ S_{i+r}^- + e^{-r \theta} S_i^- S_{i+r}^+ \right)
$$

Here, we have restricted ourselves to the case of odd $N$ for simplicity. The case for even $N$ can be constructed analogously. Moreover, the above form can be obtained for open as well as periodic boundaries in the presence of long range interactions. However, a boundary phase twist would not necessarily be well-defined with open boundaries in the context of superfluidity. 

Note that the dimensions of the associated matrix are $2^N \times 2^N$. However, for a given product state (row of the matrix) in the $S_z$ basis (which is the representation we will use), only $\mathcal{O}(N^2)$ terms have a non-zero overlap. We can use this to prevent iterating over $2^N$ bitstrings twice (once for columns and once for rows, as would be required for a dense matrix). 

In practice, this can be done by iterating over all the $2^N$ columns (hence the exponential scaling with system size). Each column number $k$ admits a binary representation with $N$ bits. Therefore, we can identify $k \rightarrow |\sigma_1,...,\sigma_N \rangle$. Then, because our off-diagonal terms correspond to flipping a pair of opposite spins, we only get a non-zero matrix element if the row $l \rightarrow |\gamma_1,...,\gamma_N \rangle$ is such that $\gamma_t = \sigma_t$ for all $t \in \{1,...,N\}$ except for a single pair of indices, $t_1, t_2 \in \{1,...,N\}$ with $\sigma_{t_1} \neq \sigma_{t_2}$. 

To isolate all rows with a non-zero overlap, given a column, we can use bit-wise operations. We proceed by iterating over all sites $i$ and all distances $r$ in the lattice. This allows us to construct all pairs of sites $(i,j)$. Then, representing the column number $k$ in binary, we can isolate the bits in $k$ corresponding to sites $i$ and $j$ via a bitwise arithmetic shift to the right by $i$ and $j$ respectively. If the resulting bits $\sigma_i, \sigma_j$ are opposite relative to one another, the XY term contributes a non-zero element. The following code block exhibits this construction: 

```
σi = (k >> i) & 1
σj = (k >> j) & 1
mask = (σi == σj) ? 0 : 1
non_zero_element = -0.5 * J * mask * exp((σj - σi) * r * theta * 1im)
```

The first two lines in the above block isolate the bits $\sigma_i, \sigma_j$ by performing a bitwise AND operation with $1$ after the shift discussed above. The `mask` determines if the spins are flipped relative to each other. The row number associated with this matrix element can also be determined via bitwise operations similarly. Representing the column number $k$ as a bitstring, the corresponding row number bitstring can be obtained by performing a bitwise XOR operation with the bitstrings $[2^i]$ and $[2^j]$. These strings have a $1$ at location $i$ and $j$ respectively, and $0$ at all other $N-1$ locations. The row number $b$ can be obtained in Julia via the following code using this approach: 

```
b = k ⊻ 2^i ⊻ 2^j
```

The bitstring associated with $b$ now is the same bitstring as $k$ except for the $i$ and $j$ positions which are flipped relative to the bitstring corresponding to $k$. To construct the diagonal entries, we can similarly use bitwise operations to map bits to the eigenvalues of $S_z$. Since those are diagonal, no row index needs to be constructed to record those in the Hamiltonian. Diagonalization of the resulting Hamiltonian is done via the native LinearAlgebra library in Julia. 

