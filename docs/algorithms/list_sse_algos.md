# List of SSE Algorithms

The general anisotropic Heisenberg model with a longitudinal field admits a variety of special cases such as the XY point ($\Delta = h = 0$) and the ferromagnetic isotropic point ($\Delta = 1, h = 0$). It can be shown that the vanilla SSE algorithm allows for tremendous optimizations for such special cases due to the increased symmetry of the resulting Hamiltonians. These are referred to as deterministic algorithms and are implemented separately in the Heisenberg Ion package. Moreover, for long range interactions, we first need to construct probability tables that can be used in the sampling process. For nearest-neighbors, no probability tables need to be constructed. Therefore, the deterministic points for nearest-neighbor Hamiltonians can be optimized further. Heisenberg Ion offers a separate implementation for these nearest neighbor edge cases as well. The complete list of SSE algorithms in the Heisenberg Ion package is as follows: 

### Long Range XXZh  
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j + \Delta S_z^i S_z^j \right) + h\sum_{i} S_z^i
\label{long_range_xxzh}
$$

* Probability Tables: Heatbath or directed loops  
* Parameters: Ferromagnetic interaction matrix $J_{ij} > 0$, $\Delta, h \in \mathbb{R}$, $h \geq 0$.  
* Remarks: We restrict to $h \geq 0$ for SSE since the physics of this system is invariant under the transformation $h \rightarrow -h$

### Long Range XXZ  
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j + \Delta S_z^i S_z^j \right)
\label{long_range_xxz}
$$

* Probability Tables: Heatbath or directed loops  
* Parameters: Ferromagnetic interaction matrix: $J_{ij} > 0$, $\Delta \in \mathbb{R}$  
* Remarks: This allows for a flipping of all the spins in the lattice because of the $\mathbb{Z}_2$ symmetry in the absence of the longitudinal field

### Long Range Isotropic (Ferromagnetic Diagonals)  
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j + S_z^i S_z^j \right)
\label{long_range_iso_fm}
$$

* Probability Tables: Deterministic  
* Parameters: Ferromagnetic interaction matrix: $J_{ij} > 0$
* Remarks: This model admits a constant SSE weight representation which allows for deterministic SSE loop construction

### Long Range Isotropic (Anti-ferromagnetic Diagonals)  
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j - S_z^i S_z^j \right)
\label{long_range_iso_afm}
$$

* Probability Tables: Deterministic  
* Parameters: Ferromagnetic interaction matrix: $J_{ij} > 0$
* Remarks: This model admits a constant SSE weight representation which allows for deterministic SSE loop construction. Note that the $S_z^i S_z^j$ term in this model has a negative coefficient, as opposed to the positive coefficient in Eq. $\eqref{long_range_iso_afm}$ 

### Long Range XY
Hamiltonian: 

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j\right)
\label{long_range_xy}
$$

* Probability Tables: Deterministic  
* Parameters: Ferromagnetic interaction matrix: $J_{ij} > 0$
* Remarks: This model admits a constant SSE weight representation which allows for deterministic SSE loop construction

### Nearest Neighbor XY
Hamiltonian:

$$
H = - J \sum_{i} \left( S_x^i S_x^{i+1} + S_y^i S_y^{i+1}\right)
\label{nn_xy}
$$

* Probability Tables: None
* Parameters: Ferromagnetic interaction strength: $J > 0$
* Remarks: Sampling is deterministic. On a bipartite lattice, this Hamiltonian is equivalent to its anti-ferromagnetic ($J < 0$) variant using a sub-lattice rotation. 

### Nearest Neighbor Ferromagnetic Isotropic
Hamiltonian:

$$
H = - J \sum_{i} \left( S_x^i S_x^{i+1} + S_y^i S_y^{i+1} + S_z^i S_z^{i+1}\right)
\label{nn_iso_fm}
$$

* Probability Tables: None
* Parameters: Ferromagnetic interaction strength: $J > 0$
* Remarks: Sampling is deterministic. This model, with $J > 0, \Delta = 1$, is equivalent to the anti-ferromagnetic variant ($J < 0$) with $\Delta = -1$ on a bipartite lattice.

### Nearest Neighbor Anti-ferromagnetic Isotropic
Hamiltonian:

$$
H = - J \sum_{i} \left( S_x^i S_x^{i+1} + S_y^i S_y^{i+1} + S_z^i S_z^{i+1}\right)
\label{nn_iso_afm}
$$

* Probability Tables: None
* Parameters: Anti-ferromagnetic interaction strength: $J < 0$
* Remarks: Sampling is deterministic. This model, with $J < 0, \Delta = 1$, is equivalent to the ferromagnetic variant ($J > 0$) with $\Delta = -1$ on a bipartite lattice.

