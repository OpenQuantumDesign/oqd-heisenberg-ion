# Hamiltonian Types

The general anisotropic Heisenberg model admits a variety of special cases such as the XY point ($\Delta = h = 0$) and the ferromagnetic isotropic point ($\Delta = 1, h = 0$). It can be shown that the vanilla SSE algorithm allows for tremendous optimizations for such special cases due to the increased symmetry of the resulting Hamiltonians. These are referred to as deterministic algorithms and are implemented separately in the Heisenberg Ion package. Moreover, for long range interactions, we first need to construct probability tables that can be used in the sampling process. For nearest-neighbors, no probability tables need to be constructed. Therefore, the deterministic points for nearest-neighbor Hamiltonians can be optimized further. Heisenberg Ion offers a separate SSE implementation for these nearest neighbor edge cases as well. The complete list of Hamiltonian types in the Heisenberg Ion package is as follows: 

### XXZ with Longitudinal & Transverse Fields
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j + \Delta S_z^i S_z^j \right) + h\sum_{i} S_z^i
\label{long_range_xxzhb}
$$

* Supported Sampling Types: None. This can only be used with ED.  
* Parameters: $J_{ij}, \Delta, h, B \in \mathbb{R}$.  
* Remarks: Because of the transverse field, the $U(1)$ symmetry is broken. Therefore, SSE is not supported for this Hamiltonian. 
* Hamiltonian Name: XXZhB

### XXZ with a Longitudinal Field
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j + \Delta S_z^i S_z^j \right) + h\sum_{i} S_z^i
\label{long_range_xxzh}
$$

* Supported Sampling Types: Heatbath or directed loops  
* Parameters: $J_{ij}, \Delta, h \in \mathbb{R}$, $h \geq 0$.  
* Remarks: We restrict to $h \geq 0$ for SSE since the physics of this system is invariant under the transformation $h \rightarrow -h$. Moreover, for SSE, we also restrict to ferromagnetic interactions $J_{ij} > 0$ to avoid the sign problem. 
* Hamiltonian Name: XXZh

### XXZ
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j + \Delta S_z^i S_z^j \right)
\label{long_range_xxz}
$$

* Supported Sampling Types: Heatbath or directed loops   
* Parameters: $J_{ij}, \Delta \in \mathbb{R}$  
* Remarks: This allows for a flipping of all the spins in the lattice in SSE because of the $\mathbb{Z}_2$ symmetry in the absence of the longitudinal field. For SSE with this model, we restrict to ferromagnetic interactions $J_{ij} > 0$ to avoid the sign problem.
* Hamiltonian Name: XXZ 

### Ferromagnetic Isotropic Heisenberg  
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j + S_z^i S_z^j \right)
\label{long_range_iso_fm}
$$

* Supported Sampling Types: Deterministic 
* Parameters: $J_{ij} > 0 \in \mathbb{R}$
* Remarks: This model admits a constant SSE weight representation which allows for deterministic SSE loop construction. We restrict to ferromagnetic interactions in general: $J_{ij} > 0$. For ED with anti-ferromagnetic interactions, this can be realized using the general XXZ Hamiltonian and appropriately setting $J_ij$ and $\Delta$. On a bipartite lattice, this model, with $J > 0, \Delta = 1$, is equivalent to the anti-ferromagnetic variant ($J < 0$) with $\Delta = -1$.
* Hamiltonian Name: fm_heisenberg_fm_Z

### Ferromagnetic XXZ with Anti-ferromagnetic Diagonals  
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j - S_z^i S_z^j \right)
\label{long_range_iso_afm}
$$

* Supported Sampling Types: Deterministic  
* Parameters: $J_{ij} > 0 \in \mathbb{R}$
* Remarks: This model admits a constant SSE weight representation which allows for deterministic SSE loop construction. Note that the $S_z^i S_z^j$ term in this model has a negative coefficient, as opposed to the positive coefficient in Eq. $\eqref{long_range_iso_afm}$. For SSE, we restrict to ferromagnetic interactions: $J_{ij} > 0$. For ED with anti-ferromagnetic interactions, this can be realized using the general XXZ Hamiltonian and appropriately setting $J_ij$ and $\Delta$. On a bipartite lattice, this model, with $J > 0, \Delta = -1$, is equivalent to the anti-ferromagnetic variant ($J < 0$) with $\Delta = 1$.
* Hamiltonian Name: fm_heisenberg_afm_Z

### XY
Hamiltonian: 

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j\right)
\label{long_range_xy}
$$

* Supported Sampling Types: Deterministic  
* Parameters: $J_{ij} \in \mathbb{R}$
* Remarks: This model admits a constant SSE weight representation which allows for deterministic SSE loop construction. For SSE, we restrict to ferromagnetic interactions: $J_{ij} > 0$. On a bipartite lattice, this Hamiltonian is equivalent to its anti-ferromagnetic (corresponding to $J_{ij} < 0$) variant via a sub-lattice rotation. 
* Hamiltonian Name: XY

### Anti-ferromagnetic Isotropic
Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^{j} + S_y^i S_y^{j} + S_z^i S_z^{j}\right)
\label{nn_iso_afm}
$$

* Supported Sampling Types: Deterministic
* Parameters:  $J_{ij} < 0 \in \mathbb{R}$
* Remarks: In SSE, this model can only be used with nearest-neighbor interactions to avoid the sign problem because of the anti-ferromagnetic couplings. Sampling is deterministic. This model, with $J < 0, \Delta = 1$, is equivalent to the ferromagnetic variant ($J > 0$) with $\Delta = -1$ on a bipartite lattice. For ED, again the more general $XXZ$ Hamiltonian can be configured to simulate both ferromagnetic and anti-ferromagnetic interactions with this model. 
* Hamiltonian Name: afm_heisenberg_fm_Z

