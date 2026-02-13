## Hamiltonian

Heisenberg Ion is a Quantum Monte Carlo (QMC) Stochastic Series Expansion (SSE) and Exact Diagonalization (ED) simulator for the anisotropic Heisenberg model with a longitudinal field, defined by the Hamiltonian:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j + \Delta S_z^i S_z^j \right) + h\sum_{i} S_z^i
\label{heisenberg_z_field}
$$

Here, $S_k^i = \frac{1}{2} \sigma_k^i$ is the spin-1/2 operator corresponding to site $i$, for each $k \in \{x,y,z\}$. The Pauli matrices $\sigma_k$ are defined as follows: 

$$
\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \quad \sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix} \quad \sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$

For QMC, if the interactions are long-range, the coupling matrix must be ferromagnetic: $J_{ij} > 0 \ \forall \ 1 \leq i,j \leq N$, to avoid the sign problem. Our ED engine allows for entirely arbitrary couplings and also supports the inclusion of a transverse field term. The general ED Hamiltonian is:

$$
H = - \sum_{i < j} J_{ij} \left( S_x^i S_x^j + S_y^i S_y^j + \Delta S_z^i S_z^j \right) + h\sum_{i} S_z^i - B\sum_{i} S_x^i
\label{heisenberg_z_x_fields}
$$

## Exact Diagonalization

ED corresponds to constructing the entire $2^N \times 2^N$ ($N$ is the total number of sites on the lattice) dimensional Hamiltonian matrix in the $S_z$ basis, and diagonalizing it numerically. Since the Hamiltonian dimensions scale exponentially with system size, ED becomes prohibitively expensive for large $N$. In fact, it is easy to show that ED exhibits exponential scaling in both time and memory. However, since it is exact and involves a relatively transparent implementation, it is instrumental for validating more complex algorithmic approaches for many-body physics. Technical details about how the Hamiltonian is constructed in Heisenberg Ion can be found in the Exact Diagonalization section of the Algorithms tab. 

## Quantum Monte Carlo

QMC is a statistical approach for studying the equilibrium properties of a many-body system. For systems without frustration, such as the isotropic Heisenberg model with nearest-neighbor interactions, it is widely believed that the SSE variant of QMC scales polynomially. Technical details about the Heisenberg Ion SSE approach can be found in the Stochastic Series Expansion section of the Algorithms section. For a complete list of the SSE algorithms implemented in Heisenberg Ion, see the List of SSE Algorithms section of the Algorithms tab.  


## General Workflow

The general approach for using any of the simulators in the Heisenberg Ion is the same. The workflow can be divided into the following categories: 

### Input Specification  
This includes defining the inputs and passing them to the package. These are then parsed to determine the workflow settings. Validation is owned by the components being used by the program. This allows for lightweight pre-processing and a modular, workflow agnostic design. Details about specifying the inputs can be found in the Input Format section. 

### Preprocessing  
This involves preparing the input files for the simulator specified by the user. Since the simulators are considerably different in the inputs they require (e.g. we need to produce probability tables for QMC with long range interaction but not for nearest neighbor interactions), this is implemented polymorphically. Details about the different preprocessors can be found in the Pre-Process section of the Development Guide. 

### Drivers  
Our QMC implementation uses C++ for optimal performance, and our ED calculator uses a Julia backend. The driver layer uses a command line subprocess to interface with these engines as needed. For nearest-neighbor QMC, the driver simply calls the relevant implementation class in Python. As with the preprocessors, each simulator carries its own driver implemented as a subclass of the base driver class. Details can be found in the Driver section of the Development Guide. 

### Engines  
These modules contain the algorithmic simulation logic. Currently, we support three simulator engines: long-range QMC, nearest-neighbor QMC and ED. The output data files in each case are produced by the engine layer. For implementation details, see the Engines section Development Guide and for technical details, see the relevant sections of the Algorithms tab. Details about the output format can be found in the Output Format section of the User Guide. 

### Postprocessing  
Since the post-processing can differ significantly between workflows, we offer some lightweight post-processing utilities. This includes statistical utilities for the QMC outputs, and simple equilibrium property calculators for the ED data. Some post-processing examples can be found in the Examples section of the User Guide. Details can be found in the Postprocess section of the Development Guide.  
