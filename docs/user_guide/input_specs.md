Heisenberg Ion offers multiple entry points. The canonical approaches are to either use a tab delimited input file to call the Python package from the command line/python script, or specify the parameters in the main Python script directly. However, each of the simulator engines can also be called directly from the command line using a tab-delimited input file prepared specifically for those simulators. This section details some of the natural entry points and also describes all the possible inputs that can to be provided, depending on user requirments.  

## List of Parameters
This section lists all of the parameters that can be used to configure a simulation. 

### Categorical Inputs
The following categorical inputs always need to be specified:

| Parameter | Description | Allowed Values |
|-----------|-----------|-----------|
| simulator | Simulator name |long_range_qmc , nearest_neighbor_qmc , exact_diagonalization |
| hamiltonian_name | Hamiltonian name | XY , XXZ , XXZh , XXZhB , fm_heisenberg_fm_Z , fm_heisenberg_afm_Z , afm_heisenberg_fm_Z |
| boundary | Spatial boundary conditions | open , periodic |
| interaction_range | Range of interactions | long_range , nearest_neighbors |
| interaction_type | Type of interactions | power_law , matrix_input |
| spatial_dimension | Lattice dimension | 1d |
| lattice_type | Type of geometry | rectangular |

All of the above parameters can be specified as strings.

### Numerical Hamiltonian Specifications
Depending on the Hamiltonian specified, more parameters might be required to configure the system. The following table lists all the numerical inputs that can be used to define the Hamiltonian:  

| Parameter     | Type  | Description |
|---------------|-------|-------------|
| N  | Integer | Number of sites in the lattice. Always required |
| Delta  | Float | Coefficient of the $S_z^i S_z^j$ term in the Hamiltonian. Required for XXZ, XXZh and XXZhB hamiltonian_names |
| alpha  | Float | Defines $\alpha$ for power law interactions $\left(J_{ij} = 1/r_{ij}^\alpha\right)$. Required if interaction_range is set to long_range and the interaction_type is power_law |
| h  | Float | Longitudinal field strength. Required for XXZh and XXZhB hamiltonian_names |
| J  | Float | Energy scale of the XXZ Hamiltonian. Always required |
| B  | Float | Transverse field strength. Required for the XXZhB hamiltonian_name |
| theta  | Float | Boundary phase twist angle (required for computing spin stiffness in ED). Only used if simulator is set to exact_diagonalization |

### Sampling & Monte Carlo Parameters
If one of the QMC simulators are selected (long_range_qmc or nearest_neighbor_qmc), we need to provide sampling parameters as well. The following table lists all of the sampling parameters that can be specified:

| Parameter     | Type  | Description |  
|---------------|-------|-------------|
| T  | Float | Temperature needed to define the propagator $\left(\rho = \frac{1}{Z} e^{H/T}\right)$ |  
| loop_type  | Categorical | QMC update type. Can be deterministic, heatbath or directed_loop |
| equilibration_steps  | Int | Number of equilibration steps |
| simulation_steps  | Int | Number of steps for observable estimation |
| operator_list_update_multiplier  | Float | Rate of operator list growth in QMC. Should always be greater than 1 |
| gamma  | Float | Additional shift used in heatbath and directed_loop probability tables. Should always be greater than 0 |
| ksi  | Float | Yes | Additional shift to prevent zero transition weights in the presence of large fields. Only used for  directed_loop probabilities. Should always be greater than 0 |
| distance_dependent_offset  | Bool | Determines whether an additional bond dependent offset should be used in directed loop probabilites |
| track_spin_configurations  | Bool | Flag to determine whether QMC shot data should be recorded |
| write_final_spin_configuration  | Bool | Flag to determine whether the final QMC configuration should be recorded |
| initial_configuration_index  | Int | Initial configuration of spins. 1 or -1 starts the simulation with all spins up or down respectively. 0 picks each spin randomly. -2 uses the configuration provided. Positive integers k > 1 correspond to lattice site i starting with spin up if i mod k = 0 |
| initial_configuration_file_path  | String | Path to file with initial spin configurations for QMC. Used if initial_configuration_index is -2 |
| initial_operator_list_size  | Int | Starting size of operator list (M) in QMC. When set to -1, uses default value of M=50 |
| number_of_threads  | Int | Number of threads to use if multiple long range QMC calculations are requested. If set to -1, the number of threads is determined to be equal to the number of parameter sets |

A complete list of QMC algorithms and their target Hamiltonians can be found in the List of SSE Algorithms section of the Algorithns tab.  

### Directory Specification
In addition to simulation parameters, we also need to provide the output directories (which can be specified as strings). The following tables lists all directory parameters: 

| Parameter | Description |  
|-----------|-------------|
| root_folder | Root folder in which the outputs will be stored. Always required |
| bin_folder | Directory containing the C++ binaries. If not provided, the Python driver will attempt to compile the C++ program itself. Only required for long_range_qmc |
| cpp_source_folder | C++ source files for the long range QMC implementation. Either the bin_folder or this must be provided if using long_range_qmc. The Python driver uses this to locate the source files if compiling itself |
| julia_path | Julia source files for the ED simulator. Only needed if the simulator is exact_diagonalization |
| uuid | The UUID corresponding to each simulation. If not provided, a unique ID for each simulation is produced by the preprocessor |
| output_folder_name | The top-level output folder for the entire program execution. If not provided, the datetime stamp of execution is used to create this directory |

## Multiple Parameter Sets
Heisenberg Ion also supports the specification of multiple parameter sets in a single input file. For long range QMC, this allows for the multi-threaded execution of multiple parameter sets simultaneously. For other simulators, this simply iterates over all the parameter sets. To use this functionality, simply replace the value of the parameter(s) that need to iterated over by a comma delimited list of values in the followsing manner: 

```
simulator	long_range_qmc
hamiltonian_name	XY
interaction_range	long_range
interaction_type	power_law
N	4,5,6,7
alpha	1.0,2.0,3.0,4.0
loop_type	deterministic
```

The above example specifies 4 long range QMC runs for the long range XY model with power law interactions using the deterministic loop SSE algorithm. The 4 parameter sets are: $(N, \alpha) = (4,1.0), (5,2.0), (6,3.0), (7,4.0)$. It should be noted that for certain parameters, multple values would not be meaningful. These parameters are as follows:  

```
simulator
root_folder
simulation_folder
number_of_threads
output_folder_name
```


## Example Parameter Settings

### Exact Diagonalization

### Long Range Quantum Monte Carlo

### Nearest Neighbor Quantum Monte Carlo

