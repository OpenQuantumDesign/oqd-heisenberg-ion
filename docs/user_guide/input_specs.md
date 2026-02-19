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

### Monte Carlo Seeds
We can also optionally specify seeds to initialize the random number generator needed for QMC. The following table lists the different seeds that can be provided to the SSE |simulators:

| Seed Name   | Description |  
|-------------|-------------|
initial_config_seed | Used if the initial_config_index is set to 0 to randomly pick the starting configuration of each spin | 
disconnected_spin_flip_seed | In SSE, disconnected spins can be flipped independently. This seed initializes the generator used for those spin flips |
off_diagonal_update_seed | Needed to seed the off-diagonal update random number generators |
metropolis_insert_seed | Needed to seed the generator used in the Metropolis step of diagonal operator insertions |
metropolis_remove_seed | Needed to seed the generator used in the Metropolis step of diagonal operator insertions |
diagonal_update_seed | Only used for long range QMC. Needed to sample the bond index for operator insertions in diagonal updates |
exit_leg_seed | Only needed for XY, XXZ and XXZh hamiltonian types. Used to initialize generator for selecting exit legs while constructing loops |
metropolis_bond_generator_seed | Only needed for XXZ and XXZh hamiltonian types. Used to probabilistically insert operators at a given bond in diagonal update step in the two-step scheme |

### Conflicting Inputs
Note that if contradictory inputs are specified for the Hamiltonian, (for instance if the input `hamiltonian_name` is `XY` and `Delta` is specified to be `1.0`), the system specification defaults to properties defined by the `hamiltonian_name` (the provided value of `Delta` is ignored). The same is true for specifying the interactions. If the `interaction_range` is set to `nearest_neighbor` and the `interaction_type` is specified to be `power_law`, the latter is ignored. Similarly, conflicting inputs pertaining to QMC sampler settings are resolved via the `loop_type` parameter. 

### Directory Specification
In addition to simulation parameters, we also need to provide the output directories (which can be specified as strings). The following tables lists all directory parameters: 

| Parameter | Description |  
|-----------|-------------|
| root_folder | Root folder in which the outputs will be stored. Always required |
| bin_folder | Directory containing the C++ binaries. If not provided, the Python driver will attempt to compile the C++ program itself. Only required for long_range_qmc |
| cpp_source_folder | C++ source files for the long range QMC implementation. If neither the bin_folder nor the cpp_source_folder is provided, the preprocessor attempts to locate the Cpp source files |
| julia_path | Path to Julia binaries for the ED simulator. Only needed if the simulator is exact_diagonalization |
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

The above example specifies 4 long range QMC runs for the long range XY model with power law interactions using the deterministic loop SSE algorithm. The 4 parameter sets are: $(N, \alpha) = (4,1.0), (5,2.0), (6,3.0), (7,4.0)$. It should be noted that for certain parameters, multple values would not be meaningful and so providing multiple inputs for those fields results in an exception. These parameters are as follows:  

```
simulator
root_folder
simulation_folder
number_of_threads
output_folder_name
```


## Example Parameter Settings
In this section, we provide example input specifications for ED and QMC. 

### Exact Diagonalization
```
simulator	exact_diagonalization

# System Parameters
hamiltonian_name	XY
boundary	periodic
interaction_range	long_range
interaction_type	power_law
alpha	1.0
spatial_dimension	1d
lattice_type	rectangular
N	4,5,6
J	1.0
theta	0.0

# Output Parameters
root_folder	./results
```

### Long Range Quantum Monte Carlo
```
simulator	long_range_qmc

# System Parameters
hamiltonian_name	XY
boundary	periodic
interaction_range	long_range
interaction_type	power_law
alpha	1.0
spatial_dimension	1d
lattice_type	rectangular
N	11,12,13
J	1.0

# Simulation Parameters
T	0.05
equilibration_steps	100000
simulation_steps	100000
operator_list_update_multiplier	1.05
loop_type	deterministic
number_of_threads	1

# Output Parameters
root_folder	./results
track_spin_configurations	True
write_final_spin_configuration	False

# Initial Configuration Parameters
initial_configuration_index	0
initial_operator_list_size	-1

# Seeds
initial_config_seed	845886
disconnected_spin_flip_seed	255995
off_diagonal_update_seed	893908
metropolis_insert_seed	222288
metropolis_remove_seed	625433
diagonal_update_seed	419814
exit_leg_seed	525338
```

### Nearest Neighbor Quantum Monte Carlo
```
simulator	nearest_neighbor_qmc

# System Parameters
hamiltonian_name	XY
boundary	periodic
interaction_range	nearest_neighbor
spatial_dimension	1d
lattice_type	rectangular
N	11,12,13
J	1.0

# Simulation Parameters
T	0.1
equilibration_steps	100000
simulation_steps	100000
operator_list_update_multiplier	1.1
loop_type	deterministic

# Output Parameters
root_folder	./results
track_spin_configurations	True
write_final_spin_configuration	False

# Initial Configuration Parameters
initial_configuration_index	0
initial_operator_list_size	-1

# Seeds
initial_config_seed	845886
disconnected_spin_flip_seed	255995
off_diagonal_update_seed	893908
metropolis_insert_seed	222288
metropolis_remove_seed	625433
diagonal_update_seed	419814
exit_leg_seed	525338
```

More examples of input file can be found in the repository root under ```examples/input_files```. Other approaches for specifying inputs are exhibited in the example notebook ```examples/long_range_XY.ipynb```.
