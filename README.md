## Overview
This repository contains Stochastic Series Expansion (SSE) Quantum Monte Carlo (QMC) and Exact Diagonalization (ED) simulators targeting the long range XXZ model in the presence of external fields. 

## Dependencies

### Runtime Dependencies
- Python 3.10+  
- Julia 1.12+   
- CMake 3.12+  
- Clang 17+ or GCC 14+ (C++ Compiler) 
CMake and a C++ compiler are required for building the long range QMC backend. Julia is needed for our ED implementation

### Python Packages
This package also requires a number of Python libraries:  
- Numpy  
- Scipy  
- Matplotlib  
These are installed automatically by pip as part of the Heisenberg Ion package install

### C++ Dependencies
- spdlog 
This is fetched and compiled automatically by CMake while building the C++ source if required

## Getting Started
To use this package, first clone the Heisenberg Ion repository:  
```
git clone https://github.com/OpenQuantumDesign/oqd-heisenberg-ion.git
```

Then, install the package:  
```
pip install heisenberg-ion
```

### Documentation 
Documentation is developed using MkDocs. To install documentation-related dependencies, use:

```
pip install -e ".[docs]"
```

The documentation server can be deployed locally using: 

```
cp -r examples/ docs/examples/
mkdocs serve
```