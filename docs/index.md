---
hide:
    - navigation
---

# Docs

!!! Note
    Heisenberg Ion is still under active development so breaking changes are possible.
    See [Features in Development](feat_dev.md) for details.

## Introduction
Welcome to Heisenberg Ion, our open source quantum many-body physics simulator targeting lattice systems that can be natively realized in trapped ion architectures. We currently support Exact Diagonalization (ED) and Quantum Monte Carlo (QMC) Stochastic Series Expansion (SSE) engines for long range anisotropic Heisenberg models in the presence of external fields. 

While Heisenberg Ion is primarily a python package, it uses a C++ engine for large-scale QMC simulations. 
Our ED calculator uses a Julia backend. 
Both of these simulators are wrapped in a lightweight Python preprocessor responsible for driving the required engine.  

## Dependencies

### Runtime Dependencies
- Python 3.10+  
- Julia 1.12+   
- CMake 3.12+  
- Clang 17+ or GCC 14+ (C++ Compiler) 
CMake and a C++ compiler are required for building the long range QMC source code 

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
pip install -e ".[docs]
```

See the [User Guide](user_guide/overview.md) and the example files in the repository for user-level documentation. A description of the SSE algorithm used in this package can be found in the [Algorithms](algorithms/sse.md) section. 