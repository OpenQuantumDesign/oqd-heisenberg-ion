---
hide:
    - navigation
---

# Docs

!!! Note
    Heisenberg Ion is still under active development so breaking changes are possible.
    See [Features in development](feat_dev.md) for details.

## Introduction
Welcome to Heisenberg Ion, our open source quantum many-body physics simulator targeting lattice systems that can be natively realized in trapped ion architectures. We currently support Exact Diagonalization (ED) and Quantum Monte Carlo (QMC) Stochastic Series Expansion (SSE) engines for long range anisotropic Heisenberg models in the presence of external fields. 

While Heisenberg Ion is primarily a python package, it uses a C++ engine for large-scale QMC simulations. 
Our ED calculator uses a Julia backend. 
Both of these simulators are wrapped in a lightweight Python preprocessor responsible for driving the required engine, as specified by the user.  

## Language Dependencies
The Heisenberg Ion library has the following language dependencies:  
- Python 3.10+  
- Clang 17+ or GCC 14+ (C++ Compiler)  
- CMake 3.12+  
- Julia 1.12+  

## Python Dependencies
This package also requires a number of Python libraries:  
- Numpy  
- Scipy  
- Matplotlib  
These are installed automatically by pip as part of the Heisenberg Ion package install.  

## Getting Started
To use this package, first clone the Heisenberg Ion repository which can be done via the following terminal command:  
```
git clone https://github.com/OpenQuantumDesign/oqd-heisenberg-ion.git
```

Then, install the package using pip:  
```
pip install heisenberg-ion
```

See the User Guide and the example files in the repository for user-level documentation. A description of the SSE algorithm used in this package can be found in the [Algorithms](algorithms/sse.md) section. 