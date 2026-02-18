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

For general user-level documentation, see the User Guide. For more detailed dev docs, see the Developer Guide. A description of the SSE algorithmic approaches used in this package can be found in the Algorithms section. To contribute, please visit the Community page and the Github repository for a list of open issues. 