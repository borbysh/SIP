# SIP Solver - Strongly Implicit Procedure for PDEs
This is a modified version of the numerical heat equation solver project by Jussi Enkovaara: https://repository.prace-ri.eu/git/CodeVault/training-material/parallel-programming/OpenMP/-/tree/master/heat-equation/c/solution?ref_type=heads
A C implementation of the Strongly Implicit Procedure solver for solving partial differential equations, specifically the 2D heat equation.

## Features

- **SIP Solver**: Implementation of Stone's Strongly Implicit Procedure
- **Heat Equation Solver**: Finite difference discretization with Dirichlet boundary conditions
- **Analytical Solution**: Green's function approach for validation
- **Error Metrics**: Comprehensive error analysis tools (L2, L-infinity norms, RMSE, etc.)
- **VTK Output**: Visualization support for numerical results

## Building

```bash
make clean
make all
