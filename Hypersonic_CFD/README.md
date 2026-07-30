Numerical Methods for PDEs

This directory contains Python implementations of numerical methods for solving and analysing partial differential equations.

The projects focus on implementing numerical schemes from first principles and assessing their accuracy, convergence, stability, and mesh sensitivity.

Topics

Finite Element Methods (FEM)

Finite Difference Methods (FDM)

Advection–diffusion equations

Heat/diffusion equations

Time integration

Method of Manufactured Solutions (MMS)

Convergence analysis

Mesh refinement

Files

Advection_problem.py

Finite Element Method implementation for one-dimensional steady and transient advection–diffusion problems.

The implementation includes:

Finite element discretisation

Global matrix assembly

Dirichlet boundary conditions

Steady-state solution

Transient solution

Backward Euler time integration

Time-step convergence analysis

Mesh refinement

Graded mesh refinement near boundary layers

Backwards_euler.py

Implementation and analysis of the Backward Euler time-integration scheme for transient numerical problems.

Convection_problem.py

Finite-difference solution of the one-dimensional heat/diffusion equation.

The implementation includes:

Fourth-order spatial discretisation

RK4 time integration

Manufactured-solution verification

Spatial convergence analysis

Spatially varying thermal conductivity

Long-time convergence towards the steady-state solution

Verification & Analysis

The projects use analytical or manufactured solutions where applicable to assess numerical accuracy and convergence.

The analysis focuses on:

Spatial convergence

Temporal convergence

Mesh sensitivity

Numerical stability

Comparison between numerical and reference solutions

Requirements

The Python implementations use standard scientific-computing packages including:

Python 3

NumPy

SciPy

SymPy

Matplotlib

Author

Mihnea Stefan Martin, VKI, Politecnico di Milano