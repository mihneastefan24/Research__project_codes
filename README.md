Aerospace & Computational Dynamics Projects

A collection of numerical simulation and analysis codes spanning spacecraft guidance and navigation, granular asteroid dynamics, and numerical methods for PDEs. Developed as part of MSc coursework and research work in astrodynamics and applied computational mechanics.


1. Spacecraft Guidance and Navigation (SGAN)

MATLAB assignments covering orbital dynamics, ephemeris-based propagation, and low-thrust trajectory optimization, using the SPICE/MICE toolkit for planetary ephemerides.

FileDescriptionMARTIN223840_Assign1_Ex1.mCircular Restricted Three-Body Problem (CR3BP): computes the collinear Lagrange points (L1–L3), propagates trajectories with the State Transition Matrix, and generates a halo orbit family around L2 via pseudo-Newton differential correction and numerical continuation.MARTIN223840_Assign1_Ex2.mEphemeris-based n-body propagation of near-Earth asteroid (99942) Apophis relative to Earth, Moon, and Sun; computes close-approach geometry and ground track. Also formulates and solves an interplanetary trajectory optimization problem (Earth → deep-space maneuver → Apophis impact) using fmincon.MARTIN223840_Assign1_Ex3.mIndirect optimal control (Pontryagin's Minimum Principle) for a continuous low-thrust Earth–Venus transfer: solves the two-point boundary value problem for the costates via fsolve, verifies the Hamiltonian, and performs numerical continuation over decreasing thrust levels.

Dependencies: MATLAB, SPICE/MICE toolkit with naif0012.tls, de432s.bsp, gm_de432.tpc, pck00010.tpc kernels (not included), Optimization Toolbox (fsolve, fmincon).


2. Granular Asteroid Dynamics: DEM Simulation & Network Analysis

A two-stage pipeline for studying the internal contact structure of rubble-pile asteroids under self-gravity and rotation.

FileDescriptionParticle_reintroduction.cppDiscrete Element Method (DEM) simulation built on Project Chrono (ChSystemMulticore, NSC contact formulation) with Irrlicht visualization. Loads a granular asteroid configuration (positions, radii, densities) from file, applies mutual gravitational attraction and spin-up, integrates contact dynamics over long timescales, and exports particle positions, velocities, contact forces/pairs, and periodic renders.Network_Code.mPost-processing pipeline that converts the raw contact-pair/force output from the DEM simulation into time-evolving contact networks. Computes graph-theoretic metrics (degree, clustering, closeness/other centralities), information-theoretic measures (Shannon, Von Neumann, and random-walker entropy), and percolation thresholds via generating functions and Karrer message-passing, across different asteroid shapes, densities, and spin rates.

Dependencies: Project Chrono (chrono_multicore, chrono_irrlicht, optionally chrono_vsg), CUDA toolkit, C++17; MATLAB Graph/Network functions.


Note: file paths in these scripts are currently hardcoded to a local Windows environment (C:/Users/mihne/...) and should be parameterized before reuse elsewhere.




3. Numerical Methods for Heat/Advection-Diffusion Equations

Standalone Python scripts implementing and verifying numerical PDE solvers using the Method of Manufactured Solutions (MMS).

FileDescriptionAdvection_problem.py1D Finite Element Method (FEM) solver for the steady and transient advection-diffusion equation: steady-state solve with Dirichlet BCs, transient integration via Backward Euler, a time-step convergence study, and an adaptive (graded) mesh refinement near a boundary layer.Convection_problem.pyHigh-order (4th-order accurate) finite-difference solver for the 1D heat/diffusion equation using RK4 time integration: includes a physical boundary-value case, an MMS-based spatial convergence study, an extension to spatially variable conductivity k(x), and a check of long-time convergence to steady state.

Dependencies: Python 3, numpy, sympy, matplotlib.
