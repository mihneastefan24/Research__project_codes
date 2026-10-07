# Numerical Methods for PDEs & Computational Fluid Dynamics

This repository contains custom numerical PDE solvers implemented in Python from first principles (FEM, FDM, MMS), along with mesh generation scripts and CFD simulation cases using **Gmsh** and **SU2**.

---

## 📌 Project Overview

The primary objective of this coursework and research project is to implement, verify, and analyze high-order numerical schemes for Partial Differential Equations (PDEs) and conduct CFD flow field analysis:

1. **First-Principles Numerical Solvers (Python):** Custom implementations of Finite Element Methods (FEM) and Finite Difference Methods (FDM) to model 1D advection–diffusion, variable-conductivity heat equation dynamics, and implicit time integration.
2. **Verification & Convergence Analysis:** Spatial/temporal error verification using the **Method of Manufactured Solutions (MMS)**, truncation error derivations, $L_2$ norm tracking, and CFL stability analysis.
3. **Mesh Generation & CFD Dynamics (Gmsh & SU2):** Domain generation, boundary layer mesh grading, and external viscous flow simulation over cylinders at low Reynolds numbers ($Re = 20$).

---

## 📂 Repository Structure

```text
.
├── Advection_problem.py             # 1D FEM solver for steady/transient advection-diffusion with graded meshes
├── Backwards_euler.py               # Implicit Backward Euler time-integration module and stability analysis
├── Convection_problem.py            # 4th-order FDM heat solver with RK4, MMS verification & variable conductivity k(x)
├── Create_geometry.m                # MATLAB script for parametric domain definitions and geometry creation
├── Cylinder_re20.cfg                # SU2 configuration file for laminar cylinder flow simulation (Re = 20)
├── NLAB3_Geometry.msh               # Baseline Gmsh mesh file for asymmetric channel/cylinder geometry
├── NLAB3_Geometry_14000.msh         # Fine Gmsh mesh (~14,000 elements) for grid refinement studies
├── NLAB3_Geometry_corse.msh         # Coarse Gmsh mesh for baseline grid testing
├── NLAB3_Geometry_corse.su2         # Converted coarse SU2 format mesh file
├── NLAB3_Geometry_medium.su2        # Converted medium-density SU2 format mesh file
├── NLAB3_Geometry_refined.su2       # Converted refined SU2 mesh file
├── cylinder_re20_BL.su2             # Boundary-layer resolved SU2 grid for cylinder flow analysis
└── README.md                        # Repository documentation