# Modelling and Simulation of Aerospace Systems (MSAS)

This repository contains the computational implementations, numerical solvers, and Simulink multi-domain physical system models for the **Modelling and Simulation of Aerospace Systems (MSAS)** course at Politecnico di Milano (AY 2024–2025).

---

## 📌 Overall Project Overview

The primary objective of this project is to develop, validate, and compare advanced numerical algorithms and physical system modeling paradigms applied to aerospace engineering dynamics. 

The project spans two major domains:
1. **Custom Numerical Algorithms:** Hand-crafted solvers developed in MATLAB to handle non-linear root-finding, fixed and adaptive step ODE integration, stability region mapping, stiff differential equations, and zero-crossing event detection for hybrid dynamics.
2. **Multi-Domain Physical Systems:** Signal-based (causal) and physical network (acausal) modeling implemented in Simulink/Simscape to evaluate closed-loop aerospace dynamics.

---

## 📂 Repository Structure

```text
.
├── Electric_circuit.slx       # Simulink model for electrical circuit and dynamic transient analysis
├── Modelling.slx              # Simulink model for dynamic physical systems and control loops
├── Numerical_Methods_1.pdf    # Problem statement & theoretical report for Assignment 1
├── Numerical_Methods_2.pdf    # Problem statement & theoretical report for Assignment 2
├── Numerical_model_1.m        # MATLAB solver suite: Root-finding & Kinematic linkage equilibrium
├── Numerical_model_2.m        # MATLAB solver suite: ODE integration, stability & variable-mass systems
├── Numerical_model_3.m        # MATLAB solver suite: Stiff systems (IEX4) & hybrid bouncing dynamics
└── README.md                  # Main repository documentation