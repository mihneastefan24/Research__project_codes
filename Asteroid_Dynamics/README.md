Asteroid Dynamics & Network Analysis

This directory contains computational tools developed for the study of rubble-pile asteroid dynamics and internal structural evolution.

The work combines Discrete Element Method (DEM) simulations, gravitational and contact dynamics, and graph/network analysis to investigate the evolution of granular asteroid systems.

Overview

The computational workflow consists of two main stages:

Generate and simulate granular asteroid configurations using C++ and Project Chrono.

Convert particle-contact information into time-dependent networks and analyse their structural properties using MATLAB.

Files

Creation_Initial_Aggregate_v1.cpp

C++ code for generating initial granular configurations used in asteroid dynamics simulations.

Particle_reintroduction.cpp

DEM simulation implemented in C++ using Project Chrono.

The simulation includes:

Granular asteroid configurations

Particle contact dynamics

Mutual gravitational attraction

Rotational spin-up

Time integration of the particle system

Particle positions and velocities

Contact forces and particle-pair information

Data export for post-processing

Network_Code.m

MATLAB post-processing pipeline for transforming DEM contact information into time-dependent contact networks.

The analysis includes:

Degree distributions

Clustering coefficients

Network centrality measures

Shannon entropy

Von Neumann entropy

Random-walker entropy

Percolation thresholds

Generating-function methods

Message-passing approaches

These quantities are used to investigate how the internal contact structure of rubble-pile asteroids evolves during their dynamical evolution and disaggregation.

Methods & Technologies

C++

MATLAB

Project Chrono

Discrete Element Method (DEM)

CUDA

Graph theory

Network science

Numerical simulation

High-performance computing

Research Context

The codes form part of research into the dynamical and structural behaviour of rubble-pile asteroids, including the application of network-based representations to granular systems.

Author

Mihnea Stefan MartinMSc Space Engineering, Politecnico di Milano