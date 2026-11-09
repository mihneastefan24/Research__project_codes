# Asteroid Dynamics & Network Analysis

This directory contains computational tools developed for the study of rubble-pile asteroid dynamics and internal structural evolution. The work combines Discrete Element Method (DEM) simulations, gravitational and contact dynamics, and graph/network analysis to investigate the evolution of granular asteroid systems.

## 📌 Overview

The computational workflow consists of two main stages:
1. Generate and simulate granular asteroid configurations using C++ and Project Chrono.
2. Convert particle-contact information into time-dependent networks and analyse their structural properties using MATLAB.

## 📂 Key Files

- **`Creation_Initial_Aggregate_v1.cpp`**  
  C++ code for generating initial granular configurations used in asteroid dynamics simulations.

- **`Particle_reintroduction.cpp`**  
  DEM simulation implemented in C++ using Project Chrono. Features granular configurations, particle contact dynamics, mutual gravitational attraction, rotational spin-up, and data export for post-processing.

- **`Network_Code.m`**  
  MATLAB post-processing pipeline for transforming DEM contact information into time-dependent contact networks. Analyzes degree distributions, clustering coefficients, centrality measures, entropies, percolation thresholds, and network evolution.

## 🛠️ Methods & Technologies

- C++, MATLAB, Project Chrono
- Discrete Element Method (DEM)
- Graph Theory & Network Science
- Numerical Simulation & High-Performance Computing

## 🔬 Research Context

These codes form part of research into the dynamical and structural behavior of rubble-pile asteroids, including the application of network-based representations to granular systems.

---
**Author:** Mihnea Stefan Martin  
*M.Sc. Space Engineering, Politecnico di Milano*
