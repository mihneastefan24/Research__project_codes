# Space Dynamics & Guidance Navigation (SGAN) Suite

This repository contains MATLAB computational implementations and analytical coursework developed for the MSc in Space Engineering at Politecnico di Milano, covering core disciplines in spacecraft guidance, navigation, orbital mechanics, and astrodynamics.

---

## 📌 Project Overview

The primary objective of this suite is to model, optimize, and estimate complex spacecraft dynamics across Earth-bound orbits, multi-body regimes, and interplanetary transfer trajectories:

1. **Spacecraft Dynamics & Optimization:** High-fidelity orbital propagation ($J_2$ perturbations, $N$-body ephemerides), periodic Halo orbit continuation in the Circular Restricted Three-Body Problem (CRTBP), and continuous low-thrust trajectory optimization via Pontryagin’s Maximum Principle (PMP).
2. **Impulsive Guidance & Planetary Defense:** Non-linear multiple-shooting trajectory design for asteroid deflection missions under impulse bounds and precise time-window constraints.
3. **Orbit Determination & State Estimation:** Comparative uncertainty propagation (LinCov, UT, Monte Carlo), Batch Weighted Non-Linear Least Squares OD using ground radar tracking, and Unscented Kalman Filtering (UKF) for relative spacecraft formation navigation.

---

## 📂 Repository Structure

```text
.
├── SGAN_Direction.m            # Frame transformations (ECI, TEME, Topocentric, LVLH, CRTBP), Jacobians & SPICE wrappers
├── SGAN_Earth_Venus_Transfer.m  # PMP-based continuous low-thrust time-optimal interplanetary transfer solver
├── SGAN_Halo_Orbit.m            # CRTBP L1/L2 Lagrange point computation, Halo orbit differential correction & continuation
├── SGAN_Kalman_Filter.m         # LinCov/UT/Monte Carlo uncertainty propagation, Batch OD & UKF relative navigation filter
├── SGAN_Propagator.m            # Multi-body, J2 perturbed, and CRTBP numerical orbit integration engine
├── SGAN_Report_1.pdf            # Detailed technical report for CRTBP Halo orbits, Apophis deflection & low-thrust transfer
├── SGAN_Report_2.pdf            # Detailed technical report for uncertainty propagation, Batch OD & UKF state estimation
├── SGAN_Station_keeping.m       # Multiple-shooting impulsive guidance solver for 99942 Apophis asteroid deflection
└── README.md                    # Repository documentation