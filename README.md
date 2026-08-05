# Chapter 3: Benchmark Methods & Numerical Simulations

This branch contains the numerical implementation and plotting scripts for the benchmark examples presented in **Chapter 3** of the PhD thesis.

---

## 1. Reproducing Thesis Figures

**Figure 3.1: Minimum Action Method (MAM) Results**
* **Computation:** `four_dimensionMAM.m` (MATLAB)
* **Computed Data:** Saved in the `MAM/` folder.
* **Plotting:** Run `MAM.py` (Python) to generate **Figure 3.1** directly using the computed data in `MAM/`.

**Figure 3.2: Stochastic Differential Equation (SDE) Simulations**
* **Computation:** `four_dimensionSDE.m` (MATLAB)
* **Note:** Since this is a stochastic simulation, each run produces different results under random fluctuations. Therefore, the output data is not stored in the branch. Run `four_dimensionSDE.m` directly in MATLAB to observe the stochastic transitions.

**Figure 3.3 & 3.4: Comparison between MAM and Simulated Trajectories**
* **Main Script:** `main.m` (MATLAB)
* **Functionality:** Executes the comparative study between the Minimum Action Method (MAM) and the Monte Carlo simulated results for different noise strength \epsilon.
* **Helper Functions Called:
  * `MAM_4D.m`
  * `SDE_epsilon.m`
  * `SDE_trajectory.m`

---

## 2. File Structure Overview

├── MAM/                   # Folder containing computed data for Figure 3.1
├── MAM.py                 # Python script to plot Figure 3.1
├── four_dimensionMAM.m    # MATLAB script computing MAM for Figures 3.1
├── four_dimensionSDE.m    # MATLAB script for stochastic simulation (Figure 3.2)
├── main.m                 # Main MATLAB script to generate data/plots for Figure 3.4
├── MAM_4D.m               # Helper function: 4D MAM solver
├── SDE_epsilon.m          # Helper function: SDE solver under different noise parameter epsilon
├── SDE_trajectory.m       # Helper function: Stochastic simulation
├── LICENSE                # License file
└── README.md              # Documentation file
