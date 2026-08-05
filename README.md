# Chapter 3: Hamiltonian Optimal Control Method

This branch contains the numerical implementation and plotting script for the **Hamiltonian Optimal Control Method** presented in **Chapter 3** of the PhD thesis.

---

## 1. Reproducing Thesis Figures

* **Computation (MATLAB):** Run `corelatednoise.m` to generate simulation data under different noise conditions. The script allows adjusting the noise correlation prior to the degenerate noise structure, producing both **correlated** and **uncorrelated** datasets.
* **Data Storage:** Output data is saved in the `Hamiltonian/` folder.
* **Plotting (Python):** Run `MAM.py` to generate the figures:
  * **Figure 3.6:** Comparison of transitions between correlated and uncorrelated noise cases.
  * **Figure 3.7:** Comparison of transition paths across different ship natural frequencies under correlated noise.

---

## 2. File Structure Overview

* **`Hamiltonian/`**: Folder containing output data for Figures 3.6 and 3.7 (both correlated and uncorrelated cases).
* **`corelatednoise.m`**: MATLAB script performing Hamiltonian Optimal Control computations with adjustable noise correlation parameters.
* **`MAM.py`**: Python script to visualise and plot Figures 3.6 and 3.7.
* **`LICENSE`**: License file.
* **`README.md`**: Documentation file.

