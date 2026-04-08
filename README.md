# PFAS Removal via Ion Exchange (HSDM & PSDM Models)

## Overview

This repository provides a framework for modeling **PFAS removal** using ion exchange and adsorption-based processes. It combines **mechanistic models** with **inverse parameter estimation** to match experimental breakthrough data and extract key transport and equilibrium parameters.

---

## Objectives

The main goals of this repository are:

* Fit experimental breakthrough data (**Ct/C0 vs time**) from laboratory ion exchange systems
* Estimate and optimize unknown model parameters, including:

  * Surface diffusion coefficient (**Ds**)
  * Ion exchange selectivity coefficient (**Kxc**)
* Quantify PFAS adsorption behavior using equilibrium models such as the **Freundlich isotherm**
* Bridge experimental data with mechanistic modeling for improved system understanding

---

## Model Framework

This work builds on the **EPA Water Treatment Models**, which provide the foundation for simulating transport and adsorption processes in fixed-bed systems.

### Included Model Groups

* **Ion Exchange Models**

  * Homogeneous Surface Diffusion Model (**HSDM**)
  * Pore Surface Diffusion Model (**PSDM**) for ion exchange

* **Adsorption Models**

  * PSDM / AdDesign_GAC for **granular activated carbon (GAC)** systems

---

## HSDM (Ion Exchange)

The HSDM implementation is extended for **multispecies inverse modeling**:

* Simulates ion exchange in fixed-bed columns
* Uses experimental breakthrough data for calibration
* Optimizes:

  * **Ds** (surface diffusion coefficient)
  * **Kxc** (selectivity coefficient)
* Minimizes fitting error (typically **RMSE**) between experimental and simulated breakthrough curves

---

## PSDM (Adsorption)

The PSDM implementation is used for **single-species adsorption systems**:

* Applicable to PFAS removal using GAC
* Models pore diffusion, surface diffusion, and equilibrium behavior
* Supports parameter estimation for adsorption and transport processes

---

## Experimental Data

Experimental datasets are included and formatted to match the **EPA model input structure**, allowing direct use without preprocessing.

These datasets typically include:

* Parameter sheets
* Ion/compound properties
* Influent concentration profiles
* Breakthrough data (**Ct/C0 vs time**)

---

## Inverse Modeling Workflow

The inverse fitting approach follows two main steps:

### 1. Breakthrough Curve Fitting

* The model is evaluated at the same time points as the experimental data
* Parameters (**Ds**, **Kxc**) are optimized to minimize error between:

  * Experimental: Ct/C0
  * Model prediction

### 2. Post-Processing

Once optimal parameters are obtained:

* **Ce** (equilibrium liquid concentration) is calculated
* **qe** (adsorbed concentration, mg/g) is determined

This separation improves computational efficiency and fitting stability.

---

## Repository Structure

```
.
├── HSDM_fit/              # Inverse modeling scripts for ion exchange
├── PSDM_fit/              # Inverse modeling scripts for adsorption
├── EPA_models/            # Original EPA model implementation
├── data/                  # Experimental datasets (Excel format)
├── notebooks/             # Jupyter notebooks for analysis and fitting
└── README.md
```

---

## Requirements

Typical dependencies include:

* Python 3.x
* numpy
* pandas
* scipy
* matplotlib
* openpyxl

---

## Usage

### 1. Prepare input data

Ensure the Excel file follows the EPA format:

* `params`
* `ions`
* `Cin`
* breakthrough data sheet

### 2. Run inverse fitting

Example (in Jupyter):

```python
r = fit_pfas("PFPeA")
print(r)
```

### 3. Plot results

```python
plot_fit(r)
```

---

## Notes

* The models use **orthogonal collocation** for solving transport equations
* Runtime depends on discretization parameters (e.g., `nr`, `nz`)
* For faster fitting:

  * reduce collocation points
  * downsample experimental data
  * limit optimizer iterations

## License

Specify your license here (e.g., MIT, GPL, etc.)

