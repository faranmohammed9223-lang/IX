# PFAS Removal via Ion Exchange (HSDM & PSDM Models)

The purpose of this repository is to:
1. Fit breakthrough data obtained from laboratory ion exchange experiments using inverse modeling with mechanistic models
2. Optimize unknown model parameters to achieve the best agreement between experimental and simulated results
3. Understand and quantify PFAS adsorption equilibrium using the Freundlich isotherm or other isotherms 

## Repository Overview

This repository contains code for both forward simulation and inverse fitting of ion exchange and diffusion-based adsorption models.

### EPA Code
The EPA code provides the original model framework for direct simulation, where the user supplies all required input parameters to generate breakthrough curves. It includes two major model groups:

- **Ion Exchange Models**: HSDM and PSDM
- **PSDM / AdDesign_GAC**: Pore Surface Diffusion Model for granular activated carbon systems

### HSDM
The HSDM folder is used for **multispecies ion exchange inverse fitting** based on experimental breakthrough data. The purpose is to estimate unknown model parameters, such as **Ds** and **Kxc**, by minimizing the fitting error between simulated and observed breakthrough curves.

### PSDM
The PSDM section is used for **inverse fitting of the pore surface diffusion model** for a **single-species system**, mainly for adsorption and equilibrium analysis.

### Experimental Data
Some experimental datasets are included in the inverse fitting files and are already organized in a format that satisfies the input layout required by the EPA code. This makes them directly readable by the original model structure, including the expected sheets, parameter names, and data arrangement.

### Relationship Between the Codes
The inverse fitting scripts are developed on top of the EPA forward models. The EPA code handles the mechanistic simulation, while the inverse fitting workflow adjusts unknown parameters to obtain the best agreement with experimental breakthrough data.
