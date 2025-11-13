---
title: 'pedon: A Python package for unsaturated soil properties'
tags:
  - hydrology
  - unsaturated zone
  - Python
authors:
  - name: Martin A. Vonk
    orcid: 0009-0007-3528-2991
    affiliation: "1, 2"
  - name: Aaron Peche
    orcid: 0000-0003-3761-2441
    affiliation: "3"

affiliations:
 - name: Artesia B.V., Schoonhoven, South Holland, The Netherlands
   index: 1
 - name: Department of Water Management, Faculty of Civil Engineering and Geosciences, Delft University of Technology, Delft, South Holland, The Netherlands
   index: 2
 - name: Federal Institute for Geosciences and Natural Resources, Hannover, Lower Saxony, Germany
   index: 3

date: 12 november 2025
bibliography: paper.bib

---

# Summary
Pedon is a Python package designed to describe and analyse unsaturated soil properties. The package offers an object-oriented modelling framework, complemented by tools for parameter retrieval from soil databases, implementation of pedotransfer functions, and optimisation routines for parameter fitting.

# Statement of need
Researchers and engineers working with unsaturated soils often need estimations of their soil parameters for their groundwater models. Pedon solves this by providing a modern, Python-native toolkit that brings together commonly used soil models, parameter databases, pedotransfer functions, and fitting routines. This makes soil analysis faster, more reproducible, and easier to integrate into existing modelling pipelines.

# Soil Models
Pedon can be installed via `pypi` using `pip install pedon` and imported using `import pedon as pe`. Different soil models are available in Pedon:
- Mualem-van Genuchten: `Genuchten` [@genuchten_mualem_1980]
- Brooks-Corey: `pe.Brooks` [@brooks_corey_1964]
- Gardner: `pe.Gardner` []
- Fredlund-Xing `pe.Fredlund` []

## Available datasets
Small description of HYDRUS, Staring and VS2D datasets

## Estimation of these curves via pedotranfser functions
<!-- List of pedotransfer functions and their effects on the estimation of the soil water retention curves -->

## Estimation from databases via HYPAGS and Rosetta
<!-- HYPAGS and Rosetta description -->

## Estimation from sample measurements
<!-- RETC replacement -->

### Going from one soil model to the other
<!-- ... -->

# Acknowledgements

# References