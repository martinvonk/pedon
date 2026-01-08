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
 - name: Federal Institute for Geosciences and Natural Resources (BGR), Hannover, Lower Saxony, Germany
   index: 3

date: 8 January 2025
bibliography: paper.bib

---

# Summary
Pedon is a Python package designed to describe and analyse unsaturated soil properties. The package offers an object-oriented modelling framework, complemented by tools for parameter retrieval from soil databases, implementation of pedotransfer functions, and optimisation routines for parameter fitting.

# Statement of need
Researchers and engineers working with unsaturated soils often need estimations of their soil parameters for their groundwater models. Pedon solves this by providing a modern Python toolkit that brings together commonly used soil models, parameter databases, pedotransfer functions, and fitting routines. This makes soil analysis faster, more reproducible, and easier to integrate into existing modelling pipelines.

# Soil Models
Pedon can be installed via `pypi` using `pip install pedon` and imported using `import pedon as pe`. Different soil models are available in Pedon. Each soil model offers a mathematical description of the soil water retention curve (SWRC) and the hydraulic conductivity function (HCF). By default the following soil models are available:
- Mualem-van Genuchten [@genuchten_mualem_1980]: `pe.Genuchten`
- Brooks-Corey [@brooks_corey_1964]: `pe.Brooks`
- Gardner-Kozeny  [@gardner_params_1970;brutsaert_kozeny_1967;bakker_gardner_2009;mathias_gardner_2006]: `pe.Gardner`
- Fredlund-Xing  [@fredlund_xing_1994]: `pe.Fredlund`
- Combination of the van Genuchten SWRC and Brooks-Corey HCF [@fuentes_burdine_1992;@panday_mfusgt_2025]: `pe.Panday`
- Gardner-Rucker [@rucker_gardner_2005]: `pe.GardnerRucker`

## Available datasets
There are a few datasets available with Brooks-Corey and Mualem-van Genuchten parameters for different soils. These parameters are obtained from a few databases:
  - Average values for selected soil water retention and hydraulic conductivity parameters for 12 major soil
  textural groups according to @carsel_dataset_1988. This dataset is also used in the popular software HYDRUS [@simunek_hydrus1d_2009] that simulates water, heat, and solute movement in one-, two- and three-dimensional variably saturated media.
  - The Staring series (Staringreeks in Dutch) is a database of soil water retention curves and hydraulic conductivity functions in the Netherlands. It contains both a description of top soils and bottom soils based of 100's of samples. These samples were processed to obtain the Mualem-van Genuchten soil models [@genuchten_mualem_1980] curves.
  - Dataset obtained from the VS2D software [@healy_vs2d_1990] containing both Brooks-Corey and Mualem-vanGenuchten parameters.

## Soil model parameter estimation
### Pedotransfer functions
Pedotranfer functions are equations that predict the soil model parameters based on easy to measure soil properties. [@bouma_pedotransfer_1989]. Different pedotransfer functions are avalaible from @wosten_pedotranfer_1999, @wosten_staringreeks_2001, @cosby_pedotransfer_1984, @cooper_pedotransfer_2021.
### Databases
Databases that are based on statistical relations such as Rosetta @schaap_rosetta_2001 and HYPAGS @peche_hypags_2024

## Estimation from sample measurements
Same routine as RETC @genuchten_retc_1991.

### Going from one soil model to the other
<!-- ... -->

# Acknowledgements

# References