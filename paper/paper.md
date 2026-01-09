---
title: 'pedon: A Python package for analyzing unsaturated soil properties'
tags:
  - hydrology
  - unsaturated zone
  - vadose zone
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
Pedon is a Python package designed to describe and analyse unsaturated soil properties. The package offers an object-oriented modelling framework, complemented by tools for parameter retrieval from soil databases, implementation of pedotransfer functions, and optimisation routines for parameter fitting. It leverages Python’s object-oriented strengths and its well-maintained scientific ecosystem, including NumPy [@numpy_article_2020], SciPy [@scipy_paper_2020], Matplotlib [@matplotlib_paper_2007], and Pandas [@pandas_software_2020; @pandas_paper_2010].

# Statement of need
Researchers and engineers working with unsaturated soils often need estimations of their soil parameters for their groundwater models. Pedon addresses this need by providing a modern Python toolkit that brings together commonly used soil models, parameter databases, pedotransfer functions, and fitting routines. This makes soil analysis faster, more reproducible, and easier to integrate into existing modelling pipelines.

# Soil models
Pedon can be installed via `pypi` using `pip install pedon` and imported using `import pedon as pe`. Different soil models are available in Pedon. A soil model is a parametric description of the soil water retention curve (SWRC) and the hydraulic conductivity function (HCF), linking soil water content and flow to pressure head or saturation for use in unsaturated flow simulations. By default the following soil models are available:
- Mualem-van Genuchten [@genuchten_mualem_1980]: `pe.Genuchten`
- Brooks-Corey [@brooks_corey_1964]: `pe.Brooks`
- Combination of the van Genuchten SWRC and Brooks-Corey HCF [@fuentes_burdine_1992; @panday_mfusgt_2025]: `pe.Panday`
- Fredlund-Xing [@fredlund_xing_1994]: `pe.Fredlund`
- Gardner-Kozeny [@gardner_params_1970; @brutsaert_kozeny_1967; @bakker_gardner_2009; @mathias_gardner_2006]: `pe.Gardner`
- Gardner-Rucker [@rucker_gardner_2005]: `pe.Rucker`

The soil models are implemented as Python classes, providing a clear structure in which model-specific methods can be consistently defined and extended.
```python
import numpy as np
import pedon as pe

mg = pe.Genuchten(
    k_s=106.1,  # saturated conductivity (cm/d)
    theta_r=0.065,  # residual water content (-)
    theta_s=0.41,  # saturated water content (-)
    alpha=0.075,  # shape parameter (1/cm)
    n=1.89,  # shape parameter (-)
)

h = np.logspace(-2, 6, 9)  # pressure head (cm)
theta = mg.theta(h)  # water content (-) at pressure head values
k = mg.k(h)  # hydraulic conductivity (cm/d) at pressure head values
```

## Parameter datasets
In Pedon there is a dataset available with Brooks-Corey and Mualem-van Genuchten parameters for different soils. These parameters are obtained from a few databases:
- Average values for selected soil water retention and hydraulic conductivity parameters for 12 major soil textural groups as defined by @carsel_dataset_1988. This dataset is also used in the popular software HYDRUS [@simunek_hydrus1d_2009] that simulates water, heat, and solute movement in one-, two- and three-dimensional variably saturated media.
- The Staring series (Staringreeks in Dutch) is a database of soil water retention curves and hydraulic conductivity functions in the Netherlands [@wosten_staringreeks_2001; @heinen_staringreeks_2020]. It contains both a description of top soils and bottom soils based on hundreds of samples. These samples were processed to obtain the Mualem-van Genuchten soil models [@genuchten_mualem_1980].
- Dataset obtained from the VS2D software [@healy_vs2d_1990] containing both Brooks-Corey and Mualem-van Genuchten parameters.

## Parameter estimation
Estimates of unsaturated soil hydraulic parameters are required for modeling water flow in the unsaturated zone, yet direct measurements are often scarce, expensive, or incomplete. Pedon therefore provides multiple, complementary approaches to obtain soil model parameters from available measurements.

### Databases and pedotransfer functions
When direct measurements are unavailable, soil hydraulic parameters can be estimated using pedotransfer functions, which relate easily measured soil properties (e.g. texture, bulk density, organic matter) to soil model parameters [@bouma_pedotransfer_1989]. Pedon implements some pedotransfer functions from the literature, including those of @wosten_pedotransfer_1999, @wosten_staringreeks_2001, @cosby_pedotransfer_1984, and @cooper_pedotransfer_2021. In addition, Pedon provides access to soil model parameter databases such as Rosetta [@schaap_rosetta_2001] and HYPAGS [@peche_hypags_2024], the latter of which enables parameter estimation based solely on saturated hydraulic conductivity.

## Estimation from sample measurements
When laboratory measurements of soil water retention and/or unsaturated hydraulic conductivity are available, Pedon supports direct parameter estimation through inverse modeling. Soil model parameters are obtained by fitting analytical retention and conductivity functions to observed data using nonlinear least-squares optimization available through SciPy [@scipy_paper_2020]. This approach minimizes the mismatch between measured and simulated values and follows the existing methodology implemented in the RETC software [@genuchten_retc_1991].

### Soil model conversion
The same fitting framework can be used to translate between different soil hydraulic models. Retention and conductivity curves generated from one model can be sampled over a range of pressure heads and refitted using another model formulation, facilitating model comparison and integration with external simulation tools.
```python
bc = pe.SoilSample(h=h, theta=theta, k=k).fit(pe.Brooks)
```
![Resulting fit of the Brooks-Corey soil model on the Mualem-van Genuchten soil model \label{fig:swrc_fit}](figures/swrc_fit.png)


# References