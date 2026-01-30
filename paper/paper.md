---
title: 'pedon: A Python package for analyzing unsaturated soil hydraulic properties'
tags:
  - hydrology
  - unsaturated zone
  - vadose zone
  - groundwater modeling
  - soil hydraulic properties
  - Python
authors:
  - name: Martin A. Vonk
    orcid: 0009-0007-3528-2991
    affiliation: "1, 2"
  - name: Aaron Peche
    orcid: 0000-0003-3761-2441
    affiliation: "3"

affiliations:
 - name: Department of Water Management, Faculty of Civil Engineering and Geosciences, Delft University of Technology, Delft, South Holland, The Netherlands
   index: 1
 - name: Artesia B.V., Schoonhoven, South Holland, The Netherlands
   index: 2
 - name: Federal Institute for Geosciences and Natural Resources (BGR), Hannover, Lower Saxony, Germany
   index: 3

date: 12 January 2026
bibliography: paper.bib

---

# Summary
`pedon` is a Python package for describing and analyzing unsaturated soil hydraulic properties. It provides a framework for soil hydraulic models, along with tools for retrieving parameters from soil databases, applying pedotransfer functions, and fitting soil hydraulic parameters to measurements.

# Statement of need
Researchers and engineers working with unsaturated soils need estimates of soil parameters for variably saturated groundwater flow models. `pedon` provides a Python toolkit that brings together soil hydraulic models, parameter databases, pedotransfer functions, and fitting routines, making soil analysis faster, more reproducible, and easier to integrate into existing groundwater modeling workflows.

## State of the field
pedon is a contribution to the groundwater modeling field that integrates, combining and extending on  existing tools, for example by incorporating parameter estimation algorithms from HYPAGS [@peche_hypags_2024]. Its primary novelty lies in the object-oriented representation of soil hydraulic models, which enables the definition of custom soil models and facilitates efficient comparison across formulations using the integrated parameter datasets and pedotransfer functions. The main functional overlap with existing software lies in the fitting of soil hydraulic functions using least-squares optimization, as implemented in tools such as `unsatfit` [@seki_unsatfit_2023], `PySWR` [@memari_pyswr_2021], and `RETC` [@genuchten_retc_1991]. However, these tools are primarily standalone fitting utilities and do not provide an integrated framework that combines (custom) soil models, parameter datasets, and estimation methods within reproducible workflows.

# Soil hydraulic models
A soil hydraulic model is a parametric description of soil hydraulic functions: the soil water retention curve (SWRC) and the unsaturated hydraulic conductivity function (HCF). These relate soil water content and flow to pressure head and vice versa for use in variably saturated groundwater flow models. At this time, `pedon` provides the following soil models:

- `pedon.Genuchten`: Mualem-van Genuchten [@genuchten_mualem_1980]
- `pedon.Brooks`: Brooks-Corey [@brooks_corey_1964]
- `pedon.Panday`: Mualem-van Genuchten SWRC and Brooks-Corey HCF [@fuentes_burdine_1992; @panday_mfusgt_2025]
- `pedon.Fredlund`: Fredlund-Xing [@fredlund_xing_1994]
- `pedon.Haverkamp`: Haverkamp [@haverkamp_model_1977]
- `pedon.Gardner`: Gardner(-Kozeny) [@gardner_model_1958; @brutsaert_kozeny_1967; @bakker_gardner_2009; @mathias_gardner_2006]
- `pedon.Rucker`: Gardner-Rucker [@rucker_gardner_2005]
- `pedon.GenuchtenGardner`: Mualem-van Genuchten SWRC and Gardner HCF [@genuchten_mualem_1980; @gardner_model_1958]

## Software design
The soil models are implemented as Python classes with model-specific methods for evaluating the SWRC and HCF. For example, the Mualem–van Genuchten model can be used as follows:

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

The object-oriented design and duck typing provides a clear and consistent structure in which users can define custom soil model classes. Additionally, `pedon` only depends on well-maintained packages in the Python scientific ecosystem such as NumPy [@numpy_article_2020], SciPy [@scipy_paper_2020], Matplotlib [@matplotlib_paper_2007], and Pandas [@pandas_software_2020; @pandas_paper_2010].

# Soil hydraulic parameters
Soil hydraulic parameters depend on soil type and determine the shape of a soil model’s SWRC and HCF. They are rarely measured directly and are usually derived from reference datasets, empirical relationships, or laboratory measurements. `pedon` links these parameters to soil models and provides a framework to obtain them from existing datasets, easily measured soil properties, and direct measurements of soil water content and hydraulic conductivity.

## Parameter datasets
`pedon` includes a dataset of Brooks–Corey and Mualem–van Genuchten parameters for a wide range of soils. At present, this dataset is compiled from three established soil hydraulic parameter databases:

- Average parameter values for twelve major soil textural groups defined by @carsel_dataset_1988, also used in the HYDRUS software for variably saturated flow modeling [@simunek_hydrus1d_2009];
- A dataset from the VS2D software [@healy_vs2d_1990] containing both Brooks–Corey and Mualem–van Genuchten parameters;
- The Staring series from the Netherlands [@wosten_staringreeks_2001; @heinen_staringreeks_2020; @heinen_bofek_2022], which describes soils using the Mualem–van Genuchten model based on hundreds of processed samples [@genuchten_mualem_1980; @wosten_texture_1988].

The databases can be called via the following code:
```python
hydrus = pe.Soil("Sand").from_name(pe.Genuchten, source="HYDRUS")
vs2d = pe.Soil("Sand").from_name(pe.Brooks, source="VS2D")
staring = pe.Soil("B01").from_name(pe.Genuchten, source="Staring_2018")
```

## Parameter estimation
`pedon` provides two approaches for obtaining soil hydraulic parameters from soil data. The first uses pedotransfer functions based on easily measured soil properties. The second relies on direct measurements of soil water content and hydraulic conductivity.

### Pedotransfer functions
Pedotransfer functions relate easily measured soil properties (e.g. sand, silt, clay or organic matter content and bulk density) to soil hydraulic parameters [@bouma_pedotransfer_1989]. `pedon` implements functions from the literature, including those of @wosten_pedotransfer_1999, @wosten_staringreeks_2001, @cosby_pedotransfer_1984, and @cooper_pedotransfer_2021. It also provides access to parameter databases such as Rosetta [@schaap_rosetta_2001] and HYPAGS [@peche_hypags_2024], the latter enabling estimation from a single value of saturated hydraulic conductivity or representative grain diameters.

```python
# Estimate parameters using Cosby's pedotransfer function
sand_p = 40.0  # sand (%)
clay_p = 10.0  # clay (%)
cosby: pe.Brooks = pe.SoilSample(sand_p=sand_p, clay_p=clay_p).cosby()

# Estimate parameters from saturated conductivity via HYPAGS
ks = 1e-4  # saturated hydraulic conductivity (m/s)
hypags: pe.Genuchten = pe.SoilSample(k=ks).hypags()
```

### Soil hydraulic measurements
`pedon` can estimate parameters directly when measurements of soil water retention and/or unsaturated hydraulic conductivity are available. A soil model, together with its SWRC and HCF, is fitted to the data by minimizing the difference between measured and simulated values. This uses nonlinear least-squares algorithm from SciPy [@scipy_paper_2020] and follows the well-established methodology of the `RETC` software [@genuchten_retc_1991].

### Soil model conversion
The same fitting procedure can translate between soil models. The SWRC and HCF generated by one model are sampled over a range of pressure heads and refitted using another formulation. This enables direct model comparison (Figure \ref{fig:swrc_fit}) and facilitates integration with external tools when a different model is required.

```python
# Fitting a Brooks-Corey soil model to existing Mualem-van Genuchten soil model
bc = pe.SoilSample(h=h, theta=theta, k=k).fit(pe.Brooks)
```

![Resulting Brooks-Corey SWRC after fitting on the Mualem-van Genuchten soil model \label{fig:swrc_fit}](figures/swrc_fit.png){height=7.5cm}

# Research impact statement
Soil hydraulic functions and their parameters are essential for simulating variably saturated groundwater flow. Determining these parameters experimentally is difficult, time-consuming, and uncertain [@genuchten_retc_1991; @brandhorst_uncertainty_2017]. Therefore, parameters are often approximated or estimated from reference databases. `pedon` bundles soil hydraulic models and parameter sources in a single framework, enabling efficient parameter derivation without extensive literature searches or ad hoc reimplementation. `pedon` is already used in scientific workflows for variably saturated groundwater flow modeling, including published studies by @vonk_nonlinear_2024 and @collenteur_signatures_2025. It is also a dependency of the Python package [`dutchsoils`](https://github.com/markvdbrink/dutchsoils), which is used in a academic context to process Dutch soil datasets [@heinen_bofek_2022].

# AI usage disclosure
GitHub Copilot was used during development for reviewing pull requests, writing unit tests, providing code completion, and sanity-checking proposed bug fixes. ChatGPT was used for this manuscript to review references, identify linguistic and grammatical errors, and verify compliance with the Journal of Open Source Software requirements. All AI-generated outputs were reviewed by the authors, who take full responsibility for the accuracy and originality of the works.

# References
