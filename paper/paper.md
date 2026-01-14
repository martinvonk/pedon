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
`pedon` is a Python package for describing and analyzing unsaturated soil hydraulic properties. It provides a framework for soil hydraulic models, along with tools for retrieving parameters from soil databases, applying pedotransfer functions, and fitting model parameters to measurements.

# Statement of need
Researchers and engineers working with unsaturated soils often need estimates of soil parameters for their variably saturated groundwater flow models. `pedon` addresses this need by providing a modern Python toolkit that brings together commonly used soil hydraulic models, parameter databases, pedotransfer functions, and fitting routines. This makes soil analysis faster, more reproducible, and easier to integrate into existing groundwater modeling workflows.

# Soil hydraulic models
Different soil hydraulic models (hereafter referred to as soil models) are available in `pedon`. A soil model is a parametric description of soil hydraulic functions, namely the soil water retention curve (SWRC) and the (unsaturated) hydraulic conductivity function (HCF). These link the soil water content and flow to pressure head or saturation for use in variably saturated groundwater flow models. At this time, the following soil models are available:

- Mualem-van Genuchten [@genuchten_mualem_1980]: `pedon.Genuchten`
- Brooks-Corey [@brooks_corey_1964]: `pedon.Brooks`
- Combination of the van Genuchten SWRC and Brooks-Corey HCF [@fuentes_burdine_1992; @panday_mfusgt_2025]: `pedon.Panday`
- Fredlund-Xing [@fredlund_xing_1994]: `pedon.Fredlund`
- Haverkamp [@haverkamp_model_1977]: `pedon.Haverkamp`
- Gardner(-Kozeny) [@gardner_model_1958; @brutsaert_kozeny_1967; @bakker_gardner_2009; @mathias_gardner_2006]: `pedon.Gardner`
- Gardner-Rucker [@rucker_gardner_2005]: `pedon.Rucker`
- Combination of the van Genuchten SWRC and Gardner HCF [@genuchten_mualem_1980; @gardner_model_1958]: `pedon.GenuchtenGardner`

# Software design
The soil models are implemented as Python classes, providing a clear and consistent structure in which model-specific methods, such as those for evaluating the SWRC and HCF, are defined. For example, the Mualem–van Genuchten model can be instantiated and used as follows:

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

Thanks to its object-oriented design and duck typing, `pedon` allows users to define their own soil model classes, in which they can implement custom or literature-based SWRC and HCF. `pedon` depends on Python's well-maintained scientific ecosystem, including NumPy [@numpy_article_2020], SciPy [@scipy_paper_2020], Matplotlib [@matplotlib_paper_2007], and Pandas [@pandas_software_2020; @pandas_paper_2010].

# Soil hydraulic parameters
Soil hydraulic parameters define the behavior of a soil model by determining the shape of its soil water retention curve and hydraulic conductivity function. These parameters are therefore essential inputs for variably saturated groundwater flow models. In practice, the parameters are rarely measured directly and are often derived from reference datasets, empirical relationships or laboratory measurements. `pedon` provides a framework for working with soil hydraulic parameters by linking them directly to soil models and offering tools to obtain and fit these parameters from existing datasets, easily measured soil properties, and direct measurements of soil water content and hydraulic conductivity.

## Parameter datasets
`pedon` includes a dataset of Brooks–Corey and Mualem–van Genuchten parameters for a wide range of soils. At present, this dataset is compiled from three established soil hydraulic parameter databases:

- Average values of selected soil hydraulic parameters for 12 major soil textural groups, as defined by @carsel_dataset_1988. This dataset is also used in the popular HYDRUS software for variably saturated groundwater flow modeling [@simunek_hydrus1d_2009].
- Dataset obtained from the VS2D software [@healy_vs2d_1990] containing both Brooks–Corey and Mualem–van Genuchten parameters.
- The Staring series (Staringreeks in Dutch) is a database of soil hydraulic functions in the Netherlands [@wosten_staringreeks_2001; @heinen_staringreeks_2020; @heinen_bofek_2022]. It contains descriptions of both topsoils (`B_`) and subsoils (`O_`) based on hundreds of samples. These samples were processed to obtain parameters for the Mualem-van Genuchten soil model [@genuchten_mualem_1980; @wosten_texture_1988].

The databases can be called via the following code:
```python
hydrus = pe.Soil("Sand").from_name(pe.Genuchten, source="HYDRUS")
vs2d = pe.Soil("Sand").from_name(pe.Brooks, source="VS2D")
staring = pe.Soil("B01").from_name(pe.Genuchten, source="Staring_2018")
```

## Parameter estimation
`pedon` provides two approaches for obtaining soil hydraulic parameters from available soil measurements. The first approach uses pedotransfer functions based on easily measured soil properties. The second approach relies on direct measurements of soil water content and hydraulic conductivity. Both methods are described below.

### Pedotransfer functions
Soil hydraulic parameters can be estimated using pedotransfer functions. Pedotransfer functions relate easily measured soil properties (e.g. sand, silt, or clay percentage, bulk density and organic matter content) to soil hydraulic parameters [@bouma_pedotransfer_1989]. `pedon` implements some pedotransfer functions from the literature, including those of @wosten_pedotransfer_1999, @wosten_staringreeks_2001, @cosby_pedotransfer_1984, and @cooper_pedotransfer_2021. In addition, `pedon` provides access to soil hydraulic parameter databases such as Rosetta [@schaap_rosetta_2001] and HYPAGS [@peche_hypags_2024], the latter of which enables parameter estimation based solely on single values of saturated hydraulic conductivity or representative grain diameters.

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
`pedon` can estimate soil hydraulic parameters directly when measurements of the soil water retention and/or unsaturated hydraulic conductivity are available. A soil model, together with its SWRC and HCF, is fitted to the measured values by minimizing the difference between measured and simulated data. This fitting procedure uses the nonlinear least-squares optimization provided by SciPy [@scipy_paper_2020] and follows the well-established methodology implemented in the RETC software [@genuchten_retc_1991].

### Soil model conversion
The same fitting framework for sample measurements can also be used to translate between different soil hydraulic models. The SWRC and HCF generated by one model can be sampled over a range of pressure heads and then refitted using a different soil model formulation. This enables direct model comparison (Figure \ref{fig:swrc_fit}) and facilitates integration with external modeling tools when another soil model formulation is required.

```python
# Fitting a Brooks-Corey soil model to existing Mualem-van Genuchten soil model
bc = pe.SoilSample(h=h, theta=theta, k=k).fit(pe.Brooks)
```

![Resulting fit of the Brooks-Corey soil model on the Mualem-van Genuchten soil model \label{fig:swrc_fit}](figures/swrc_fit.png)

# Research impact statement
Soil hydraulic functions and their parameters are essential for simulating variably saturated groundwater flow. Determining these parameters experimentally is difficult, time-consuming, expensive, and uncertain [@genuchten_retc_1991; @brandhorst_uncertain_2017]. As a result, parameters are often approximated or taken from reference databases. `pedon` bundles soil hydraulic models and parameter sources in a single framework, enabling efficient parameter derivation without extensive literature searches or ad hoc reimplementation.

As a result, `pedon` is already used in scientific workflows for variably saturated groundwater flow modeling, including in published research workflows from @vonk_nonlinear_2024 and @collenteur_signatures_2025. Additionally, `pedon` is a dependency of the Python package [`dutchsoils`](https://github.com/markvdbrink/dutchsoils), which is used in a scientific context to access and process Dutch soil datasets.

# AI usage disclosure
Generative AI tools (GitHub Copilot) were used during the development of this project for reviewing pull requests, writing unit tests, providing code completion suggestions, and sanity-checking proposed bug solutions.

For the preparation of this manuscript, generative AI (OpenAI ChatGPT) was used to review the reference list, identify linguistic, grammatical, and textual errors, and verify that the manuscript complies with the requirements of the Journal of Open Source Software.

All outputs generated by AI tools were thoroughly reviewed and validated by the authors. Human authors made all primary architectural, design, and scientific decisions. The authors take full responsibility for the accuracy, originality, licensing compliance, and ethical and legal standards of all submitted materials.

# References
