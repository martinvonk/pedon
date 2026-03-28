# pedon documentation

Pedon (*from Greek: πέδον, pedon -> soil*) is a Python package for working with unsaturated soil hydraulic properties. It provides an object-oriented framework for soils, soil samples, and soil hydraulic models, making it easy to describe, analyze, and parameterize soil water retention and hydraulic conductivity behavior. At its core, Pedon treats each soil model as a Python class that defines a soil water retention curve and an unsaturated hydraulic conductivity function. The package currently includes implementations of several widely used models:
    - Mualem-Van Genuchten
    - Brooks-Corey
    - Fredlund & Xing
    - Haverkamp
    - Gardner
    - Panday (MODFLOW USG-Transport)

Pedon allows these models to be evaluated directly and fitted to data. When laboratory measurements of soil water retention and/or unsaturated hydraulic conductivity are available, model parameters can be estimated by fitting analytical curves to the data using nonlinear least squares, following the same methodology as the widely used RETC.

To support parameterization without direct measurements, Pedon includes access to several established soil hydraulic parameter datasets:
    - HYDRUS
    - VS2D
    - Staring Series (2001 & 2018)

In addition, Pedon implements multiple pedotransfer functions that estimate soil hydraulic parameters from easily measured soil properties:
    - Van Genuchten: Wösten
    - Van Genuchten: Staring Series
    - Van Genuchten: Rosetta v1, 2 & 3 (Schaap et al. 2001)
    - Brooks-Corey: Cosby
    - van Genuchten: HYPAGS

By combining soil hydraulic models, reference datasets, pedotransfer functions, and fitting routines in a single, consistent framework, Pedon makes it straightforward to move from soil information—whether coarse texture data or detailed laboratory measurements—to parameterized soil models ready for use in variably saturated flow simulations.

The pedon package is open-source and hosted on [GitHub](https://github.com/martinvonk/pedon), where you can find more information about the ongoing developments. The package is published on [PyPi](https://pypi.org/project/pedon/) from which it can be installed using `pip install pedon`. If you use this software in your research or analyses, please cite the package via its Zenodo archive:
> Vonk, M. A. & Peche, A. (XXXX). Pedon: A Python package for analyzing unsaturated soil hydraulic properties (vX.X.X). Zenodo. [doi.org/10.5281/zenodo.18222514](https://doi.org/10.5281/zenodo.18222514).



```{toctree}
:maxdepth: 2

examples/index.md
_api/modules.rst
```