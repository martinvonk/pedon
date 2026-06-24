# pedon

[![PyPI](https://img.shields.io/pypi/v/pedon?style=flat-square&color=007ec6)](https://pypi.org/project/pedon/)
[![PyPi Supported Python Versions](https://img.shields.io/pypi/pyversions/pedon?style=flat-square&color=007ec6)](https://pypi.org/project/pedon/)
[![Code Size](https://img.shields.io/github/languages/code-size/martinvonk/pedon?style=flat-square&color=007ec6)](https://pypi.org/project/pedon/)
[![PyPi Downloads](https://img.shields.io/pypi/dm/pedon?style=flat-square&color=0a3d62)](https://pypi.org/project/pedon/)
[![License](https://img.shields.io/pypi/l/pedon?style=flat-square&color=0a3d62&logo=open-source-initiative&logoColor=white)](https://pypi.org/project/pedon/)

[![JOSS](https://img.shields.io/badge/JOSS-10.21105/joss.10409-ff6600?style=flat-square)](https://doi.org/10.21105/joss.10409)
[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281/zenodo.18222514-ff6600?style=flat-square)](https://doi.org/10.5281/zenodo.18222514)

[![Tests](https://img.shields.io/github/actions/workflow/status/martinvonk/pedon/tests.yml?style=flat-square&color=green)](https://github.com/martinvonk/pedon/actions/workflows/tests.yml)
[![Coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/martinvonk/a891dd883f258c46826ba657cf9a89cf/raw/coverage.json&style=flat-square&color=green)](https://github.com/martinvonk/pedon/actions/workflows/tests.yml)
[![Typed: Mypy](https://img.shields.io/badge/typing-mypy-darkgreen?style=flat-square)](https://mypy-lang.org/)
[![Typed: Ty](https://img.shields.io/badge/typing-ty-darkgreen?style=flat-square)](https://docs.astral.sh/ty/)
[![Formatter and Linter: Ruff](https://img.shields.io/badge/linter-ruff-darkgreen?style=flat-square)](https://docs.astral.sh/ruff/)

`pedon` (*from Greek: πέδον, pedon -> soil*) is a Python package for working with unsaturated soil hydraulic properties. It provides an object-oriented framework for soils, soil samples, and soil hydraulic models, making it easy to describe, analyze, and parameterize soil water retention and hydraulic conductivity behavior. At its core, `pedon` treats each soil model as a Python class that defines a soil water retention curve and an unsaturated hydraulic conductivity function. The package currently includes implementations of several widely used models:
  - [Mualem-Van Genuchten](https://doi.org/10.2136/sssaj1980.03615995004400050002x)
  - [Brooks-Corey](https://mountainscholar.org/items/3c7b98df-13e3-486c-9d1e-949a7a869f76)
  - [Gardner](https://doi.org/10.1097/00010694-195804000-00006)
  - [Campbell](https://doi.org/10.1097/00010694-197406000-00001)
  - [Haverkamp](https://doi.org/10.2136/sssaj1977.03615995004100020024x)
  - [Kosugi](https://doi.org/10.1029/96WR01776)
  - [Fredlund & Xing](https://doi.org/10.1139/t94-061)
  - [Brunswick](https://doi.org/10.1029/2018WR024584)
  - Hysteresis model from [Kool & Parker](https://doi.org/10.1029/WR023i001p00105)
  - Dual-porosity model from [Gerke & van Genuchten](https://doi.org/10.1029/92WR02339)
  - van Genuchten-like SWRC and Gardner HCF [Rucker](https://doi.org/10.1016/j.advwatres.2005.01.004)
  - van Genuchten SWRC and Gardner HCF [Genuchten-Gardner](https://doi.org/10.2136/sssaj1980.03615995004400050002x)
  - van Genuchten SWRC and Brooks-Corey HCF [Panday (MODFLOW USG-Transport)](https://www.gsienv.com/software/modflow-usg/modflow-usg/)

`pedon` allows these models to be evaluated directly and fitted to data. When laboratory measurements of soil water retention and/or unsaturated hydraulic conductivity are available, model parameters can be estimated by fitting analytical curves to the data using nonlinear least squares, following the same methodology as the widely used [RETC](https://www.pc-progress.com/Documents/programs/retc.pdf).

To support parameterization without direct measurements, `pedon` includes access to several established soil hydraulic parameter datasets:
  - [HYDRUS](https://www.pc-progress.com/downloads/pgm_hydrus1d/hydrus1d-4.08.pdf)
  - [VS2D](https://doi.org/10.3133/wri904025)
  - Staring Series: [2001](https://edepot.wur.nl/43272) and [2018](https://doi.org/10.18174/512761)
  - [Clapp & Hornberger (1978)](https://doi.org/10.1029/WR014i004p00601)

In addition, `pedon` implements multiple pedotransfer functions that estimate soil hydraulic parameters from easily measured soil properties:
  - Van Genuchten: [Wösten](https://doi.org/10.1016/S0016-7061(98)00132-3)
  - Van Genuchten: [Staring Series](https://edepot.wur.nl/43272)
  - Van Genuchten: [Rosetta v1, 2 & 3 (Schaap et al. 2001)](https://doi.org/10.1016/S0022-1694(01)00466-8)
  - Van Genuchten: [Vereecken](https://doi.org/10.1097/00010694-198912000-00001)
  - Van Genuchten: [Weynants](https://doi.org/10.2136/vzj2008.0062)
  - Van Genuchten: [Tóth](https://doi.org/10.1111/ejss.12192)
  - Van Genuchten: [Hodnett & Tomasella](https://doi.org/10.1016/S0016-7061(02)00105-2)
  - Brooks-Corey: [Cosby](https://doi.org/10.1029/WR020i006p00682) (see also [Cooper et al. 2021](https://doi.org/10.5194/hess-25-2445-2021))
  - Brooks-Corey: [Saxton](https://doi.org/10.2136/sssaj1986.03615995005000040039x)
  - van Genuchten: HYPAGS ([Peche et al. 2024](https://doi.org/10.1111/gwat.13365); [Peche & Houben 2023](https://doi.org/10.1111/gwat.13266))

By combining soil hydraulic models, reference datasets, pedotransfer functions, and fitting routines in a single, consistent framework, `pedon` makes it straightforward to move from soil information—whether coarse texture data or detailed laboratory measurements—to parameterized soil models ready for use in variably saturated flow simulations.

If you use this software in your research or analyses, you can cite the package the Journal of Open Source Software
> Vonk, M. A. & Peche, A. (2026). `pedon`: A Python package for analyzing unsaturated soil hydraulic properties. Journal of Open Source Software, 11(122), 10409. [doi.org/10.21105/joss.10409](https://doi.org/10.21105/joss.10409).

You can also cite a specific version of the package via its Zenodo archive:
> Vonk, M. A. & Peche, A. (XXXX). `pedon`: A Python package for analyzing unsaturated soil hydraulic properties (vX.X.X). Zenodo. [doi.org/10.5281/zenodo.18222514](https://doi.org/10.5281/zenodo.18222514).

## Installation
To get the latest stable version install using:

`pip install pedon`

To get the development version download the GitHub code to your computer. Use cd to get to the download directory and install using:

`pip install -e .`

## Documentation

The documentation for `pedon` is available at [https://martinvonk.github.io/pedon/](https://martinvonk.github.io/pedon/). It includes an extensive user guide with Jupyter Notebooks as examples and an API reference.