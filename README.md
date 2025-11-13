# pedon

*from Greek: πέδον, pedon -> soil*

[![PyPI](https://img.shields.io/pypi/v/pedon?style=flat-square)](https://pypi.org/project/pedon/)
[![PyPi Supported Python Versions](https://img.shields.io/pypi/pyversions/pedon?style=flat-square)](https://pypi.org/project/pedon/)
[![Code Size](https://img.shields.io/github/languages/code-size/martinvonk/pedon?style=flat-square)](https://pypi.org/project/pedon/)
[![PyPi Downloads](https://img.shields.io/pypi/dm/pedon?style=flat-square) ![License](https://img.shields.io/pypi/l/pedon?style=flat-square)](https://pypi.org/project/pedon/)

[![Tests](https://img.shields.io/github/actions/workflow/status/martinvonk/pedon/tests.yml?style=flat-square)](https://github.com/martinvonk/pedon/actions/workflows/tests.yml)
[![MyPy](https://img.shields.io/badge/type_checker-mypy-2A6DB2?style=flat-square)](https://mypy-lang.org/)
[![ruff](https://img.shields.io/badge/linter-ruff-red?style=flat-square)](https://github.com/charliermarsh/ruff)


Python package for (unsaturated) soil properties including pedotransfer functions. This package takes an object-oriented approach to soils, soil samples and soil models. Soil models that are available in this package are:
  - [Mualem-Van Genuchten](https://www.soilphysics.okstate.edu/teaching/soil-6583/references-folder/van%20Genuchten%201980.pdf)
  - [Brooks-Corey](https://www.wipp.energy.gov/library/cra/2009_cra/references/Others/Brooks_Corey_1964_Hydraulic_Properties_ERMS241117.pdf)
  - [Fredlund & Xing](https://projects.mans.edu.eg/heepf/ilppp/cources/12/pdf%20course/1/pressure/osmotic%20soilsalinity22.pdf)
  - Gardner
  - Panday ([Modflow-USG](https://www.gsienv.com/product/modflow-usg/))
  - [Haverkamp](https://doi.org/10.2136/sssaj1977.03615995004100020024x)

This package can fit soil water retention curves and  using least squares (same as [RETC](https://www.pc-progress.com/Documents/programs/retc.pdf)) to measurements.

Measurements for different soil properties and parameters are available from these datasets:
  - [HYDRUS](https://www2.pc-progress.com/downloads/Pgm_Hydrus3D5/HYDRUS_user_Manual_V5.pdf)
  - [VS2D](https://pubs.usgs.gov/wri/1983/4099/report.pdf)
  - Staring Series ([2001](https://edepot.wur.nl/43272) & [2018](https://edepot.wur.nl/512761))

Additionaly, there are pedotransfer functions implemented such as:
  - Van Genuchten: [Wösten](https://www.sciencedirect.com/science/article/pii/S0016706198001323/pdfft?md5=6844f89c07deb81001c2a6eea6fc9e32&pid=1-s2.0-S0016706198001323-main.pdf)
  - Van Genuchten: [Staring Series](https://edepot.wur.nl/43272)
  - Van Genuchten: [Rosetta v1, 2 & 3 (Schaap et al. 2001)](https://doi.org/10.1016/S0022-1694(01)00466-8)
  - Brooks-Corey: [Cosby](https://hess.copernicus.org/articles/25/2445/2021/hess-25-2445-2021.pdf)

## Installation
To get the latest stable version install using:

`pip install pedon`

To get the development version download the GitHub code to your computer. Use cd to get to the download directory and install using:

`pip install -e .`

