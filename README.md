# pedon

[![PyPI](https://img.shields.io/pypi/v/pedon?style=flat-square&color=007ec6)](https://pypi.org/projectpedonspei/)
[![PyPi Supported Python Versions](https://img.shields.io/pypi/pyversions/pedon?style=flat-square&color=007ec6)](https://pypi.org/project/pedon/)
[![Code Size](https://img.shields.io/github/languages/code-size/martinvonk/pedon?style=flat-square&color=007ec6)](https://pypi.org/project/pedon/)
[![PyPi Downloads](https://img.shields.io/pypi/dm/pedon?style=flat-square&color=0a3d62)](https://pypi.org/project/pedon/)
[![License](https://img.shields.io/pypi/l/pedon?style=flat-square&color=0a3d62&logo=open-source-initiative&logoColor=white)](https://pypi.org/project/pedon/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.18222514-ff6600?style=flat-square)](https://doi.org/10.5281/zenodo.18222514)

[![Tests](https://img.shields.io/github/actions/workflow/status/martinvonk/pedon/tests.yml?style=flat-square&color=green)](https://github.com/martinvonk/pedon/actions/workflows/tests.yml)
[![Typed: MyPy](https://img.shields.io/badge/type_checker-mypy-darkgreen?style=flat-square)](https://mypy-lang.org/)
[![Formatter and Linter: ruff](https://img.shields.io/badge/linter-ruff-darkgreen?style=flat-square)](https://github.com/charliermarsh/ruff)

Pedon (*from Greek: πέδον, pedon -> soil*) is a Python package for working with unsaturated soil hydraulic properties. It provides an object-oriented framework for soils, soil samples, and soil hydraulic models, making it easy to describe, analyze, and parameterize soil water retention and hydraulic conductivity behavior. At its core, Pedon treats each soil model as a Python class that defines a soil water retention curve and an unsaturated hydraulic conductivity function. The package currently includes implementations of several widely used models:
  - [Mualem-Van Genuchten](https://www.soilphysics.okstate.edu/teaching/soil-6583/references-folder/van%20Genuchten%201980.pdf)
  - [Brooks-Corey](https://www.wipp.energy.gov/library/cra/2009_cra/references/Others/Brooks_Corey_1964_Hydraulic_Properties_ERMS241117.pdf)
  - [Fredlund & Xing](https://projects.mans.edu.eg/heepf/ilppp/cources/12/pdf%20course/1/pressure/osmotic%20soilsalinity22.pdf)
  - [Haverkamp](https://doi.org/10.2136/sssaj1977.03615995004100020024x)
  - [Gardner](https://doi.org/10.1097/00010694-195804000-00006)
  - Panday ([MODFLOW USG-Transport](https://www.gsienv.com/product/modflow-usg/))

Pedon allows these models to be evaluated directly and fitted to data. When laboratory measurements of soil water retention and/or unsaturated hydraulic conductivity are available, model parameters can be estimated by fitting analytical curves to the data using nonlinear least squares, following the same methodology as the widely used [RETC](https://www.pc-progress.com/Documents/programs/retc.pdf).

To support parameterization without direct measurements, Pedon includes access to several established soil hydraulic parameter datasets:
  - [HYDRUS](https://www2.pc-progress.com/downloads/Pgm_Hydrus3D5/HYDRUS_user_Manual_V5.pdf)
  - [VS2D](https://pubs.usgs.gov/wri/1983/4099/report.pdf)
  - Staring Series ([2001](https://edepot.wur.nl/43272) & [2018](https://edepot.wur.nl/512761))

In addition, Pedon implements multiple pedotransfer functions that estimate soil hydraulic parameters from easily measured soil properties:
  - Van Genuchten: [Wösten](https://www.sciencedirect.com/science/article/pii/S0016706198001323/pdfft?md5=6844f89c07deb81001c2a6eea6fc9e32&pid=1-s2.0-S0016706198001323-main.pdf)
  - Van Genuchten: [Staring Series](https://edepot.wur.nl/43272)
  - Van Genuchten: [Rosetta v1, 2 & 3 (Schaap et al. 2001)](https://doi.org/10.1016/S0022-1694(01)00466-8)
  - Brooks-Corey: [Cosby](https://hess.copernicus.org/articles/25/2445/2021/hess-25-2445-2021.pdf)
  - van Genuchten: [HYPAGS](https://doi.org/10.1111/gwat.13266)

By combining soil hydraulic models, reference datasets, pedotransfer functions, and fitting routines in a single, consistent framework, Pedon makes it straightforward to move from soil information—whether coarse texture data or detailed laboratory measurements—to parameterized soil models ready for use in variably saturated flow simulations.

If you use this software in your research or analyses, please cite the package via its Zenodo archive:
> Vonk, M. A. & Peche, A. (XXXX). Pedon: A Python package for analyzing unsaturated soil hydraulic properties (vX.X.X). Zenodo. [doi.org/10.5281/zenodo.18222514](https://doi.org/10.5281/zenodo.18222514).

## Installation
To get the latest stable version install using:

`pip install pedon`

To get the development version download the GitHub code to your computer. Use cd to get to the download directory and install using:

`pip install -e .`

