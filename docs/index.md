# pedon documentation

Python package for (unsaturated) soil properties including pedotransfer functions. This package takes an object-oriented approach to soils, soil samples and soil models. Soil models that are available in this package are:

- Mualem-Van Genuchten
- Brooks-Corey
- Fredlund & Xing
- Gardner
- Panday (Modflow-USG)
- Haverkamp

This package can fit soil water retention curves and using least squares (same as RETC) to measurements. Measurements for different soil properties and parameters are available from these datasets:

- HYDRUS
- VS2D
- Staring Series (2001 & 2018)

Additionaly, there are pedotransfer functions implemented such as:

- Van Genuchten: Wösten
- Van Genuchten: Staring Series
- Van Genuchten: Rosetta v1, 2 & 3 (Schaap et al. 2001)
- Van Genuchten: HYPAGS
- Brooks-Corey: Cosby

The pedon package is open-source and hosted on [GitHub]((https://github.com/martinvonk/pedon)), where you can find more information about the available drought indices and ongoing development. The package is published on [PyPi](https://pypi.org/project/pedon/) from which it can be installed using `pip install pedon`.

```{toctree}
:maxdepth: 2

examples/index.md
_api/modules.rst
```