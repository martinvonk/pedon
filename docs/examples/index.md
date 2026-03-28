# Pedon

In variably saturated groundwater flow modeling, soil hydraulic functions and their parameters are essential for accurate simulations. However, obtaining these parameters experimentally is difficult, time-consuming, and uncertain. `pedon` integrates soil hydraulic models, parameter databases, and parameter estimation methods into a single reproducible framework.

## Key Concepts

**Soil Hydraulic Models:** Parametric descriptions of the soil water retention curve (SWRC) and unsaturated hydraulic conductivity function (HCF). `pedon` provides different formulations from the literature. Additionally, users can easily create their own custom soil models.

**Parameter Datasets:** Pre-determined parameter values for common soils from established databases (HYDRUS, VS2D, Staring series).

**Pedotransfer Functions:** Empirical relationships that estimate soil parameters from easily measured properties (sand/silt/clay content, bulk density, organic matter).

**Parameter Fitting:** Direct estimation of model parameters by fitting to laboratory measurements using nonlinear least-squares optimization (RETC methodology).

**Model Conversion:** Translate between different soil model formulations to enable comparison and integration with external tools.

# Examples

This collection of Jupyter Notebooks demonstrates how to use `pedon` for analyzing and modeling unsaturated soil hydraulic properties. These notebooks show practical workflows for obtaining soil hydraulic parameters from various data sources and integrating them into groundwater flow modeling studies.


```{toctree}
:maxdepth: 1

00_pedon_structure
01_soil_models
02_datasets
03_pedotransfer_functions
04_hypags
05_curve_fitting
06_hypags_curvefitting
```