"""Image regression tests for plotting functions.

Baseline files can be regenerated with tests/baseline_images/generate_plot_baselines.py.
"""

import matplotlib
import pytest

import pedon as pe
from pedon.plot import curves, hcf, swrc

matplotlib.use("Agg")

savefig_kwargs = {"dpi": 100, "bbox_inches": "tight"}


@pytest.mark.mpl_image_compare(
    baseline_dir="baseline_images/test_plots",
    filename="swrc_genuchten.png",
    tolerance=0.0,
    style="default",
    savefig_kwargs=savefig_kwargs,
)
def test_swrc_plot(gen: pe.SoilModel):
    """Compare the soil water retention curve against a baseline image."""
    ax = swrc(gen, color="tab:blue")
    return ax.figure


@pytest.mark.mpl_image_compare(
    baseline_dir="baseline_images/test_plots",
    filename="hcf_genuchten.png",
    tolerance=0.0,
    style="default",
    savefig_kwargs=savefig_kwargs,
)
def test_hcf_plot(gen: pe.SoilModel):
    """Compare the hydraulic conductivity curve against a baseline image."""
    ax = hcf(gen, color="tab:green")
    return ax.figure


@pytest.mark.mpl_image_compare(
    baseline_dir="baseline_images/test_plots",
    filename="curves_genuchten.png",
    tolerance=0.0,
    style="default",
    savefig_kwargs=savefig_kwargs,
)
def test_curves_plot(gen: pe.SoilModel):
    """Compare the combined curves plot against a baseline image."""
    axes = curves(gen, color="tab:orange")
    return axes[0].figure
