"""Plotting modules for pedon."""

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from numpy import logspace, nan

if TYPE_CHECKING:
    from .soilmodel import SoilModel


def swrc(
    sm: "SoilModel",
    saturation: bool = False,
    ax: Axes | None = None,
    **kwargs,
) -> Axes:
    """Plot soil water retention curve."""
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(3.0, 6.0))

    h = logspace(-6, 10, num=1000)

    if saturation:
        sw = sm.s(h=h)
        ax.set_xlim(-0.01, 1.01)
    else:
        sw = sm.theta(h=h)
        right = getattr(sm, "theta_s", nan) + 0.01
        ax.set_xlim(left=0.0, right=right)

    label = kwargs.pop("label", sm.__class__.__name__)

    ax.plot(sw, h, label=label, **kwargs)
    ax.set(
        ylim=(min(h), max(h)),
        yscale="log",
        xlabel="\N{GREEK SMALL LETTER THETA}",
        ylabel="|\N{GREEK SMALL LETTER PSI}|",
    )
    ax.grid(True)
    return ax


def hcf(
    sm: "SoilModel",
    ax: Axes | None = None,
    **kwargs,
) -> Axes:
    """Plot the hydraulic conductivity function."""
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(3.0, 6.0))

    h = logspace(-6, 10, num=1000)
    k = sm.k(h=h)

    label = kwargs.pop("label", sm.__class__.__name__)

    ax.plot(k, h, label=label, **kwargs)
    ax.set(
        ylim=(min(h), max(h)),
        yscale="log",
        xscale="log",
        xlabel="$K$",
        ylabel="|\N{GREEK SMALL LETTER PSI}|",
    )
    ax.grid(True)
    return ax


def curves(sm: "SoilModel", axes: list[Axes] | None = None, **kwargs) -> list[Axes]:
    """Plot both the soil water retention curve and the hydraulic conductivity function."""
    if axes is None:
        _, axs = plt.subplots(1, 2, figsize=(6.0, 6.0), sharey=True, layout="tight")
    else:
        axs: list[Axes] = axes

    axs[0] = swrc(sm, ax=axs[0], **kwargs)
    axs[1] = hcf(sm, ax=axs[1], **kwargs)
    return axs
