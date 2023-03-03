import matplotlib.pyplot as plt
import numpy as np

from ._typing import SoilModel


def swrc(
    sm: SoilModel,
    saturation: bool = False,
    ax: plt.Axes | None = None,
    **kwargs: dict,
) -> plt.Axes:
    """Plot soil water retention curve"""

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(3, 6))
        ax.set_yscale("log")

    h = -np.logspace(-6, 10, num=1000)

    if saturation:
        sw = sm.s(h=h)
        ax.set_xlim(-0.01, 1.01)
    else:
        sw = sm.theta(h=h)

    ax.plot(sw, -h, label=sm.__class__.__name__, **kwargs)
    ax.set_ylim(1e-3, 1e6)
    ax.grid(True)
    return ax


def hcf(
    sm: SoilModel,
    ax: plt.Axes | None = None,
    **kwargs: dict,
) -> plt.Axes:
    """Plot the hydraulic conductivity function"""

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(3, 6))
        ax.set_yscale("log")

    h = np.logspace(-6, 10, num=1000)

    k = sm.k(h=h)

    ax.plot(k, h, label=sm.__class__.__name__, **kwargs)
    ax.set_ylim(1e-3, 1e6)
    ax.grid(True)
    return ax
