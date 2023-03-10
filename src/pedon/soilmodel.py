from dataclasses import dataclass
from typing import Protocol

import matplotlib.pyplot as plt
from numpy import abs as npabs
from numpy import exp, full, linspace, log, logspace

from ._typing import FloatArray


class SoilModel(Protocol):
    def theta(self, h: FloatArray) -> FloatArray:
        """Method to calculate the soil moisture content from the pressure head h"""
        ...

    def s(self, h: FloatArray) -> FloatArray:
        """Method to calculate the effective saturation from the pressure head h"""
        ...

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        """Method to calcualte the permeability from the pressure head h"""
        ...

    def plot(self):
        """Method to plot the soil water retention curve"""
        ...


@dataclass
class Genuchten:
    """Mualem- van Genuchten Soil Model

    van Genuchten, M. Th. (1970) - A Closed-form Equation for Predicting the
    Hydraulic Conductivity of Unsaturated Soil
    """

    k_s: float
    theta_r: float
    theta_s: float
    alpha: float
    n: float
    l: float = 0.5  # noqa: E741

    def __post_init__(self):
        self.m = 1 - 1 / self.n

    def theta(self, h: FloatArray) -> FloatArray:
        theta = (
            self.theta_r
            + (self.theta_s - self.theta_r)
            / (1 + npabs(self.alpha * h) ** self.n) ** self.m
        )
        return theta

    def s(self, h: FloatArray) -> FloatArray:
        return (self.theta(h) - self.theta_r) / (self.theta_s - self.theta_r)

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return self.k_s * s**self.l * (1 - (1 - s ** (1 / self.m)) ** self.m) ** 2

    def plot(self):
        return plot_swrc(self)


@dataclass
class Brooks:
    """Brooks and Corey Soil Model

    Brooks, R.H. and Corey, A.T. (1964) - Hydraulic Properties of Porous Media
    """

    k_s: float
    theta_r: float
    theta_s: float
    h_b: float
    l: float  # noqa: E741

    def theta(self, h: FloatArray) -> FloatArray:
        h = npabs(h)
        if isinstance(h, float):
            if h >= self.h_b:
                return self.theta_r + self.s(h) * (self.theta_s - self.theta_r)
            else:
                return self.theta_s
        else:
            theta = full(h.shape, self.theta_s)
            theta[h >= self.h_b] = self.theta_r + self.s(h[h > self.h_b]) * (
                self.theta_s - self.theta_r
            )
            return theta

    def s(self, h: FloatArray) -> FloatArray:
        h = npabs(h)
        if isinstance(h, float):
            if h >= self.h_b:
                return (h / self.h_b) ** -self.l
            else:
                return 1
        else:
            s = full(h.shape, 1.0)
            s[h >= self.h_b] = (h[h >= self.h_b] / self.h_b) ** -self.l
            return s

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return self.k_s * s ** (3 + 2 / self.l)

    def plot(self):
        return plot_swrc(self)


@dataclass
class Gardner:
    """Gardner Soil Model

    Gardner et al (1970) - Post-irrigation movement of soil water
    """

    k_s: float
    theta_r: float
    theta_s: float
    a: float
    b: float
    m: float

    def theta(self, h: FloatArray) -> FloatArray:
        return self.a * npabs(h) ** -self.b

    def s(self, h: FloatArray) -> FloatArray:
        return (self.theta(h) - self.theta_r) / (self.theta_s - self.theta_r)

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is not None:
            theta = s * (self.theta_s - self.theta_r) + self.theta_r
            return self.k_s * self.a * theta**self.m
        return self.k_s * (self.a / (self.b + npabs(h) ** self.m))

    def plot(self):
        return plot_swrc(self)


@dataclass
class Panday:
    """Panday Soil Model (MODFLOW-USG)

    Panday, S. - USG-Transport: Transport and other Enhancements to MODFLOW-USG
    """

    k_s: float
    theta_r: float
    theta_s: float
    alpha: float  # alpha
    beta: float  # n
    brook: float  # brooks-corey l

    def __post_init__(self):
        self.gamma = 1 - 1 / self.beta  # m
        self.sr = self.theta_r / self.theta_s  # theta_r / theta_s

    def theta(self, h: FloatArray) -> FloatArray:
        return (self.sr + self.s(h) * (1 - self.sr)) * self.theta_s

    def s(self, h: FloatArray) -> FloatArray:
        return (1 + npabs(self.alpha * h) ** self.beta) ** -self.gamma

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return self.k_s * s**self.brook

    def plot(self):
        return plot_swrc(self)


@dataclass
class Fredlund:
    """Fredlund and Xing Soil Model

    Fredlund, D.G. and Xing, A. (1994) - Equations for the soil-water
    characteristic curve
    """

    k_s: float
    theta_s: float
    a: float
    n: float
    m: float

    def theta(self, h: FloatArray) -> FloatArray:
        return self.theta_s / (log(exp(1) + npabs(h / self.a) ** self.n)) ** self.m

    def s(self, h: FloatArray) -> FloatArray:
        return self.theta(h) / self.theta_s

    def k(self, h: FloatArray, s: FloatArray | None = None):
        if s is not None:
            raise NotImplementedError(
                "Can only calculate the hydraulic conductivity"
                "using the pressure head, not the saturation"
            )

        def theta_d(
            h: FloatArray, a: float, n: float, m: float, theta_s: float
        ) -> FloatArray:
            """Derivative of theta (according to WolframAlpha)"""
            return -(
                theta_s
                * m
                * n
                * (h / a) ** n
                * (log((h / a) ** n + exp(1))) ** (-m - 1)
            ) / (h * (h / a) ** n + exp(1) * h)

        h_b = 0.03
        n_steps = 100

        teller = 0
        for y in linspace(log(h), log(1e6), n_steps):
            teller += (
                (self.theta(exp(y)) - self.theta(h))
                / exp(y)
                * theta_d(exp(y), self.a, self.n, self.m, self.theta_s)
            )

        noemer = 0
        for y in linspace(log(h_b), log(1e6), n_steps):
            noemer += (
                (self.theta(exp(y)) - self.theta_s)
                / exp(y)
                * theta_d(exp(y), self.a, self.n, self.m, self.theta_s)
            )

        return self.k_s * (teller / noemer)

    def plot(self):
        return plot_swrc(self)


def plot_swrc(
    sm: SoilModel,
    saturation: bool = False,
    ax: plt.Axes | None = None,
    **kwargs: dict,
) -> plt.Axes:
    """Plot soil water retention curve"""

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(3, 6))
        ax.set_yscale("log")

    h = -logspace(-6, 10, num=1000)

    if saturation:
        sw = sm.s(h=h)
        ax.set_xlim(-0.01, 1.01)
    else:
        sw = sm.theta(h=h)

    ax.plot(sw, -h, label=sm.__class__.__name__, **kwargs)
    ax.set_ylim(1e-3, 1e6)
    ax.grid(True)
    return ax


def plot_hcf(
    sm: SoilModel,
    ax: plt.Axes | None = None,
    **kwargs: dict,
) -> plt.Axes:
    """Plot the hydraulic conductivity function"""

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(3, 6))
        ax.set_yscale("log")
        ax.set_xscale("log")

    h = logspace(-6, 10, num=1000)

    k = sm.k(h=h)

    ax.plot(k, h, label=sm.__class__.__name__, **kwargs)
    ax.set_ylim(1e-3, 1e6)
    ax.set_xlim()
    ax.grid(True)
    return ax
