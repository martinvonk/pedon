from dataclasses import dataclass, field
from typing import Protocol, Type, runtime_checkable

import matplotlib.pyplot as plt
from numpy import abs as npabs
from numpy import exp, full, linspace, log, log10, logspace

from ._typing import FloatArray, SoilModelNames


@runtime_checkable
class SoilModel(Protocol):
    def theta(self, h: FloatArray) -> FloatArray:
        """Method to calculate the soil moisture content from the pressure head h"""
        ...

    def s(self, h: FloatArray) -> FloatArray:
        """Method to calculate the effective saturation from the pressure head h"""
        ...

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        """Method to calcualte the relative permeability from the pressure head h"""
        ...

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        """Method to calcualte the permeability from the pressure head h"""
        ...

    def h(self, theta: FloatArray) -> FloatArray:
        """Method to calcualte the pressure head h from the water content"""
        ...

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        """Method to plot the soil water retention curve"""
        ...


@dataclass
class Genuchten:
    """Mualem-van Genuchten Soil Model

    van Genuchten, M. Th. (1970) - A Closed-form Equation for Predicting the
    Hydraulic Conductivity of Unsaturated Soil
    """

    k_s: float
    theta_r: float
    theta_s: float
    alpha: float
    n: float
    l: float = 0.5  # noqa: E741
    m: float = field(init=False, repr=False)

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

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return s**self.l * (1 - (1 - s ** (1 / self.m)) ** self.m) ** 2

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        se = (theta - self.theta_r) / (self.theta_s - self.theta_r)
        h = 1 / self.alpha * ((1 / se) ** (1 / self.m) - 1) ** (1 / self.n)
        return h

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


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
            theta[h >= self.h_b] = self.theta_r + self.s(h[h >= self.h_b]) * (
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

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return s ** (3 + 2 / self.l)

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        if isinstance(theta, float):
            if theta >= self.theta_r:
                return self.h_b * ((theta - self.theta_r) / (self.s(theta))) ** (
                    -1 / self.l
                )
            else:
                return self.h_b
        else:
            h = full(theta.shape, self.h_b)
            mask = theta >= self.theta_r
            h[mask] = self.h_b * (
                (theta[mask] - self.theta_r) / (self.s(theta[mask]))
            ) ** (-1 / self.l)
            return h

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class Haverkamp:
    """Haverkamp Soil Model

    Haverkamp, R., Vauclin, M., Touma, J., Wierenga, P. J., & Vachaud, G. (1977).
    A comparison of numerical simulation models for one-dimensional infiltration.
    """

    k_s: float
    theta_r: float
    theta_s: float
    alpha: float
    beta: float
    a: float

    def theta(self, h: FloatArray) -> FloatArray:
        return (
            self.alpha
            * (self.theta_s - self.theta_r)
            / (self.alpha + npabs(h) ** self.beta)
            + self.theta_r
        )

    def s(self, h: FloatArray) -> FloatArray:
        return (self.theta(h) - self.theta_r) / (self.theta_s - self.theta_r)

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is not None:
            return (self.a * s) / (self.a * s + self.alpha * (1.0 - s))
        return self.a / (self.a + npabs(h) ** self.beta)

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        s = (theta - self.theta_r) / (self.theta_s - self.theta_r)
        return (self.alpha * ((1.0 / s) - 1.0)) ** (1.0 / self.beta)

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class Gardner:
    """Gardner(-Kozeny) Soil Model

    Gardner, W.H. (1958) - Some steady-state solutions of the unsaturated
    moisture flow equation with application to evaporation from soils
    Bakker and Nieber (2009) - Damping of Sinusoidal Surface Flux Fluctuations
    with Soil Depth
    """

    k_s: float
    theta_s: float
    c: float
    m: float

    def theta(self, h: FloatArray) -> FloatArray:
        return self.theta_s * exp(-self.m * npabs(h))

    def s(self, h: FloatArray) -> FloatArray:
        return self.theta(h) / self.theta_s

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return exp(-self.c * npabs(h))

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        return -(1.0 / self.m) * log(theta / self.theta_s)

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class Panday:
    """Panday Soil Model (MODFLOW-USG)

    Panday, S. - USG-Transport: Transport and other Enhancements to MODFLOW-USG
    """

    k_s: float
    theta_r: float = field(repr=False)
    theta_s: float = field(repr=False)
    alpha: float  # alpha
    beta: float  # n
    brook: float  # brooks-corey l
    h_b: float = field(default=0.0, repr=False)
    sr: float = field(init=False, repr=True)
    gamma: float = field(init=False, repr=False)  # 1 - 1 / beta
    sy: float = field(init=False, repr=False)
    ss: float = field(default=1e-6, repr=False)

    def __post_init__(self):
        self.sr = self.theta_r / self.theta_s  # theta_r / theta_s
        self.gamma = 1 - 1 / self.beta  # m
        theta_fc = (
            self.beta ** -(0.60 * (2 + log10(self.k_s))) * (self.theta_s - self.theta_r)
            + self.theta_r
        )  # assumes k_s is in [cm]
        self.sy = self.theta_s - theta_fc

    def theta(self, h: FloatArray) -> FloatArray:
        return (self.sr + self.s(h) * (1 - self.sr)) * self.theta_s

    def s(self, h: FloatArray) -> FloatArray:
        return (1 + npabs(self.alpha * (h - self.h_b)) ** self.beta) ** -self.gamma

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return s**self.brook

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        se = (theta - self.theta_r) / (self.theta_s - self.theta_r)
        h = 1 / self.alpha * ((1 / se) ** (1 / self.gamma) - 1) ** (1 / self.beta)
        return h

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


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

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
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
        return teller / noemer

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        return self.a * (exp(self.theta_s / theta) ** (1 / self.m) - exp(1)) ** (
            1 / self.n
        )

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


def get_soilmodel(
    soilmodel_name: SoilModelNames,
) -> Type[SoilModel]:
    sms = {
        "Genuchten": Genuchten,
        "Brooks": Brooks,
        "Gardner": Gardner,
        "Panday": Panday,
        "Fredlund": Fredlund,
        "Haverkamp": Haverkamp,
    }
    return sms[soilmodel_name]


def plot_swrc(
    sm: SoilModel, saturation: bool = False, ax: plt.Axes | None = None, **kwargs
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

    if "label" in kwargs:
        label = kwargs.pop("label")
    else:
        label = getattr(getattr(sm, "__class__"), "__name__")

    ax.plot(sw, -h, label=label, **kwargs)
    ax.set_ylim(1e-3, 1e6)
    ax.grid(True)
    return ax


def plot_hcf(
    sm: SoilModel,
    ax: plt.Axes | None = None,
    **kwargs,
) -> plt.Axes:
    """Plot the hydraulic conductivity function"""

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(3, 6))
        ax.set_yscale("log")
        ax.set_xscale("log")

    h = logspace(-6, 10, num=1000)
    k = sm.k(h=h)

    if "label" in kwargs:
        label = kwargs.pop("label")
    else:
        label = getattr(getattr(sm, "__class__"), "__name__")

    ax.plot(k, h, label=label, **kwargs)
    ax.set_ylim(1e-3, 1e6)
    ax.set_xlim()
    ax.grid(True)
    return ax
