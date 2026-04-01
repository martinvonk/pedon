from dataclasses import dataclass, field
from typing import Protocol, Type, runtime_checkable

import matplotlib.pyplot as plt
from numpy import abs as npabs
from numpy import exp, full, linspace, log, log10, logspace, maximum

from ._typing import FloatArray, MatplotlibAxes, SoilModelNames


@runtime_checkable
class SoilModel(Protocol):
    """Protocol for soil models

    This protocol defines the interface for custom soil models. To create a custom
    soil model, implement all methods defined below. Your class will automatically
    conform to this protocol.

    Example:
        >>> @dataclass
        ... class CustomModel:
        ...     k_s: float
        ...     theta_s: float
        ...     theta_r: float
        ...
        ...     def theta(self, h: FloatArray) -> FloatArray:
        ...         # Calculate soil moisture content from pressure head
        ...         return ...
        ...
        ...     def s(self, h: FloatArray) -> FloatArray:
        ...         # Calculate effective saturation from pressure head
        ...         return (self.theta(h) - self.theta_r) / (self.theta_s - self.theta_r)
        ...
        ...     def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        ...         # Calculate relative permeability from pressure head or saturation
        ...         return ...
        ...
        ...     def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        ...         # Calculate hydraulic conductivity from pressure head or saturation
        ...         return self.k_s * self.k_r(h=h, s=s)
        ...
        ...     def h(self, theta: FloatArray) -> FloatArray:
        ...         # Inverse of theta method
        ...         return ...
        ...
        ...     def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        ...         # Plot the soil water retention curve by calling `plot_swrc`
        ...         return plot_swrc(sm=self, ax=ax)
    """

    def theta(self, h: FloatArray) -> FloatArray:
        """Calculate soil moisture content (water content) from
        pressure head h."""
        ...

    def s(self, h: FloatArray) -> FloatArray:
        """Calculate effective saturation from pressure head h.
        Effective saturation is normalized between 0 (dry) and 1 (saturated)."""
        ...

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        """Calculate relative permeability (or relative hydraulic conductivity).
        Relative permeability is normalized between 0 (dry) and 1 (saturated).
        Can be calculated from either pressure head h or saturation s."""
        ...

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        """Calculate hydraulic conductivity from pressure head h or saturation s."""
        ...

    def h(self, theta: FloatArray) -> FloatArray:
        """Calculate pressure head h from soil moisture content (inverse of theta)."""
        ...

    def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        """Plot the soil water retention curve by calling `plot_swrc`."""
        ...

    def convert(self, method: str | None = None) -> "list[str] | SoilModel":
        """Convert this soil model to another model type.

        If `method` is None, returns a list of available conversion methods for this
        soil model. Otherwise, performs the specified conversion and returns the
        converted soil model instance.
        """
        ...

    def fit(self, sm: "Type[SoilModel]", h: FloatArray, **kwargs) -> "SoilModel":
        """Fit this soil model to another soil model type using the provided pressure
        head data `h` and optional keyword arguments for fitting."""
        ...


@dataclass
class Genuchten:
    """Mualem-van Genuchten Soil Model

    van Genuchten, M. Th. (1980) - A Closed-form Equation for Predicting the
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

    def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        return plot_swrc(sm=self, ax=ax)

    def convert(self, method: str | None = None) -> list[str] | SoilModel:

        if method is None:
            return SoilModelConverter.list_methods(sm=self)
        return getattr(SoilModelConverter, method)(self)

    def fit(self, sm: Type[SoilModel], h: FloatArray, **kwargs):
        from .soil import SoilSample

        ss = SoilSample(h=h, theta=self.theta(h), k=self.k(h))
        return ss.fit(sm=sm, **kwargs)


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
        if isinstance(theta, (float, int)):
            if theta >= self.theta_r:
                return (
                    self.h_b
                    * (theta - self.theta_r)
                    / (self.theta_s - self.theta_r) ** (-1 / self.l)
                )
            else:
                return self.h_b
        else:
            h = full(theta.shape, self.h_b, dtype=float)
            mask = theta >= self.theta_r
            h[mask] = self.h_b * (
                (theta[mask] - self.theta_r) / (self.s(theta[mask]))
            ) ** (-1 / self.l)
            return h

    def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        return plot_swrc(sm=self, ax=ax)

    def convert(self, method: str | None = None) -> list[str] | SoilModel:

        if method is None:
            return SoilModelConverter.list_methods(sm=self)
        return getattr(SoilModelConverter, method)(self)

    def fit(self, sm: Type[SoilModel], h: FloatArray, **kwargs):
        from .soil import SoilSample

        ss = SoilSample(h=h, theta=self.theta(h), k=self.k(h))
        return ss.fit(sm=sm, **kwargs)


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

    def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        return plot_swrc(sm=self, ax=ax)

    def convert(self, method: str | None = None) -> list[str] | SoilModel:

        if method is None:
            return SoilModelConverter.list_methods(sm=self)
        return getattr(SoilModelConverter, method)(self)

    def fit(self, sm: Type[SoilModel], h: FloatArray, **kwargs):
        from .soil import SoilSample

        ss = SoilSample(h=h, theta=self.theta(h), k=self.k(h))
        return ss.fit(sm=sm, **kwargs)


@dataclass
class Gardner:
    """Gardner(-Kozeny) Soil Model with optional air entry pressure.

    Gardner, W.H. (1958) - Some steady-state solutions of the unsaturated
    moisture flow equation with application to evaporation from soils
    Bakker and Nieber (2009) - Damping of Sinusoidal Surface Flux Fluctuations
    with Soil Depth
    """

    k_s: float
    theta_s: float
    m: float
    c: float
    h_b: float = 0.0

    def theta(self, h: FloatArray) -> FloatArray:
        h_abs = maximum(npabs(h) - self.h_b, 0.0)
        return self.theta_s * exp(-self.m * h_abs)

    def s(self, h: FloatArray) -> FloatArray:
        return self.theta(h) / self.theta_s

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is not None:
            theta = s * self.theta_s
            h = self.h(theta)
        h_abs = maximum(npabs(h) - self.h_b, 0.0)
        return exp(-self.c * h_abs)

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        return self.h_b - (1.0 / self.m) * log(theta / self.theta_s)

    def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        return plot_swrc(sm=self, ax=ax)

    def convert(self, method: str | None = None) -> list[str] | SoilModel:

        if method is None:
            return SoilModelConverter.list_methods(sm=self)
        return getattr(SoilModelConverter, method)(self)

    def fit(self, sm: Type[SoilModel], h: FloatArray, **kwargs):
        from .soil import SoilSample

        ss = SoilSample(h=h, theta=self.theta(h), k=self.k(h))
        return ss.fit(sm=sm, **kwargs)


@dataclass
class Rucker:
    """Gardner(-Rucker) Soil Model

    Gardner, W.H. (1958) - Some steady-state solutions of the unsaturated
    moisture flow equation with application to evaporation from soils
    Rucker, D. F., Warrick, A. W., & Ferré, T. P. (2005). Parameter equivalence
    for the Gardner and van Genuchten soil hydraulic conductivity functions
    for steady vertical flow with inclusions.
    """

    k_s: float
    theta_r: float
    theta_s: float
    m: float
    c: float

    def theta(self, h: FloatArray) -> FloatArray:
        return self.theta_r + (self.theta_s - self.theta_r) * (
            ((1 + 0.5 * self.c * npabs(h)) * exp(-0.5 * self.c * npabs(h)))
            ** (2 / (self.m + 2))
        )

    def s(self, h: FloatArray) -> FloatArray:
        return (self.theta(h) - self.theta_r) / (self.theta_s - self.theta_r)

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is not None:
            theta = s * (self.theta_s - self.theta_r) + self.theta_r
            h = self.h(theta)
        return exp(-self.c * npabs(h))

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        return -(1.0 / self.m) * log(
            (theta - self.theta_r) / (self.theta_s - self.theta_r)
        )

    def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        return plot_swrc(sm=self, ax=ax)

    def convert(self, method: str | None = None) -> list[str] | SoilModel:

        if method is None:
            return SoilModelConverter.list_methods(sm=self)
        return getattr(SoilModelConverter, method)(self)

    def fit(self, sm: Type[SoilModel], h: FloatArray, **kwargs):
        from .soil import SoilSample

        ss = SoilSample(h=h, theta=self.theta(h), k=self.k(h))
        return ss.fit(sm=sm, **kwargs)


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

    def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        return plot_swrc(sm=self, ax=ax)

    def convert(self, method: str | None = None) -> list[str] | SoilModel:

        if method is None:
            return SoilModelConverter.list_methods(sm=self)
        return getattr(SoilModelConverter, method)(self)

    def fit(self, sm: Type[SoilModel], h: FloatArray, **kwargs):
        from .soil import SoilSample

        ss = SoilSample(h=h, theta=self.theta(h), k=self.k(h))
        return ss.fit(sm=sm, **kwargs)


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

    def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        return plot_swrc(sm=self, ax=ax)

    def convert(self, method: str | None = None) -> list[str] | SoilModel:

        if method is None:
            return SoilModelConverter.list_methods(sm=self)
        return getattr(SoilModelConverter, method)(self)

    def fit(self, sm: Type[SoilModel], h: FloatArray, **kwargs):
        from .soil import SoilSample

        ss = SoilSample(h=h, theta=self.theta(h), k=self.k(h))
        return ss.fit(sm=sm, **kwargs)


@dataclass
class GenuchtenGardner:
    """Combination soil model using the van Genuchten soil water retention
    curve and the Gardner hydraulic conductivity function.

    Gardner, W.H. (1958) - Some steady-state solutions of the unsaturated
    moisture flow equation with application to evaporation from soils
    van Genuchten, M. Th. (1970) - A Closed-form Equation for Predicting the
    Hydraulic Conductivity of Unsaturated Soil
    """

    k_s: float
    theta_r: float
    theta_s: float
    alpha: float
    n: float
    c: float
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
        if s is not None:
            theta = s * (self.theta_s - self.theta_r) + self.theta_r
            h = self.h(theta)
        return exp(-self.c * npabs(h))

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        se = (theta - self.theta_r) / (self.theta_s - self.theta_r)
        h = 1 / self.alpha * ((1 / se) ** (1 / self.m) - 1) ** (1 / self.n)
        return h

    def plot(self, ax: MatplotlibAxes | None = None) -> MatplotlibAxes:
        return plot_swrc(sm=self, ax=ax)

    def convert(self, method: str | None = None) -> list[str] | SoilModel:

        if method is None:
            return SoilModelConverter.list_methods(sm=self)
        return getattr(SoilModelConverter, method)(self)

    def fit(self, sm: Type[SoilModel], h: FloatArray, **kwargs):
        from .soil import SoilSample

        ss = SoilSample(h=h, theta=self.theta(h), k=self.k(h))
        return ss.fit(sm=sm, **kwargs)


def get_soilmodel(
    soilmodel_name: SoilModelNames,
) -> Type[SoilModel]:
    sms = {
        "Genuchten": Genuchten,
        "Brooks": Brooks,
        "Haverkamp": Haverkamp,
        "Gardner": Gardner,
        "Rucker": Rucker,
        "Panday": Panday,
        "Fredlund": Fredlund,
        "GenuchtenGardner": GenuchtenGardner,
    }
    return sms[soilmodel_name]


def plot_swrc(
    sm: SoilModel, saturation: bool = False, ax: MatplotlibAxes | None = None, **kwargs
) -> MatplotlibAxes:
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
    ax: MatplotlibAxes | None = None,
    **kwargs,
) -> MatplotlibAxes:
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


class SoilModelConverter:
    """
    Namespace class for soil hydraulic model conversion methods.

    All conversion methods are static and take a SoilModel instance as their
    first argument ``sm``. Use :meth:`list_methods` to discover which conversions
    are available for a given model type — it reads the ``sm`` type annotation of
    each static method directly, so no extra registration is needed.

    Example
    --------
    >>> vg = Genuchten(k_s=1e-5, alpha=2.0, n=2.5, theta_s=0.45, theta_r=0.05)
    >>> SoilModelConverter.list_methods(vg)
    ['ghezzehei', 'morel', 'peche']
    >>> gardner = SoilModelConverter.peche(vg)
    >>> print(gardner.c)  # Output: 4.4

    """

    @classmethod
    def list_methods(cls, sm: SoilModel) -> list[str]:
        """
        Get all available conversion methods for the given SoilModel instance.

        Inspects the ``sm`` parameter type annotation of every static method on
        this class and returns those where ``isinstance(sm, annotation)`` is True.

        Parameters
        ----------
        sm : SoilModel
            The soil model instance to check conversions for.

        Returns
        -------
        list[str]
            Sorted list of method names applicable for the given model type.

        Example
        --------
        >>> vg = Genuchten(k_s=1e-5, alpha=2.0, n=2.5, theta_s=0.45, theta_r=0.05)
        >>> SoilModelConverter.list_methods(vg)
        ['ghezzehei', 'morel', 'peche']
        """
        from inspect import getmembers_static
        from typing import get_type_hints

        return sorted(
            name
            for name, obj in getmembers_static(cls)
            if isinstance(obj, staticmethod)
            and not name.startswith("_")
            and isinstance(sm, get_type_hints(obj.__func__).get("sm", ()))
        )

    @staticmethod
    def morel(sm: Genuchten | Brooks | Panday) -> Panday | Brooks | Genuchten:
        """
        Convert between van Genuchten and Brooks-Corey models using the Morel-Seytoux
        equivalence method, which preserves the effective capillary drive.

        This method implements the parameter equivalence from Morel-Seytoux et al. (1996),
        providing conversion between Brooks-Corey (BC) and van Genuchten (vG) parameters
        while preserving the maximum value of the effective capillary drive Hc.M,
        defined as the integral of relative permeability over capillary pressure.

        The conversion is based on two criteria:
        1. Preserve the effective capillary drive (primary criterion)
        2. Preserve asymptotic behavior at low water contents (secondary criterion)

        Parameters
        ----------
        sm : Genuchten or Brooks
            The soil model instance to convert.

        Returns
        -------
        Brooks or Genuchten
            Converted soil model instance:
            - If input is Genuchten, returns Brooks model
            - If input is Brooks, returns Genuchten model

        References
        ----------
        Morel-Seytoux, H. J., Meyer, P. D., Nachabe, M., Touma, J., van Genuchten, M. T.,
        & Lenhard, R. J. (1996). Parameter equivalence for the Brooks-Corey and van
        Genuchten soil characteristics: Preserving the effective capillary drive.
        Water Resources Research, 32(5), 1251-1258.
        https://doi.org/10.1029/96WR00069

        """
        if isinstance(sm, Genuchten):
            # Convert van Genuchten to Brooks-Corey
            # Calculate p from m (relationship 16a in Morel-Seytoux et al. 1996)
            # p = 1 + (2/m)
            m = sm.m
            p = 1 + (2 / m)

            # Calculate h_ce (entry/bubble pressure) using formula (17)
            # hce = (1/a) * 2p(p-1)/(p+3) * (147.8 + 8.1p + 0.092p^2) / (55.6 + 7.4p + p^2)
            a_inv = (
                (1 / sm.alpha)
                * (2 * p * (p - 1))
                / (p + 3)
                * (147.8 + 8.1 * p + 0.092 * p**2)
                / (55.6 + 7.4 * p + p**2)
            )

            # Calculate lambda (l) from p using Corey relationship (8b)
            # M = (p - 3) / 2, where M is lambda
            lamb = (p - 3) / 2

            return Brooks(
                k_s=sm.k_s,
                theta_r=sm.theta_r,
                theta_s=sm.theta_s,
                h_b=a_inv,  # h_ce maps to h_b
                l=lamb,  # lambda maps to l
            )
        elif isinstance(sm, Brooks):  # noqa: RET505
            # Convert Brooks-Corey to van Genuchten
            # Calculate p from l using Corey relationship (8a)
            # p = 3 + 2*M, where M is lambda (l)
            p = 3 + 2 * sm.l

            # Calculate m from p using relationship (16a)
            # m = 2 / (p - 1)
            m = 2 / (p - 1)

            # Calculate n from m using relationship (10b)
            # n = 1 / (1 - m)
            n = 1 / (1 - m)

            # Calculate alpha (1/a) from h_b using formula (18), rearranged
            # 1/a = hce * [2p(p-1)/(p+3)] * [(55.6 + 7.4p + p^2) / (147.8 + 8.1p + 0.092p^2)]
            # Therefore: a = hce / {[2p(p-1)/(p+3)] * [(55.6 + 7.4p + p^2) / (147.8 + 8.1p + 0.092p^2)]}
            # Or: alpha = 1/a = 1 / {hce * [2p(p-1)/(p+3)] * [(55.6 + 7.4p + p^2) / (147.8 + 8.1p + 0.092p^2)]}
            alpha = 1 / (
                sm.h_b
                * (2 * p * (p - 1))
                / (p + 3)
                * (55.6 + 7.4 * p + p**2)
                / (147.8 + 8.1 * p + 0.092 * p**2)
            )

            return Genuchten(
                k_s=sm.k_s,
                theta_r=sm.theta_r,
                theta_s=sm.theta_s,
                alpha=alpha,
                n=n,
            )

        else:
            raise TypeError(
                f"Unsupported model type: {type(sm).__name__}. "
                "Only Genuchten and Brooks soil models are supported by `morel` method."
            )

    @staticmethod
    def ghezzehei(sm: Genuchten) -> Gardner:
        """
        Converts van Genuchten model to Gardner model using Ghezzehei et al. (2007).

        Estimates the Gardner sorptive number c (= alpha_G) from van Genuchten's n
        and alpha parameters using eq. (17): c = 1.3 * n * alpha. The air entry
        pressure h_b is computed from the point of maximum downward concavity of
        the vGM retention curve using eq. (13).

        Parameters
        ----------
        sm : Genuchten
            The van Genuchten soil model instance to convert.

        Returns
        -------
        Gardner
            Gardner model instance with converted parameters, including h_b.

        Raises
        ------
        TypeError
            If sm is not a Genuchten model instance.
        ValueError
            If n < 2 (outside model validity range; eq. 13 requires m > 0.5).

        Notes
        -----
        Units: unit-agnostic. alpha and h_b share the same length unit;
        c has units of 1/[length]. k_s retains its original units.

        References
        ----------
        Ghezzehei, T. A., Kneafsey, T. J., & Su, G. W. (2007). Correspondence of
        the Gardner and van Genuchten-Mualem relative permeability function
        parameters. Water Resources Research, 43(10). https://doi.org/10.1029/2006WR005339

        Examples
        --------
        >>> vg = Genuchten(k_s=1e-5, alpha=2.0, n=2.5, theta_s=0.45, theta_r=0.05)
        >>> gardner = SoilModelConverter.ghezzehei(vg)
        >>> print(gardner.c) # Output: 6.5
        >>> print(round(gardner.h_b, 3)) # Output: 0.668
        """
        if not isinstance(sm, Genuchten):
            raise TypeError(
                f"Unsupported model type: {type(sm).__name__}. "
                "Only Genuchten soil models are supported by `ghezzehei` method."
            )

        # Validate n parameter — eq. (13) requires m > 0.5, i.e. n > 2
        if sm.n < 2:
            raise ValueError(
                f"n = {sm.n} is outside model validity range for Ghezzehei et al. (2007). "
                f"Method requires n >= 2 such that m > 0.5 (m = 1-1/n)."
            )

        m = sm.m  # m = 1 - 1/n

        # Calculate Gardner sorptive number c = alpha_G using eq. (17)
        c = 1.3 * sm.n * sm.alpha

        # Calculate air entry pressure h_b using eq. (13)
        # Derived from d^3 Theta / d psi^3 = 0 on the vGM retention curve
        inner = (5 * m - m**2 + (8 * m + 5 * m**2 - 2 * m**3 + m**4) ** 0.5) / (
            4 * m**2 - 2 * m
        )
        h_b = (1 / sm.alpha) * inner ** (m - 1)

        return Gardner(k_s=sm.k_s, theta_s=sm.theta_s, m=c, c=c, h_b=h_b)

    @staticmethod
    def peche(sm: Genuchten) -> Gardner:
        """
        Converts van Genuchten model to Gardner model using Peche et al. (in prep.).

        Estimates the Gardner c parameter from van Genuchten's alpha parameter
        using the empirical relationship c = 2.2 * alpha.

        Parameters
        ----------
        sm : Genuchten
            The van Genuchten soil model instance to convert.

        Returns
        -------
        Gardner
            Gardner model instance with converted parameters.

        Raises
        ------
        TypeError
            If sm is not a Genuchten model instance.
        ValueError
            If n < 1.8 (outside model validity range).

        References
        ----------
        Peche, A., Vonk, M.A., Altfelder, S., Houben, G. & Bakker, M. (in preparation).
        A new model for the approximation of the Gardner relative
        conductivity curve parameter based on van Genuchten's alpha.

        Examples
        --------
        >>> vg = Genuchten(k_s=1e-5, alpha=2.0, n=2.5, theta_s=0.45, theta_r=0.05)
        >>> gardner = SoilModelConverter.peche(vg)
        >>> print(gardner.c)
        4.4
        """
        if not isinstance(sm, Genuchten):
            raise TypeError(
                f"Unsupported model type: {type(sm).__name__}. "
                "Only Genuchten soil models are supported by `peche` method."
            )

        if sm.n < 1.8:
            raise ValueError(
                f"n = {sm.n} is outside model validity range for Peche et al. (in prep.). "
                f"Method requires n >= 1.8."
            )

        c = 2.2 * sm.alpha

        return Gardner(k_s=sm.k_s, theta_s=sm.theta_s, m=c, c=c)
