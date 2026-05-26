"""Soil models for soil water retention and hydraulic conductivity."""

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

import matplotlib.pyplot as plt
from numpy import abs as npabs
from numpy import asarray, exp, full, linspace, log, log10, maximum, sqrt
from scipy.integrate import trapezoid
from scipy.special import erfc, erfcinv, lambertw

from ._typing import FloatArray, SoilModelNames
from .plot import swrc as plot_swrc


@runtime_checkable
class SoilModel(Protocol):
    """Protocol for soil models.

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
        ...     def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        ...         # Plot the soil water retention curve by calling `pe.plot.swrc`
        ...         return pe.plot.swrc(self, ax=ax)

    """

    def theta(self, h: FloatArray) -> FloatArray:
        """Calculate soil moisture content (water content) from pressure head h."""
        ...

    def s(self, h: FloatArray) -> FloatArray:
        """Calculate effective saturation from pressure head h.

        Effective saturation is normalized between 0 (dry) and 1 (saturated).
        """
        ...

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        """Calculate relative permeability (or relative hydraulic conductivity).

        Relative permeability is normalized between 0 (dry) and 1 (saturated).
        Can be calculated from either pressure head h or saturation s.
        """
        ...

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        """Calculate hydraulic conductivity from pressure head h or saturation s."""
        ...

    def h(self, theta: FloatArray) -> FloatArray:
        """Calculate pressure head h from soil moisture content (inverse of theta)."""
        ...

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        """Plot the soil water retention curve by calling `plot_swrc`."""
        ...


@dataclass
class Genuchten:
    """Mualem-van Genuchten Soil Model.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_r: float
        Residual soil water content [-]
    theta_s: float
        Saturated soil water content [-]
    n: float
        Dimensionless measure of the pore-size distribution
    alpha: float
        Empirical shape parameter that is physically related
        to the inverse of the air-entry pressure [1/L]
    l: float, optional
        Pore-connectivity or tortuosity parameter [-]
        Default value 0.5 is commonly used and is based on the Mualem model.

    Attributes
    ----------
    m: float, optional
        Calculated from n as m = 1 - 1/n [-] (Mualem restriction)

    References
    ----------
    van Genuchten, M. Th. (1980) - A Closed-form Equation for Predicting the
    Hydraulic Conductivity of Unsaturated Soil.

    Mualem, Y. (1976). A new model for predicting the hydraulic conductivity
    of unsaturated porous media.

    """

    k_s: float
    theta_r: float
    theta_s: float
    alpha: float
    n: float
    l: float = 0.5  # noqa: E741
    m: float = field(init=False, repr=False)

    def __post_init__(self):
        """Calculate m from n after initialization."""
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
    """Brooks and Corey Soil Model.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_r: float
        Residual soil water content [-]
    theta_s: float
        Saturated soil water content [-]
    h_b: float
        Bubbling pressure or air-entry suction [L]
    l: float
        Pore-size distribution index [-]

    References
    ----------
    Brooks, R.H. and Corey, A.T. (1964) - Hydraulic Properties of Porous Media.

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
                (theta[mask] - self.theta_r) / (self.theta_s - self.theta_r)
            ) ** (-1 / self.l)
            return h

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class Haverkamp:
    """Haverkamp Soil Model.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_r: float
        Residual soil water content [-]
    theta_s: float
        Saturated soil water content [-]
    alpha: float
        Empirical constant for the soil water retention curve [-]
    beta: float
        Empirical exponent for the soil water retention curve.
        In this implementation, it is also used as the exponent
        for the hydraulic conductivity curve (assuming B = beta) [-]
    a: float
        Empirical constant for the hydraulic conductivity curve
        (Often denoted as 'A' in literature)

    References
    ----------
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
    """Gardner(-Kozeny) Soil Model.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_s: float
        Saturated soil water content [-]
    m: float
        Empirical exponent parameter for the soil water retention curve [1/L]
    c: float
        Empirical exponent parameter for the hydraulic conductivity curve [1/L]
        (Often denoted as alpha in the classic Gardner exponential model)

    References
    ----------
    Gardner, W. H. (1958). Some steady-state solutions of the unsaturated
    moisture flow equation with application to evaporation from soils.

    Bakker, M., & Nieber, J. L. (2009). Damping of sinusoidal surface flux
    fluctuations with soil depth.

    """

    k_s: float
    theta_s: float
    m: float
    c: float

    def theta(self, h: FloatArray) -> FloatArray:
        return self.theta_s * exp(-self.m * npabs(h))

    def s(self, h: FloatArray) -> FloatArray:
        return self.theta(h) / self.theta_s

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is not None:
            theta = s * self.theta_s
            h = self.h(theta)
        return exp(-self.c * npabs(h))

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        return -(1.0 / self.m) * log(theta / self.theta_s)

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class Rucker:
    """Gardner(-Rucker) Soil Model.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_r: float
        Residual soil water content [-]
    theta_s: float
        Saturated soil water content [-]
    m: float
        Empirical shape parameter for the soil water retention curve [-]
    c: float
        Empirical parameter for the hydraulic conductivity curve [1/L]

    References
    ----------
    Gardner, W. H. (1958). Some steady-state solutions of the unsaturated
    moisture flow equation with application to evaporation from soils.

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
        se = (theta - self.theta_r) / (self.theta_s - self.theta_r)
        y = se ** ((self.m + 2.0) / 2.0)
        # Invert (1 + x) * exp(-x) = y via x = -W_{-1}(-y / e) - 1.
        x = -lambertw(-y / exp(1.0), k=-1).real - 1.0
        return (2.0 / self.c) * x

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class Panday:
    """Panday Soil Model as used in MODFLOW-USG Transport and MODFLOW UZF.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
        Note: The empirical specific yield calculation assumes this is provided
        in cm-based units (e.g., cm/day or cm/s).
    theta_r: float
        Residual soil water content [-]
    theta_s: float
        Saturated soil water content [-]
    alpha: float
        Inverse of the air-entry suction [1/L]
        (Represents the van Genuchten alpha parameter)
    beta: float
        Pore-size distribution parameter [-]
        (Represents the van Genuchten n parameter)
    brook: float
        Exponent for the relative permeability function [-]
        (Represents the Brooks-Corey pore-connectivity parameter)
    h_b: float, optional
        Bubbling pressure or air-entry offset [L] (default: 0.0)
    ss: float, optional
        Specific storage [1/L] (default: 1e-6)

    Attributes
    ----------
    sr: float
        Residual saturation calculated as theta_r / theta_s [-]
    gamma: float
        Calculated from beta as 1 - 1 / beta [-]
        (Represents the Mualem m restriction)
    sy: float
        Specific yield, calculated empirically based on k_s and beta [-]

    References
    ----------
    Panday, S. (2026). USG-Transport: Transport and other Enhancements to
    MODFLOW-USG. Documentation and User's Guide.

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
        """Calculate additional parameters after initialization."""
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
        h = self.h_b + (1.0 / self.alpha) * (
            (1.0 / se) ** (1.0 / self.gamma) - 1.0
        ) ** (1.0 / self.beta)
        return h

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class Fredlund:
    """Fredlund and Xing Soil Model.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_s: float
        Saturated soil water content [-]
    a: float
        Soil parameter related to the air-entry value [L]
    n: float
        Empirical soil parameter controlling the slope of the SWCC [-]
    m: float
        Empirical soil parameter related to the residual water content [-]

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
                "Can only calculate the hydraulic conductivity "
                "using the pressure head, not the saturation"
            )

        def theta_d(
            h_val: FloatArray, a: float, n: float, m: float, theta_s: float
        ) -> FloatArray:
            """Compute the derivative of theta."""
            return -(
                theta_s
                * m
                * n
                * (h_val / a) ** n
                * (log((h_val / a) ** n + exp(1))) ** (-m - 1)
            ) / (h_val * (h_val / a) ** n + exp(1) * h_val)

        h_b = 0.03
        n_steps = 100

        h = asarray(npabs(h), dtype=float)
        # Integrate along the integration steps (axis=0) using trapezoidal rule
        y_num = linspace(log(h), log(1e6), n_steps)  # shape: (n_steps, len(h_arr))
        exp_y_num = exp(y_num)
        integrand_num = (
            (self.theta(exp_y_num) - self.theta(h))
            / exp_y_num
            * theta_d(exp_y_num, self.a, self.n, self.m, self.theta_s)
        )  # shape: (n_steps, len(h_arr))
        numerator = trapezoid(integrand_num, x=y_num, axis=0)

        y_den = linspace(log(h_b), log(1e6), n_steps)  # shape: (n_steps,)
        exp_y_den = exp(y_den)
        integrand_den = (
            (self.theta(exp_y_den) - self.theta_s)
            / exp_y_den
            * theta_d(exp_y_den, self.a, self.n, self.m, self.theta_s)
        )  # shape: (n_steps,)
        denominator = trapezoid(integrand_den, x=y_den)

        result = numerator / denominator

        # Return a clean scalar if the user originally passed a scalar instead of an array
        return result

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        return self.a * (exp((self.theta_s / theta) ** (1.0 / self.m)) - exp(1)) ** (
            1.0 / self.n
        )

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class Kosugi:
    """Kosugi 2-Parameter Lognormal Soil Model.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_r: float
        Residual soil water content [-]
    theta_s: float
        Saturated soil water content [-]
    h_m: float
        Median capillary pressure head [L]
    sigma: float
        Standard deviation of log-transformed soil pore radii [-]

    References
    ----------
    Kosugi, K. (1996). Lognormal distribution model for unsaturated soil
    hydraulic properties.

    """

    k_s: float
    theta_r: float
    theta_s: float
    h_m: float
    sigma: float

    def theta(self, h: FloatArray) -> FloatArray:
        se = 0.5 * erfc(log(h / self.h_m) / (sqrt(2.0) * self.sigma))
        return self.theta_r + se * (self.theta_s - self.theta_r)

    def s(self, h: FloatArray) -> FloatArray:
        return (self.theta(h) - self.theta_r) / (self.theta_s - self.theta_r)

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)

        return s**0.5 * (0.5 * erfc(erfcinv(2.0 * s) + self.sigma / sqrt(2.0))) ** 2

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        se = (theta - self.theta_r) / (self.theta_s - self.theta_r)
        return self.h_m * exp(sqrt(2.0) * self.sigma * erfcinv(2.0 * se))

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class Campbell:
    """Campbell Soil Model.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_s: float
        Saturated soil water content [-]
    h_e: float
        Air-entry potential or bubbling pressure [L]
    b: float
        Empirical shape parameter controlling the slope of the curve [-]

    References
    ----------
    Campbell, G. S. (1974). A simple method for determining unsaturated
    conductivity from moisture retention data. Soil Science, 117(6), 311-314.

    """

    k_s: float
    theta_s: float
    h_e: float
    b: float

    def theta(self, h: FloatArray) -> FloatArray:
        # Use maximum to bound h to the air-entry pressure (h_e).
        h = maximum(npabs(h), self.h_e)
        return self.theta_s * (self.h_e / h) ** (1.0 / self.b)

    def s(self, h: FloatArray) -> FloatArray:
        return self.theta(h) / self.theta_s

    def k_r(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)

        return s ** (2.0 * self.b + 3.0)

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        return self.k_s * self.k_r(h=h, s=s)

    def h(self, theta: FloatArray) -> FloatArray:
        return self.h_e * (theta / self.theta_s) ** -self.b

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class GenuchtenGardner:
    """Combination soil model.

    Uses the van Genuchten soil water retention curve and the Gardner
    hydraulic conductivity function.

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_r: float
        Residual soil water content [-]
    theta_s: float
        Saturated soil water content [-]
    alpha: float
        Inverse of the air-entry pressure [1/L]
        (van Genuchten parameter)
    n: float
        Pore-size distribution parameter [-]
        (van Genuchten parameter)
    c: float
        Empirical parameter for the hydraulic conductivity curve [1/L]
        (Gardner parameter)

    Attributes
    ----------
    m: float
        Calculated from n as m = 1 - 1/n [-] (Mualem restriction)

    References
    ----------
    Gardner, W. H. (1958). Some steady-state solutions of the unsaturated
    moisture flow equation with application to evaporation from soils.

    van Genuchten, M. Th. (1980). A closed-form equation for predicting the
    hydraulic conductivity of unsaturated soils.

    """

    k_s: float
    theta_r: float
    theta_s: float
    alpha: float
    n: float
    c: float
    m: float = field(init=False, repr=False)

    def __post_init__(self):
        """Calculate m from n after initialization."""
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

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


@dataclass
class GenuchtenKool:
    """Soil model with a scaling factor for hysteresis.

    Uses the van Genuchten soil water retention curve and hydraulic
    conductivity function with a scaling factor xi to approximate hysteresis
    behavior as described in Kool and Parker (1987).

    Parameters
    ----------
    k_s: float
        Saturated hydraulic conductivity [L/T]
    theta_r: float
        Residual soil water content [-]
    theta_s: float
        Saturated soil water content [-]
    alpha: float
        Inverse of the air-entry pressure for the main drying curve [1/L]
    n: float
        Pore-size distribution parameter [-]
    l: float, optional
        Pore-connectivity or tortuosity parameter [-] (default: 0.5)
    xi: float, optional
        Hysteresis scaling factor, typically defined as alpha_wetting / alpha_drying [-]
        (default: 2.0)

    Attributes
    ----------
    m: float
        Calculated from n as m = 1 - 1/n [-] (Mualem restriction)
    alpha_w: float
        Scaled alpha parameter for the water retention curve [1/L]

    References
    ----------
    Kool, J. B., & Parker, J. C. (1987). Development and evaluation of
    closed-form expressions for hysteretic soil hydraulic properties.

    van Genuchten, M. Th. (1980). A closed-form equation for predicting the
    hydraulic conductivity of unsaturated soils.

    """

    k_s: float
    theta_r: float
    theta_s: float
    alpha: float
    n: float
    l: float = 0.5  # noqa: E741
    m: float = field(init=False, repr=False)
    xi: float = field(default=2.0, repr=True)

    def __post_init__(self):
        """Calculate m from n after initialization."""
        self.m = 1 - 1 / self.n

    @property
    def alpha_w(self) -> float:
        """Alpha parameter for the water retention curve, which is scaled by xi."""
        return self.alpha * self.xi

    def theta(self, h: FloatArray) -> FloatArray:
        theta = (
            self.theta_r
            + (self.theta_s - self.theta_r)
            / (1 + npabs(self.alpha_w * h) ** self.n) ** self.m
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
        h = 1 / self.alpha_w * ((1 / se) ** (1 / self.m) - 1) ** (1 / self.n)
        return h

    def plot(self, ax: plt.Axes | None = None) -> plt.Axes:
        return plot_swrc(self, ax=ax)


def get_soilmodel(
    soilmodel_name: SoilModelNames,
) -> type[SoilModel]:
    """Get soil model class by name."""
    sms = {
        "Genuchten": Genuchten,
        "Brooks": Brooks,
        "Haverkamp": Haverkamp,
        "Gardner": Gardner,
        "Rucker": Rucker,
        "Panday": Panday,
        "Fredlund": Fredlund,
        "Kosugi": Kosugi,
        "Campbell": Campbell,
        "GenuchtenGardner": GenuchtenGardner,
        "GenuchtenKool": GenuchtenKool,
    }
    return sms[soilmodel_name]
