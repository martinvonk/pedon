# type: ignore
import logging
from bisect import bisect_right
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Self, Type

from numpy import abs as npabs
from numpy import (
    append,
    array,
    array2string,
    cos,
    divide,
    exp,
    full,
    log,
    log10,
    multiply,
    ndarray,
)
from pandas import DataFrame, isna, read_csv
from scipy.optimize import fixed_point, least_squares

from ._params import get_params
from ._typing import FloatArray
from .soilmodel import Brooks, Genuchten, SoilModel, get_soilmodel


@dataclass
class SoilSample:
    """
    A container for measured soil properties and in-situ soil hydraulic data, with
    convenience routines to predict hydraulic parameters using a range of pedo-
    transfer functions and empirical models.

    Attributes
    ----------
    sand_p : float | None
        Sand fraction in percent (%).
    silt_p : float | None
        Silt fraction in percent (%).
    clay_p : float | None
        Clay fraction in percent (%).
    rho : float | None
        Bulk density (g cm^-3).
    th33 : float | None
        Water content at -33 kPa (cm).
    th1500 : float | None
        Water content at -1500 kPa (cm).
    om_p : float | None
        Organic matter content in percent (%).
    m50 : float | None
        Median sand fraction (μm).
    d10 : float | None
        Representative grain diameter from sieve curve (m).
    d20 : float | None
        Representative grain diameter from sieve curve (m).
    h : FloatArray | None
        Pressure head measurements (cm).
    k : FloatArray | None
        Hydraulic conductivity measurements (m s^-1 or cm d^-1).
    theta : FloatArray | None
        Measured volumetric water content (dimensionless).

    Notes
    -----
        - Verify units before use as different methods expect specific conventions.
        - Methods return soil hydraulic model instances (e.g., Genuchten, Brooks).
        - Some methods make external API calls (rosetta) or use empirical relationships.
    """

    sand_p: float | None = None  # sand %
    silt_p: float | None = None  # silt %
    clay_p: float | None = None  # clay %
    rho: float | None = None  # soil density g/cm3
    th33: float | None = None  # cm
    th1500: float | None = None  # cm
    om_p: float | None = None  # organic matter %
    m50: float | None = None  # median sand fraction
    d10: float | None = (
        None  # d10 representative grain diameter (e.g. from sieve curve) [m]
    )
    d20: float | None = (
        None  # d20 representative grain diameter (e.g. from sieve curve) [m]
    )
    h: FloatArray | None = field(default=None, repr=False)  # pressure head measurement
    k: FloatArray | None = field(
        default=None, repr=False
    )  # hydraulic conductivity measurement
    theta: FloatArray | None = field(
        default=None, repr=False
    )  # moisture content measurement

    def from_staring(self, name: str, year: str = "2018") -> "SoilSample":
        """Get properties and measurements from Staring series"""
        if year not in ("2001", "2018"):
            raise ValueError(
                f"No Staring series available for year '{year}'"
                "please use either '2001' or '2018'"
            )
        path = Path(__file__).parent / "datasets/soilsamples.csv"
        properties = read_csv(path, delimiter=";")
        staring_properties = properties[
            properties["source"] == f"Staring_{year}"
        ].set_index("name")

        self.silt_p = staring_properties.loc[name, "silt_p"]
        self.clay_p = staring_properties.loc[name, "clay_p"]
        self.om_p = staring_properties.loc[name, "om_p"]
        self.m50 = staring_properties.loc[name, "m50"]
        if year == "2001":
            self.rho = staring_properties.loc[name, "rho"]
        return self

    def fit(
        self,
        sm: Type[SoilModel],
        pbounds: DataFrame | None = None,
        weights: FloatArray | float = 1.0,
        W1: float = 0.1,
        W2: float | None = None,
        k_s: float | None = None,
        silent: bool = True,
    ) -> SoilModel:
        """
        Fit the provided SoilModel (e.g., van Genuchten, Brooks-Corey class) to the
        stored measurements (theta, k, h) using nonlinear least squares. If
        pbounds is not provided, default parameter bounds are retrieved for the
        requested model name. The objective combines water retention and log10(k)
        errors; weighting terms W1 and W2 control the relative contribution of k.
        Returns a model instance with optimized parameters.

        Notes
        -----
            - Requires theta and k (and optionally h) to be set.
            - If k_s is provided it will be fixed during optimization.
            - Input/outputs and bounds are expected in the units used by the soil model.
        """
        theta = self.theta
        N = len(theta)
        k = self.k
        M = N + len(k)

        if pbounds is None:
            pbounds = get_params(sm.__name__)
            if k_s is not None:
                pbounds = pbounds.drop("k_s")
            else:
                pbounds.loc["k_s", "p_ini"] = max(k)
                pbounds.loc["k_s", "p_max"] = max(k) * 10
            pbounds.loc["theta_s", "p_ini"] = max(theta)
            pbounds.loc["theta_s", "p_max"] = max(theta) + 0.02

        if isinstance(weights, float):
            weights = full(M, weights)

        if W2 is None:
            W2 = (
                (M - N)
                * sum(weights[0:N] * theta)
                / (N * sum(weights[N:M] * npabs(log10(k))))
            )

        def get_diff(p: FloatArray) -> FloatArray:
            est_pars = dict(zip(pbounds.index, p))
            if k_s is not None:
                est_pars["k_s"] = k_s
            sml = sm(**est_pars)
            theta_diff = sml.theta(h=self.h) - theta
            k_diff = log10(sml.k(h=self.h)) - log10(k)
            diff = append(weights[0:N] * theta_diff, weights[N:M] * W1 * W2 * k_diff)
            return diff

        res = least_squares(
            get_diff,
            x0=pbounds.loc[:, "p_ini"],
            bounds=(
                pbounds.loc[:, "p_min"],
                pbounds.loc[:, "p_max"],
            ),
        )
        opt_pars = dict(zip(pbounds.index, res.x))
        if k_s is not None:
            opt_pars["k_s"] = k_s

        if not silent:
            print("SciPy Optimization Result\n", res)

        return sm(**opt_pars)

    def wosten(self, ts: bool = False) -> Genuchten:
        """Wosten et al (1999) - Development and use of a database of hydraulic
        properties of European soils"""

        topsoil = 1.0 * ts

        theta_s = (
            0.7919
            + 0.001691 * self.clay_p
            - 0.29619 * self.rho
            - 0.000001419 * self.silt_p**2
            + 0.0000821 * self.om_p**2
            + 0.02427 * self.clay_p**-1
            + 0.01113 * self.silt_p**-1
            + 0.01472 * log(self.silt_p)
            - 0.0000733 * self.om_p * self.clay_p
            - 0.000619 * self.rho * self.clay_p
            - 0.001183 * self.rho * self.om_p
            - 0.0001664 * topsoil * self.silt_p
        )
        alpha_ = (
            -14.96
            + 0.03135 * self.clay_p
            + 0.0351 * self.silt_p
            + 0.646 * self.om_p
            + 15.29 * self.rho
            - 0.192 * topsoil
            - 4.671 * self.rho**2
            - 0.000781 * self.clay_p**2
            - 0.00687 * self.om_p**2
            + 0.0449 * self.om_p**-1
            + 0.0663 * log(self.silt_p)
            + 0.1482 * log(self.om_p)
            - 0.04546 * self.rho * self.silt_p
            - 0.4852 * self.rho * self.om_p
            + 0.00673 * topsoil * self.clay_p
        )
        n_ = (
            -25.23
            - 0.02195 * self.clay_p
            + 0.0074 * self.silt_p
            - 0.1940 * self.om_p
            + 45.5 * self.rho
            - 7.24 * self.rho**2
            + 0.0003658 * self.clay_p**2
            + 0.002885 * self.om_p**2
            - 12.81 * self.rho**-1
            - 0.1524 * self.silt_p**-1
            - 0.01958 * self.om_p**-1
            - 0.2876 * log(self.silt_p)
            - 0.0709 * log(self.om_p)
            - 44.6 * log(self.rho)
            - 0.02264 * self.rho * self.clay_p
            + 0.0896 * self.rho * self.om_p
            + 0.00718 * topsoil * self.clay_p
        )
        l_ = (
            0.0202
            + 0.0006193 * self.clay_p**2
            - 0.001136 * self.om_p**2
            - 0.2316 * log(self.om_p)
            - 0.03544 * self.rho * self.clay_p
            + 0.00283 * self.rho * self.silt_p
            + 0.0488 * self.rho * self.om_p
        )
        ks_ = (
            7.755
            + 0.0352 * self.silt_p
            + 0.93 * topsoil
            - 0.967 * self.rho**2
            - 0.000484 * self.clay_p**2
            - 0.000322 * self.silt_p**2
            + 0.001 * self.silt_p**-1
            - 0.0748 * self.om_p**-1
            - 0.643 * log(self.silt_p)
            - 0.01398 * self.rho * self.clay_p
            - 0.1673 * self.rho * self.om_p
            + 0.02986 * topsoil * self.clay_p
            - 0.03305 * topsoil * self.silt_p
        )
        theta_r = 0.01
        return Genuchten(
            k_s=max(exp(ks_), 0),
            theta_r=theta_r,
            theta_s=theta_s,
            alpha=exp(alpha_),
            n=exp(n_) + 1,
            l=(10 * exp(l_) - 10) / (1 + exp(l_)),
        )

    def wosten_sand(self, ts: bool = False) -> Genuchten:
        """Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
        van boven- en ondergronden in Nederland: de Staringreeks"""

        topsoil = 1.0 * ts
        theta_s = (
            -35.7
            - 0.1843 * self.rho
            - 0.03576 * self.m50
            + 0.0000261 * self.m50**2
            - 0.0564 * self.silt_p**-1
            + 0.008 * self.om_p**-1
            + 496.0 * self.m50**-1
            + 0.02244 * log(self.om_p)
            + 7.56 * log(self.m50)
        )

        k_s = exp(
            45.8
            - 14.34 * self.rho
            + 0.001481 * self.silt_p**2
            - 27.5 * self.rho**-1
            - 0.891 * log(self.silt_p)
            - 0.34 * log(self.om_p)
        )

        alpha = exp(
            13.66
            - 5.91 * self.rho
            - 0.172 * topsoil
            + 0.003248 * self.m50
            - 11.89 * self.rho**-1
            - 2.121 * self.silt_p**-1
            - 0.3742 * log(self.silt_p)
        )

        n = max(
            exp(
                -1.057
                + 0.1003 * self.om_p
                + 1.119 * self.rho
                + 0.000764 * self.silt_p**2
                - 0.1397 * self.om_p**-1
                - 57.2 * self.m50**-1
                - 0.557 * log(self.om_p)
                - 0.02997 * self.rho * self.silt_p
            )
            + 1,
            1.0,
        )

        l_ = (
            -76.4
            - 0.097 * self.silt_p
            + 59.6 * self.rho
            + 0.0332 * self.m50
            - 13.45 * self.rho**2
            + 0.001127 * self.silt_p**2
            + 0.00431 * self.om_p**2
            - 0.0000399 * self.m50**2
            + 40.8 * self.rho**-1
            + 2.364 * self.silt_p**-1
            + 1.014 * log(self.silt_p)
        )
        l = min(max(2 * (exp(l_) - 1) / (exp(l_) + 1), -2.0), 2.0)  # noqa: E741

        return Genuchten(
            k_s=round(k_s, 4),
            theta_r=0.01,
            theta_s=round(theta_s, 4),
            alpha=round(alpha, 7),
            n=round(n, 4),
            l=round(l, 4),
        )

    def wosten_clay(self) -> Brooks:
        """Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
        van boven- en ondergronden in Nederland: de Staringreeks"""

        theta_s = (
            0.6311
            + 0.003383 * self.clay_p
            - 0.09699 * self.rho**2
            - 0.00204 * self.rho * self.clay_p
        )

        k_s = exp(
            -42.6
            + 8.71 * self.om_p
            + 61.9 * self.rho
            - 20.79 * self.rho**2
            - 0.2107 * self.om_p**2
            - 0.01622 * self.clay_p * self.om_p
            - 5.382 * self.rho * self.om_p
        )

        alpha = exp(
            -19.13
            + 0.812 * self.om_p
            + 23.4 * self.rho
            - 8.16 * self.rho**2
            + 0.423 * self.om_p**-1
            + 2.388 * log(self.om_p)
            - 1.338 * self.rho * self.om_p
        )

        n = max(
            exp(
                -0.235
                + 0.972 * self.rho**-1
                - 0.7743 * log(self.clay_p)
                - 0.3154 * log(self.om_p)
                + 0.0678 * self.rho * self.om_p
            )
            + 1.0,
            1.0,
        )

        l_ = 0.102 + 0.0222 * self.clay_p - 0.043 * self.rho * self.clay_p
        l = min(max(10.0 * (exp(l_) - 1) / (exp(l_) + 1), -10.0), 10.0)  # noqa: E741
        return Genuchten(
            k_s=round(k_s, 4),
            theta_r=0.01,
            theta_s=round(theta_s, 4),
            alpha=round(alpha, 7),
            n=round(n, 4),
            l=round(l, 4),
        )

    def cosby(self) -> Brooks:
        """Cooper (2021) - Using data assimilation to optimize pedotransfer
        functions using field-scale in situ soil moisture observations"""
        c = self.clay_p / 100
        s = self.sand_p / 100
        b = 3.100 + 15.70 * c - 0.300 * s
        theta_s = 0.505 - 0.037 * c - 0.142 * s
        psi_s = 0.01 * 10 ** (2.170 - 0.630 * c - 1.580 * s)
        k_s = 10 ** (-0.600 - 0.640 * c + 1.260 * s) * 25.2 / 3600
        labda = 1 / b
        k_s = k_s * 8640000 / 1000  # kg/m2/s to cm/d
        psi_s = psi_s * 100  # m to cm
        return Brooks(
            k_s=round(k_s, 4),
            theta_r=0.0,
            theta_s=round(theta_s, 4),
            h_b=round(psi_s, 5),
            l=round(labda, 5),
        )

    def rosetta(self, version: Literal[1, 2, 3] = 3) -> Genuchten:
        """Rosetta (Schaap et al., 2001) - Predicting soil water retention from soil"""
        try:
            from httpx import post as httpx_post
        except ImportError:
            raise ImportError(
                "httpx is required for the rosetta method to make api calls. "
                "Please install it with 'pip install httpx'."
            )

        soildata = [
            self.sand_p,  # %
            self.silt_p,  # %
            self.clay_p,  # %
            -9.9 if self.rho is None else self.rho,  # g/cm3
            -9.9 if self.th33 is None else self.th33,  # [-]
            -9.9 if self.th1500 is None else self.th1500,  # [-]
        ]

        data = {"soildata": [soildata]}
        r = httpx_post(
            f"https://www.handbook60.org/api/v1/rosetta/{version}",
            json=data,
        )
        if r.is_error:
            raise ValueError(f"Rosetta API error: {r}")

        rjson = r.json()
        vgpar = rjson["van_genuchten_params"][0]
        # stdev = rjson["stdev"][0] # TODO: return extra Genuchten classes with stdev or report/print
        if None in vgpar:
            raise ValueError(f"Rosetta API returned None values: {vgpar}")
        return Genuchten(
            k_s=10 ** vgpar[4],
            theta_r=vgpar[0],
            theta_s=vgpar[1],
            alpha=10 ** vgpar[2],
            n=10 ** vgpar[3],
        )

    def hypags(self) -> Genuchten:
        """
        Estimate van Genuchten parameters using the HYPAGS method.

        Implements the Kozeny-Carman-based parameterization from Peche & Houben
        (2023, 2024) to derive van Genuchten parameters and characteristic grain-size
        metrics from hydraulic conductivity and/or grain size inputs.

        The routine:
        - Accepts k (hydraulic conductivity), d10, or d20 as input (at least one required)
        - Iteratively estimates effective porosity (ne), d50, and d60
        - Computes van Genuchten alpha and n parameters with bootstrap-derived errors
        - Returns a Genuchten model with k_s, theta_r, theta_s, alpha, and n

        Notes
        -----
        - Grain diameters should be in meters; k in m s^-1
        - Input values outside empirically supported ranges trigger warnings
        - theta_r is estimated from k; theta_s equals effective porosity (ne)

        Returns
        -------
        Genuchten
            Van Genuchten soil model with estimated parameters.

        References
        ----------
        Peche, A., Houben, G., & Altfelder, S. (2024). Approximation of van Genuchten
        Parameter Ranges from Hydraulic Conductivity Data. Groundwater, 62(3), 469-479.

        Peche, A., & Houben, G. J. (2023). Estimating characteristic grain sizes and
        effective porosity from hydraulic conductivity data. Groundwater, 61(4), 574-585.
        """
        # constants and coefficients
        rho_f = 999.7  # fluid density [kg/m³] (assumed 20°C)
        g = 9.81  # gravitational acceleration [m/s²]
        mu = 1.1306e-3  # dynamic viscosity [kg/(m·s)]
        c = mu / (rho_f * g)
        a1, a2, a3 = 1.00401066e00, 1.50991411e-04, 5.78888587e-03  # alpha-coefficients
        gamma = 0.0728  # surface tension [J/m²]
        theta = 0.0  # contact angle, parameter for Young-Laplace eq.
        factor = (
            2 * gamma * cos(theta) / (rho_f * g)
        )  # factor Young-Laplace eq. (fluid properties, contact angle, surface tension)
        Pi, c1, c2 = (
            0.0009,
            1.2,
            1.13,
        )  # dimensionless number Pi and c-coefficients --> [Kozeny-Carman-based parameterization]

        # hypags functions
        def solve_kozeny_carman(
            k: float, ne_i: float, d50_i: float, const: float
        ) -> tuple[float, float]:
            """
            Solve for d50 and effective porosity based on Kozeny-Carman equation
            using scipy's fixed_point for robust successive substitution iteration.

            Parameters (Parametrization-dependent)
            ----------
            k : input hydraulic conductivity
            ne_i : initial guess for effective porosity
            d50_i : initial guess for d50
            const : fluid properties and gravitational acceleration constant

            Returns
            -------
            ne : effective porosity
            d50 : d50
            """
            c = 180 * const

            def update_func(vars):
                """
                Fixed-point iteration function.
                Computes the next iteration values for [n, d] given current values.
                """
                n, d = vars
                n_next = (k / (d**2) * c * (1 - n) ** 2) ** (1 / 3)
                d_next = (c * k * ((1 - n_next) ** 2 / (n_next**3))) ** (1 / 2)
                return array([n_next, d_next])

            solution = fixed_point(
                update_func, array([ne_i, d50_i]), xtol=1e-10, maxiter=100
            )
            ne, d50 = solution

            return ne, d50

        def get_alpha_errors(k: float) -> tuple[float, float]:
            """
            Return bootstrap confidence interval errors for a given k value.

            Uses binary search to quickly find the corresponding [+Δα, -Δα] pair
            based on predefined k thresholds.

            Parameters
            ----------
            k : float
                Hydraulic conductivity.

            Returns
            -------
            tuple[float, float]
                (added_error, subtracted_error)
            """
            bounds = [
                0.00000025,
                0.0000005,
                0.00000075,
                0.000001,
                0.0000025,
                0.000005,
                0.0000075,
                0.00001,
                0.000025,
                0.00005,
                0.000075,
                0.0001,
                0.0005,
            ]
            errors = [
                [3.002, 1.144],
                [5.364, 1.258],
                [7.031, 1.366],
                [6.728, 1.519],
                [6.758, 2.072],
                [6.607, 2.393],
                [5.91, 2.985],
                [5.643, 3.229],
                [4.84, 4.0],
                [4.082, 4.233],
                [3.182, 4.486],
                [3.051, 5.03],
                [3.442, 5.483],
            ]
            default_error = (1.811, 2.858)

            i = bisect_right(bounds, k)
            return tuple(errors[i - 1]) if i <= len(errors) else default_error

        def get_n_errors(alpha: float) -> tuple[float, float]:
            """
            Return bootstrap confidence interval errors for a given alpha value.

            Uses binary search to find the corresponding [+Δn, -Δn] pair
            based on predefined alpha thresholds.

            Parameters
            ----------
            alpha : float
            Alpha parameter.

            Returns
            -------
            tuple[float, float]
            (added_error, subtracted_error)
            """

            bounds = [1, 2, 3, 4, 5, 6, 7, 8, 9]
            errors = [
                [2.836, 0.358],
                [3.061, 0.826],
                [3.481, 0.828],
                [3.182, 0.605],
                [1.905, 0.509],
                [1.184, 0.297],
                [4.032, 0.352],
                [1.385, 0.241],
                [1.017, 0.155],
            ]
            default_error = (3.006, 0.377)

            i = bisect_right(bounds, alpha)
            return tuple(errors[i - 1]) if i <= len(errors) else default_error

        def get_res_water_content(k: float) -> tuple[float, float, float]:
            """
            Return residual water content (θr) with lower and upper bounds for a given k value.

            Uses binary search to find the corresponding [θr_lower, θr, θr_upper]
            triple based on predefined k thresholds.

            Parameters
            ----------
            k : float
                Hydraulic conductivity.

            Returns
            -------
            tuple[float, float, float]
                [θr_lower, θr, θr_upper]
            """
            bounds = [
                0.00000025,
                0.0000005,
                0.00000075,
                0.000001,
                0.0000025,
                0.000005,
                0.0000075,
                0.00001,
                0.000025,
                0.00005,
                0.000075,
                0.0001,
                0.0005,
            ]
            values = [
                [0.089, 0.097, 0.198],
                [0.097, 0.107, 0.207],
                [0.09, 0.099, 0.227],
                [0.095, 0.102, 0.221],
                [0.091, 0.098, 0.216],
                [0.087, 0.097, 0.227],
                [0.104, 0.111, 0.212],
                [0.099, 0.108, 0.219],
                [0.079, 0.085, 0.233],
                [0.082, 0.087, 0.226],
                [0.067, 0.07, 0.223],
                [0.105, 0.11, 0.23],
                [0.096, 0.105, 0.23],
            ]
            default_value = (0.068, 0.111, 0.096)

            i = bisect_right(bounds, k)
            return tuple(values[i - 1]) if i <= len(values) else default_value

        # mathematical model of hypags calculation of k, d10, d20:
        # CONDITION: in HYPAGS, either k, d10 or d20 have to be given
        if self.k is not None:
            if isinstance(self.k, float):
                k = self.k
            elif isinstance(self.k, ndarray):
                if len(self.k) > 1:
                    logging.warning(
                        "HYPAGS routine only accepts single k value, choosing the first k value in the array."
                    )
                k = float(self.k[0])
            else:
                k = float(self.k)

            # case 0: mathematical model where k is given
            # check for non-valid input
            if k > 2.6e-2 or k < 2.87e-7:
                logging.error("k out of hypags model limits.")
            logging.debug("Using case 0 of hypags model (k given).")
            self.d10 = (k / Pi * c) ** (0.5)  # calculation of d10
            self.d20 = c1 * self.d10  # calculation of d20
        elif self.d10 is not None:
            # case 1: mathematical model where d10 is given
            if not 5.35e-5 <= self.d10 <= 8.3e-4:
                logging.error(
                    f"d10 ({self.d10:.3e}) out of hypags model limits: 5.35e-5 to 8.3e-4."
                )
            logging.debug("Using case 1 of hypags model (d10 given).")
            self.d20 = c1 * self.d10  # calculation of d20
            k = (Pi * rho_f * g * self.d10**2) / mu
            self.k = array([k], dtype=float)  # calculation of k
        elif self.d20 is not None:
            # case 2: mathematical model where d20 is given
            if not 6.25e-5 <= self.d20 <= 1.2e-3:
                logging.error(
                    f"d20 ({self.d20:.3e}) out of hypags model limits: 6.25e-5 to 1.2e-3."
                )
            logging.debug("Using case 2 of hypags model (d20 given).")
            self.d10 = self.d20 / c1  # calculation of d10
            k = (Pi * rho_f * g * self.d10**2) / mu
            self.k = array([k], dtype=float)  # calculation of k
        else:
            raise ValueError(
                "No parameter (k, d10, or d20) was provided for hypags routine."
            )

        # calculation of d50 and ne
        d50_0 = a3 * self.d20  # starting value for iterative approximation of d50
        ne_0 = a1 * k**a2
        ne, d50 = solve_kozeny_carman(k=k, ne_i=ne_0, d50_i=d50_0, const=c)

        # calculation of d60
        d60 = c2 * d50

        # calculation of van Genuchten alpha
        drep = d60 if k <= 5e-5 else 0.0001803 + k * 0.10
        h = factor * divide(
            1, multiply(0.5, drep ** (1))
        )  # Note the conversion of diameter to radius
        alpha = divide(1, h)

        # calculation of van Genuchten n
        nt = (
            1.12411782 + alpha * 0.55750592
            if alpha <= 1.9
            else 1.67561574 / (alpha - 0.30307062) + 1.16193138
        )

        # error range of van Genuchten parameters.
        # TODO: implement if needed #https://github.com/martinvonk/pedon/issues/17
        # van Genuchten alpha
        alpha_errors = get_alpha_errors(k)
        alpha_plus_dalpha = alpha + alpha_errors[0]
        alpha_minus_dalpha = max(
            alpha - alpha_errors[1], 0.1
        )  # ensure a positive lower bound for alpha
        _, _ = (
            alpha_plus_dalpha,
            alpha_minus_dalpha,
        )  # don't use for now, see TODO above
        # van Genuchten n
        n_errors = get_n_errors(alpha)
        n_plus_dnt = nt + n_errors[0]
        n_minus_dnt = max(nt - n_errors[1], 0.1)  # ensure a positive lower bound for n
        _, _ = n_plus_dnt, n_minus_dnt  # don't use for now, see TODO above
        # residual water content theta_r
        thetar_minus_dthetar, thetar, thetar_plus_dthetar = get_res_water_content(k)
        _, _ = (
            thetar_minus_dthetar,
            thetar_plus_dthetar,
        )  # don't use for now, see TODO above

        return Genuchten(
            k_s=k,
            theta_r=thetar,
            theta_s=ne,
            alpha=alpha,
            n=nt,
        )


@dataclass
class Soil:
    """
    A class representing soil properties and models.

    Attributes
    ----------
    name : str
        The name identifier for the soil.
    model : SoilModel | None, optional
        The soil model instance containing hydraulic parameters, by default None.
    sample : SoilSample | None, optional
        The soil sample data associated with this soil, by default None.
    source : str | None, optional
        The data source for the soil parameters (e.g., 'HYDRUS', 'Staring_2018'), by default None.
    description : str | None, optional
        A text description of the soil type, by default None.
    """

    name: str
    model: SoilModel | None = None
    sample: SoilSample | None = None
    source: str | None = None
    description: str | None = None

    def from_name(
        self, sm: Type[SoilModel] | SoilModel | str, source: str | None = None
    ) -> Self:
        """Load soil parameters from a CSV database by soil name and model type."""
        if isinstance(sm, SoilModel):
            if hasattr(sm, "__name__"):
                smn = sm.__name__
            else:
                smn = sm.__class__.__name__
                sm = type(sm)
        elif isinstance(sm, str):
            smn = sm
            sm = get_soilmodel(smn)
        else:
            raise ValueError(
                f"Argument must either be Type[SoilModel] | SoilModel | str,"
                f"not {type(sm)}"
            )

        path = Path(__file__).parent / "datasets/soilsamples.csv"
        ser = read_csv(path, delimiter=";", index_col=0)
        if "HYDRUS_" in self.name:
            self.name = self.name.replace("HYDRUS_", "")
            source = "HYDRUS"
            logging.warning(
                "Removed 'HYDRUS_' from soil name. For future use"
                "please provide source='HYDRUS' argument"
            )
        sersm = ser[ser["soilmodel"] == smn].loc[[self.name], :]
        if source is None and len(sersm) > 1:
            raise Exception(
                f"Multiple sources for soil {self.name}: "
                f"{array2string(sersm.loc[:, 'source'].values)}. "
                f"Please provide the source using the source argument"
            )
        elif (source is not None) and len(sersm) > 1:
            sersm = sersm[sersm["source"] == source]

        serd = sersm.squeeze().to_dict()
        if isna(serd["description"]):
            serd["description"] = serd["soil type"]
        self.__setattr__("source", serd.pop("source"))
        self.__setattr__("description", serd.pop("description"))
        smserd = {
            x: serd[x]
            for x in sm.__dataclass_fields__.keys()
            if sm.__dataclass_fields__[x].init
        }
        self.__setattr__("model", sm(**smserd))
        return self

    @staticmethod
    def list_names(sm: Type[SoilModel] | SoilModel | str) -> list[str]:
        """Return a list of available soil names for a given soil model."""
        if isinstance(sm, SoilModel):
            if hasattr(sm, "__name__"):
                smn = sm.__name__
            else:
                smn = sm.__class__.__name__
                sm = type(sm)
        elif isinstance(sm, str):
            smn = sm
            sm = get_soilmodel(smn)

        else:
            raise ValueError(
                f"Argument must either be Type[SoilModel] | SoilModel | str,"
                f"not {type(sm)}"
            )

        path = Path(__file__).parent / "datasets/soilsamples.csv"
        names = read_csv(path, delimiter=";")

        return names[names["soilmodel"] == smn].loc[:, "name"].unique().tolist()

    def from_staring(self, year: str = "2018") -> Self:
        """Load soil parameters from the Staring series database."""
        if year not in ("2001", 2001, "2018", 2018):
            raise ValueError(f"Year must either be '2001' or '2018', not {year}")

        self.from_name(sm=Genuchten, source=f"Staring_{year}")
        ss = SoilSample().from_staring(name=self.name, year=year)
        self.__setattr__("sample", ss)
        return self
