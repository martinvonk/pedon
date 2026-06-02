"""Soil sample class and pedotransfer functions."""

from bisect import bisect_right
from dataclasses import dataclass, field, fields
from logging import getLogger
from pathlib import Path
from typing import Any, Literal, Self, cast

from numpy import abs as npabs
from numpy import (
    append,
    array,
    asarray,
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
from ._typing import FloatArray, SoilModelNames, SourceNames
from .soilmodel import Brooks, Genuchten, SoilModel, resolve_soilmodel

logger = getLogger(__name__)


@dataclass
class SoilSample:
    """A container for measured soil properties and in-situ soil hydraulic data.

    Class provides convenience routines to predict hydraulic parameters using
    a range of pedotransfer functions and empirical models.

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
        """Get properties and measurements from Staring series.

        References
        ----------
        Wösten, J. H. M. and Veerman, G. J. and de Groot, W. J. M. and Stolte, J. (2001).
        Waterretentie- en Doorlatendheidskarakteristieken van Boven- en Ondergronden in
        Nederland: De Staringreeks. url: https://edepot.wur.nl/43272

        Heinen, M. and Bakker, G. and Wösten, J. H. M. (2020). Waterretentie- en
        Doorlatendheidskarakteristieken van Boven- en Ondergronden in Nederland:
        De Staringreeks; (Update 2018). doi: 10.18174/512761

        """
        if year not in ("2001", "2018"):
            raise ValueError(
                f"No Staring series available for year '{year}'"
                "please use either '2001' or '2018'"
            )
        path = Path(__file__).parent / "datasets/soilsamples.csv"
        properties = read_csv(path, delimiter=";")
        staring_properties = properties.query(f"source == 'Staring_{year}'").set_index(
            "name"
        )

        self.silt_p = cast(float, staring_properties.at[name, "silt_p"])
        self.clay_p = cast(float, staring_properties.at[name, "clay_p"])
        self.om_p = cast(float, staring_properties.at[name, "om_p"])
        self.m50 = cast(float, staring_properties.at[name, "m50"])
        if year == "2001":
            self.rho = cast(float, staring_properties.at[name, "rho"])
        return self

    def fit(
        self,
        sm: type[SoilModel],
        pbounds: DataFrame | None = None,
        weights: FloatArray | float = 1.0,
        W1: float = 0.1,
        W2: float | None = None,
        k_s: float | None = None,
        silent: bool = True,
        **kwargs,
    ) -> SoilModel:
        """Fit the provided SoilModel to the measurements.

        Fit the provided SoilModel (e.g., van Genuchten, Brooks-Corey class) to the
        stored measurements (theta, k, h) using nonlinear least squares. If
        pbounds is not provided, default parameter bounds are retrieved for the
        requested model name. The objective combines water retention and log10(k)
        errors; weighting terms W1 and W2 control the relative contribution of k.
        Returns a model instance with optimized parameters.

        Parameters
        ----------
        sm : Type[SoilModel]
            The soil model class to fit (e.g., Genuchten, Brooks).
        pbounds : DataFrame | None, optional
            DataFrame with parameter bounds and initial values. If None, defaults are
            used based on the model name. Expected columns: 'p_ini', 'p_min', 'p_max'.
        weights : FloatArray | float, optional
            Weights for the objective function. Can be a single float (applied to all)
            or an array of length N+M (N for theta, M-N for k). Default is 1.0.
        W1 : float, optional
            Scaling factor for the k error in the objective function. Default is 0.1.
        W2 : float | None, optional
            Additional scaling factor for the k error. If None, it is computed to balance
            the contributions of theta and k errors based on their magnitudes and weights.
        k_s : float | None, optional
            If provided, this value of saturated hydraulic conductivity will be fixed during
            optimization. This means that the relative hydraulic conductivity curve will be
            estimated.
        silent : bool, optional
            If False, prints the optimization result. Default is True.
        kwargs : dict, optional
            Additional keyword arguments to pass to the scipy.optimize.least_squares function.

        Notes
        -----
            - Requires theta and k (and optionally h) to be set.
            - If k_s is provided it will be fixed during optimization.
            - Input/outputs and bounds are expected in the units used by the soil model.

        References
        ----------
        van Genuchten, M. Th. and Leij, F. J. and Yates, S. R. (1991). The RETC Code for
        Quantifying the Hydraulic Functions.
        url: https://cfpub.epa.gov/si/si_public_record_report.cfm?dirEntryId=130162

        """
        if self.theta is None or self.k is None or self.h is None:
            raise ValueError("theta, k, and h measurements are required for fitting")

        theta = asarray(self.theta, dtype=float)
        h = asarray(self.h, dtype=float)
        N = len(theta)
        k = asarray(self.k, dtype=float)
        M = N + len(k)

        if pbounds is None:
            pbounds = get_params(sm)
            if k_s is not None:
                pbounds = pbounds.drop("k_s")
            else:
                pbounds.loc["k_s", "p_ini"] = max(k)
                pbounds.loc["k_s", "p_max"] = max(k) * 10
            pbounds.loc["theta_s", "p_ini"] = max(theta)
            pbounds.loc["theta_s", "p_max"] = max(theta) + 0.02

        weights = (
            full(M, weights)
            if isinstance(weights, float)
            else asarray(weights, dtype=float)
        )

        if W2 is None:
            W2 = (
                (M - N)
                * sum(weights[0:N] * theta)
                / (N * sum(weights[N:M] * npabs(log10(k))))
            )
            logger.debug(f"Computed W2: {W2}")

        def get_diff(p: FloatArray) -> FloatArray:
            """Objective function for least squares optimization.

            Computes the difference between measured and model-predicted
            theta and log10(k) values, applying the specified weights and
            scaling factors.
            """
            p = asarray(p, dtype=float)
            est_pars = dict(zip(pbounds.index, p))
            if k_s is not None:
                est_pars["k_s"] = k_s
            sml = sm(**est_pars)
            theta_diff = sml.theta(h=h) - theta
            k_diff = log10(sml.k(h=h)) - log10(k)
            diff = append(weights[0:N] * theta_diff, weights[N:M] * W1 * W2 * k_diff)
            return diff

        kwargs_merged: dict = {"jac": "3-point", "x_scale": "jac"} | (kwargs or {})
        res = least_squares(
            get_diff,
            x0=pbounds.loc[:, "p_ini"],
            bounds=(
                pbounds.loc[:, "p_min"],
                pbounds.loc[:, "p_max"],
            ),
            **kwargs_merged,
        )
        opt_pars = dict(zip(pbounds.index, res.x))
        if k_s is not None:
            opt_pars["k_s"] = k_s

        if not silent:
            print("SciPy Optimization Result\n", res)

        return sm(**opt_pars)

    def wosten(self, topsoil: bool = False) -> Genuchten:
        """Pedotransfer function for general soils.

        Sometimes also referred to as HYPRES: HYdraulic PRoperties of European Soils

        Parameters
        ----------
        topsoil : bool, optional
            If True, applies the topsoil adjustment to the pedotransfer
            function. Default is False.

        References
        ----------
        Wosten et al (1999) - Development and use of a database of hydraulic
        properties of European soils. doi: 10.1016/S0016-7061(98)00132-3

        """
        assert (
            self.clay_p is not None
            and self.rho is not None
            and self.silt_p is not None
            and self.om_p is not None
        )
        ts = 1.0 * topsoil

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
            - 0.0001664 * ts * self.silt_p
        )
        alpha_ = (
            -14.96
            + 0.03135 * self.clay_p
            + 0.0351 * self.silt_p
            + 0.646 * self.om_p
            + 15.29 * self.rho
            - 0.192 * ts
            - 4.671 * self.rho**2
            - 0.000781 * self.clay_p**2
            - 0.00687 * self.om_p**2
            + 0.0449 * self.om_p**-1
            + 0.0663 * log(self.silt_p)
            + 0.1482 * log(self.om_p)
            - 0.04546 * self.rho * self.silt_p
            - 0.4852 * self.rho * self.om_p
            + 0.00673 * ts * self.clay_p
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
            + 0.00718 * ts * self.clay_p
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
            + 0.93 * ts
            - 0.967 * self.rho**2
            - 0.000484 * self.clay_p**2
            - 0.000322 * self.silt_p**2
            + 0.001 * self.silt_p**-1
            - 0.0748 * self.om_p**-1
            - 0.643 * log(self.silt_p)
            - 0.01398 * self.rho * self.clay_p
            - 0.1673 * self.rho * self.om_p
            + 0.02986 * ts * self.clay_p
            - 0.03305 * ts * self.silt_p
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

    def wosten_sand(self, topsoil: bool = False) -> Genuchten:
        """Pedotransfer function for sandy soils.

        Parameters
        ----------
        topsoil : bool, optional
            If True, applies the topsoil adjustment to the pedotransfer
            function. Default is False.

        Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
        van boven- en ondergronden in Nederland: de Staringreeks.
        url: https://edepot.wur.nl/43272

        """
        msg = "Wosten sand pedotransfer function requires 'rho', 'm50', 'silt_p' and 'om_p' to be set."
        assert self.rho is not None, msg
        assert self.m50 is not None, msg
        assert self.silt_p is not None, msg
        assert self.om_p is not None, msg

        ts = 1.0 * topsoil
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
            - 0.172 * ts
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

    def wosten_clay(self) -> Genuchten:
        """Pedotransfer function for clay soils.

        References
        ----------
        Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
        van boven- en ondergronden in Nederland: de Staringreeks.
        url: https://edepot.wur.nl/43272

        """
        msg = "Wosten clay pedotransfer function requires 'clay_p', 'rho', and 'om_p' to be set."
        assert self.clay_p is not None, msg
        assert self.rho is not None, msg
        assert self.om_p is not None, msg
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
        """Pedotransfer function returning Brooks-Corey parameters.

        Cooper (2021) - Using data assimilation to optimize pedotransfer
        functions using field-scale in situ soil moisture observations.
        doi: 10.5194/hess-25-2445-2021

        Cosby et al. (1984) - A statistical exploration of the relationships of soil moisture
        characteristics to the physical properties of soils. doi: 10.1029/WR020i006p00682

        """
        msg = "Cosby pedotransfer function requires 'clay_p' and 'sand_p' to be set."
        assert self.clay_p is not None, msg
        assert self.sand_p is not None, msg
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

    def saxton(self, df: float = 1.0) -> Brooks:
        """Pedotransfer function returning Brooks-Corey parameters.

        Implements the equations from Saxton and Rawls (2006) which estimate
        soil water characteristics from soil texture and organic matter.

        Parameters
        ----------
        df : float, optional
            Density adjustment factor (normally between 0.9 and 1.3) to account
            for compaction or loosening. Default is 1.0 (normal density).

        References
        ----------
        Saxton, K. E., Rawls, W. J., Romberger, J. S., & Papendick, R. I. (1986).
        Estimating generalized soil-water characteristics from texture.
        doi: 10.2136/sssaj1986.03615995005000040039x

        Saxton, K. E., & Rawls, W. J. (2006). Soil water characteristic estimates
        by texture and organic matter for hydrologic solutions.
        doi: 10.2136/sssaj2005.0117

        """
        msg = "Saxton-Rawls pedotransfer function requires 'sand_p', 'clay_p', and 'om_p' to be set."
        assert self.sand_p is not None, msg
        assert self.clay_p is not None, msg
        assert self.om_p is not None, msg

        # Sand and Clay are calculated as fractions, OM is kept as a percentage
        s_f = self.sand_p / 100.0
        c_f = self.clay_p / 100.0

        # 1500 kPa moisture (Wilting Point)
        th1500t = (
            -0.024 * s_f
            + 0.487 * c_f
            + 0.006 * self.om_p
            + 0.005 * (s_f * self.om_p)
            - 0.013 * (c_f * self.om_p)
            + 0.068 * (s_f * c_f)
            + 0.031
        )
        th1500 = th1500t + (0.14 * th1500t - 0.02)

        # 33 kPa moisture (Field capacity)
        th33t = (
            -0.251 * s_f
            + 0.195 * c_f
            + 0.011 * self.om_p
            + 0.006 * (s_f * self.om_p)
            - 0.027 * (c_f * self.om_p)
            + 0.452 * (s_f * c_f)
            + 0.299
        )
        th33 = th33t + (1.283 * th33t**2 - 0.374 * th33t - 0.015)

        # Saturation - 33 kPa moisture
        ths_33t = (
            0.278 * s_f
            + 0.034 * c_f
            + 0.022 * self.om_p
            - 0.018 * (s_f * self.om_p)
            - 0.027 * (c_f * self.om_p)
            - 0.584 * (s_f * c_f)
            + 0.078
        )
        ths_33 = ths_33t + (0.636 * ths_33t - 0.107)

        # Saturation moisture & density
        ths = th33 + ths_33 - 0.097 * s_f + 0.043
        rho_n = (1 - ths) * 2.65

        # Apply density factor adjustment
        if df != 1.0:
            rho_df = rho_n * df
            ths_df = 1 - (rho_df / 2.65)
            th33_df = th33 - 0.2 * (ths - ths_df)
            ths_33_df = ths_df - th33_df

            # Update variables with adjusted density
            ths = ths_df
            th33 = th33_df
            ths_33 = ths_33_df

        # Air entry tension (bubbling pressure) in kPa
        psi_et = (
            -21.67 * s_f
            - 27.93 * c_f
            - 81.97 * ths_33
            + 71.12 * (s_f * ths_33)
            + 8.29 * (c_f * ths_33)
            + 14.05 * (s_f * c_f)
            + 27.16
        )
        # Prevent mathematically negative air-entry tensions in edge-case soil bounds
        psi_e = max(psi_et + (0.02 * psi_et**2 - 0.113 * psi_et - 0.70), 0.1)
        # Convert air entry tension from kPa to cm water column (~10.1972)
        h_b = psi_e * 10.1972

        # Brooks-Corey parameters (Corrected inversion)
        lp = (log(th33) - log(th1500)) / (log(1500) - log(33))

        # Saturated hydraulic conductivity (cm/d = 2.4 mm/h)
        k_s = 1930 * (ths - th33) ** (3 - lp) * 2.4

        return Brooks(
            k_s=round(k_s, 4),
            theta_r=0.0,
            theta_s=round(ths, 4),
            h_b=round(h_b, 5),
            l=round(lp, 5),
        )

    def weynants(self) -> Genuchten:
        """Pedotransfer function returning Mualem-van Genuchten parameters.

        The Weynants PTF relies on Organic Carbon (OC), not Organic Matter (OM).
        This method converts the `om_p` to OC using a factor 1.724.

        References
        ----------
        Weynants, M., Vereecken, H., & Javaux, M. (2009). Revisiting Vereecken
        Pedotransfer Functions: Introducing a Closed-Form Hydraulic Model.
        doi: 10.2136/vzj2008.0062

        Weihermüller, L., Herbst, M., Javaux, M., & Weynants, M. (2017). Erratum to
        "Revisiting Vereecken Pedotransfer Functions: Introducing a Closed-Form
        Hydraulic Model". doi: 10.2136/vzj2008.0062er

        Vereecken, H., Maes, J., Feyen, J., & Darius, P. (1989). Estimating the
        soil moisture retention characteristic from texture, bulk density, and
        carbon content. doi: 10.1097/00010694-198912000-00001

        """
        msg = "Weynants pedotransfer function requires 'sand_p', 'clay_p', 'rho', and 'om_p' to be set."
        assert self.sand_p is not None, msg
        assert self.clay_p is not None, msg
        assert self.rho is not None, msg
        assert self.om_p is not None, msg

        # Convert organic matter percentage to organic carbon percentage
        oc_p = self.om_p / 1.724

        theta_s = 0.6355 + 0.0013 * self.clay_p - 0.1631 * self.rho
        alpha = exp(
            -4.3003 - 0.0097 * self.clay_p + 0.0138 * self.sand_p - 0.0992 * oc_p
        )
        n = (
            exp(
                -1.0846
                - 0.0236 * self.clay_p
                - 0.0085 * self.sand_p
                + 0.0001 * (self.sand_p**2)
            )
            + 1.0
        )

        k0 = exp(1.9582 + 0.0308 * self.sand_p - 0.6142 * self.rho - 0.1566 * oc_p)
        lt = -1.8642 - 0.1317 * self.clay_p + 0.0067 * self.sand_p

        return Genuchten(
            k_s=round(k0, 4),
            theta_r=0.0,
            theta_s=round(theta_s, 4),
            alpha=round(alpha, 7),
            n=round(n, 4),
            l=round(lt, 4),
        )

    def toth(self, topsoil: bool = False) -> Genuchten:
        """Pedotransfer function returning Mualem-van Genuchten parameters.

        Implements the continuous pedotransfer functions from Tóth et al. (2015)
        based on the EU-HYDI database. Uses Eq (21) for the Moisture Retention
        Characteristic and Eq (16) for Saturated Hydraulic Conductivity.

        Parameters
        ----------
        topsoil : bool, optional
            If True, applies the topsoil adjustment to the pedotransfer
            function. Default is False.

        Returns
        -------
        Genuchten
            Van Genuchten soil model with estimated parameters.

        References
        ----------
        Tóth, B., Weynants, M., Nemes, A., Makó, A., Bilas, G., & Tóth, G. (2015).
        New generation of hydraulic pedotransfer functions for Europe.
        doi: 10.1111/ejss.12192

        """
        msg = "Tóth 2015 PTF requires 'sand_p', 'silt_p', 'clay_p', 'rho', and 'om_p'."
        assert self.sand_p is not None, msg
        assert self.silt_p is not None, msg
        assert self.clay_p is not None, msg
        assert self.rho is not None, msg
        assert self.om_p is not None, msg

        oc_p = self.om_p / 1.724  # Converts OM to Organic Carbon
        ts = 1.0 * topsoil

        theta_r = 0.041 if self.sand_p >= 2.00 else 0.179
        theta_s = (
            0.83080
            - 0.28217 * self.rho
            + 0.0002728 * self.clay_p
            + 0.000187 * self.silt_p
        )
        alpha = 10 ** (
            -0.43348
            - 0.41729 * self.rho
            - 0.04762 * oc_p
            + 0.21810 * ts
            - 0.01581 * self.clay_p
            - 0.01207 * self.silt_p
        )
        n = (
            10
            ** (
                0.22236
                - 0.30189 * self.rho
                - 0.05558 * ts
                - 0.005306 * self.clay_p
                - 0.003084 * self.silt_p
                - 0.01072 * oc_p
            )
        ) + 1.0

        if oc_p < 0.07:
            log10_ks = 0.55
        elif 0.07 <= oc_p < 0.40:
            if self.sand_p < 5.77:
                log10_ks = -0.11
            elif 5.77 <= self.sand_p < 69.72:
                log10_ks = 1.28
            else:  # self.sand_p >= 69.72
                log10_ks = 1.96
        elif 0.40 <= oc_p < 0.41:
            if self.silt_p >= 32.11:
                log10_ks = -1.81
            else:
                log10_ks = -0.40
        elif 0.41 <= oc_p < 0.96:
            if self.clay_p >= 37.4:
                log10_ks = 0.67
            else:
                log10_ks = 1.53
        else:  # oc_p >= 0.96
            if topsoil:  # Subsoil logic
                if 0.96 <= oc_p < 2.09:
                    if self.silt_p < 10.85:
                        log10_ks = 0.01
                    else:  # self.silt_p >= 10.85
                        if 0.96 <= oc_p < 1.52:
                            log10_ks = 1.82
                        elif 1.52 <= oc_p < 1.54:
                            log10_ks = -0.46
                        else:  # 1.54 <= oc_p < 2.09
                            log10_ks = 1.72
                elif 2.09 <= oc_p < 2.10:
                    log10_ks = -0.87
                else:  # oc_p >= 2.10
                    if self.sand_p >= 38.95:
                        log10_ks = 1.44
                    else:  # self.sand_p < 38.95
                        if 2.10 <= oc_p < 2.42:
                            log10_ks = 1.82
                        else:  # oc_p >= 2.42
                            log10_ks = -0.22
            else:
                if 0.96 <= oc_p < 0.97:
                    log10_ks = -0.95
                elif 0.97 <= oc_p < 1.52:
                    log10_ks = 1.13
                elif 1.52 <= oc_p < 1.54:
                    log10_ks = -0.75
                elif 1.54 <= oc_p < 2.04:
                    log10_ks = 1.33
                else:  # oc_p >= 2.04
                    log10_ks = -0.25

        k_s = 10**log10_ks

        return Genuchten(
            k_s=round(k_s, 4),
            theta_r=round(theta_r, 4),
            theta_s=round(theta_s, 4),
            alpha=round(alpha, 7),
            n=round(n, 4),
            l=0.5,
        )

    def rosetta(self, version: Literal[1, 2, 3] = 3) -> Genuchten:
        """Pedotransfer function using the Rosetta API.

        References
        ----------
        Schaap et al., (2001) - Predicting soil water retention from soil.
        doi: 10.1016/S0022-1694(01)00466-8

        Zhang Y. and Schaap, M. G. (2017) - Weighted recalibration of the
        Rosetta pedotransfer model with improved estimates of hydraulic
        parameter distributions and summary statistics (Rosetta3).
        doi: 10.1016/j.jhydrol.2017.01.004

        """
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
        """Estimate van Genuchten parameters using the HYPAGS method.

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
        doi: 10.1111/gwat.13365

        Peche, A., & Houben, G. J. (2023). Estimating characteristic grain sizes and
        effective porosity from hydraulic conductivity data. Groundwater, 61(4), 574-585.
        doi: 10.1111/gwat.13266

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
            """Solve for d50 and effective porosity.

            Based on Kozeny-Carman equation using scipy's fixed_point
            for robust successive substitution iteration.

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
                """Fixed-point iteration function.

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
            """Return bootstrap confidence interval errors for a given k value.

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
                (3.002, 1.144),
                (5.364, 1.258),
                (7.031, 1.366),
                (6.728, 1.519),
                (6.758, 2.072),
                (6.607, 2.393),
                (5.91, 2.985),
                (5.643, 3.229),
                (4.84, 4.0),
                (4.082, 4.233),
                (3.182, 4.486),
                (3.051, 5.03),
                (3.442, 5.483),
            ]
            default_error = (1.811, 2.858)

            i = bisect_right(bounds, k)
            sel_error = errors[i - 1] if i <= len(errors) else default_error

            return sel_error

        def get_n_errors(alpha: float) -> tuple[float, float]:
            """Return bootstrap confidence interval errors for a given alpha value.

            Uses binary search to find the corresponding [+Δn, -Δn] pair
            based on predefined alpha thresholds.

            Parameters
            ----------
            alpha : float
                Van Genuchten alpha parameter.

            Returns
            -------
            tuple[float, float]
            (added_error, subtracted_error)

            """
            bounds = [1, 2, 3, 4, 5, 6, 7, 8, 9]
            errors = [
                (2.836, 0.358),
                (3.061, 0.826),
                (3.481, 0.828),
                (3.182, 0.605),
                (1.905, 0.509),
                (1.184, 0.297),
                (4.032, 0.352),
                (1.385, 0.241),
                (1.017, 0.155),
            ]
            default_error = (3.006, 0.377)

            i = bisect_right(bounds, alpha)
            sel_error = errors[i - 1] if i <= len(errors) else default_error
            return sel_error

        def get_res_water_content(k: float) -> tuple[float, float, float]:
            """Return residual water content (θr) with lower and upper bounds for a given k value.

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
                (0.089, 0.097, 0.198),
                (0.097, 0.107, 0.207),
                (0.09, 0.099, 0.227),
                (0.095, 0.102, 0.221),
                (0.091, 0.098, 0.216),
                (0.087, 0.097, 0.227),
                (0.104, 0.111, 0.212),
                (0.099, 0.108, 0.219),
                (0.079, 0.085, 0.233),
                (0.082, 0.087, 0.226),
                (0.067, 0.07, 0.223),
                (0.105, 0.11, 0.23),
                (0.096, 0.105, 0.23),
            ]
            default_water_content = (0.068, 0.111, 0.096)

            i = bisect_right(bounds, k)
            res = values[i - 1] if i <= len(values) else default_water_content
            return res

        # mathematical model of hypags calculation of k, d10, d20:
        # CONDITION: in HYPAGS, either k, d10 or d20 have to be given
        if self.k is not None:
            if isinstance(self.k, float):
                k = self.k
            elif isinstance(self.k, ndarray):
                if len(self.k) > 1:
                    logger.warning(
                        "HYPAGS routine only accepts single k value, choosing the first k value in the array."
                    )
                k = float(self.k.item(0))
            else:
                k = float(self.k)

            # case 0: mathematical model where k is given
            # check for non-valid input
            if k > 2.6e-2 or k < 2.87e-7:
                logger.error("k out of hypags model limits.")
            logger.debug("Using case 0 of hypags model (k given).")
            d10 = (k / Pi * c) ** (0.5)  # calculation of d10
            self.d20 = c1 * d10  # calculation of d20
            self.d10 = d10
        elif self.d10 is not None:
            # case 1: mathematical model where d10 is given
            if not 5.35e-5 <= self.d10 <= 8.3e-4:
                logger.error(
                    f"d10 ({self.d10:.3e}) out of hypags model limits: 5.35e-5 to 8.3e-4."
                )
            logger.debug("Using case 1 of hypags model (d10 given).")
            self.d20 = c1 * self.d10  # calculation of d20
            k = (Pi * rho_f * g * self.d10**2) / mu
            self.k = array([k], dtype=float)  # calculation of k
        elif self.d20 is not None:
            # case 2: mathematical model where d20 is given
            if not 6.25e-5 <= self.d20 <= 1.2e-3:
                logger.error(
                    f"d20 ({self.d20:.3e}) out of hypags model limits: 6.25e-5 to 1.2e-3."
                )
            logger.debug("Using case 2 of hypags model (d20 given).")
            self.d10 = self.d20 / c1  # calculation of d10
            k = (Pi * rho_f * g * self.d10**2) / mu
            self.k = array([k], dtype=float)  # calculation of k
        else:
            raise ValueError(
                "No parameter (k, d10, or d20) was provided for hypags routine."
            )

        # calculation of d50 and ne
        assert self.d20 is not None, (
            "d20 should have been calculated in the previous steps."
        )
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
    """A class representing soil properties and models.

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
    source: SourceNames | None = None
    description: str | None = None

    def from_name(
        self,
        sm: type[SoilModel] | SoilModel | SoilModelNames,
        source: SourceNames | None = None,
    ) -> Self:
        """Load soil parameters from a CSV database by soil name and model type.

        Available sources include HYDRUS, VS2D, Staring_2001, Staring_2018.

        References
        ----------
        Clapp, R. B., & Hornberger, G. M. (1978). Empirical equations for some soil
        hydraulic properties. doi: 10.1029/WR014i004p00601

        Carsel, R. F. and Parrish, R. S. (1988). Developing Joint Probability
        Distributions of Soil Water Retention Characteristics.
        doi: 10.1029/WR024i005p00755

        Simunek, J. and Sejna, M. and Saito, H. and Sakai, M. and van Genuchten, M. Th. (2009).
        The HYDRUS-1D software package for simulating the one-dimensional movement of water,
        heat, and multiple solutes in variably-saturated media. Version 4.08.
        url: https://www.pc-progress.com/downloads/pgm_hydrus1d/hydrus1d-4.08.pdf}

        Healy, R. W. (1990). Simulation of Solute Transport in Variably Saturated Porous Media
        with Supplemental Information on Modifications to the U.S. Geological Survey Computer
        Program VS2D. doi: 10.3133/wri904025

        Wösten, J. H. M. and Veerman, G. J. and de Groot, W. J. M. and Stolte, J. (2001).
        Waterretentie- en Doorlatendheidskarakteristieken van Boven- en Ondergronden in
        Nederland: De Staringreeks. url: https://edepot.wur.nl/43272

        Heinen, M. and Bakker, G. and Wösten, J. H. M. (2020). Waterretentie- en
        Doorlatendheidskarakteristieken van Boven- en Ondergronden in Nederland:
        De Staringreeks; (Update 2018). doi: 10.18174/512761

        """
        smn, sm_cls = resolve_soilmodel(sm)

        path = Path(__file__).parent / "datasets/soilsamples.csv"
        ser = read_csv(path, delimiter=";", index_col=0)
        if "HYDRUS_" in self.name:
            self.name = self.name.replace("HYDRUS_", "")
            source = "HYDRUS"
            logger.warning(
                "Removed 'HYDRUS_' from soil name. For future use"
                "please provide source='HYDRUS' argument"
            )
        sersm = ser.query(f"soilmodel == '{smn}'").loc[[self.name], :]
        if source is None and len(sersm) > 1:
            raise ValueError(
                f"Multiple sources for soil {self.name}: "
                f"{sersm.loc[:, 'source'].to_numpy(dtype=str)}. "
                f"Please provide the source using the source argument"
            )
        elif (source is not None) and len(sersm) > 1:
            sersm = sersm.query(f"source == '{source}'")

        sers: Any = sersm.iloc[0].copy()
        if isna(sers.at["description"]):
            sers.loc["description"] = sers.at["soil type"]
        setattr(self, "source", sers.pop("source"))
        setattr(self, "description", sers.pop("description"))
        model_fields = [f.name for f in fields(cast(Any, sm_cls)) if f.init]
        smserd = {x: sers.at[x] for x in model_fields}
        setattr(self, "model", cast(Any, sm_cls)(**smserd))
        return self

    @staticmethod
    def list_names(sm: type[SoilModel] | SoilModel | SoilModelNames) -> list[str]:
        """Return a list of available soil names for a given soil model."""
        smn, _ = resolve_soilmodel(sm)

        path = Path(__file__).parent / "datasets/soilsamples.csv"
        names = read_csv(path, delimiter=";")

        return names.query(f"soilmodel == '{smn}'").loc[:, "name"].unique().tolist()

    def from_staring(self, year: Literal["2001", "2018"] = "2018") -> Self:
        """Load soil parameters from the Staring series database."""
        if year not in ("2001", 2001, "2018", 2018):
            raise ValueError(f"Year must either be '2001' or '2018', not {year}")

        year_str = str(year)
        source: SourceNames = "Staring_2001" if year_str == "2001" else "Staring_2018"
        self.from_name(sm=Genuchten, source=source)
        ss = SoilSample().from_staring(name=self.name, year=year_str)
        setattr(self, "sample", ss)
        return self
