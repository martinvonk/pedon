from dataclasses import dataclass

from numpy import exp, log
from scipy.optimize import least_squares

from ._typing import FloatArray
from .soilmodel import Brooks, Genuchten, SoilModel


@dataclass
class SoilSample:
    sand_p: float | None = None  # sand %
    silt_p: float | None = None  # silt %
    clay_p: float | None = None  # clay %
    rho: float | None = None  # soil density g/cm3
    th33: float | None = None  # cm
    th1500: float | None = None  # cm
    om_p: float | None = None  # organic matter %
    m50: float | None = None  # median sand fraction
    h: FloatArray | None = None  # pressure head measurement
    k: FloatArray | None = None  # hydraulic conductivity measurement
    theta: FloatArray | None = None  # moisture content measurement


def fit(ss: SoilSample, sm: SoilModel) -> SoilModel:
    h = ss.h
    k = ss.k
    theta = ss.k


def wosten(ss: SoilSample, ts: bool = False) -> Genuchten:
    """Wosten et al (1999) - Development and use of a database of hydraulic
    properties of European soils"""

    topsoil = 1.0 * ts

    theta_s = (
        0.7919
        + 0.001691 * ss.clay_p_p
        - 0.29619 * ss.rho
        - 0.000001419 * ss.silt_p_p**2
        + 0.0000821 * ss.om_p**2
        + 0.02427 * ss.clay_p**-1
        + 0.01113 * ss.silt_p**-1
        + 0.01472 * log(ss.silt_p)
        - 0.0000733 * ss.om_p * ss.clay_p
        - 0.000619 * ss.rho * ss.clay_p
        - 0.001183 * ss.rho * ss.om_p
        - 0.0001664 * topsoil * ss.silt_p
    )
    alpha_ = (
        -14.96
        + 0.03135 * ss.clay_p
        + 0.0351 * ss.silt_p
        + 0.646 * ss.om_p
        + 15.29 * ss.rho
        - 0.192 * topsoil
        - 4.671 * ss.rho**2
        - 0.000781 * ss.clay_p**2
        - 0.00687 * ss.om_p**2
        + 0.0449 * ss.om_p**-1
        + 0.0663 * log(ss.silt_p)
        + 0.1482 * log(ss.om_p)
        - 0.04546 * ss.rho * ss.silt_p
        - 0.4852 * ss.rho * ss.om_p
        + 0.00673 * topsoil * ss.clay_p
    )
    n_ = (
        -25.23
        - 0.02195 * ss.clay_p
        + 0.0074 * ss.silt_p
        - 0.1940 * ss.om_p
        + 45.5 * ss.rho
        - 7.24 * ss.rho**2
        - 0.0003658 * ss.clay_p**2
        + 0.002885 * ss.om_p**2
        - 12.81 * ss.rho**-1
        - 0.1524 * ss.silt_p**-1
        - 0.01958 * ss.om_p**-1
        - 0.2876 * log(ss.silt_p)
        - 0.0709 * log(ss.om_p)
        - 44.6 * log(ss.rho)
        - 0.02264 * ss.rho * ss.clay_p
        + 0.0896 * ss.rho * ss.om_p
        + 0.00718 * topsoil * ss.clay_p
    )
    l_ = (
        0.0202
        + 0.0006193 * ss.clay_p**2
        - 0.001136 * ss.om_p**2
        - 0.2316 * log(ss.om_p)
        - 0.03544 * ss.rho * ss.clay_p
        + 0.00283 * ss.rho * ss.silt_p
        + 0.0488 * ss.rho * ss.om_p
    )
    ks_ = (
        7.755
        + 0.0352 * ss.silt_p
        + 0.93 * topsoil
        - 0.967 * ss.rho**2
        - 0.000484 * ss.clay_p**2
        - 0.000322 * ss.silt_p**2
        + 0.001 * ss.silt_p**-1
        - 0.0748 * ss.om_p**-1
        - 0.643 * log(ss.silt_p)
        - 0.01398 * ss.rho * ss.clay_p
        - 0.1673 * ss.rho * ss.om_p
        + 0.02986 * topsoil * ss.clay_p
        - 0.03305 * topsoil * ss.silt_p
    )
    theta_r = 0.01
    return Genuchten(
        k_s=exp(ks_),
        theta_r=theta_r,
        theta_s=theta_s,
        alpha=exp(alpha_),
        n=exp(n_),
        l=exp(l_),
    )


def wosten_sand(ss: SoilSample, ts: bool = False) -> Genuchten:
    """Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
    van boven- en ondergronden in Nederland: de Staringreeks"""

    topsoil = 1.0 * ts
    theta_s = (
        -35.7
        - 0.1843 * ss.rho
        - 0.03576 * ss.m50
        + 0.0000261 * ss.m50**2
        - 0.0564 * ss.silt_p**-1
        + 0.008 * ss.om_p**-1
        + 496.0 * ss.m50**-1
        + 0.02244 * log(ss.om_p)
        + 7.56 * log(ss.m50)
    )

    k_s = (
        45.8
        - 14.34 * ss.rho
        + 0.001481 * ss.silt_p**2
        - 27.5 * ss.rho**-1
        - 0.891 * log(ss.silt_p)
        - 0.34 * log(ss.om_p)
    )

    alpha = (
        13.66
        - 5.91 * ss.rho
        - 0.172 * topsoil
        + 0.003248 * ss.m50
        - 11.89 * ss.rho**-1
        - 2.121 * ss.silt_p**-1
        - 0.3742 * log(ss.silt_p)
    )

    n = (
        -1.057
        + 0.1003 * ss.om_p
        + 1.119 * ss.rho
        + 0.000764 * ss.silt_p**2
        - 0.1397 * ss.om_p**-1
        - 57.2 * ss.m50**-1
        - 0.557 * log(ss.om_p)
        - 0.02997 * ss.rho * ss.silt_p
    )

    l = (
        -76.4
        - 0.097 * ss.silt_p
        + 59.6 * ss.rho
        + 0.0332 * ss.m50
        - 13.45 * ss.rho**2
        + 0.001127 * ss.silt_p**2
        + 0.00431 * ss.om_p**2
        - 0.0000399 * ss.m50**2
        + 40.8 * ss.rho**-1
        + 2.364 * ss.silt_p**-1
        + 1.014 * log(ss.silt_p)
    )

    return Genuchten(k_s=k_s, theta_r=0.01, theta_s=theta_s, alpha=alpha, n=n, l=l)


def wosten_clay(ss: SoilSample) -> Brooks:
    """Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
    van boven- en ondergronden in Nederland: de Staringreeks"""

    theta_s = (
        0.6311
        + 0.003383 * ss.clay_p
        - 0.09699 * ss.rho**2
        - 0.00204 * ss.rho * ss.clay_p
    )

    k_s = (
        -42.6
        + 8.71 * ss.om_p
        + 61.9 * ss.rho
        - 20.79 * ss.rho**2
        - 0.2107 * ss.om_p**2
        - 0.01622 * ss.clay_p * ss.om_p
        - 5.382 * ss.rho * ss.om_p
    )

    alpha = (
        -19.13
        + 0.812 * ss.om_p
        + 23.4 * ss.rho
        - 8.16 * ss.rho**2
        + 0.423 * ss.om_p**-1
        + 2.388 * log(ss.om_p)
        - 1.338 * ss.rho * ss.om_p
    )

    n = (
        -0.235
        + 0.972 * ss.rho**-1
        - 0.7743 * log(ss.clay_p)
        - 0.3154 * log(ss.om_p)
        + 0.0678 * ss.rho_p * ss.om_p
    )

    l = 0.102 + 0.0222 * ss.clay_p - 0.043 * ss.rho * ss.clay_p

    return Genuchten(k_s=k_s, theta_r=0.01, theta_s=theta_s, alpha=alpha, n=n, l=l)


def cosby(ss: SoilSample) -> Brooks:
    """Cooper (2021) - Using data assimilation to optimize pedotransfer
    functions using field-scale in situ soil moisture observations"""
    c = ss.clay_p / 100
    s = ss.sand_p / 100
    b = 3.100 + 15.70 * c - 0.300 * s
    theta_s = 0.505 - 0.037 * c - 0.142 * s
    psi_s = 0.01 * 10 ** (2.170 - 0.630 * c - 1.580 * s)
    k_s = 10 ** (-0.600 - 0.640 * c + 1.260 * s) * 25.2 / 3600
    theta_r = 0.0
    labda = 1 / b
    k_s = k_s * 8640000 / 1000  # kg/m2/s to cm/d
    psi_s = psi_s * 100  # m to cm
    return Brooks(k_s=k_s, theta_r=theta_r, theta_s=theta_s, h_b=psi_s, l=labda)
