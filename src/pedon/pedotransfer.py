from dataclasses import dataclass

from numpy import exp, log

from .soilmodel import Brooks, Genuchten


@dataclass
class SoilSample:
    sand_p: float  # sand %
    silt_p: float  # silt %
    clay_p: float  # clay %
    rho: float  # soil density g/cm3
    th33: float  # cm
    th1500: float  # cm
    om_p: float  # organic matter %
    m50: float  # median sand fraction


def wosten(sm: SoilSample, ts: bool = False) -> Genuchten:
    """Wosten et al (1999) - Development and use of a database of hydraulic
    properties of European soils"""

    topsoil = 1.0 * ts

    theta_s = (
        0.7919
        + 0.001691 * sm.clay_p_p
        - 0.29619 * sm.rho
        - 0.000001419 * sm.silt_p_p**2
        + 0.0000821 * sm.om_p**2
        + 0.02427 * sm.clay_p**-1
        + 0.01113 * sm.silt_p**-1
        + 0.01472 * log(sm.silt_p)
        - 0.0000733 * sm.om_p * sm.clay_p
        - 0.000619 * sm.rho * sm.clay_p
        - 0.001183 * sm.rho * sm.om_p
        - 0.0001664 * topsoil * sm.silt_p
    )
    alpha_ = (
        -14.96
        + 0.03135 * sm.clay_p
        + 0.0351 * sm.silt_p
        + 0.646 * sm.om_p
        + 15.29 * sm.rho
        - 0.192 * topsoil
        - 4.671 * sm.rho**2
        - 0.000781 * sm.clay_p**2
        - 0.00687 * sm.om_p**2
        + 0.0449 * sm.om_p**-1
        + 0.0663 * log(sm.silt_p)
        + 0.1482 * log(sm.om_p)
        - 0.04546 * sm.rho * sm.silt_p
        - 0.4852 * sm.rho * sm.om_p
        + 0.00673 * topsoil * sm.clay_p
    )
    n_ = (
        -25.23
        - 0.02195 * sm.clay_p
        + 0.0074 * sm.silt_p
        - 0.1940 * sm.om_p
        + 45.5 * sm.rho
        - 7.24 * sm.rho**2
        - 0.0003658 * sm.clay_p**2
        + 0.002885 * sm.om_p**2
        - 12.81 * sm.rho**-1
        - 0.1524 * sm.silt_p**-1
        - 0.01958 * sm.om_p**-1
        - 0.2876 * log(sm.silt_p)
        - 0.0709 * log(sm.om_p)
        - 44.6 * log(sm.rho)
        - 0.02264 * sm.rho * sm.clay_p
        + 0.0896 * sm.rho * sm.om_p
        + 0.00718 * topsoil * sm.clay_p
    )
    l_ = (
        0.0202
        + 0.0006193 * sm.clay_p**2
        - 0.001136 * sm.om_p**2
        - 0.2316 * log(sm.om_p)
        - 0.03544 * sm.rho * sm.clay_p
        + 0.00283 * sm.rho * sm.silt_p
        + 0.0488 * sm.rho * sm.om_p
    )
    ks_ = (
        7.755
        + 0.0352 * sm.silt_p
        + 0.93 * topsoil
        - 0.967 * sm.rho**2
        - 0.000484 * sm.clay_p**2
        - 0.000322 * sm.silt_p**2
        + 0.001 * sm.silt_p**-1
        - 0.0748 * sm.om_p**-1
        - 0.643 * log(sm.silt_p)
        - 0.01398 * sm.rho * sm.clay_p
        - 0.1673 * sm.rho * sm.om_p
        + 0.02986 * topsoil * sm.clay_p
        - 0.03305 * topsoil * sm.silt_p
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


def wosten_sand(sm: SoilSample, ts: bool = False) -> Genuchten:
    """Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
    van boven- en ondergronden in Nederland: de Staringreeks"""

    topsoil = 1.0 * ts
    theta_s = (
        -35.7
        - 0.1843 * sm.rho
        - 0.03576 * sm.m50
        + 0.0000261 * sm.m50**2
        - 0.0564 * sm.silt_p**-1
        + 0.008 * sm.om_p**-1
        + 496.0 * sm.m50**-1
        + 0.02244 * log(sm.om_p)
        + 7.56 * log(sm.m50)
    )

    k_s = (
        45.8
        - 14.34 * sm.rho
        + 0.001481 * sm.silt_p**2
        - 27.5 * sm.rho**-1
        - 0.891 * log(sm.silt_p)
        - 0.34 * log(sm.om_p)
    )

    alpha = (
        13.66
        - 5.91 * sm.rho
        - 0.172 * topsoil
        + 0.003248 * sm.m50
        - 11.89 * sm.rho**-1
        - 2.121 * sm.silt_p**-1
        - 0.3742 * log(sm.silt_p)
    )

    n = (
        -1.057
        + 0.1003 * sm.om_p
        + 1.119 * sm.rho
        + 0.000764 * sm.silt_p**2
        - 0.1397 * sm.om_p**-1
        - 57.2 * sm.m50**-1
        - 0.557 * log(sm.om_p)
        - 0.02997 * sm.rho * sm.silt_p
    )

    l = (
        -76.4
        - 0.097 * sm.silt_p
        + 59.6 * sm.rho
        + 0.0332 * sm.m50
        - 13.45 * sm.rho**2
        + 0.001127 * sm.silt_p**2
        + 0.00431 * sm.om_p**2
        - 0.0000399 * sm.m50**2
        + 40.8 * sm.rho**-1
        + 2.364 * sm.silt_p**-1
        + 1.014 * log(sm.silt_p)
    )

    return Genuchten(k_s=k_s, theta_r=0.01, theta_s=theta_s, alpha=alpha, n=n, l=l)


def wosten_clay(sm: SoilSample) -> Brooks:
    """Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
    van boven- en ondergronden in Nederland: de Staringreeks"""

    theta_s = (
        0.6311
        + 0.003383 * sm.clay_p
        - 0.09699 * sm.rho**2
        - 0.00204 * sm.rho * sm.clay_p
    )

    k_s = (
        -42.6
        + 8.71 * sm.om_p
        + 61.9 * sm.rho
        - 20.79 * sm.rho**2
        - 0.2107 * sm.om_p**2
        - 0.01622 * sm.clay_p * sm.om_p
        - 5.382 * sm.rho * sm.om_p
    )

    alpha = (
        -19.13
        + 0.812 * sm.om_p
        + 23.4 * sm.rho
        - 8.16 * sm.rho**2
        + 0.423 * sm.om_p**-1
        + 2.388 * log(sm.om_p)
        - 1.338 * sm.rho * sm.om_p
    )

    n = (
        -0.235
        + 0.972 * sm.rho**-1
        - 0.7743 * log(sm.clay_p)
        - 0.3154 * log(sm.om_p)
        + 0.0678 * sm.rho_p * sm.om_p
    )

    l = 0.102 + 0.0222 * sm.clay_p - 0.043 * sm.rho * sm.clay_p

    return Genuchten(k_s=k_s, theta_r=0.01, theta_s=theta_s, alpha=alpha, n=n, l=l)


def cosby(sm: SoilSample) -> Brooks:
    """Cooper (2021) - Using data assimilation to optimize pedotransfer
    functions using field-scale in situ soil moisture observations"""
    c = sm.clay_p / 100
    s = sm.sand_p / 100
    b = 3.100 + 15.70 * c - 0.300 * s
    theta_s = 0.505 - 0.037 * c - 0.142 * s
    psi_s = 0.01 * 10 ** (2.170 - 0.630 * c - 1.580 * s)
    k_s = 10 ** (-0.600 - 0.640 * c + 1.260 * s) * 25.2 / 3600
    theta_r = 0.0
    labda = 1 / b
    k_s = k_s * 8640000 / 1000  # kg/m2/s to cm/d
    psi_s = psi_s * 100  # m to cm
    return Brooks(k_s=k_s, theta_r=theta_r, theta_s=theta_s, h_b=psi_s, l=labda)
