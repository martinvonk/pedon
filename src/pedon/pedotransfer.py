from dataclasses import dataclass
from numpy import log, exp

from .soilmodel import Genuchten, Brooks


@dataclass
class SoilSample:
    sand_p: float  # sand %
    silt_p: float  # silt %
    clay_p: float  # clay %
    rho: float  # soil density g/cm3
    th33: float  # cm
    th1500: float  # cm
    om_p: float  # organic matter %


def wosten(sm: SoilSample, ts: bool = False) -> Genuchten:
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


def cosby(sm: SoilSample) -> Brooks:
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
