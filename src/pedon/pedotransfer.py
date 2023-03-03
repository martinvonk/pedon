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


def wosten(sm: SoilSample, ts: bool = False):
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
    return theta_r, theta_s, exp(ks_), exp(alpha_), exp(n_), exp(l_)


def cosby(
    sm: SoilSample,
):
    k01 = 3.100
    k02 = 15.70
    k03 = 0.300
    k04 = 0.505
    k05 = 0.037
    k06 = 0.142
    k07 = 2.170
    k08 = 0.630
    k09 = 1.580
    k10 = 0.600
    k11 = 0.640
    k12 = 1.260
    c = sm.clay_p / 100
    s = sm.sand_p / 100
    b = k01 + k02 * c - k03 * s
    theta_s = k04 - k05 * c - k06 * s
    psi_s = 0.01 * 10 ** (k07 - k08 * c - k09 * s)
    ks = 10 ** (-k10 - k11 * c + k12 * s) * 25.2 / 3600
    theta_r = 0.0
    labda = 1 / b
    ks = ks * 8640000 / 1000  # kg/m2/s to cm/d
    psi_s = psi_s * 100  # m to cm
    return theta_r, theta_s, ks, psi_s, labda
