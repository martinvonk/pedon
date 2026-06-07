"""Fixtures for testing."""

from pytest import fixture

import pedon as pe


@fixture
def genuchten() -> pe.SoilModel:
    """Fixture for a van Genuchten soil model with specific parameters for testing."""
    return pe.Genuchten(k_s=10.0, theta_r=0.01, theta_s=0.43, alpha=0.02, n=1.1, l=0.5)


@fixture
def brooks() -> pe.SoilModel:
    """Fixture for a Brooks-Corey soil model with specific parameters for testing."""
    return pe.Brooks(k_s=10.0, theta_r=0.01, theta_s=0.43, h_b=10, l=2)


@fixture
def panday() -> pe.SoilModel:
    """Fixture for a Panday soil model with specific parameters for testing."""
    return pe.Panday(
        k_s=10.0, theta_r=0.01, theta_s=0.43, alpha=0.02, beta=1.1, brook=3
    )


@fixture
def gardner() -> pe.SoilModel:
    """Fixture for a Gardner soil model with specific parameters for testing."""
    return pe.Gardner(k_s=10.0, theta_s=0.43, c=0.02, m=1.1)


@fixture
def rucker() -> pe.SoilModel:
    """Fixture for a Rucker soil model with specific parameters for testing."""
    return pe.Rucker(k_s=10.0, theta_r=0.01, theta_s=0.43, c=0.02, m=1.1)


@fixture
def genuchtengardner() -> pe.SoilModel:
    """Fixture for a Genuchten-Gardner soil model with specific parameters for testing."""
    return pe.GenuchtenGardner(
        k_s=10.0, theta_s=0.43, theta_r=0.01, c=0.02, alpha=0.04, n=1.4
    )


@fixture
def fredlund() -> pe.SoilModel:
    """Fixture for a Fredlund soil model with specific parameters for testing."""
    return pe.Fredlund(k_s=10.0, theta_s=0.43, a=0.5, n=1.4, m=1.1)


@fixture
def kosugi() -> pe.Kosugi:
    """Fixture for a Kosugi soil model with specific parameters for testing."""
    return pe.Kosugi(k_s=10.0, theta_r=0.01, theta_s=0.43, h_m=100.0, sigma=1.2)


@fixture
def campbell() -> pe.Campbell:
    """Fixture for a Campbell soil model with specific parameters for testing."""
    return pe.Campbell(k_s=10.0, theta_s=0.43, h_b=10.0, b=4.0)


@fixture
def haverkamp() -> pe.Haverkamp:
    """Fixture for a Haverkamp soil model with specific parameters for testing."""
    # Example parameters similar in style to other model tests
    return pe.Haverkamp(
        k_s=10.0,
        theta_r=0.01,
        theta_s=0.43,
        alpha=0.02,
        beta=1.2,
        a=0.5,
    )


@fixture
def kool() -> pe.Kool:
    """Fixture for a Kool soil model with specific parameters for testing."""
    return pe.Kool(
        k_s=10.0,
        theta_r=0.01,
        theta_s=0.43,
        alpha=0.02,
        n=1.1,
        l=0.5,
        xi=2.5,
    )


@fixture
def brunswick() -> pe.Brunswick:
    """Fixture for a Brunswick soil model with specific parameters for testing."""
    return pe.Brunswick(
        theta_snc=0.08,
        theta_sc=0.42,
        alpha=0.05,
        n=1.6,
        k_sc=1e2,
        k_snc=1e-2,
        l=0.5,
    )
