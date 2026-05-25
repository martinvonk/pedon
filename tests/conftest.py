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
    return pe.Brooks(k_s=10, theta_r=0.01, theta_s=0.43, h_b=10, l=2)


@fixture
def panday() -> pe.SoilModel:
    """Fixture for a Panday soil model with specific parameters for testing."""
    return pe.Panday(k_s=10, theta_r=0.01, theta_s=0.43, alpha=0.02, beta=1.1, brook=3)


@fixture
def gardner() -> pe.SoilModel:
    """Fixture for a Gardner soil model with specific parameters for testing."""
    return pe.Gardner(k_s=10, theta_s=0.43, c=0.02, m=1.1)


@fixture
def rucker() -> pe.SoilModel:
    """Fixture for a Rucker soil model with specific parameters for testing."""
    return pe.Rucker(k_s=10, theta_r=0.01, theta_s=0.43, c=0.02, m=1.1)


@fixture
def genuchtengardner() -> pe.SoilModel:
    """Fixture for a Genuchten-Gardner soil model with specific parameters for testing."""
    return pe.GenuchtenGardner(
        k_s=10, theta_s=0.43, theta_r=0.01, c=0.02, alpha=0.04, n=1.4
    )


@fixture
def haverkamp() -> pe.soilmodel.Haverkamp:
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
def genuchtenkool() -> pe.SoilModel:
    """Fixture for a Genuchten-Kool soil model with specific parameters for testing."""
    return pe.GenuchtenKool(
        k_s=10.0,
        theta_r=0.01,
        theta_s=0.43,
        alpha=0.02,
        n=1.1,
        l=0.5,
        xi=2.5,
    )