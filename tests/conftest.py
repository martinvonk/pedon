"""Fixtures for testing."""

from pytest import fixture

from pedon import Genuchten


@fixture
def gen() -> Genuchten:
    """Fixture for a van Genuchten soil model with specific parameters for testing."""
    return Genuchten(k_s=10.0, theta_r=0.01, theta_s=0.43, alpha=0.02, n=1.1, l=0.5)
