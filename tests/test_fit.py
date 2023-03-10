import pytest
from numpy import array

import pedon as pe


@pytest.fixture
def sample() -> pe.soil.SoilSample:
    h = array(
        [
            -1.0e00,
            -1.0e01,
            -2.0e01,
            -3.1e01,
            -5.0e01,
            -1.0e02,
            -2.5e02,
            -5.0e02,
            -1.0e03,
            -2.5e03,
            -5.0e03,
            -1.0e04,
            -1.6e04,
        ]
    )
    theta = array(
        [
            0.43,
            0.417,
            0.391,
            0.356,
            0.302,
            0.21,
            0.118,
            0.077,
            0.053,
            0.036,
            0.029,
            0.025,
            0.024,
        ]
    )

    k = array(
        [
            2.341e01,
            1.138e01,
            6.040e00,
            3.130e00,
            1.140e00,
            1.600e-01,
            7.500e-03,
            6.500e-04,
            5.400e-05,
            2.000e-06,
            1.600e-07,
            1.400e-08,
            2.600e-09,
        ]
    )

    return pe.soil.SoilSample(h=h, theta=theta, k=k)


def test_fit(sample: pe.soil.SoilSample) -> None:
    sample.fit(pe.soilmodel.Genuchten)


def test_fit_seperate(sample: pe.soil.SoilSample) -> None:
    sample.fit(pe.soilmodel.Brooks)
