import pytest

import pedon as pe


@pytest.fixture
def ss() -> pe.soil.SoilSample:
    return pe.soil.SoilSample(
        sand_p=40, silt_p=10, clay_p=30, rho=1.5, om_p=20, m50=10000
    )


def test_wosten(ss: pe.soil.SoilSample) -> None:
    ss.wosten(ts=True)


def test_wosten_sand(ss: pe.soil.SoilSample) -> None:
    ss.wosten_sand(ts=True)


def test_wosten_clay(ss: pe.soil.SoilSample) -> None:
    ss.wosten_clay()


def test_cosby(ss: pe.soil.SoilSample) -> None:
    ss.cosby()
