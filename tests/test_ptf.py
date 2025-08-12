import pytest

import pedon as pe


@pytest.fixture
def ss() -> pe.soil.SoilSample:
    return pe.soil.SoilSample(
        sand_p=50, silt_p=10, clay_p=40, rho=1.5, om_p=20, m50=150
    )


def test_wosten(ss: pe.soil.SoilSample) -> None:
    sm = ss.wosten(ts=True)
    assert isinstance(sm, pe.Genuchten)


def test_wosten_sand(ss: pe.soil.SoilSample) -> None:
    sm = ss.wosten_sand(ts=True)
    assert isinstance(sm, pe.Genuchten)


def test_wosten_clay(ss: pe.soil.SoilSample) -> None:
    sm = ss.wosten_clay()
    assert isinstance(sm, pe.Genuchten)


def test_cosby(ss: pe.soil.SoilSample) -> None:
    sm = ss.cosby()
    assert isinstance(sm, pe.Brooks)


def test_rosetta(ss: pe.soil.SoilSample) -> None:
    sm = ss.rosetta()
    assert isinstance(sm, pe.Genuchten)


def test_rosetta_invalidpercentage(ss: pe.soil.SoilSample) -> None:
    ss.sand_p = 10
    with pytest.raises(ValueError):
        _ = ss.rosetta()
