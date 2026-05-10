"""Tests for pedotransfer functions."""

import pytest

import pedon as pe

REL = 1e-5


@pytest.fixture
def ss() -> pe.soil.SoilSample:
    """Fixutre for a soil sample with specific properties for testing pedotransfer functions."""
    return pe.soil.SoilSample(
        sand_p=50.0, silt_p=10.0, clay_p=40.0, rho=1.5, om_p=20.0, m50=150.0
    )


def test_wosten(ss: pe.soil.SoilSample) -> None:
    """Test Wösten pedotransfer function with texture-dependent parameters."""
    sm = ss.wosten(ts=True)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(0.6547824638330241, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.3506329025688723, rel=REL)
    assert sm.alpha == pytest.approx(0.0014029528494249415, rel=REL)
    assert sm.n == pytest.approx(1.1209348010913704, rel=REL)
    assert sm.l == pytest.approx(-3.6143956260446446, rel=REL)
    assert sm.m == pytest.approx(0.10788745337697181, rel=REL)


def test_wosten_sand(ss: pe.soil.SoilSample) -> None:
    """Test Wösten sand-specific pedotransfer function."""
    sm = ss.wosten_sand(ts=True)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(20.793, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.4959, rel=REL)
    assert sm.alpha == pytest.approx(0.0204414, rel=REL)
    assert sm.n == pytest.approx(2.2182, rel=REL)
    assert sm.l == pytest.approx(2.0, rel=REL)
    assert sm.m == pytest.approx(0.5491840230817779, rel=REL)


def test_wosten_clay(ss: pe.soil.SoilSample) -> None:
    """Test Wösten clay-specific pedotransfer function."""
    sm = ss.wosten_clay()
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(0.0, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.4258, rel=REL)
    assert sm.alpha == pytest.approx(0.0, rel=REL)
    assert sm.n == pytest.approx(1.2582, rel=REL)
    assert sm.l == pytest.approx(-6.6123, rel=REL)
    assert sm.m == pytest.approx(0.20521379748847557, rel=REL)


def test_cosby(ss: pe.soil.SoilSample) -> None:
    """Test Cosby pedotransfer function."""
    sm = ss.cosby()
    assert isinstance(sm, pe.Brooks)
    assert sm.k_s == pytest.approx(35.9428, rel=REL)
    assert sm.theta_r == pytest.approx(0.0, rel=REL)
    assert sm.theta_s == pytest.approx(0.4192, rel=REL)
    assert sm.h_b == pytest.approx(13.42765, rel=REL)
    assert sm.l == pytest.approx(0.10834, rel=REL)


def test_rosetta(ss: pe.soil.SoilSample) -> None:
    """Test ROSETTA pedotransfer function."""
    sm = ss.rosetta()
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(13.7019747459841, rel=REL)
    assert sm.theta_r == pytest.approx(0.1142278836409842, rel=REL)
    assert sm.theta_s == pytest.approx(0.42963731868993743, rel=REL)
    assert sm.alpha == pytest.approx(0.013518839874177326, rel=REL)
    assert sm.n == pytest.approx(1.274931762885784, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(0.2156442963374613, rel=REL)


def test_rosetta_invalidpercentage(ss: pe.soil.SoilSample) -> None:
    """Test that ROSETTA raises ValueError for invalid soil percentages."""
    ss.sand_p = 10
    with pytest.raises(ValueError):
        _ = ss.rosetta()
