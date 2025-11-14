import pytest
from numpy import allclose, array, logspace

import pedon as pe
from pedon._typing import FloatArray

h = -logspace(-2, 6, num=50)

theta = array([0.1, 0.2, 0.3, 0.4])


@pytest.fixture
def gen() -> pe.soilmodel.SoilModel:
    return pe.Genuchten(k_s=10, theta_r=0.01, theta_s=0.43, alpha=0.02, n=1.1, l=0.5)


@pytest.fixture
def bro() -> pe.soilmodel.SoilModel:
    return pe.Brooks(k_s=10, theta_r=0.01, theta_s=0.43, h_b=10, l=2)


@pytest.fixture
def sor() -> pe.soilmodel.SoilModel:
    return pe.Panday(k_s=10, theta_r=0.01, theta_s=0.43, alpha=0.02, beta=1.1, brook=3)


@pytest.fixture
def gar() -> pe.soilmodel.SoilModel:
    return pe.Gardner(k_s=10, theta_s=0.43, c=0.02, m=1.1)


@pytest.fixture
def hav() -> pe.soilmodel.Haverkamp:
    # Example parameters similar in style to other model tests
    return pe.Haverkamp(
        k_s=10.0,
        theta_r=0.01,
        theta_s=0.43,
        alpha=0.02,
        beta=1.2,
        a=0.5,
    )


def test_theta_genuchten(gen: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    gen.theta(h=h)


def test_s_genuchten(gen: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    gen.s(h=h)


def test_k_genuchten(gen: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    gen.k(h=h)


def test_h_genuchten(gen: pe.soilmodel.SoilModel, theta: FloatArray = theta) -> None:
    gen.h(theta=theta)


def test_theta_brooks(bro: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    bro.theta(h=h)


def test_s_brooks(bro: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    bro.s(h=h)


def test_k_brooks(bro: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    bro.k(h=h)


def test_h_brooks(bro: pe.soilmodel.SoilModel, theta: FloatArray = theta) -> None:
    bro.h(theta=theta)


def test_theta_panday(sor: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    sor.theta(h=h)


def test_s_panday(sor: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    sor.s(h=h)


def test_k_panday(sor: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    sor.k(h=h)


def test_h_panday(sor: pe.soilmodel.SoilModel, theta: FloatArray = theta) -> None:
    sor.h(theta=theta)


def test_theta_gardner(gar: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    gar.theta(h=h)


def test_s_gardner(gar: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    gar.s(h=h)


def test_k_gardner(gar: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    gar.k(h=h)


def test_h_gardner(gar: pe.soilmodel.SoilModel, theta: FloatArray = theta) -> None:
    gar.h(theta=theta)


def test_theta_haverkamp(hav: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    hav.theta(h=h)


def test_s_haverkamp(hav: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    hav.s(h=h)


def test_k_haverkamp(hav: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    hav.k(h=h)


def test_k_r_haverkamp_accepts_s_kwarg(
    hav: pe.soilmodel.Haverkamp,
) -> None:
    s = array([0.1, 0.2, 0.3])
    # Expect k_r(s) = (a*s)/(a*s + alpha*(1-s))
    expected = (hav.a * s) / (hav.a * s + hav.alpha * (1.0 - s))
    got = hav.k_r(h=array([1.0, 2.0, 3.0]), s=s)
    assert allclose(got, expected)


def test_h_haverkamp_inverse(
    hav: pe.soilmodel.Haverkamp, theta: FloatArray = theta
) -> None:
    # Check theta(h(h(theta))) ≈ theta
    h_out = hav.h(theta=theta)
    theta_back = hav.theta(h=h_out)
    assert allclose(theta_back, theta)
