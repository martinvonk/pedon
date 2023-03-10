import pytest
from numpy import logspace

import pedon as pe
from pedon._typing import FloatArray

h = -logspace(-2, 6, num=50)


@pytest.fixture
def gen() -> pe.soilmodel.SoilModel:
    return pe.Genuchten(k_s=10, theta_r=0.01, theta_s=0.43, alpha=0.02, n=1.1, l=0.5)


@pytest.fixture
def bro() -> pe.soilmodel.SoilModel:
    return pe.Brooks(k_s=10, theta_r=0.01, theta_s=0.43, h_b=10, l=2)


@pytest.fixture
def sor() -> pe.soilmodel.SoilModel:
    return pe.Panday(k_s=10, theta_r=0.01, theta_s=0.43, alpha=0.02, beta=1.1, brook=3)


def test_theta_genuchten(gen: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    gen.theta(h=h)


def test_s_genuchten(gen: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    gen.s(h=h)


def test_k_genuchten(gen: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    gen.k(h=h)


def test_theta_brooks(bro: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    bro.theta(h=h)


def test_s_brooks(bro: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    bro.s(h=h)


def test_k_brooks(bro: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    bro.k(h=h)


def test_theta_panday(sor: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    sor.theta(h=h)


def test_s_panday(sor: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    sor.s(h=h)


def test_k_panday(sor: pe.soilmodel.SoilModel, h: FloatArray = h) -> None:
    sor.k(h=h)
