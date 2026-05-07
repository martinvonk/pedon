import pytest
from numpy import logspace

import pedon as pe
from pedon.soilmodel import SoilModelConverter

h = -logspace(-2, 6, num=50)


@pytest.fixture
def gen() -> pe.Genuchten:
    return pe.Genuchten(k_s=10, theta_r=0.01, theta_s=0.43, alpha=0.02, n=2.1)


@pytest.fixture
def gen_low_n() -> pe.Genuchten:
    # n < 1.8, invalid for peche; n < 2, invalid for ghezzehei
    return pe.Genuchten(k_s=10, theta_r=0.01, theta_s=0.43, alpha=0.02, n=1.5)


@pytest.fixture
def bro() -> pe.Brooks:
    return pe.Brooks(k_s=10, theta_r=0.01, theta_s=0.43, h_b=10, l=2)


# --- list_methods ---


def test_list_methods_genuchten(gen: pe.Genuchten) -> None:
    methods = SoilModelConverter.list_methods(gen)
    assert isinstance(methods, list)
    assert methods == sorted(methods)
    assert "ghezzehei" in methods
    assert "peche" in methods
    assert "morel" in methods


def test_list_methods_brooks(bro: pe.Brooks) -> None:
    methods = SoilModelConverter.list_methods(bro)
    assert "morel" in methods
    assert "ghezzehei" not in methods
    assert "peche" not in methods


def test_list_methods_gardner(
    gar: pe.SoilModel = pe.Gardner(k_s=10, theta_s=0.43, c=0.02, m=1.1),
) -> None:
    assert SoilModelConverter.list_methods(gar) == []


# --- morel ---


def test_morel_genuchten_to_brooks(gen: pe.Genuchten) -> None:
    result = SoilModelConverter.morel(gen)
    assert isinstance(result, pe.Brooks)
    assert result.k_s == gen.k_s
    assert result.theta_r == gen.theta_r
    assert result.theta_s == gen.theta_s


def test_morel_brooks_to_genuchten(bro: pe.Brooks) -> None:
    result = SoilModelConverter.morel(bro)
    assert isinstance(result, pe.Genuchten)
    assert result.k_s == bro.k_s
    assert result.theta_r == bro.theta_r
    assert result.theta_s == bro.theta_s


def test_morel_roundtrip(gen: pe.Genuchten) -> None:
    bro = SoilModelConverter.morel(gen)
    gen2 = SoilModelConverter.morel(bro)
    assert isinstance(gen2, pe.Genuchten)


def test_morel_invalid_type_raises() -> None:
    gar = pe.Gardner(k_s=10, theta_s=0.43, c=0.02, m=1.1)
    with pytest.raises(TypeError):
        SoilModelConverter.morel(gar)  # type: ignore[arg-type]


# --- ghezzehei ---


def test_ghezzehei_returns_gardner(gen: pe.Genuchten) -> None:
    result = SoilModelConverter.ghezzehei(gen)
    assert isinstance(result, pe.Gardner)
    assert result.k_s == gen.k_s


def test_ghezzehei_invalid_type_raises() -> None:
    bro = pe.Brooks(k_s=10, theta_r=0.01, theta_s=0.43, h_b=10, l=2)
    with pytest.raises(TypeError):
        SoilModelConverter.ghezzehei(bro)  # type: ignore[arg-type]


def test_ghezzehei_low_n_raises(gen_low_n: pe.Genuchten) -> None:
    with pytest.raises(ValueError):
        SoilModelConverter.ghezzehei(gen_low_n)


def test_ghezzehei_c_value(gen: pe.Genuchten) -> None:
    result = SoilModelConverter.ghezzehei(gen)
    assert result.c == pytest.approx(1.3 * gen.n * gen.alpha)


# --- peche ---


def test_peche_returns_gardner(gen: pe.Genuchten) -> None:
    result = SoilModelConverter.peche(gen)
    assert isinstance(result, pe.Gardner)
    assert result.k_s == gen.k_s


def test_peche_invalid_type_raises() -> None:
    bro = pe.Brooks(k_s=10, theta_r=0.01, theta_s=0.43, h_b=10, l=2)
    with pytest.raises(TypeError):
        SoilModelConverter.peche(bro)  # type: ignore[arg-type]


def test_peche_low_n_raises(gen_low_n: pe.Genuchten) -> None:
    with pytest.raises(ValueError):
        SoilModelConverter.peche(gen_low_n)


def test_peche_c_value(gen: pe.Genuchten) -> None:
    result = SoilModelConverter.peche(gen)
    assert result.c == pytest.approx(2.2 * gen.alpha)
