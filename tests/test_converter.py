import pytest

import pedon as pe
from pedon.soilmodel import SoilModelConverter


def test_list_methods_genuchten(genuchten: pe.Genuchten) -> None:
    """Test listing available conversion methods for a Genuchten soil model."""
    methods = SoilModelConverter.list_methods(genuchten)
    assert isinstance(methods, list)
    assert methods == sorted(methods)
    assert "ghezzehei" in methods
    assert "peche" in methods
    assert "morel" in methods


def test_list_methods_brooks(brooks: pe.Brooks) -> None:
    """Test listing available conversion methods for a Brooks-Corey soil model."""
    methods = SoilModelConverter.list_methods(brooks)
    assert "morel" in methods
    assert "ghezzehei" not in methods
    assert "peche" not in methods


def test_list_methods_gardner(
    gardner: pe.Gardner,
) -> None:
    """Test listing available conversion methods for a Gardner soil model."""
    assert SoilModelConverter.list_methods(gardner) == []


# --- morel ---


def test_morel_genuchten_to_brooks(genuchten: pe.Genuchten) -> None:
    """Test converting a Genuchten soil model to a Brooks-Corey soil model using the Morel method."""
    result = SoilModelConverter.morel(genuchten)
    assert isinstance(result, pe.Brooks)
    assert result.k_s == genuchten.k_s
    assert result.theta_r == genuchten.theta_r
    assert result.theta_s == genuchten.theta_s


def test_morel_brooks_to_genuchten(bro: pe.Brooks) -> None:
    """Test converting a Brooks-Corey soil model to a Genuchten soil model using the Morel method."""
    result = SoilModelConverter.morel(bro)
    assert isinstance(result, pe.Genuchten)
    assert result.k_s == bro.k_s
    assert result.theta_r == bro.theta_r
    assert result.theta_s == bro.theta_s


def test_morel_roundtrip(genuchten: pe.Genuchten) -> None:
    """Test round-trip conversion from Genuchten to Brooks-Corey and back to Genuchten using the Morel method."""
    bro = SoilModelConverter.morel(genuchten)
    gen2 = SoilModelConverter.morel(bro)
    assert isinstance(gen2, pe.Genuchten)


def test_morel_invalid_type_raises() -> None:
    """Test that passing an unsupported soil model type to the Morel method raises a TypeError."""
    gar = pe.Gardner(k_s=10, theta_s=0.43, c=0.02, m=1.1)
    with pytest.raises(TypeError):
        SoilModelConverter.morel(gar)  # type: ignore[arg-type]


# --- ghezzehei ---


def test_ghezzehei_returns_gardner(genuchten: pe.Genuchten) -> None:
    """Test converting a Genuchten soil model to a Gardner soil model using the Ghezzehei method."""
    genuchten.n = 2.1
    result = SoilModelConverter.ghezzehei(genuchten)
    assert isinstance(result, pe.Gardner)
    assert result.k_s == genuchten.k_s


def test_ghezzehei_invalid_type_raises() -> None:
    """Test that passing an unsupported soil model type to the Ghezzehei method raises a TypeError."""
    bro = pe.Brooks(k_s=10, theta_r=0.01, theta_s=0.43, h_b=10, l=2)
    with pytest.raises(TypeError):
        SoilModelConverter.ghezzehei(bro)  # type: ignore[arg-type]


def test_ghezzehei_low_n_raises(genuchten: pe.Genuchten) -> None:
    """Test that passing a Genuchten model with n < 1 to the Ghezzehei method raises a ValueError."""
    with pytest.raises(ValueError):
        SoilModelConverter.ghezzehei(genuchten)


def test_ghezzehei_c_value(genuchten: pe.Genuchten) -> None:
    """Test that the c value of the Gardner model returned by the Ghezzehei method is approximately 1.3 * n * alpha."""
    genuchten.n = 2.1
    result = SoilModelConverter.ghezzehei(genuchten)
    assert result.c == pytest.approx(1.3 * genuchten.n * genuchten.alpha)


# --- peche ---


def test_peche_returns_gardner(genuchten: pe.Genuchten) -> None:
    """Test converting a Genuchten soil model to a Gardner soil model using the Peche method."""
    result = SoilModelConverter.peche(genuchten)
    assert isinstance(result, pe.Gardner)
    assert result.k_s == genuchten.k_s


def test_peche_invalid_type_raises() -> None:
    """Test that passing an unsupported soil model type to the Peche method raises a TypeError."""
    bro = pe.Brooks(k_s=10, theta_r=0.01, theta_s=0.43, h_b=10, l=2)
    with pytest.raises(TypeError):
        SoilModelConverter.peche(bro)  # type: ignore[arg-type]


def test_peche_low_n_raises(genuchten: pe.Genuchten) -> None:
    """Test that passing a Genuchten model with n < 1 to the Peche method raises a ValueError."""
    with pytest.raises(ValueError):
        SoilModelConverter.peche(genuchten)


def test_peche_c_value(genuchten: pe.Genuchten) -> None:
    """Test that the c value of the Gardner model returned by the Peche method is approximately 2.2 * alpha."""
    result = SoilModelConverter.peche(genuchten)
    assert result.c == pytest.approx(2.2 * genuchten.alpha)
