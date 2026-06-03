"""Tests for pedotransfer functions."""

import numpy as np
import pytest

import pedon as pe

REL = 1e-3


@pytest.fixture
def ss() -> pe.soil.SoilSample:
    """Fixture for a soil sample with specific properties for testing pedotransfer functions."""
    return pe.soil.SoilSample(
        sand_p=50.0, silt_p=10.0, clay_p=38.0, rho=1.5, om_p=2.0, m50=150.0
    )


def test_wosten(ss: pe.soil.SoilSample) -> None:
    """Test Wösten pedotransfer function with texture-dependent parameters."""
    sm = ss.wosten(topsoil=True)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(61.41854253249469, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.4016384367793987, rel=REL)
    assert sm.alpha == pytest.approx(0.07053697798996453, rel=REL)
    assert sm.n == pytest.approx(1.137061406702455, rel=REL)
    assert sm.l == pytest.approx(-4.936831819765616, rel=REL)
    assert sm.m == pytest.approx(0.12054002175655676, rel=REL)


def test_wosten_no_topsoil(ss: pe.soil.SoilSample) -> None:
    """Test Wösten pedotransfer function without topsoil adjustment."""
    sm = ss.wosten(topsoil=False)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(10.843130926922838, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.4033024367793987, rel=REL)
    assert sm.alpha == pytest.approx(0.06618124289279576, rel=REL)
    assert sm.n == pytest.approx(1.104333140070533, rel=REL)
    assert sm.l == pytest.approx(-4.936831819765616, rel=REL)
    assert sm.m == pytest.approx(0.09447614699298912, rel=REL)


def test_wosten_sand(ss: pe.soil.SoilSample) -> None:
    """Test Wösten sand-specific pedotransfer function."""
    sm = ss.wosten_sand(topsoil=True)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(45.49, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.4478, rel=REL)
    assert sm.alpha == pytest.approx(0.0204414, rel=REL)
    assert sm.n == pytest.approx(1.6782, rel=REL)
    assert sm.l == pytest.approx(2.0, rel=REL)
    assert sm.m == pytest.approx(0.4041234656179239, rel=REL)


def test_wosten_sand_no_topsoil(ss: pe.soil.SoilSample) -> None:
    """Test Wösten sand-specific pedotransfer function without topsoil adjustment."""
    sm = ss.wosten_sand(topsoil=False)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(45.49, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.4478, rel=REL)
    assert sm.alpha == pytest.approx(0.0242778, rel=REL)
    assert sm.n == pytest.approx(1.6782, rel=REL)
    assert sm.l == pytest.approx(2.0, rel=REL)
    assert sm.m == pytest.approx(0.4041234656179239, rel=REL)


def test_wosten_clay(ss: pe.soil.SoilSample) -> None:
    """Test Wösten clay-specific pedotransfer function."""
    sm = ss.wosten_clay()
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(14.4541, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.4251, rel=REL)
    assert sm.alpha == pytest.approx(0.0542982, rel=REL)
    assert sm.n == pytest.approx(1.089, rel=REL)
    assert sm.l == pytest.approx(-6.3676, rel=REL)
    assert sm.m == pytest.approx(0.08172635445362719, rel=REL)


def test_cosby(ss: pe.soil.SoilSample) -> None:
    """Test Cosby pedotransfer function."""
    sm = ss.cosby()
    assert isinstance(sm, pe.Brooks)
    assert sm.k_s == pytest.approx(37.0179, rel=REL)
    assert sm.theta_r == pytest.approx(0.0, rel=REL)
    assert sm.theta_s == pytest.approx(0.4199, rel=REL)
    assert sm.h_b == pytest.approx(13.82293, rel=REL)
    assert sm.l == pytest.approx(0.11216, rel=REL)


def test_saxton(ss: pe.soil.SoilSample) -> None:
    """Test Saxton-Rawls pedotransfer function."""
    sm = ss.saxton()
    assert isinstance(sm, pe.Brooks)
    assert sm.k_s == pytest.approx(4.4696, rel=REL)
    assert sm.theta_r == pytest.approx(0.0, rel=REL)
    assert sm.theta_s == pytest.approx(0.4387, rel=REL)
    assert sm.h_b == pytest.approx(34.47941, rel=REL)
    assert sm.l == pytest.approx(0.10203, rel=REL)


def test_saxton_density_factor(ss: pe.soil.SoilSample) -> None:
    """Test Saxton-Rawls pedotransfer function with a density adjustment."""
    sm = ss.saxton(df=1.1)
    assert isinstance(sm, pe.Brooks)
    assert sm.k_s == pytest.approx(0.6079, rel=REL)
    assert sm.theta_r == pytest.approx(0.0, rel=REL)
    assert sm.theta_s == pytest.approx(0.3825, rel=REL)
    assert sm.h_b == pytest.approx(58.91006, rel=REL)
    assert sm.l == pytest.approx(0.09343, rel=REL)


def test_rawls(ss: pe.soil.SoilSample) -> None:
    """Test Rawls-Braakensiek pedotransfer function with measured bulk density."""
    sm = ss.rawls()
    assert isinstance(sm, pe.Brooks)
    assert sm.k_s == pytest.approx(6.9961132379760835, rel=REL)
    assert sm.theta_r == pytest.approx(0.1134085532217871, rel=REL)
    assert sm.theta_s == pytest.approx(0.43396226415094336, rel=REL)
    assert sm.h_b == pytest.approx(19.767218479230795, rel=REL)
    assert sm.l == pytest.approx(0.18238934593362763, rel=REL)


def test_rawls_with_cecc(ss: pe.soil.SoilSample) -> None:
    """Test Rawls-Braakensiek pedotransfer function when bulk density is estimated."""
    ss.rho = None
    sm = ss.rawls(cecc=20.0)

    assert isinstance(sm, pe.Brooks)
    assert sm.k_s == pytest.approx(23.56273325252222, rel=REL)
    assert sm.theta_r == pytest.approx(0.10248398945322691, rel=REL)
    assert sm.theta_s == pytest.approx(0.48550943396226426, rel=REL)
    assert sm.h_b == pytest.approx(13.370760200694354, rel=REL)
    assert sm.l == pytest.approx(0.20448715068488366, rel=REL)


def test_rawls_requires_cecc_when_rho_missing(ss: pe.soil.SoilSample) -> None:
    """Test Rawls-Braakensiek raises when rho is missing and cecc is not provided."""
    ss.rho = None

    with pytest.raises(AssertionError, match="requires 'cecc' to be provided"):
        ss.rawls()


def test_vereecken(ss: pe.soil.SoilSample) -> None:
    """Test Vereecken pedotransfer function."""
    sm = ss.vereecken()
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(6.2506, rel=REL)
    assert sm.theta_r == pytest.approx(0.2212, rel=REL)
    assert sm.theta_s == pytest.approx(0.4235, rel=REL)
    assert sm.alpha == pytest.approx(0.0015921, rel=REL)
    assert sm.n == pytest.approx(0.5969, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(1.0, rel=REL)


def test_weynants(ss: pe.soil.SoilSample) -> None:
    """Test Weynants pedotransfer function."""
    sm = ss.weynants()
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(10.9707, rel=REL)
    assert sm.theta_r == pytest.approx(0.0, rel=REL)
    assert sm.theta_s == pytest.approx(0.4402, rel=REL)
    assert sm.alpha == pytest.approx(0.0166727, rel=REL)
    assert sm.n == pytest.approx(1.1157, rel=REL)
    assert sm.l == pytest.approx(-6.5338, rel=REL)
    assert sm.m == pytest.approx(0.1037017119297301, rel=REL)


def test_toth(ss: pe.soil.SoilSample) -> None:
    """Test Tóth pedotransfer function."""
    sm = ss.toth()
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(13.4896, rel=REL)
    assert sm.theta_r == pytest.approx(0.041, rel=REL)
    assert sm.theta_s == pytest.approx(0.4198, rel=REL)
    assert sm.alpha == pytest.approx(0.0145835, rel=REL)
    assert sm.n == pytest.approx(1.3347, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(0.250767962838091, rel=REL)


def test_toth_topsoil(ss: pe.soil.SoilSample) -> None:
    """Test Tóth pedotransfer function with topsoil adjustment."""
    sm = ss.toth(topsoil=True)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(1.0233, rel=REL)
    assert sm.theta_r == pytest.approx(0.041, rel=REL)
    assert sm.theta_s == pytest.approx(0.4198, rel=REL)
    assert sm.alpha == pytest.approx(0.0240969, rel=REL)
    assert sm.n == pytest.approx(1.2945, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(0.22750096562379296, rel=REL)


def test_toth_no_topsoil(ss: pe.soil.SoilSample) -> None:
    """Test Tóth pedotransfer function without topsoil adjustment."""
    sm = ss.toth(topsoil=False)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(13.4896, rel=REL)
    assert sm.theta_r == pytest.approx(0.041, rel=REL)
    assert sm.theta_s == pytest.approx(0.4198, rel=REL)
    assert sm.alpha == pytest.approx(0.0145835, rel=REL)
    assert sm.n == pytest.approx(1.3347, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(0.250767962838091, rel=REL)


def test_hodnett(ss: pe.soil.SoilSample) -> None:
    """Test Hodnett and Tomasella pedotransfer function."""
    sm = ss.hodnett()
    assert isinstance(sm, pe.Genuchten)
    assert np.isnan(sm.k_s)
    assert sm.theta_r == pytest.approx(0.1904, rel=REL)
    assert sm.theta_s == pytest.approx(0.4069, rel=REL)
    assert sm.alpha == pytest.approx(0.036191, rel=REL)
    assert sm.n == pytest.approx(1.4516, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(0.3111049875998898, rel=REL)


def test_rosetta(ss: pe.soil.SoilSample) -> None:
    """Test ROSETTA pedotransfer function with mocked HTTP response."""
    from unittest.mock import Mock, patch

    # Mock response data matching the structure returned by the ROSETTA API
    mock_response = Mock()
    mock_response.is_error = False
    # van_genuchten_params: [theta_r, theta_s, log10(alpha), log10(n), log10(k_s)]
    mock_response.json.return_value = {
        "van_genuchten_params": [
            [
                0.1142278836409842,  # theta_r
                0.42963731868993743,  # theta_s
                -1.8689880666646325,  # log10(alpha)
                0.10546903898903125,  # log10(n)
                1.1367547932846525,  # log10(k_s)
            ]
        ]
    }

    with patch("httpx.post", return_value=mock_response) as mock_post:
        sm = ss.rosetta()

        # Verify the API was called with correct parameters
        assert mock_post.called
        call_args = mock_post.call_args
        assert "rosetta/3" in call_args[0][0]  # Default version is 3

        # Verify response is parsed and mapped correctly
        assert isinstance(sm, pe.Genuchten)
        assert sm.k_s == pytest.approx(13.7019747459841, rel=REL)
        assert sm.theta_r == pytest.approx(0.1142278836409842, rel=REL)
        assert sm.theta_s == pytest.approx(0.42963731868993743, rel=REL)
        assert sm.alpha == pytest.approx(0.013518839874177326, rel=REL)
        assert sm.n == pytest.approx(1.274931762885784, rel=REL)
        assert sm.l == pytest.approx(0.5, rel=REL)
        assert sm.m == pytest.approx(0.2156442963374613, rel=REL)


def test_rosetta_invalidpercentage(ss: pe.soil.SoilSample) -> None:
    """Test that ROSETTA raises ValueError when API returns None values."""
    from unittest.mock import Mock, patch

    mock_response = Mock()
    mock_response.is_error = False
    # Response with None values in van_genuchten_params
    mock_response.json.return_value = {
        "van_genuchten_params": [[None, None, None, None, None]]
    }

    with patch("httpx.post", return_value=mock_response):
        with pytest.raises(ValueError, match="Rosetta API returned None values"):
            _ = ss.rosetta()


def test_rosetta_api_error(ss: pe.soil.SoilSample) -> None:
    """Test that ROSETTA raises ValueError on API error responses."""
    from unittest.mock import Mock, patch

    mock_response = Mock()
    mock_response.is_error = True

    with patch("httpx.post", return_value=mock_response):
        with pytest.raises(ValueError, match="Rosetta API error"):
            _ = ss.rosetta()
