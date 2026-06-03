"""Tests for pedotransfer functions."""

import numpy as np
import pytest

import pedon as pe

REL = 1e-3


@pytest.fixture
def ss() -> pe.soil.SoilSample:
    """Fixture for a soil sample with specific properties for testing pedotransfer functions."""
    return pe.soil.SoilSample(
        sand_p=50.0, silt_p=10.0, clay_p=40.0, rho=1.5, om_p=20.0, m50=150.0
    )


def test_wosten(ss: pe.soil.SoilSample) -> None:
    """Test Wösten pedotransfer function with texture-dependent parameters."""
    sm = ss.wosten(topsoil=True)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(0.6547824638330241, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.3506329025688723, rel=REL)
    assert sm.alpha == pytest.approx(0.0014029528494249415, rel=REL)
    assert sm.n == pytest.approx(1.1209348010913704, rel=REL)
    assert sm.l == pytest.approx(-3.6143956260446446, rel=REL)
    assert sm.m == pytest.approx(0.10788745337697181, rel=REL)


def test_wosten_no_topsoil(ss: pe.soil.SoilSample) -> None:
    """Test Wösten pedotransfer function without topsoil adjustment."""
    sm = ss.wosten(topsoil=False)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(0.10889706346920575, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.3522969025688723, rel=REL)
    assert sm.alpha == pytest.approx(0.001298720038387063, rel=REL)
    assert sm.n == pytest.approx(1.0907448358614613, rel=REL)
    assert sm.l == pytest.approx(-3.6143956260446446, rel=REL)
    assert sm.m == pytest.approx(0.08319529268253811, rel=REL)


def test_wosten_sand(ss: pe.soil.SoilSample) -> None:
    """Test Wösten sand-specific pedotransfer function."""
    sm = ss.wosten_sand(topsoil=True)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(20.793, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.4959, rel=REL)
    assert sm.alpha == pytest.approx(0.0204414, rel=REL)
    assert sm.n == pytest.approx(2.2182, rel=REL)
    assert sm.l == pytest.approx(2.0, rel=REL)
    assert sm.m == pytest.approx(0.5491840230817779, rel=REL)


def test_wosten_sand_no_topsoil(ss: pe.soil.SoilSample) -> None:
    """Test Wösten sand-specific pedotransfer function without topsoil adjustment."""
    sm = ss.wosten_sand(topsoil=False)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(20.793, rel=REL)
    assert sm.theta_r == pytest.approx(0.01, rel=REL)
    assert sm.theta_s == pytest.approx(0.4959, rel=REL)
    assert sm.alpha == pytest.approx(0.0242778, rel=REL)
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


def test_saxton(ss: pe.soil.SoilSample) -> None:
    """Test Saxton-Rawls pedotransfer function."""
    sm = ss.saxton()
    assert isinstance(sm, pe.Brooks)
    assert sm.k_s == pytest.approx(17.1624, rel=REL)
    assert sm.theta_r == pytest.approx(0.0, rel=REL)
    assert sm.theta_s == pytest.approx(0.5961, rel=REL)
    assert sm.h_b == pytest.approx(6.39505, rel=REL)
    assert sm.l == pytest.approx(0.09416, rel=REL)


def test_saxton_density_factor(ss: pe.soil.SoilSample) -> None:
    """Test Saxton-Rawls pedotransfer function with a density adjustment."""
    sm = ss.saxton(df=1.1)
    assert isinstance(sm, pe.Brooks)
    assert sm.k_s == pytest.approx(8.1961, rel=REL)
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
    assert sm.k_s == pytest.approx(2.0632, rel=REL)
    assert sm.theta_r == pytest.approx(0.3774, rel=REL)
    assert sm.theta_s == pytest.approx(0.4255, rel=REL)
    assert sm.alpha == pytest.approx(3.89e-05, rel=REL)
    assert sm.n == pytest.approx(0.5816, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(1.0, rel=REL)


def test_weynants(ss: pe.soil.SoilSample) -> None:
    """Test Weynants pedotransfer function."""
    sm = ss.weynants()
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(2.1387, rel=REL)
    assert sm.theta_r == pytest.approx(0.0, rel=REL)
    assert sm.theta_s == pytest.approx(0.4429, rel=REL)
    assert sm.alpha == pytest.approx(0.0058046, rel=REL)
    assert sm.n == pytest.approx(1.1104, rel=REL)
    assert sm.l == pytest.approx(-6.7972, rel=REL)
    assert sm.m == pytest.approx(0.09942363112391939, rel=REL)


def test_toth(ss: pe.soil.SoilSample) -> None:
    """Test Tóth pedotransfer function."""
    sm = ss.toth()
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(0.5623, rel=REL)
    assert sm.theta_r == pytest.approx(0.041, rel=REL)
    assert sm.theta_s == pytest.approx(0.4203, rel=REL)
    assert sm.alpha == pytest.approx(0.0043157, rel=REL)
    assert sm.n == pytest.approx(1.2524, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(0.2015330565314596, rel=REL)


def test_toth_topsoil(ss: pe.soil.SoilSample) -> None:
    """Test Tóth pedotransfer function with topsoil adjustment."""
    sm = ss.toth(topsoil=True)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(27.5423, rel=REL)
    assert sm.theta_r == pytest.approx(0.041, rel=REL)
    assert sm.theta_s == pytest.approx(0.4203, rel=REL)
    assert sm.alpha == pytest.approx(0.007131, rel=REL)
    assert sm.n == pytest.approx(1.2221, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(0.18173635545372713, rel=REL)


def test_toth_no_topsoil(ss: pe.soil.SoilSample) -> None:
    """Test Tóth pedotransfer function without topsoil adjustment."""
    sm = ss.toth(topsoil=False)
    assert isinstance(sm, pe.Genuchten)
    assert sm.k_s == pytest.approx(0.5623, rel=REL)
    assert sm.theta_r == pytest.approx(0.041, rel=REL)
    assert sm.theta_s == pytest.approx(0.4203, rel=REL)
    assert sm.alpha == pytest.approx(0.0043157, rel=REL)
    assert sm.n == pytest.approx(1.2524, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(0.2015330565314596, rel=REL)


def test_hodnett(ss: pe.soil.SoilSample) -> None:
    """Test Hodnett and Tomasella pedotransfer function."""
    sm = ss.hodnett()
    assert isinstance(sm, pe.Genuchten)
    assert np.isnan(sm.k_s)
    assert sm.theta_r == pytest.approx(0.1878, rel=REL)
    assert sm.theta_s == pytest.approx(0.4121, rel=REL)
    assert sm.alpha == pytest.approx(0.0466917, rel=REL)
    assert sm.n == pytest.approx(1.3657, rel=REL)
    assert sm.l == pytest.approx(0.5, rel=REL)
    assert sm.m == pytest.approx(0.2677747675184886, rel=REL)


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
