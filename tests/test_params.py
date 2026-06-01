"""Tests for parameter helpers."""

import pedon as pe
from pedon._params import get_params


def test_get_params_warns_for_model_without_defaults(caplog) -> None:
    """Unknown default bounds should log a warning and still return defaults."""
    with caplog.at_level("WARNING"):
        params = get_params(pe.Haverkamp)

    assert "No default parameter bounds for SoilModel type Haverkamp" in caplog.text
    assert "p_ini" in params.columns
    assert "p_min" in params.columns
    assert "p_max" in params.columns
    assert "k_s" in params.index
