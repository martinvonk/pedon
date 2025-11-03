import numpy as np
import pytest

from pedon.soil import SoilSample  # adjust import path if different
from pedon.soilmodel import Genuchten


@pytest.mark.parametrize("k_value", [1e-5, 1e-6])
def test_hypags_with_k(k_value):
    """Test hypags with valid k input."""
    k = np.array(
        [
            k_value,
        ],
        dtype=float,
    )
    sample = SoilSample(k=k)

    result = sample.hypags()
    print(result)
    # basic type check
    assert isinstance(result, Genuchten)
    # check that d10/d20 were generated
    assert sample.d10 is not None
    assert sample.d20 is not None
    # check that some output params are positive
    # TODO: assert results with values from HYPAGS gui
    assert result.k_s > 0
    assert result.alpha > 0
    assert result.n > 0


def test_hypags_with_d10():
    """Test hypags when only d10 is given."""
    sample = SoilSample(d10=1e-4)
    result = sample.hypags()

    assert isinstance(result, Genuchten)
    assert isinstance(sample.k, np.ndarray)
    # TODO: assert results with values from HYPAGS gui
    assert sample.d20 > sample.d10  # by c1 * d10
    assert result.alpha > 0
    assert result.n > 0


def test_hypags_with_d20():
    """Test hypags when only d20 is given."""
    sample = SoilSample(d20=1.2e-4)
    result = sample.hypags()

    assert isinstance(result, Genuchten)
    assert isinstance(sample.k, np.ndarray)
    # TODO: assert results with values from HYPAGS gui
    assert np.isclose(sample.d10 * 1.2, sample.d20, rtol=1e-2)
    assert result.alpha > 0
    assert result.n > 0


def test_hypags_raises_value_error():
    """Test that hypags raises ValueError when no k, d10, or d20 are provided."""
    sample = SoilSample()
    with pytest.raises(ValueError, match="No parameter"):
        sample.hypags()


def test_hypags_logs_out_of_range_k(caplog):
    """Test that hypags logs warning when k is out of limits."""
    sample = SoilSample(
        k=np.array(
            [
                1e-1,
            ],
            dtype=float,
        )
    )  # much higher than supported range
    sample.hypags()

    assert any("k out of hypags model limits" in rec.message for rec in caplog.records)


def test_hypags_logs_out_of_range_d10(caplog):
    """Test that hypags logs warning when d10 is out of limits."""
    sample = SoilSample(d10=1e-3)
    sample.hypags()

    assert any(
        "d10 out of hypags model limits" in rec.message for rec in caplog.records
    )


def test_hypags_internal_values_stable():
    """Smoke test: check internal relationships."""
    sample = SoilSample(
        k=np.array(
            [
                1e-6,
            ],
            dtype=float,
        )
    )
    result = sample.hypags()

    # TODO: assert results with values from HYPAGS gui
    # d50, d60 relations
    assert np.isclose(sample.d60 / sample.d50, 1.13, atol=0.05)
    # ne must be between 0 and 1
    assert 0 < sample.ne < 1
    # ensure Genuchten parameters have expected ranges
    assert 0 < result.alpha < 10
    assert 1 < result.n < 10
