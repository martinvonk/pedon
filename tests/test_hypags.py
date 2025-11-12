import numpy as np
import pandas as pd
import pytest

from pedon.soil import SoilSample
from pedon.soilmodel import Genuchten

# Reference data from HYPAGS GUI for validation
HYPAGS_K_REFERENCE = pd.DataFrame(
    {
        "ne": [0.308, 0.283, 0.258, 0.23, 0.201],
        "alpha": [9.44, 6.41, 3.11, 1.21, 0.48],
        "alpha_upper": [11.25, 9.46, 8.77, 7.94, 3.49],
        "alpha_lower": [6.58, 1.38, 0.1, 0.1, 0.1],
        "n": [1.35, 1.44, 1.76, 1.8, 1.4],
        "n_upper": [4.35, 5.57, 4.94, 4.86, 4.23],
        "n_lower": [0.97, 1.08, 1.15, 0.97, 1.04],
    },
    index=pd.Index([1e-3, 1e-4, 1e-5, 1e-6, 1e-7], name="k"),
    dtype=float,
)

HYPAGS_D10_REFERENCE = pd.DataFrame(
    {
        "ne": [0.269, 0.278, 0.296, 0.311],
        "alpha": [4.81, 6.28, 7.12, 10.25],
        "alpha_upper": [8.89, 9.47, 10.56, 12.06],
        "alpha_lower": [0.57, 1.8, 1.63, 7.4],
        "n": [1.53, 1.44, 1.41, 1.33],
        "n_upper": [3.44, 5.47, 2.79, 4.34],
        "n_lower": [1.03, 1.09, 1.17, 0.95],
        "k": [2.7934e-5, 6.2851e-5, 3.1038e-4, 1.2415e-3],
    },
    index=pd.Index([6e-5, 9e-5, 2e-4, 4e-4], name="d10"),
    dtype=float,
)

HYPAGS_D20_REFERENCE = pd.DataFrame(
    {
        "ne": [0.272, 0.301, 0.315, 0.324],
        "alpha": [5.25, 7.71, 12.61, 20.77],
        "alpha_upper": [9.33, 11.15, 14.42, 22.58],
        "alpha_lower": [1.02, 2.22, 9.75, 17.91],
        "n": [1.5, 1.39, 1.3, 1.24],
        "n_upper": [2.68, 2.77, 4.3, 4.25],
        "n_lower": [1.2, 1.15, 0.92, 0.87],
        "k": [3.4486e-5, 4.8496e-4, 1.9399e-3, 4.3647e-3],
    },
    index=pd.Index([8e-5, 3e-4, 6e-4, 9e-4], name="d20"),
    dtype=float,
)

RTOL = 1e-2


@pytest.mark.parametrize("k", HYPAGS_K_REFERENCE.index.tolist())
def test_hypags_with_k(k: float):
    """Test hypags with valid k input."""
    sample = SoilSample(k=np.array([k], dtype=float))
    result = sample.hypags()

    # basic type check
    assert isinstance(result, Genuchten)
    # check that d10/d20 were generated
    assert isinstance(sample.d10, float)
    assert isinstance(sample.d20, float)
    # check that some output params same as HYPAGS GUI
    assert np.isclose(result.theta_s, HYPAGS_K_REFERENCE.at[k, "ne"], rtol=RTOL)
    assert np.isclose(result.alpha, HYPAGS_K_REFERENCE.at[k, "alpha"], rtol=RTOL)
    assert np.isclose(result.n, HYPAGS_K_REFERENCE.at[k, "n"], rtol=RTOL)


@pytest.mark.parametrize("d10", HYPAGS_D10_REFERENCE.index.tolist())
def test_hypags_with_d10(d10: float):
    """Test hypags when only d10 is given."""
    sample = SoilSample(d10=d10)
    result = sample.hypags()

    assert isinstance(result, Genuchten)
    assert isinstance(sample.k, np.ndarray)
    # check that some output params same as HYPAGS GUI
    assert sample.d20 > sample.d10  # by c1 * d10
    assert np.isclose(result.theta_s, HYPAGS_D10_REFERENCE.at[d10, "ne"], rtol=RTOL)
    assert np.isclose(result.alpha, HYPAGS_D10_REFERENCE.at[d10, "alpha"], rtol=RTOL)
    assert np.isclose(result.n, HYPAGS_D10_REFERENCE.at[d10, "n"], rtol=RTOL)
    assert np.isclose(result.k_s, HYPAGS_D10_REFERENCE.at[d10, "k"], rtol=RTOL)


@pytest.mark.parametrize("d20", HYPAGS_D20_REFERENCE.index.tolist())
def test_hypags_with_d20(d20: float):
    """Test hypags when only d20 is given."""
    sample = SoilSample(d20=d20)
    result = sample.hypags()

    assert isinstance(result, Genuchten)
    assert isinstance(sample.k, np.ndarray)
    # check that some output params same as HYPAGS GUI
    assert sample.d20 > sample.d10  # by c1 * d10
    assert np.isclose(
        result.theta_s, HYPAGS_D20_REFERENCE.at[sample.d20, "ne"], rtol=RTOL
    )
    assert np.isclose(
        result.alpha, HYPAGS_D20_REFERENCE.at[sample.d20, "alpha"], rtol=RTOL
    )
    assert np.isclose(result.n, HYPAGS_D20_REFERENCE.at[sample.d20, "n"], rtol=RTOL)
    assert np.isclose(result.k_s, HYPAGS_D20_REFERENCE.at[sample.d20, "k"], rtol=RTOL)


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

    assert any("out of hypags model limits" in rec.message for rec in caplog.records)


def test_hypags_logs_out_of_range_d10(caplog):
    """Test that hypags logs warning when d10 is out of limits."""
    sample = SoilSample(d10=1e-3)
    sample.hypags()

    assert any("out of hypags model limits" in rec.message for rec in caplog.records)


def test_hypags_logs_out_of_range_d20(caplog):
    """Test that hypags logs warning when d20 is out of limits."""
    sample = SoilSample(d20=1e-6)
    sample.hypags()

    assert any("out of hypags model limits" in rec.message for rec in caplog.records)
