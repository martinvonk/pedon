from typing import Literal

from pandas import DataFrame

pGenuchten = DataFrame(
    data={
        "p_ini": {
            "k_s": 10.0,
            "theta_r": 0.01,
            "theta_s": 0.40,
            "alpha": 0.02,
            "n": 2.0,
            "l": 0.5,
        },
        "p_min": {
            "k_s": 0.001,
            "theta_r": 0.0,
            "theta_s": 0.2,
            "alpha": 0.001,
            "n": 1.000001,
            "l": -7.0,
        },
        "p_max": {
            "k_s": 100000.0,
            "theta_r": 0.2,
            "theta_s": 0.9,
            "alpha": 0.20,
            "n": 12.0,
            "l": 8.0,
        },
    },
    dtype=float,
)

pBrooks = DataFrame(
    data={
        "p_ini": {"k_s": 50.0, "theta_r": 0.02, "theta_s": 0.4, "h_b": 0.003, "l": 1.5},
        "p_min": {
            "k_s": 0.001,
            "theta_r": 0.0,
            "theta_s": 0.2,
            "h_b": 0.0001,
            "l": 0.1,
        },
        "p_max": {
            "k_s": 100000.0,
            "theta_r": 0.2,
            "theta_s": 0.5,
            "h_b": 100.0,
            "l": 5.0,
        },
    },
    dtype=float,
)
pPanday = DataFrame(
    data={
        "p_ini": {
            "k_s": 50.0,
            "theta_r": 0.02,
            "theta_s": 0.4,
            "alpha": 0.02,
            "beta": 2.3,
            "brook": 10.0,
        },
        "p_min": {
            "k_s": 0.001,
            "theta_r": 1e-5,
            "theta_s": 0.2,
            "alpha": 0.001,
            "beta": 1.0,
            "brook": 1.0,
        },
        "p_max": {
            "k_s": 100000.0,
            "theta_r": 0.2,
            "theta_s": 0.5,
            "alpha": 0.30,
            "beta": 12.0,
            "brook": 50.0,
        },
    },
    dtype=float,
)


def get_params(sm_name: Literal["Genuchten", "Brooks", "Panday"]) -> DataFrame:
    """Get the parameter bounds for a specific soil model."""
    params = {"Genuchten": pGenuchten, "Brooks": pBrooks, "Panday": pPanday}
    return params[sm_name].copy()
