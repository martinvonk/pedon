import logging
from dataclasses import MISSING, fields
from typing import Type

from numpy import inf, nan
from pandas import DataFrame

from pedon._typing import SoilModelNames
from pedon.soilmodel import SoilModel, get_soilmodel


def _get_default_params(sm: Type[SoilModel]) -> DataFrame:
    """Return an empty DataFrame with the same structure as parameter DataFrames."""
    index = [
        f.name
        for f in fields(sm)
        if f.init and f.default is MISSING and f.default_factory is MISSING
    ]
    df = DataFrame(
        data={"p_ini": nan, "p_min": -inf, "p_max": inf},
        index=index,
        columns=["p_ini", "p_min", "p_max"],
        dtype=float,
    )
    return df


def get_params(
    sm: Type[SoilModel] | SoilModel | SoilModelNames,
) -> DataFrame:
    """Get the parameter bounds for a specific soil model."""
    if isinstance(sm, type) and issubclass(sm, SoilModel):
        smn = sm.__name__
    elif isinstance(sm, SoilModel):
        smn = getattr(sm, "__name__", sm.__class__.__name__)
        sm = type(sm)
    elif isinstance(sm, str):
        smn = sm
        sm = get_soilmodel(smn)
    else:
        raise ValueError(
            f"Argument must either be Type[SoilModel] | SoilModel | str, not {type(sm)}"
        )

    params = _get_default_params(sm)

    param_bounds = {
        "Genuchten": {
            "k_s": [10.0, 0.001, 100000.0],
            "theta_r": [0.01, 0.0, 0.2],
            "theta_s": [0.40, 0.2, 0.9],
            "alpha": [0.02, 0.001, 0.20],
            "n": [2.0, 1.000001, 12.0],
            "l": [0.5, -7.0, 8.0],
        },
        "Brooks": {
            "k_s": [50.0, 0.001, 100000.0],
            "theta_r": [0.02, 0.0, 0.2],
            "theta_s": [0.4, 0.2, 0.5],
            "h_b": [0.003, 0.0001, 100.0],
            "l": [1.5, 0.1, 5.0],
        },
        "Panday": {
            "k_s": [50.0, 0.001, 100000.0],
            "theta_r": [0.02, 1e-5, 0.2],
            "theta_s": [0.4, 0.2, 0.5],
            "alpha": [0.02, 0.001, 0.30],
            "beta": [2.3, 1.0, 12.0],
            "brook": [10.0, 1.0, 50.0],
        },
    }

    if smn in param_bounds:
        for param, bounds in param_bounds[smn].items():
            params.loc[param] = bounds
    else:
        logging.warning(f"No default parameter bounds for SoilModel type {smn}")

    return params
