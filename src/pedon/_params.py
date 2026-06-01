from dataclasses import MISSING, fields
from logging import getLogger
from typing import Any, Type, cast

from numpy import inf, nan
from pandas import DataFrame

from pedon._typing import SoilModelNames
from pedon.soilmodel import SoilModel, resolve_soilmodel

logger = getLogger(__name__)


def _get_default_params(sm: Type[SoilModel]) -> DataFrame:
    """Return an empty DataFrame with the same structure as parameter DataFrames."""
    index = [
        f.name
        for f in fields(cast(Any, sm))
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
    smn, sm_cls = resolve_soilmodel(sm)

    params = _get_default_params(sm_cls)

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
        logger.warning(f"No default parameter bounds for SoilModel type {smn}")

    return params
