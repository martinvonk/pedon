from typing import Literal
from numpy import nan, inf
from pandas import DataFrame
from dataclasses import fields
from pedon.soilmodel import Genuchten, Brooks, Panday
import logging

def get_default_params(model: pe.soilmodel.SoilModel) -> DataFrame:
    """Return an empty DataFrame with the same structure as parameter DataFrames."""
    index = [f.name for f in fields(model) if f.init]
    df = DataFrame(
        data={"p_ini": nan, "p_min": -inf, "p_max": inf},
        index=index,
        columns=["p_ini", "p_min", "p_max"],
        dtype=float,
    )
    return df

def get_default_bounds(model: pe.soilmodel.SoilModel) -> DataFrame:
    """Get the parameter bounds for a specific soil model."""
    params = get_default_params(type(model))
    if isinstance(model, Genuchten):
        params.loc["k_s"] = [10.0, 0.001, 100000.0]
        params.loc["theta_r"] = [0.01, 0.0, 0.2]
        params.loc["theta_s"] = [0.40, 0.2, 0.9]
        params.loc["alpha"] = [0.02, 0.001, 0.20]
        params.loc["n"] = [2.0, 1.000001, 12.0]
        params.loc["l"] = [0.5, -7.0, 8.0]
    elif isinstance(model, Brooks):
        params.loc["k_s"] = [50.0, 0.001, 100000.0]
        params.loc["theta_r"] = [0.02, 0.0, 0.2]
        params.loc["theta_s"] = [0.4, 0.2, 0.5]
        params.loc["h_b"] = [0.003, 0.0001, 100.0]
        params.loc["l"] = [1.5, 0.1, 5.0]
    elif isinstance(model, Panday):
        params.loc["k_s"] = [50.0, 0.001, 100000.0]
        params.loc["theta_r"] = [0.02, 1e-5, 0.2]
        params.loc["theta_s"] = [0.4, 0.2, 0.5]
        params.loc["alpha"] = [0.02, 0.001, 0.30]
        params.loc["beta"] = [2.3, 1.0, 12.0]
        params.loc["brook"] = [10.0, 1.0, 50.0]
    else:
        logging.warning(f"No default parameter bounds for model type {type(model)}")
    return params
