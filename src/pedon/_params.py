from pandas import DataFrame

pGenuchten = DataFrame(
    data={
        "p_i": {
            "k_s": 50.0,
            "theta_r": 0.02,
            "theta_s": 0.4,
            "alpha": 0.02,
            "n": 2.3,
            "l": 0.5,
        },
        "p_min": {
            "k_s": 0.001,
            "theta_r": 0.0,
            "theta_s": 0.2,
            "alpha": 0.001,
            "n": 1.0,
            "l": -7,
        },
        "p_max": {
            "k_s": 100000.0,
            "theta_r": 0.2,
            "theta_s": 0.9,
            "alpha": 0.15,
            "n": 12,
            "l": 8,
        },
        "swrc": {
            "k_s": False,
            "theta_r": True,
            "theta_s": True,
            "alpha": True,
            "n": True,
            "l": False,
        },
    },
)
pGenuchten.loc[:, ["p_i", "p_min", "p_max"]] = pGenuchten.loc[
    :, ["p_i", "p_min", "p_max"]
].astype(float)
pGenuchten.loc[:, "swrc"] = pGenuchten.loc[:, "swrc"].astype(bool)

pBrooks = DataFrame(
    data={
        "p_i": {"k_s": 50.0, "theta_r": 0.02, "theta_s": 0.4, "h_b": 0.003, "l": 1.5},
        "p_min": {
            "k_s": 0.001,
            "theta_r": 0.0,
            "theta_s": 0.2,
            "h_b": 0.0001,
            "l": 0.1,
        },
        "p_max": {"k_s": 100000.0, "theta_r": 0.2, "theta_s": 0.5, "h_b": 100, "l": 5},
        "swrc": {
            "k_s": False,
            "theta_r": True,
            "theta_s": True,
            "h_b": True,
            "l": True,
        },
    },
)
pBrooks.loc[:, ["p_i", "p_min", "p_max"]] = pBrooks.loc[
    :, ["p_i", "p_min", "p_max"]
].astype(float)
pBrooks.loc[:, "swrc"] = pBrooks.loc[:, "swrc"].astype(bool)


def get_params(sm_name: str) -> DataFrame:
    params = {"Genuchten": pGenuchten, "Brooks": pBrooks}
    return params[sm_name]
