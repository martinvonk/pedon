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
            "l": -7,
        },
        "p_max": {
            "k_s": 100000.0,
            "theta_r": 0.2,
            "theta_s": 0.9,
            "alpha": 0.20,
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
pGenuchten.loc[:, ["p_ini", "p_min", "p_max"]] = pGenuchten.loc[
    :, ["p_ini", "p_min", "p_max"]
].astype(float)
pGenuchten.loc[:, "swrc"] = pGenuchten.loc[:, "swrc"].astype(bool)

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
pBrooks.loc[:, ["p_ini", "p_min", "p_max"]] = pBrooks.loc[
    :, ["p_ini", "p_min", "p_max"]
].astype(float)
pBrooks.loc[:, "swrc"] = pBrooks.loc[:, "swrc"].astype(bool)

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
            "beta": 12,
            "brook": 50.0,
        },
        "swrc": {
            "k_s": False,
            "theta_r": True,
            "theta_s": True,
            "alpha": True,
            "beta": True,
            "brook": False,
        },
    },
)
pPanday.loc[:, ["p_ini", "p_min", "p_max"]] = pPanday.loc[
    :, ["p_ini", "p_min", "p_max"]
].astype(float)
pPanday.loc[:, "swrc"] = pPanday.loc[:, "swrc"].astype(bool)


def get_params(sm_name: str) -> DataFrame:
    params = {"Genuchten": pGenuchten, "Brooks": pBrooks, "Panday": pPanday}
    return params[sm_name].copy()
