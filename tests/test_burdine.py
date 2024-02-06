# %% [markdown]
# Soil Properties, Types and Fitting

# Since MODFLOW-USG uses both Brooks-Corey and van Genuchten. We need to
# transform known van Genuchten soil types to Brooks-Corey parameters

from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pedon as pe

def flip(items: list, ncol: int):
    return chain(*[items[i::ncol] for i in range(ncol)])

from pedon._params import pPanday
# %%
# Define Soil Names
soilnames = [
    'Sand',
    'Loamy Sand',
    'Sandy Loam',
    'Loam',
    'Silt',
    'Silt Loam',
    'Sandy Clay Loam',
    'Clay Loam',
    'Silty Clay Loam',
    'Sandy Clay',
    'Silty Clay',
    'Clay',
]
# %%
# fit brooks on van genuchten and plot
plotting = True

df = pd.DataFrame(
    index=pd.MultiIndex.from_product([soilnames, ["Genuchten-Mualem", "Genuchten-BC-Mualem", "Genuchten-BC-Burdine"]]), columns=np.unique(np.append(list(pe.Panday.__dict__["__annotations__"].keys()), list(pe.Genuchten.__dict__["__annotations__"].keys()))), dtype=float
)

for sn in soilnames:
    soilms = []
    soil = pe.Soil(sn).from_name(pe.Genuchten, source="HYDRUS")
    soilm = getattr(soil, "model")
    soilms.append(soilm)
    h = np.logspace(-4, 6, num=11)
    k = soilm.k(h)
    theta = soilm.theta(h)
    soil_usg = pe.SoilSample(h=h, k=k, theta=theta)
    soilm_usg = []
    for c in [1.0, 2.0]:
        pbounds = pd.DataFrame(
            data={
                "p_ini": {
                    "k_s": getattr(soilm, "k_s"),
                    # "theta_r": getattr(soilm, "theta_r"),
                    # "theta_s": getattr(soilm, "theta_s"),
                    "theta_r": pPanday.at["theta_r", "p_ini"],
                    "theta_s": pPanday.at["theta_s", "p_ini"],
                    "alpha": pPanday.at["alpha", "p_ini"],
                    "beta": pPanday.at["beta", "p_ini"],
                    "brook": 10.0,
                    "c": c,
                },
                "p_min": {
                    "k_s": getattr(soilm, "k_s") - 1e-10,
                    # "theta_r": getattr(soilm, "theta_r") - 1e-10,
                    # "theta_s": getattr(soilm, "theta_s") - 1e-10,
                    "theta_r": pPanday.at["theta_r", "p_min"],
                    "theta_s": pPanday.at["theta_s", "p_min"],
                    "alpha": pPanday.at["alpha", "p_min"],
                    "beta": pPanday.at["beta", "p_min"],
                    "brook": 1.0,
                    "c": c - 1e-10,
                },
                "p_max": {
                    "k_s": getattr(soilm, "k_s") + 1e-10,
                    # "theta_r": getattr(soilm, "theta_r") + 1e-10,
                    # "theta_s": getattr(soilm, "theta_s") + 1e-10,
                    "theta_r": pPanday.at["theta_r", "p_max"],
                    "theta_s": pPanday.at["theta_s", "p_max"],
                    "alpha": pPanday.at["alpha", "p_max"],
                    "beta": pPanday.at["beta", "p_max"],
                    "brook": 50.0,
                    "c": c + 1e-10,
                },
                "swrc": {
                    "k_s": False,
                    "theta_r": True,
                    "theta_s": True,
                    "alpha": True,
                    "beta": True,
                    "brook": False,
                    "c": True,
                },
            },
        )
        soilm_usg = soil_usg.fit(
            pe.Panday,
            pbounds=pbounds,
            k_s=soilm.k_s,
            W1=0.1
        )
        soilms.append(soilm_usg)

    for sm in soilms:
        for col in sm.__annotations__:
            if isinstance(sm, pe.Genuchten):
                mi = "Genuchten-Mualem"
            else:
                if np.isclose(sm.c, 1.0):
                    mi = "Genuchten-BC-Mualem"
                elif np.isclose(sm.c, 2.0):
                    mi = "Genuchten-BC-Burdine"

            df.loc[(sn, mi), col] = sm.__getattribute__(col)


    f, ax = plt.subplots(
        1, 2, sharey=True, figsize=(6, 4)
    )
    for sm, ls, lab in zip(soilms, ["-", "--", "-."], ["Genuchten-Mualem", "Genuchten-BC-Mualem", "Genuchten-BC-Burdine"]):
        _ = pe.soilmodel.plot_swrc(sm, ax=ax[0], linestyle=ls, label=lab)
        _ = pe.soilmodel.plot_hcf(sm, ax=ax[1], linestyle=ls, label=lab)

    ax[0].set_yscale("log")
    ax[1].set_yscale("log")

    ax[0].set_xlim(0, 0.5)
    ax[0].set_xticks(np.linspace(0, 0.5, 6))
    ax[0].set_yticks(h)

    ax[1].set_xscale("log")
    ax[1].set_xlim(1e-25, 5e3)

    ax[0].set_ylabel(r"|$\psi$| [cm]")
    ax[0].set_xlabel(r"$\theta$ [-]")
    ax[1].set_xlabel(r"$K_s$ [cm/d]")
    handles, labels = ax[0].get_legend_handles_labels()
    ncol = 3
    ax[0].legend(
        flip(handles, ncol),
        flip(labels, ncol),
        loc=(-0.02, 1),
        fontsize=7.5,
        frameon=False,
        ncol=ncol,
        columnspacing=0.8,
        handlelength=2.5,
    )
    f.savefig(
        f"Soil_Properties_{sn.replace(' ', '')}.png",
        bbox_inches="tight",
        dpi=300,
    )
    f.align_xlabels()
    plt.close(f)

df.loc[(slice(None), "Genuchten-Mualem"), "beta"] = df.loc[(slice(None), "Genuchten-Mualem"), "n"]
df.loc[(slice(None), "Genuchten-Mualem"), "gamma"] = df.loc[(slice(None), "Genuchten-Mualem"), "m"]
df = df.drop(columns=["n", "m", "c", "sr", "sy", "ss", "h_b"])
df.round(4).to_excel("Fitted_Parameters_Genuchten_BC_Mualem_Burdine.xlsx")
