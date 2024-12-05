# /// script
# dependencies = [
#   "requests",
#   "pandas",
#   "geopandas",
#   "matplotlib",
#   "openpyxl",
#   "pedon",
# ]
# ///
import requests
import pedon as pe
import pandas as pd
import geopandas as gpd
from shapely import Polygon


# %%
# voor 1 punt
# send get request
r = requests.get(
    url="https://www.soilphysics.wur.nl/soil.php",
    params={"latitude": 52, "longitude": 5.5},
)

# extract data in json format
data = r.json()

# further processing
pd.DataFrame(data["horizon"])

# %%
# voor bofek
# haal cluster en bijbehorende profielen op
data = pd.read_excel(
    "tabel_BOFEK_2020_V_1_0.xlsx", sheet_name="Data (2)", skiprows=9, index_col=0
)
profiel_cluster = dict(zip(data["Cluster"], data["Profiel"]))
# %%
# haal de geometrie van de clusters op
gdf = gpd.read_file("BOFEK2020.gdb")
# snij bij met extent van wijde wormer
extent = Polygon(
    [
        (117500.0, 497250.0),
        (125000.0, 497250.0),
        (125000.0, 503000.0),
        (117500.0, 503000.0),
    ]
)  # wijde wormer extent
gdf_clip = gpd.clip(gdf, extent).dropna(subset=["BOFEK2020"])

# %%
# pak unieke profielen op basis van clusters
unique_profiles = gdf_clip["BOFEK2020"].astype(int).map(profiel_cluster).unique()
# %%
# haal de profielen op
profiles = pd.read_csv(
    "AllProfiles_368.csv",
    index_col=0,
    na_values=99999,
    dtype={
        "iSoil1": "int64",
        "iSoil2": "int64",
        "iSoil3": "int64",
        "iSoil4": "int64",
        "iSoil5": "int64",
        "iSoil6": "int64",
        "iSoil7": "int64",
        "iSoil8": "int64",
        "iSoil9": "int64",
        "iZ1": "float64",
        "iZ2": "float64",
        "iZ3": "float64",
        "iZ4": "float64",
        "iZ5": "float64",
        "iZ6": "float64",
        "iZ7": "float64",
        "iZ8": "float64",
        "iZ9": "float64",
    },
).drop(columns=["iBOFEK2012", "iHoofd"])
# %%
# haal de bodemparameters op (kan netter)
soil_mapper = {
    1: pe.Soil("B01").from_staring("2018").model,
    2: pe.Soil("B02").from_staring("2018").model,
    3: pe.Soil("B03").from_staring("2018").model,
    4: pe.Soil("B04").from_staring("2018").model,
    5: pe.Soil("B05").from_staring("2018").model,
    6: pe.Soil("B06").from_staring("2018").model,
    7: pe.Soil("B07").from_staring("2018").model,
    8: pe.Soil("B08").from_staring("2018").model,
    9: pe.Soil("B09").from_staring("2018").model,
    10: pe.Soil("B10").from_staring("2018").model,
    11: pe.Soil("B11").from_staring("2018").model,
    12: pe.Soil("B12").from_staring("2018").model,
    13: pe.Soil("B13").from_staring("2018").model,
    14: pe.Soil("B14").from_staring("2018").model,
    15: pe.Soil("B15").from_staring("2018").model,
    16: pe.Soil("B16").from_staring("2018").model,
    17: pe.Soil("B17").from_staring("2018").model,
    18: pe.Soil("B18").from_staring("2018").model,
    19: pe.Soil("O01").from_staring("2018").model,
    20: pe.Soil("O02").from_staring("2018").model,
    21: pe.Soil("O03").from_staring("2018").model,
    22: pe.Soil("O04").from_staring("2018").model,
    23: pe.Soil("O05").from_staring("2018").model,
    24: pe.Soil("O06").from_staring("2018").model,
    25: pe.Soil("O07").from_staring("2018").model,
    26: pe.Soil("O08").from_staring("2018").model,
    27: pe.Soil("O09").from_staring("2018").model,
    28: pe.Soil("O10").from_staring("2018").model,
    29: pe.Soil("O11").from_staring("2018").model,
    30: pe.Soil("O12").from_staring("2018").model,
    31: pe.Soil("O13").from_staring("2018").model,
    32: pe.Soil("O14").from_staring("2018").model,
    33: pe.Soil("O15").from_staring("2018").model,
    34: pe.Soil("O16").from_staring("2018").model,
    35: pe.Soil("O17").from_staring("2018").model,
    36: pe.Soil("O18").from_staring("2018").model,
}

# %%
# zet Bofek profiel om te zetten naar boring met bodemparameters


def bofek_to_boring(profiles: pd.DataFrame, iprofile: int):
    profile = profiles.loc[[iprofile]]
    depth = (
        profile.loc[:, ["iZ1", "iZ2", "iZ3", "iZ4", "iZ5", "iZ6", "iZ7", "iZ8", "iZ9"]]
        .squeeze()
        .rename("Depth [cm]")
    )
    soil = (
        profile.loc[
            :,
            [
                "iSoil1",
                "iSoil2",
                "iSoil3",
                "iSoil4",
                "iSoil5",
                "iSoil6",
                "iSoil7",
                "iSoil8",
                "iSoil9",
            ],
        ]
        .squeeze()
        .rename("Soil")
    )
    soil.index = depth
    # soil = soil[~soil.duplicated(keep="last")] # remove duplicates maar dit klopt zo nog niet als je meerdere trajecten hebt in bodem met dezelfde bodem
    return soil.map(soil_mapper).dropna().rename(iprofile)


{iprofile: bofek_to_boring(profiles, iprofile) for iprofile in unique_profiles}
# %%
import pedon as pe

import numpy as np


gs = [
    pe.Genuchten(0.1, 0.05, 0.4, 0.3, 1.2),
    pe.Genuchten(0.1, 0.05, 0.4, 0.3, 1.2),
    pe.Genuchten(0.1, 0.05, 0.4, 0.3, 1.2),
]

depth = np.array([0.0, 20.0, 40.0, 70.0, 120.0])
ghg_depth = 50.0
h = depth - ghg_depth  # pressure head
# pressure head = head - z = h + depths
dh = 1.0
for hs, he, g in zip(h[:-1], h[1:], gs):
    hr = np.arange(hs, he, dh)
    print(hr)
    theta = g.theta(hr)


# %%
def get_brooks(gen: pe.Genuchten, h: np.ndarray[float] | None) -> pe.Brooks:
    if h is not None:
        k = gen.k(h)
        theta = gen.theta(h)
        soilsample = pe.SoilSample(h=h, k=k, theta=theta)
        pbounds = pe._params.pBrooks.copy()
        pbounds.loc["theta_r"] = (
            gen.theta_r,
            max(g.theta_r - 0.1, 0.0),
            gen.theta_r + 0.1,
        )
        pbounds.loc["theta_s"] = (gen.theta_s, gen.theta_s - 0.1, gen.theta_s + 0.1)
        pbounds.loc["l"] = (5.0, 0.01, 20.0)
        pbounds.loc["h_b"] = (1.0, 0.1, 100.0)
        bc = soilsample.fit(pe.Brooks, pbounds=pbounds, k_s=g.k_s)
    else:
        # Morel-Seytoux (1996) - Parameter equivalence for the Brooks-Corey and van Genuchten soil characteristics
        eps = 1 + 2 / gen.m  # eq 16b
        h_b = (
            1
            / gen.alpha
            * (eps + 3)
            / (2 * eps * (eps - 1))
            * (147.8 + 8.1 * eps + 0.0928 * eps**2)
            / (55.6 + 7.4 * eps + eps**2)
        )  # eq 17
        l = 2 / (eps - 3)  # because eps = 3 + 2 / l
        bc = pe.Brooks(
            k_s=gen.k_s, theta_r=gen.theta_r, theta_s=gen.theta_s, h_b=h_b, l=l
        )
    return bc


# %%

# name = "B05"
# g = pe.Soil(name=name).from_staring("2018").model
hb = []
for name, g in soil_mapper.items():
    h = np.logspace(-2, 6, num=11)
    bc1 = get_brooks(g, h)
    bc2 = get_brooks(g, None)

    print(name)
    print(3 + 2 / bc1.l, 3 + 2 / bc2.l)
    print(bc1.h_b, bc2.h_b)

    ax1 = pe.plot_swrc(g)
    pe.plot_swrc(bc1, ax=ax1)
    pe.plot_swrc(bc2, ax=ax1)
    ax1.set_title(name)

    ax2 = pe.plot_hcf(g)
    pe.plot_hcf(bc1, ax=ax2)
    pe.plot_hcf(bc2, ax=ax2)
    ax2.set_title(name)
    hb.append(bc1.h_b)
    # break

# %%
