from dataclasses import dataclass
from pathlib import Path

from pandas import read_csv

from ._typing import SoilModel
from .soilmodel import Brooks, Fredlund, Gardner, Genuchten, Sorab


@dataclass
class Soil:
    name: str
    type: str | None = None
    model: SoilModel | None = None
    source: str | None = None
    description: str | None = None


def get_soil(name: str, soilmodelname: str = "Genuchten") -> Soil:
    sms = {
        "Genuchten": Genuchten,
        "Brooks": Brooks,
        "Gardner": Gardner,
        "Sorab": Sorab,
        "Fredlund": Fredlund,
    }
    path = Path(f"src/pedon/datasets/{soilmodelname}.csv")
    ser = read_csv(path, index_col=["name"]).loc[name].to_dict()
    soil = Soil(
        name=name,
        type=ser.pop("soil type"),
        source=ser.pop("source"),
        description=ser.pop("description"),
    )
    soil.__setattr__("model", sms[soilmodelname](**ser))
    return soil
