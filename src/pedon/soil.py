from dataclasses import dataclass
from pathlib import Path
from typing import Type
from pandas import read_csv, read_excel

from .soilmodel import Genuchten, SoilModel
from .pedotransfer import SoilSample


@dataclass
class Soil:
    name: str
    type: str | None = None
    model: SoilModel | None = None
    sample: SoilSample | None = None
    source: str | None = None
    description: str | None = None

    def from_name(self, sm: Type[SoilModel]) -> "Soil":
        path = Path(__file__).parent / f"datasets/{sm.__name__}.csv"
        ser = read_csv(path, index_col=["name"]).loc[self.name].to_dict()
        self.__setattr__("type", ser.pop("soil type"))
        self.__setattr__("source", ser.pop("source"))
        self.__setattr__("description", ser.pop("description"))
        self.__setattr__("model", sm(**ser))
        return self

    def from_staring(self, year: str = "2018") -> "Soil":
        path = Path(__file__).parent / f"datasets/Staring_{year}.xlsx"
        parameters = read_excel(path, sheet_name="parameters", index_col=0)
        ser = parameters.loc[self.name].to_dict()
        self.__setattr__("type", ser.pop("soil type"))
        self.__setattr__("source", ser.pop("source"))
        self.__setattr__("description", ser.pop("description"))
        sm = Genuchten(**ser)
        self.__setattr__("model", sm)
        ss = SoilSample().from_staring(name=self.name, year=year)
        self.__setattr__("sample", ss)
        return self
