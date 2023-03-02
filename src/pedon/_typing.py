from numpy.typing import NDArray
from numpy import float64
from typing import Protocol

FloatArray = float | NDArray[float64]


class SoilModel(Protocol):
    def theta(self, h: FloatArray) -> FloatArray:
        """Method to calculate the soil moisture content from the pressure head h"""
        ...

    def s(self, h: FloatArray) -> FloatArray:
        """Method to calculate the effective saturation from the pressure head h"""
        ...

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        """Method to calcualte the permeability from the pressure head h"""
        ...