from typing import Literal

from matplotlib.axes import Axes
from numpy import float64
from numpy.typing import NDArray

MatplotlibAxes = Axes
FloatArray = float | NDArray[float64]
SoilModelNames = Literal[
    "Genuchten",
    "GenuchtenGardner",
    "Brooks",
    "Panday",
    "Gardner",
    "Rucker",
    "Haverkamp",
    "Fredlund",
]
