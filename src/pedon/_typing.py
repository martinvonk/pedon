from typing import Literal

from numpy import float64
from numpy.typing import NDArray

FloatArray = float | NDArray[float64]
SoilModelNames = Literal[
    "Genuchten",
    "Brooks",
    "Panday",
    "Gardner",
    "Rucker",
    "Haverkamp",
    "Fredlund",
    "Kosugi",
    "Campbell",
    "GenuchtenKool",
    "GenuchtenGardner",
]
