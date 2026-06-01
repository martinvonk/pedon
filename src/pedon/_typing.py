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

SourceNames = Literal[
    "HYDRUS",
    "VS2D",
    "Staring_2001",
    "Staring_2018",
]
