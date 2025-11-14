from typing import Literal

from numpy import float64
from numpy.typing import NDArray

FloatArray = NDArray[float64]
SoilModelNames = Literal[
    "Genuchten",
    "Brooks",
    "Panday",
    "Gardner",
    "Haverkamp",
    "Fredlund",
]
