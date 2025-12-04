from typing import Literal

from numpy import float64
from numpy.typing import NDArray

FloatArray = NDArray[float64]
SoilModelNames = Literal[
    "Genuchten",
    "Brooks",
    "Panday",
    "Gardner",
    "Mod_Gardner",
    "Haverkamp",
    "Fredlund",
]
