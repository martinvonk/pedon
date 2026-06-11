from typing import Literal

from matplotlib.axes import Axes
from numpy import float64
from numpy.typing import NDArray

MatplotlibAxes = Axes
FloatArray = float | NDArray[float64]
SoilModelNames = Literal[
    "Genuchten",
    "Brooks",
    "Brunswick",
    "Panday",
    "Gardner",
    "Rucker",
    "Haverkamp",
    "Fredlund",
    "Kosugi",
    "Campbell",
    "Kool",
    "GenuchtenGardner",
    "Gerke",
]

SourceNames = Literal[
    "HYDRUS", "VS2D", "Staring_2001", "Staring_2018", "Clapp", "Rawls"
]
