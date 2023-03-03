from dataclasses import dataclass


@dataclass
class SoilSample:
    sand_p: float  # %
    silt_p: float  # %
    clay_p: float  # %
    rho: float  # g/cm3
    th33: float  # cm
    th1500: float  # cm
