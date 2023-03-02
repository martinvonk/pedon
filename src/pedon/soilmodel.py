from dataclasses import dataclass
from numpy import abs as npabs
from numpy import full, log, e
from traitlets import Float

from ._typing import FloatArray


@dataclass
class Genuchten:
    theta_r: float
    theta_s: float
    alpha: float
    n: float
    k_s: float
    l: float = 0.5

    def __post_init__(self):
        self.m = 1 - 1 / self.n

    def theta(self, h: FloatArray) -> FloatArray:
        theta = (
            self.theta_r
            + (self.theta_s - self.theta_r)
            / (1 + npabs(self.alpha * h) ** self.n) ** self.m
        )
        return theta

    def s(self, h: FloatArray) -> FloatArray:
        return (self.theta(h) - self.theta_r) / (self.theta_s - self.theta_r)

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return self.k_s * s**self.l * (1 - (1 - s ** (1 / self.m)) ** self.m) ** 2


@dataclass
class Brooks:
    theta_r: float
    theta_s: float
    h_b: float
    l: float
    k_s: float

    def theta(self, h: FloatArray) -> FloatArray:
        h = npabs(h)
        if isinstance(h, float):
            if h >= self.h_b:
                return self.theta_r + self.s(h) * (self.theta_s - self.theta_r)
            else:
                return self.theta_s
        else:
            theta = full(h.shape, self.theta_s)
            theta[h >= self.h_b] = self.theta_r + self.s(h[h > self.h_b]) * (
                self.theta_s - self.theta_r
            )
            return theta

    def s(self, h: FloatArray) -> FloatArray:
        h = npabs(h)
        if isinstance(h, float):
            if h >= self.h_b:
                return (h / self.h_b) ** -self.l
            else:
                return 1
        else:
            s = full(h.shape, 1.0)
            s[h >= self.h_b] = (h[h >= self.h_b] / self.h_b) ** -self.l
            return s

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return self.k_s * s ** (3 + 2 / self.l)


@dataclass
class Gardner:
    theta_r: float
    theta_s: float
    a: float
    b: float
    m: float
    k_s: float

    def theta(self, h: FloatArray) -> FloatArray:
        return self.theta_r + (self.theta_s - self.theta_r) * self.s(h)

    def s(self, h: FloatArray) -> FloatArray:
        return 1 / (1 + npabs(h / self.a) ** self.b)

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            return self.a / self.b + npabs(h) ** self.m
        else:
            return self.a * self.theta_r + s * (self.theta_s - self.theta_r) ** self.m


@dataclass
class Sorab:
    sr: float  # theta_r / theta_s
    alpha: float  # alpha
    beta: float  # n
    k_s: float
    brook: float  # brooks-corey l

    def __post_init__(self):
        self.gamma = 1 - 1 / self.beta  # m

    def theta(self, h: FloatArray, theta_s: FloatArray) -> FloatArray:
        return (self.sr + self.s(h) * (1 - self.sr)) * theta_s

    def s(self, h: FloatArray) -> FloatArray:
        return (1 + self.alpha * h**self.beta) ** -self.gamma

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return self.k_s * s**self.brook


@dataclass
class Fredlund:
    theta_s: float
    a: float
    n: float
    m: float
    theta_r: float = 0

    def theta(self, h: FloatArray) -> FloatArray:
        return (
            self.theta_r
            + (self.theta_s - self.theta_r)
            / (log(e + npabs(h / self.a) ** self.n)) ** self.m
        )

    def s(self, h: FloatArray) -> FloatArray:
        return (self.theta(h) - self.theta_r) / (self.theta_s - self.theta_r)

    def k(self, h: FloatArray, s: FloatArray | None = None):
        raise NotImplementedError("See Fredlund & Xing 1994")
