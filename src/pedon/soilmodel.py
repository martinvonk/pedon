from dataclasses import dataclass
from numpy import abs as npabs
from numpy import full, log, exp, linspace
from .plot import swrc
from ._typing import FloatArray


@dataclass
class Genuchten:
    k_s: float
    theta_r: float
    theta_s: float
    alpha: float
    n: float
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

    def plot(self):
        return swrc(self)


@dataclass
class Brooks:
    k_s: float
    theta_r: float
    theta_s: float
    h_b: float
    l: float

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

    def plot(self):
        return swrc(self)


@dataclass
class Gardner:
    k_s: float
    theta_r: float
    theta_s: float
    a: float
    b: float
    m: float

    def theta(self, h: FloatArray) -> FloatArray:
        return self.theta_r + (self.theta_s - self.theta_r) * self.s(h)

    def s(self, h: FloatArray) -> FloatArray:
        return 1 / (1 + npabs(h / self.a) ** self.b)

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is not None:
            raise NotImplementedError(
                "Can only calculate the hydraulic conductivity using the pressure head, not the saturation"
            )
        return self.k_s * exp(-self.a * h)

    def plot(self):
        return swrc(self)


@dataclass
class Sorab:
    k_s: float
    sr: float  # theta_r / theta_s
    alpha: float  # alpha
    beta: float  # n
    brook: float  # brooks-corey l
    theta_s: float | None = None

    def __post_init__(self):
        self.gamma = 1 - 1 / self.beta  # m

    def theta(self, h: FloatArray) -> FloatArray:
        if self.theta_s is not None:
            return (self.sr + self.s(h) * (1 - self.sr)) * self.theta_s
        raise ValueError("theta_s must not be none")

    def s(self, h: FloatArray) -> FloatArray:
        return (1 + self.alpha * npabs(h) ** self.beta) ** -self.gamma

    def k(self, h: FloatArray, s: FloatArray | None = None) -> FloatArray:
        if s is None:
            s = self.s(h)
        return self.k_s * s**self.brook

    def plot(self):
        return swrc(self)


@dataclass
class Fredlund:
    k_s: float
    theta_s: float
    a: float
    n: float
    m: float

    def theta(self, h: FloatArray) -> FloatArray:
        return self.theta_s / (log(exp(1) + npabs(h / self.a) ** self.n)) ** self.m

    def s(self, h: FloatArray) -> FloatArray:
        return self.theta(h) / self.theta_s

    def k(self, h: FloatArray, s: FloatArray | None = None):
        if s is not None:
            raise NotImplementedError(
                "Can only calculate the hydraulic conductivity using the pressure head, not the saturation"
            )

        def theta_d(
            h: FloatArray, a: float, n: float, m: float, theta_s: float
        ) -> FloatArray:
            """Derivative of theta (according to WolframAlpha)"""
            return -(
                theta_s
                * m
                * n
                * (h / a) ** n
                * (log((h / a) ** n + exp(1))) ** (-m - 1)
            ) / (h * (h / a) ** n + exp(1) * h)

        h_b = 0.03
        n_steps = 100

        teller = 0
        for y in linspace(log(h), log(1e6), n_steps):
            teller += (
                (self.theta(exp(y)) - self.theta(h))
                / exp(y)
                * theta_d(exp(y), self.a, self.n, self.m, self.theta_s)
            )

        noemer = 0
        for y in linspace(log(h_b), log(1e6), n_steps):
            noemer += (
                (self.theta(exp(y)) - self.theta_s)
                / exp(y)
                * theta_d(exp(y), self.a, self.n, self.m, self.theta_s)
            )

        return self.k_s * (teller / noemer)

    def plot(self):
        return swrc(self)
