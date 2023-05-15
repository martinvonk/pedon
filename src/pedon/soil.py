# type: ignore
from dataclasses import dataclass, field
from pathlib import Path
from typing import Type

from numpy import append, exp, log, ones
from pandas import DataFrame, read_excel
from scipy.optimize import least_squares

from ._params import get_params
from ._typing import FloatArray
from .soilmodel import Brooks, Genuchten, SoilModel, get_soilmodel


@dataclass
class SoilSample:
    sand_p: float | None = None  # sand %
    silt_p: float | None = None  # silt %
    clay_p: float | None = None  # clay %
    rho: float | None = None  # soil density g/cm3
    th33: float | None = None  # cm
    th1500: float | None = None  # cm
    om_p: float | None = None  # organic matter %
    m50: float | None = None  # median sand fraction
    h: FloatArray | None = field(default=None, repr=False)  # pressure head measurement
    k: FloatArray | None = field(
        default=None, repr=False
    )  # hydraulic conductivity measurement
    theta: FloatArray | None = field(
        default=None, repr=False
    )  # moisture content measurement

    def from_staring(self, name: str, year: str = "2018") -> "SoilSample":
        """Get properties and measurements from Staring series"""
        if year not in ("2001", "2018"):
            raise ValueError(
                f"No Staring series available for year '{year}'"
                "please use either '2001' or '2018'"
            )
        path = Path(__file__).parent / f"datasets/Staring_{year}.xlsx"
        properties = read_excel(path, sheet_name="properties", index_col=0)
        measurements = read_excel(path, sheet_name="measurements", index_col=[0, 1])
        self.h = measurements.columns.astype(float).values
        self.k = measurements.loc[name, "k"].astype(float).values
        self.theta = measurements.loc[name, "theta"].astype(float).values

        self.silt_p = properties.loc[name, "silt_p"]
        self.clay_p = properties.loc[name, "clay_p"]
        self.om_p = properties.loc[name, "om_p"]
        self.m50 = properties.loc[name, "m50"]
        if year == "2001":
            self.rho = properties.loc[name, "rho"]
        return self

    def fit_seperate(
        self,
        sm: Type[SoilModel],
        pbounds: DataFrame | None = None,
        weights: FloatArray | None = None,
        return_res: bool = False,
    ) -> SoilModel:
        """Fit the soil water retention and conductivity seperate."""
        if pbounds is None:
            pbounds = get_params(sm.__name__)
            pbounds.loc["k_s", "p_ini"] = max(self.k)
            pbounds.loc["theta_s", "p_ini"] = max(self.theta)
            pbounds.loc["theta_s", "p_max"] = max(self.theta) + 0.01

        if weights is None:
            weights = ones(self.h.shape)

        sml = sm(**dict(zip(pbounds.index, pbounds.loc[:, "p_ini"])))

        def fit_swrc(p: FloatArray) -> FloatArray:
            for pname, pv in zip(pbounds.index[pbounds.loc[:, "swrc"]], p):
                sml.__setattr__(pname, pv)
                diff = weights * (sml.theta(h=self.h) - self.theta)
            return diff

        def fit_k(p: FloatArray) -> FloatArray:
            for pname, pv in zip(pbounds.index[~pbounds.loc[:, "swrc"]], p):
                sml.__setattr__(pname, pv)
                diff = weights * (log(sml.k(h=self.h)) - log(self.k))
            return diff

        res_swrc = least_squares(
            fit_swrc,
            x0=pbounds.loc[pbounds.swrc, "p_ini"],
            bounds=(
                pbounds.loc[pbounds.swrc, "p_min"],
                pbounds.loc[pbounds.swrc, "p_max"],
            ),
        )

        res_k = least_squares(
            fit_k,
            x0=pbounds.loc[~pbounds.swrc, "p_ini"],
            bounds=(
                pbounds.loc[~pbounds.swrc, "p_min"],
                pbounds.loc[~pbounds.swrc, "p_max"],
            ),
        )

        opt_pars = dict(zip(pbounds.index[pbounds.loc[:, "swrc"]], res_swrc.x))
        opt_pars.update(dict(zip(pbounds.index[~pbounds.loc[:, "swrc"]], res_k.x)))
        opt_sm = sm(**opt_pars)
        if return_res:
            return opt_sm, {"res_swrc": res_swrc, "res_k": res_k}
        return opt_sm

    def fit(
        self,
        sm: Type[SoilModel],
        pbounds: DataFrame | None = None,
        weights: FloatArray | None = None,
        W1: float | None = None,
        W2: float | None = None,
        return_res: bool = False,
        k_s: float | None = None,
    ) -> SoilModel:
        """Same method as RETC"""

        theta = self.theta
        k = self.k

        if pbounds is None:
            pbounds = get_params(sm.__name__)
            if k_s is not None:
                pbounds = pbounds.drop("k_s")
            else:
                pbounds.loc["k_s", "p_ini"] = max(k)
                pbounds.loc["k_s", "p_max"] = max(k) * 10
            pbounds.loc["theta_s", "p_ini"] = max(theta)
            pbounds.loc["theta_s", "p_max"] = max(theta) + 0.02

        if weights is None:
            weights = ones(self.h.shape)

        if W1 is None:
            W1 = 0.1

        if W2 is None:
            M = len(k) + len(theta)
            N = len(theta)
            W2 = (M - N) * sum(weights * theta) / (N * sum(weights * k))

        def fit_staring(p: FloatArray) -> FloatArray:
            est_pars = dict(zip(pbounds.index, p))
            if k_s is not None:
                est_pars["k_s"] = k_s
            sml = sm(**est_pars)
            theta_diff = sml.theta(h=self.h) - theta
            k_diff = log(sml.k(h=self.h)) - log(k)
            diff = append(weights * theta_diff, weights * W1 * W2 * k_diff)
            return diff

        res = least_squares(
            fit_staring,
            x0=pbounds.loc[:, "p_ini"],
            bounds=(
                pbounds.loc[:, "p_min"],
                pbounds.loc[:, "p_max"],
            ),
        )
        opt_pars = dict(zip(pbounds.index, res.x))
        if k_s is not None:
            opt_pars["k_s"] = k_s
        opt_sm = sm(**opt_pars)
        if return_res:
            return opt_sm, {"res": res}
        return opt_sm

    def wosten(self, ts: bool = False) -> Genuchten:
        """Wosten et al (1999) - Development and use of a database of hydraulic
        properties of European soils"""

        topsoil = 1.0 * ts

        theta_s = (
            0.7919
            + 0.001691 * self.clay_p
            - 0.29619 * self.rho
            - 0.000001419 * self.silt_p**2
            + 0.0000821 * self.om_p**2
            + 0.02427 * self.clay_p**-1
            + 0.01113 * self.silt_p**-1
            + 0.01472 * log(self.silt_p)
            - 0.0000733 * self.om_p * self.clay_p
            - 0.000619 * self.rho * self.clay_p
            - 0.001183 * self.rho * self.om_p
            - 0.0001664 * topsoil * self.silt_p
        )
        alpha_ = (
            -14.96
            + 0.03135 * self.clay_p
            + 0.0351 * self.silt_p
            + 0.646 * self.om_p
            + 15.29 * self.rho
            - 0.192 * topsoil
            - 4.671 * self.rho**2
            - 0.000781 * self.clay_p**2
            - 0.00687 * self.om_p**2
            + 0.0449 * self.om_p**-1
            + 0.0663 * log(self.silt_p)
            + 0.1482 * log(self.om_p)
            - 0.04546 * self.rho * self.silt_p
            - 0.4852 * self.rho * self.om_p
            + 0.00673 * topsoil * self.clay_p
        )
        n_ = (
            -25.23
            - 0.02195 * self.clay_p
            + 0.0074 * self.silt_p
            - 0.1940 * self.om_p
            + 45.5 * self.rho
            - 7.24 * self.rho**2
            - 0.0003658 * self.clay_p**2
            + 0.002885 * self.om_p**2
            - 12.81 * self.rho**-1
            - 0.1524 * self.silt_p**-1
            - 0.01958 * self.om_p**-1
            - 0.2876 * log(self.silt_p)
            - 0.0709 * log(self.om_p)
            - 44.6 * log(self.rho)
            - 0.02264 * self.rho * self.clay_p
            + 0.0896 * self.rho * self.om_p
            + 0.00718 * topsoil * self.clay_p
        )
        l_ = (
            0.0202
            + 0.0006193 * self.clay_p**2
            - 0.001136 * self.om_p**2
            - 0.2316 * log(self.om_p)
            - 0.03544 * self.rho * self.clay_p
            + 0.00283 * self.rho * self.silt_p
            + 0.0488 * self.rho * self.om_p
        )
        ks_ = (
            7.755
            + 0.0352 * self.silt_p
            + 0.93 * topsoil
            - 0.967 * self.rho**2
            - 0.000484 * self.clay_p**2
            - 0.000322 * self.silt_p**2
            + 0.001 * self.silt_p**-1
            - 0.0748 * self.om_p**-1
            - 0.643 * log(self.silt_p)
            - 0.01398 * self.rho * self.clay_p
            - 0.1673 * self.rho * self.om_p
            + 0.02986 * topsoil * self.clay_p
            - 0.03305 * topsoil * self.silt_p
        )
        theta_r = 0.01
        return Genuchten(
            k_s=exp(ks_),
            theta_r=theta_r,
            theta_s=theta_s,
            alpha=exp(alpha_),
            n=exp(n_),
            l=exp(l_),
        )

    def wosten_sand(self, ts: bool = False) -> Genuchten:
        """Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
        van boven- en ondergronden in Nederland: de Staringreeks"""

        topsoil = 1.0 * ts
        theta_s = (
            -35.7
            - 0.1843 * self.rho
            - 0.03576 * self.m50
            + 0.0000261 * self.m50**2
            - 0.0564 * self.silt_p**-1
            + 0.008 * self.om_p**-1
            + 496.0 * self.m50**-1
            + 0.02244 * log(self.om_p)
            + 7.56 * log(self.m50)
        )

        k_s = exp(
            45.8
            - 14.34 * self.rho
            + 0.001481 * self.silt_p**2
            - 27.5 * self.rho**-1
            - 0.891 * log(self.silt_p)
            - 0.34 * log(self.om_p)
        )

        alpha = exp(
            13.66
            - 5.91 * self.rho
            - 0.172 * topsoil
            + 0.003248 * self.m50
            - 11.89 * self.rho**-1
            - 2.121 * self.silt_p**-1
            - 0.3742 * log(self.silt_p)
        )

        n = max(
            exp(
                -1.057
                + 0.1003 * self.om_p
                + 1.119 * self.rho
                + 0.000764 * self.silt_p**2
                - 0.1397 * self.om_p**-1
                - 57.2 * self.m50**-1
                - 0.557 * log(self.om_p)
                - 0.02997 * self.rho * self.silt_p
            )
            + 1,
            1.0,
        )

        l_ = (
            -76.4
            - 0.097 * self.silt_p
            + 59.6 * self.rho
            + 0.0332 * self.m50
            - 13.45 * self.rho**2
            + 0.001127 * self.silt_p**2
            + 0.00431 * self.om_p**2
            - 0.0000399 * self.m50**2
            + 40.8 * self.rho**-1
            + 2.364 * self.silt_p**-1
            + 1.014 * log(self.silt_p)
        )
        l = min(max(2 * (exp(l_) - 1) / (exp(l_) + 1), -2.0), 2.0)  # noqa: E741

        return Genuchten(
            k_s=round(k_s, 4),
            theta_r=0.01,
            theta_s=round(theta_s, 4),
            alpha=round(alpha, 7),
            n=round(n, 4),
            l=round(l, 4),
        )

    def wosten_clay(self) -> Brooks:
        """Wosten et al. (2001) - Waterretentie- en doorlatendheidskarakteristieken
        van boven- en ondergronden in Nederland: de Staringreeks"""

        theta_s = (
            0.6311
            + 0.003383 * self.clay_p
            - 0.09699 * self.rho**2
            - 0.00204 * self.rho * self.clay_p
        )

        k_s = exp(
            -42.6
            + 8.71 * self.om_p
            + 61.9 * self.rho
            - 20.79 * self.rho**2
            - 0.2107 * self.om_p**2
            - 0.01622 * self.clay_p * self.om_p
            - 5.382 * self.rho * self.om_p
        )

        alpha = exp(
            -19.13
            + 0.812 * self.om_p
            + 23.4 * self.rho
            - 8.16 * self.rho**2
            + 0.423 * self.om_p**-1
            + 2.388 * log(self.om_p)
            - 1.338 * self.rho * self.om_p
        )

        n = max(
            exp(
                -0.235
                + 0.972 * self.rho**-1
                - 0.7743 * log(self.clay_p)
                - 0.3154 * log(self.om_p)
                + 0.0678 * self.rho * self.om_p
            )
            + 1.0,
            1.0,
        )

        l_ = 0.102 + 0.0222 * self.clay_p - 0.043 * self.rho * self.clay_p
        l = min(max(10.0 * (exp(l_) - 1) / (exp(l_) + 1), -10.0), 10.0)  # noqa: E741
        return Genuchten(
            k_s=round(k_s, 4),
            theta_r=0.01,
            theta_s=round(theta_s, 4),
            alpha=round(alpha, 7),
            n=round(n, 4),
            l=round(l, 4),
        )

    def cosby(self) -> Brooks:
        """Cooper (2021) - Using data assimilation to optimize pedotransfer
        functions using field-scale in situ soil moisture observations"""
        c = self.clay_p / 100
        s = self.sand_p / 100
        b = 3.100 + 15.70 * c - 0.300 * s
        theta_s = 0.505 - 0.037 * c - 0.142 * s
        psi_s = 0.01 * 10 ** (2.170 - 0.630 * c - 1.580 * s)
        k_s = 10 ** (-0.600 - 0.640 * c + 1.260 * s) * 25.2 / 3600
        labda = 1 / b
        k_s = k_s * 8640000 / 1000  # kg/m2/s to cm/d
        psi_s = psi_s * 100  # m to cm
        return Brooks(
            k_s=round(k_s, 4),
            theta_r=0.0,
            theta_s=round(theta_s, 4),
            h_b=round(psi_s, 5),
            l=round(labda, 5),
        )


@dataclass
class Soil:
    name: str
    type: str | None = None
    model: SoilModel | None = None
    sample: SoilSample | None = None
    source: str | None = None
    description: str | None = None

    def from_name(self, sm: Type[SoilModel] | SoilModel | str) -> "Soil":
        if isinstance(sm, SoilModel):
            if hasattr(sm, "__name__"):
                smn = sm.__name__
            else:
                smn = sm.__class__.__name__
                sm = type(sm)
        elif isinstance(sm, str):
            smn = sm
            sm = get_soilmodel(smn)

        else:
            raise ValueError(
                f"Argument must either be Type[SoilModel] | SoilModel | str,"
                f"not {type(sm)}"
            )

        path = Path(__file__).parent / "datasets/Soil_Parameters.xlsx"
        ser = read_excel(path, sheet_name=smn, index_col=0).loc[self.name].to_dict()
        # path = Path(__file__).parent / f"datasets/{sm.__name__}.csv"
        # ser = read_csv(path, index_col=["name"]).loc[self.name].to_dict()
        self.__setattr__("type", ser.pop("soil type"))
        self.__setattr__("source", ser.pop("source"))
        self.__setattr__("description", ser.pop("description"))
        self.__setattr__("model", sm(**ser))
        return self

    @staticmethod
    def list_names(sm: Type[SoilModel] | SoilModel | str) -> list[str]:
        if isinstance(sm, SoilModel):
            if hasattr(sm, "__name__"):
                smn = sm.__name__
            else:
                smn = sm.__class__.__name__
                sm = type(sm)
        elif isinstance(sm, str):
            smn = sm
            sm = get_soilmodel(smn)

        else:
            raise ValueError(
                f"Argument must either be Type[SoilModel] | SoilModel | str,"
                f"not {type(sm)}"
            )

        path = Path(__file__).parent / "datasets/Soil_Parameters.xlsx"
        names = read_excel(path, sheet_name=smn).loc[:, "name"].to_list()
        return names

    def from_staring(self, year: str = "2018") -> "Soil":
        if year not in ("2001", 2001, "2018", 2018):
            raise ValueError(f"Year must either be '2001' or '2018', not {year}")
        path = Path(__file__).parent / f"datasets/Staring_{year}.xlsx"
        parameters = read_excel(path, sheet_name="parameters", index_col=0)
        ser = parameters.loc[self.name].to_dict()
        self.__setattr__("type", ser.pop("soil type"))
        self.__setattr__("source", ser.pop("source"))
        self.__setattr__("description", ser.pop("description"))
        sm = Genuchten(**ser)
        self.__setattr__("model", sm)
        ss = SoilSample().from_staring(name=self.name, year=year)
        self.__setattr__("sample", ss)
        return self
