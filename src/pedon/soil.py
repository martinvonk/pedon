# type: ignore
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Type

from numpy import abs as npabs
# from numpy import append, array2string, exp, full, log, log10
#### AP1025: how do you treat operations using arrays? i continue using np.multiply etc. in my modifications to pedon
from numpy import append, array2string, exp, full, log, log10, cos, divide, multiply, zeros, abs 
####
from pandas import DataFrame, isna, read_csv
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
    #### AP1025: modification for HYPAGS
    d10: float | None = None # d10 representative grain diameter (e.g. from sieve curve) [m]
    d20: float | None = None # d20 representative grain diameter (e.g. from sieve curve) [m]
    d50: float | None = None # d50 (mean) representative grain diameter (e.g. from sieve curve) [m]
    d60: float | None = None # d60 representative grain diameter (e.g. from sieve curve) [m]
    ne: float | None = None # effective porosity [-]
    ####
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
        path = Path(__file__).parent / "datasets/soilsamples.csv"
        properties = read_csv(path, delimiter=";")
        staring_properties = properties[
            properties["source"] == f"Staring_{year}"
        ].set_index("name")

        self.silt_p = staring_properties.loc[name, "silt_p"]
        self.clay_p = staring_properties.loc[name, "clay_p"]
        self.om_p = staring_properties.loc[name, "om_p"]
        self.m50 = staring_properties.loc[name, "m50"]
        if year == "2001":
            self.rho = staring_properties.loc[name, "rho"]
        return self

    def fit(
        self,
        sm: Type[SoilModel],
        pbounds: DataFrame | None = None,
        weights: FloatArray | float = 1.0,
        W1: float = 0.1,
        W2: float | None = None,
        k_s: float | None = None,
        silent: bool = True,
    ) -> SoilModel:
        """Same method as RETC"""

        theta = self.theta
        N = len(theta)
        k = self.k
        M = N + len(k)

        if pbounds is None:
            pbounds = get_params(sm.__name__)
            if k_s is not None:
                pbounds = pbounds.drop("k_s")
            else:
                pbounds.loc["k_s", "p_ini"] = max(k)
                pbounds.loc["k_s", "p_max"] = max(k) * 10
            pbounds.loc["theta_s", "p_ini"] = max(theta)
            pbounds.loc["theta_s", "p_max"] = max(theta) + 0.02

        if isinstance(weights, float):
            weights = full(M, weights)

        if W2 is None:
            W2 = (
                (M - N)
                * sum(weights[0:N] * theta)
                / (N * sum(weights[N:M] * npabs(log10(k))))
            )

        def get_diff(p: FloatArray) -> FloatArray:
            est_pars = dict(zip(pbounds.index, p))
            if k_s is not None:
                est_pars["k_s"] = k_s
            sml = sm(**est_pars)
            theta_diff = sml.theta(h=self.h) - theta
            k_diff = log10(sml.k(h=self.h)) - log10(k)
            diff = append(weights[0:N] * theta_diff, weights[N:M] * W1 * W2 * k_diff)
            return diff

        res = least_squares(
            get_diff,
            x0=pbounds.loc[:, "p_ini"],
            bounds=(
                pbounds.loc[:, "p_min"],
                pbounds.loc[:, "p_max"],
            ),
        )
        opt_pars = dict(zip(pbounds.index, res.x))
        if k_s is not None:
            opt_pars["k_s"] = k_s

        if not silent:
            print("SciPy Optimization Result\n", res)

        return sm(**opt_pars)

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
            + 0.0003658 * self.clay_p**2
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
            k_s=max(exp(ks_), 0),
            theta_r=theta_r,
            theta_s=theta_s,
            alpha=exp(alpha_),
            n=exp(n_) + 1,
            l=(10 * exp(l_) - 10) / (1 + exp(l_)),
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

    def rosetta(self, version: Literal[1, 2, 3] = 3) -> Genuchten:
        """Rosetta (Schaap et al., 2001) - Predicting soil water retention from soil"""
        try:
            from httpx import post as httpx_post
        except ImportError:
            raise ImportError(
                "httpx is required for the rosetta method to make api calls. "
                "Please install it with 'pip install httpx'."
            )

        soildata = [
            self.sand_p,  # %
            self.silt_p,  # %
            self.clay_p,  # %
            -9.9 if self.rho is None else self.rho,  # g/cm3
            -9.9 if self.th33 is None else self.th33,  # [-]
            -9.9 if self.th1500 is None else self.th1500,  # [-]
        ]

        data = {"soildata": [soildata]}
        r = httpx_post(
            f"https://www.handbook60.org/api/v1/rosetta/{version}",
            json=data,
        )
        if r.is_error:
            raise ValueError(f"Rosetta API error: {r}")

        rjson = r.json()
        vgpar = rjson["van_genuchten_params"][0]
        # stdev = rjson["stdev"][0] # TODO: return extra Genuchten classes with stdev or report/print
        if None in vgpar:
            raise ValueError(f"Rosetta API returned None values: {vgpar}")
        return Genuchten(
            k_s=10 ** vgpar[4],
            theta_r=vgpar[0],
            theta_s=vgpar[1],
            alpha=10 ** vgpar[2],
            n=10 ** vgpar[3],
        )
        
    #### AP1025: modification for HYPAGS
    def hypags(self) -> Genuchten:
        # Peche, A., Houben, G., & Altfelder, S. (2024). Approximation of van Genuchten Parameter Ranges from Hydraulic Conductivity Data. Groundwater, 62(3), 469-479.
        # Peche, A., & Houben, G. J. (2023). Estimating characteristic grain sizes and effective porosity from hydraulic conductivity data. Groundwater, 61(4), 574-585.
        
        """function carries out HYPAGS procedure

        Parameters
        ----------
        KK : Kozeny-Carman-based K --> note that only the Kozeny-Carman-based parameterization is implemented in pedon
    
        Calculates and eturns
        -------
        d1K : Kozeny-Carman-based d10
        d2K : Kozeny-Carman-based d20
        d5K : Kozeny-Carman-based d50
        d6K : Kozeny-Carman-based effective porosity
        poroK : Kozeny-Carman-based effective porosity
        h : capillary rise height
        drep : capillary rise representative grain diameter
        alpha : approximation of mean van Genuchten alpha
        alpha_plus_dalpha : upper range limit of van Genuchten alpha
        alpha_minus_dalpha : lower range limit of van Genuchten alpha
        n : approximation of mean van Genuchten n
        n_plus_dn : upper range limit of van Genuchten n
        n_minus_dn : lower range limit of van Genuchten n
        theta_r : residual water content
        theta_r_plus_dtheta_r : upper range limit of residual water content
        theta_r_minus_dtheta_r : lower range limit of residual water content
        
        """
        ### constants and coefficients
        
        rho_f, g, mu = 999.7, 9.81, 1.1306e-3 # fluid density [kg/m³], grav. accel. [m²/s], dynamic visc. [kg/(ms)] -> assumption: fluid temperature is 20°C
        c = mu/(rho_f*g) 
        a1, a2, a3 = 1.00401066e+00, 1.50991411e-04, 5.78888587e-03 # alpha-coefficients
        gamma, theta = .0728, .0 # surface tension [J/m²], contact angle, parameters for Young-Laplace eq.
        factor = 2*gamma*cos(theta)/(rho_f*g) # factor Young-Laplace eq. (fluid properties, contact angle, surface tension)
        crgd1, crgd2 = .0001803,  0.10 # coefficients for the calculation of the capillary-rise-representative grain diameter
        vGn1, vGn2, vGn3, vGn4, vGn5 = 1.12411782, 0.55750592, 1.67561574, 0.30307062, 1.16193138 # coefficients for the calculation of van Genuchten n
        Pi, c1, c2 = .0009, 1.2, 1.13 # dimensionless number Pi and c-coefficients --> [Kozeny-Carman-based parameterization]
        
        ### hypags functions
        def iterator(ne_i, d50_i, const):
            """
            

            Parameters (Parametrization-dependent)
            ----------
            K : input hydraulic conductivity
            d2 : input d20
            a1 : input coefficient a1
            a2 : input coefficient a2
            a3 : input coefficient a3
            const : fluid properties and gravitaiotnal acceleration constant

            Returns
            -------
            n_r : effective porosity
            d5 : d50

            """
            # iterative scheme to approximate d50 and effective porosity
            # constants 
            error = 10
            c = 180 * const
            n = zeros(100)
            d = zeros(100)
            d[0] = d50_i
            n[0] = ne_i
            j = 0
            while error > 1e-10:
                n[j+1] = (self.k/(d[j]**2) * c * (1 - n[j])**2)**(1/3)
                d[j+1] = (c * self.k * ( (1-n[j+1])**2 / (n[j+1]**3) ))**(1/2)
                error = abs(n[j+1]-n[j]+d[j+1]-d[j])
                j = j + 1
            n_r = n[j-1]
            d5 = d[j-1]
            return n_r, d5
         
        def getAlphaErrors(Kt):
        # Errors as determined using bootstrap resampling confidence intervals [note that confidence interval is the 2.5th and 97.5th quantile]
        # Error = [+ Delta alpha, - Delta alpha] = [added error, subtracted error]
            a_error = [0, 0]
            if Kt < .00000025:
                a_error = [3.002, 1.144]
            elif Kt > .00000025 and Kt <= .0000005:
                a_error = [5.364, 1.258]
            elif Kt > .0000005 and Kt <= .00000075:
                a_error = [7.031, 1.366]
            elif Kt > .00000075 and Kt <= .000001:
                a_error = [6.728, 1.519]
            elif Kt > .000001 and Kt <= .0000025:
                a_error = [6.758, 2.072]
            elif Kt > .0000025 and Kt <= .000005:
                a_error = [6.607, 2.393]
            elif Kt > .000005 and Kt <= .0000075:
                a_error = [5.91,  2.985]
            elif Kt > .0000075 and Kt <= .00001:
                a_error = [5.643, 3.229]
            elif Kt > .00001 and Kt <= .000025:
                a_error = [4.84, 4.]
            elif Kt > .000025 and Kt <= .00005:
                a_error = [4.082, 4.233]
            elif Kt > .00005 and Kt <= .000075:
                a_error = [3.182, 4.486]
            elif Kt > .000075 and Kt <= .0001:
                a_error = [3.051, 5.03 ]
            elif Kt > .0001 and Kt <= .0005:
                a_error = [3.442, 5.483]
            elif Kt > .0005:
                a_error = [1.811, 2.858]
            return a_error

        def getNErrors(alpha_t): 
        # Errors as determined using bootstrap resampling confidence intervals
        # Error = [+ Delta n, - Delta n] = [added error, subtracted error]
            if alpha_t < 1:
                n_error = [2.836, 0.358]
            elif alpha_t > 1 and alpha_t <= 2:
                n_error = [3.061, 0.826]
            elif alpha_t > 2 and alpha_t <= 3:
                n_error = [3.481, 0.828]
            elif alpha_t > 3 and alpha_t <= 4:
                n_error = [3.182, 0.605]
            elif alpha_t > 4 and alpha_t <= 5:
                n_error = [1.905, 0.509]
            elif alpha_t > 5 and alpha_t <= 6:
                n_error = [1.184, 0.297]
            elif alpha_t > 6 and alpha_t <= 7:
                n_error = [4.032, 0.352]
            elif alpha_t > 7 and alpha_t <= 8:
                n_error = [1.385, 0.241]
            elif alpha_t > 8 and alpha_t <= 9:
                n_error = [1.017, 0.155]
            elif alpha_t > 9:
                n_error = [3.006, 0.377]   
            return n_error
        
        def getResWaterContent(Kt):
            if Kt < .00000025:
                theta_rt = [0.097, 0.089, 0.198]
            elif Kt >= .00000025 and Kt < .0000005:
                theta_rt = [0.107, 0.097, 0.207]
            elif Kt >= .0000005 and Kt < .00000075:
                theta_rt = [0.099, 0.09, 0.227]
            elif Kt >= .00000075 and Kt < .000001:
                theta_rt = [0.102, 0.095, 0.221]
            elif Kt >= .000001 and Kt < .0000025:
                theta_rt = [0.098, 0.091, 0.216]
            elif Kt >= .0000025 and Kt < .000005:
                theta_rt = [0.097, 0.087, 0.227]
            elif Kt >= .000005 and Kt < .0000075:
                theta_rt = [0.111, 0.104, 0.212]
            elif Kt >= .0000075 and Kt < .00001:
                theta_rt = [0.108, 0.099, 0.219]
            elif Kt >= .00001 and Kt < .000025:
                theta_rt = [0.085, 0.079, 0.233]
            elif Kt >= .000025 and Kt < .00005:
                theta_rt = [0.087, 0.082, 0.226]
            elif Kt >= .00005 and Kt < .000075:
                theta_rt = [0.07, 0.067, 0.223]
            elif Kt >= .000075 and Kt < .0001:
                theta_rt = [0.11, 0.105, 0.23]
            elif Kt >= .0001 and Kt < .0005:
                theta_rt = [0.105, 0.096, 0.23]
            elif Kt >= .0005:
                theta_rt = [0.111, 0.068, 0.096]  
            return theta_rt
        
        ### mathematical model of hypags
        
        ## calculation of k, d10, d20: 
        
        # CONDITION: in HYPAGS, either k, d10 or d20 have to be given
        # case 0: k is given
        # case 1: d10 is given
        # case 2: d20 is given
        if self.k is not None:
            # mathematical model case 0
            # check for non-valid input
            if self.k > 2.6e-2 or self.k < 2.87e-7:
                print("Warning: k out of hypags model limits.")
            self.d10 = (self.k/Pi * c)**(.5) # calculation of d10
            self.d20 = c1 * self.d10 # calculation of d20
        elif self.d10 is not None:
            # mathematical model case 1
            if self.d10 > 8.3e-4 or self.d10 < 5.35e-5:
                print("Warning: d10 out of hypags model limits.")       
            self.k = (Pi * rho_f * g * self.d10**2)/mu # calculation of k
            self.d20 = c1 * self.d10    # calculation of d20
        elif self.d20 is not None:
            # mathematical model case 2
            if self.d20 > 1.2e-3 or self.d20 < 6.25e-5:
                print("Warning: d20 out of hypags model limits.")
            self.d10 = self.d20/c1 # calculation of d10
            self.k = (Pi * rho_f * g * self.d10**2)/mu
        else:
            # error: if CONDITION is not true
            raise ImportError(
                "No parameter (k, d10, or d20) was provided for hypags routine."
            )
            
        ## calculation of d50 and ne
        
        d50_0 = a3 * self.d20 # starting value for iterative approximation of d50
        ne_0 = a1 * self.k**a2
        self.ne, self.d50 = iterator(ne_0, d50_0, c)
        
        ## calculation of d60
        
        self.d60 = c2 * self.d50
        
        ## calculation of van Genuchten alpha
        
        if self.k <= 5e-5:
            drep = self.d60
        if self.k > 5e-5:
            drep = .0001803 + self.k *  0.10
        h = factor * divide(1, multiply(.5, drep**(1))) # Note the conversion of diameter to radius
        alphat = divide(1, h)
        
        
        ## calculation of van Genuchten n
        
        if alphat <= 1.9:
            nt = 1.12411782 + (alphat) * 0.55750592
        elif alphat > 1.9:
            nt = (1.67561574 / ((alphat)-0.30307062)) + 1.16193138
            
        ## error range of van Genuchten parameters
        
        # van Genuchten alpha 
        a_error = getAlphaErrors(self.k)
        alpha_plus_dalphat = alphat + a_error[0]
        alpha_minus_dalphat = alphat - a_error[1]
        if alpha_minus_dalphat <= 0: # Condition in case negative values for alpha arise
            alpha_minus_dalphat = .1
        # van Genuchten n
        n_error = getNErrors(alphat)
        n_plus_dnt = nt + n_error[0]
        n_minus_dnt = nt - n_error[1]  
        if n_minus_dnt <= 0:
            n_minus_dnt = .1

        ## return van Genuchten parameters
        
        return Genuchten(
            k_s=self.k,
            theta_r=None,
            theta_s=None,
            alpha=alphat,
            n=nt,
#            d_alpha_upper=alpha_plus_dalphat,
#            d_alpha_lower=alpha_minus_dalphat,
#            d_n_upper=n_plus_dnt,
#            d_n_lower=n_minus_dnt,
        )
        
        ####
@dataclass
class Soil:
    name: str
    model: SoilModel | None = None
    sample: SoilSample | None = None
    source: str | None = None
    description: str | None = None

    def from_name(
        self, sm: Type[SoilModel] | SoilModel | str, source: str | None = None
    ) -> "Soil":
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

        path = Path(__file__).parent / "datasets/soilsamples.csv"
        ser = read_csv(path, delimiter=";", index_col=0)
        if "HYDRUS_" in self.name:
            self.name = self.name.replace("HYDRUS_", "")
            source = "HYDRUS"
            logging.warning(
                "Removed 'HYDRUS_' from soil name. For future use"
                "please provide source='HYDRUS' argument"
            )
        sersm = ser[ser["soilmodel"] == smn].loc[[self.name], :]
        if source is None and len(sersm) > 1:
            raise Exception(
                f"Multiple sources for soil {self.name}: "
                f"{array2string(sersm.loc[:, 'source'].values)}. "
                f"Please provide the source using the source argument"
            )
        elif (source is not None) and len(sersm) > 1:
            sersm = sersm[sersm["source"] == source]

        serd = sersm.squeeze().to_dict()
        if isna(serd["description"]):
            serd["description"] = serd["soil type"]
        self.__setattr__("source", serd.pop("source"))
        self.__setattr__("description", serd.pop("description"))
        smserd = {
            x: serd[x]
            for x in sm.__dataclass_fields__.keys()
            if sm.__dataclass_fields__[x].init
        }
        self.__setattr__("model", sm(**smserd))
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

        path = Path(__file__).parent / "datasets/soilsamples.csv"
        names = read_csv(path, delimiter=";")

        return names[names["soilmodel"] == smn].loc[:, "name"].unique().tolist()

    def from_staring(self, year: str = "2018") -> "Soil":
        if year not in ("2001", 2001, "2018", 2018):
            raise ValueError(f"Year must either be '2001' or '2018', not {year}")

        self.from_name(sm=Genuchten, source=f"Staring_{year}")
        ss = SoilSample().from_staring(name=self.name, year=year)
        self.__setattr__("sample", ss)
        return self
