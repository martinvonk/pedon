# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 09:50:47 2025

@author: Peche.A
"""

import sys
sys.path.insert(0, r"C:\Peche.a\Softwareentwicklungen\pedon_mod_gardner\pedon\src")

import pedon
print(pedon.__file__)


sample = pedon.SoilSample(k=9e-5)
result = sample.hypags()
print(result.alpha, result.n)
print(sample.k, result.theta_s)

# shared properties
k_s = 100  # saturated conductivity [cm/d]
theta_r = 0.03  # residual water content [-]
theta_s = 0.42  # saturated water content [-]



c = 0.005  # shape parameter k curve
m = 0.01  # shape parameter theta curve
gar = pedon.Gardner(k_s=k_s, theta_s=theta_s, c=c, m=m)
print(gar)

c = 0.005  # shape parameter k curve
m = 0.01  # shape parameter theta curve
mod_gar = pedon.Mod_Gardner(k_s=k_s, theta_s=theta_s, theta_r=theta_r, c=c, m=m)
print(mod_gar)


