# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:07:19 2026

@author: Peche.A
"""

import sys

sys.path.insert(
    0, r"C:\Peche.a\Softwareentwicklungen\pedon_gardner_approximator_ruff_v03\src"
)
import matplotlib.pyplot as plt

import pedon as pe


def main():
    """Main function
    this function wraps around
    testing of the new container for conversion methods

    g := gardner
    vg := van genuchten
    """
    # shared properties
    k_s = 0.00000722335 * 100 * 86400  # saturated conductivity [cm/d]
    theta_r = 0.0028  # residual water content [-]
    theta_s = 0.35  # saturated water content [-]

    # Mualem-van Genuchten
    alpha = 0.041  # shape parameter [1/cm]
    n = 2.56  # shape parameter [-]

    vg = pe.Genuchten(k_s=k_s, alpha=alpha, n=n, theta_s=theta_s, theta_r=theta_r)
    converter = pe.SoilModelConverter(sm=vg)
    g_p = converter.peche()
    print(g_p.c)
    g_g = converter.ghezzehei()
    print(g_g.c)

    # plot results
    plotFun(vg, g_p, g_g)
    # print if run successful
    printRun()


def plotFun(vg, g_p, g_g):
    f, axs = plt.subplots(1, 2, figsize=(7.0, 6.0), sharey=True, layout="tight")
    # f.suptitle("Testing converter", fontsize=16)
    # retention curve
    pe.plot_swrc(vg, ax=axs[0])
    pe.plot_swrc(g_p, ax=axs[0])
    pe.plot_swrc(g_g, ax=axs[0])
    axs[0].set_title("Soil Water Retention Curve")
    axs[0].set_xlabel("\N{GREEK SMALL LETTER THETA} [-]")
    axs[0].set_yscale("log")
    axs[0].set_ylabel("|\N{GREEK SMALL LETTER PSI}| [cm]")
    # rel. perm. fun
    pe.plot_hcf(vg, ax=axs[1])
    pe.plot_hcf(g_p, ax=axs[1])
    pe.plot_hcf(g_g, ax=axs[1])
    axs[1].set_title("Hydraulic Conductivity Function")
    axs[1].set_xlabel("Ks [cm/d]")
    axs[1].set_yscale("log")
    axs[0].legend(
        [
            "van Genuchten",
            "Gardner - converted, \n our new method",
            "Gardner - converted, \n Ghezzehei et al. (2007) method",
        ],
        loc="lower center",
        bbox_to_anchor=(0.4, 1.05),
        frameon=False,
    )


def printRun():
    """Print ascii string if succesful

    Returns
    -------
    None.

    """
    print("+-+-+-+-+-+-+-+-+-+-+ +-+-+-+")
    print("|s|u|c|c|e|s|s|f|u|l| |r|u|n|")
    print("+-+-+-+-+-+-+-+-+-+-+ +-+-+-+")


main()
