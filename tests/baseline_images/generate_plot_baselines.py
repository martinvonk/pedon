"""Generate baseline images for pytest-mpl plot tests."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

import pedon as pe
from pedon.plot import curves, hcf, swrc

BASELINE_DIR = Path(__file__).resolve().parent / "test_plots"

soilmodel = pe.Genuchten(
    k_s=10.0,
    theta_r=0.01,
    theta_s=0.43,
    alpha=0.02,
    n=1.1,
    l=0.5,
)


def generate_baselines(dpi: int, bbox_inches: str) -> None:
    """Generate baseline images for plot tests."""
    BASELINE_DIR.mkdir(parents=True, exist_ok=True)

    savefig_kwargs = {
        "dpi": dpi,
        "bbox_inches": bbox_inches,
    }

    ax = swrc(soilmodel, color="tab:blue")
    ax.figure.savefig(BASELINE_DIR / "swrc_genuchten.png", **savefig_kwargs)
    plt.close(ax.figure)

    ax = hcf(soilmodel, color="tab:green")
    ax.figure.savefig(BASELINE_DIR / "hcf_genuchten.png", **savefig_kwargs)
    plt.close(ax.figure)

    axs = curves(soilmodel, color="tab:orange")
    axs[0].figure.savefig(BASELINE_DIR / "curves_genuchten.png", **savefig_kwargs)
    plt.close(axs[0].figure)


def main() -> None:
    """CLI entrypoint."""
    kwargs = {"dpi": 100, "bbox_inches": "tight"}
    generate_baselines(**kwargs)
    print(f"Generated plot baselines in {BASELINE_DIR}.")


if __name__ == "__main__":
    main()
