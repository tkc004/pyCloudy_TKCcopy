#!/usr/bin/env python3
"""Plot H/He heating and cooling versus temperature for density slices."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable, get_cmap


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path,
        default=Path("data/photoionization_cooling_Z1_HHe.h5"),
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("data/hhe_rates_vs_temperature_density.png"),
    )
    parser.add_argument(
        "--density-stride", type=int, default=1,
        help="Plot every Nth density; default plots every density.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.density_stride < 1:
        raise SystemExit("--density-stride must be positive")

    with h5py.File(args.input, "r") as handle:
        temperature = np.asarray(handle["temperature_K"])
        density = np.asarray(handle["hydrogen_density_cm-3"])
        log_u = float(np.asarray(handle["log10_ionization_parameter"])[0])
        cooling = np.asarray(handle["cooling_erg_cm-3_s"][0, :, :, 0], dtype=float)
        heating = np.asarray(handle["heating_erg_cm-3_s"][0, :, :, 0], dtype=float)

    if cooling.shape != heating.shape:
        raise SystemExit("Heating and cooling arrays have different shapes")

    indices = np.arange(0, len(density), args.density_stride)
    colors = plt.get_cmap("viridis")(
        np.linspace(0.0, 1.0, max(len(indices), 2))
    )
    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    for color, density_index in zip(colors, indices):
        cool = cooling[:, density_index]
        heat = heating[:, density_index]
        valid_cooling = np.isfinite(cool) & (cool > 0)
        valid_heating = np.isfinite(heat) & (heat > 0)
        label = fr"$n_{{\rm H}}={density[density_index]:g}\,\mathrm{{cm}}^{{-3}}$"
        if np.any(valid_cooling):
            ax.plot(temperature[valid_cooling], cool[valid_cooling], color=color,
                    linestyle="-", linewidth=1.0)
        if np.any(valid_heating):
            ax.plot(temperature[valid_heating], heat[valid_heating], color=color,
                    linestyle="--", linewidth=1.0)

    norm = LogNorm(vmin=density[indices[0]], vmax=density[indices[-1]])
    scalar_map = ScalarMappable(norm=norm, cmap="viridis")
    scalar_map.set_array([])
    colorbar = fig.colorbar(scalar_map, ax=ax)
    colorbar.set_label(r"Hydrogen density $n_{\rm H}$ [cm$^{-3}$]")

    ax.plot([], [], color="black", linestyle="-", label="Cooling")
    ax.plot([], [], color="black", linestyle="--", label="Heating")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel(r"Rate [erg cm$^{-3}$ s$^{-1}$]")
    ax.set_title(fr"H/He heating and cooling, $\log U={log_u:g}$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print(f"Wrote {args.output}")
    print(f"Plotted {len(indices)} density slices")


if __name__ == "__main__":
    main()
