#!/usr/bin/env python3
"""Plot metal cooling versus temperature for metallicity-specific HDF5 tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


LINE_STYLES = ("-", "--", ":", "-.", (0, (5, 1)))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir", type=Path, default=Path("data/data"),
        help="Directory containing *_metals.h5 files.",
    )
    parser.add_argument(
        "--pattern", default="*_metals.h5",
        help="Input filename pattern.",
    )
    parser.add_argument(
        "--component", choices=("metals", "HHe"), default="metals",
        help="Cooling-table component to plot.",
    )
    parser.add_argument(
        "--nH", type=float, default=1.0,
        help="Hydrogen density in cm^-3 for the plotted slice.",
    )
    parser.add_argument(
        "--logU", type=float, default=-2.0,
        help="log10 ionization parameter for the plotted slice.",
    )
    parser.add_argument(
        "--output", type=Path,
        default=None,
    )
    return parser.parse_args()


def nearest_index(axis, value):
    index = int(np.argmin(np.abs(np.asarray(axis) - value)))
    return index, float(np.asarray(axis)[index])


def main():
    args = parse_args()
    if args.pattern == "*_metals.h5" and args.component == "HHe":
        args.pattern = "*_HHe.h5"
    if args.output is None:
        label = "metal" if args.component == "metals" else "hhe"
        args.output = Path(f"data/data/{label}_cooling_vs_temperature.png")
    paths = sorted(args.data_dir.glob(args.pattern))
    if not paths:
        raise SystemExit(f"No files matching {args.pattern!r} in {args.data_dir}")

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    plotted = 0
    selected_nH = None
    selected_logU = None

    for line_style, path in zip(LINE_STYLES, paths):
        with h5py.File(path, "r") as handle:
            temperatures = np.asarray(handle["temperature_K"])
            densities = np.asarray(handle["hydrogen_density_cm-3"])
            log_us = np.asarray(handle["log10_ionization_parameter"])
            metallicity = float(np.asarray(handle["metallicity_Zsun"])[0])
            density_index, selected_nH = nearest_index(densities, args.nH)
            u_index, selected_logU = nearest_index(log_us, args.logU)
            cooling = np.asarray(
                handle["cooling_erg_cm-3_s"][0, :, density_index, u_index],
                dtype=float,
            )

        valid = np.isfinite(cooling) & (cooling > 0)
        if not np.any(valid):
            continue
        ax.plot(
            temperatures[valid], cooling[valid],
            linestyle=line_style, linewidth=2,
            label=fr"$Z={metallicity:g}\,Z_\odot$",
        )
        plotted += 1

    if plotted == 0:
        raise SystemExit("No positive finite cooling values were found to plot")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Temperature [K]")
    component_label = "Metal" if args.component == "metals" else "H/He"
    ax.set_ylabel(fr"{component_label} cooling rate [erg cm$^{{-3}}$ s$^{{-1}}$]")
    ax.set_title(
        fr"{component_label} cooling: $n_{{\rm H}}={selected_nH:g}\,\mathrm{{cm}}^{{-3}}$, "
        fr"$\log U={selected_logU:g}$"
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print(f"Wrote {args.output}")
    print(f"Used nearest grid slice nH={selected_nH:g} cm^-3, logU={selected_logU:g}")


if __name__ == "__main__":
    main()
