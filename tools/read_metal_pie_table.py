#!/usr/bin/env python3
"""Read a MetalPIE HDF5 table and display its Cloudy input and rates."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


AXES = (
    "log10_temperature_K",
    "log10_hydrogen_density_cm-3",
    "log10_ionization_parameter",
    "metallicity_Zsun",
)
RATES = (
    "metal_photoheating_erg_cm3_s",
    "metal_cooling_erg_cm3_s",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input", type=Path,
        help="MetalPIE HDF5 file, for example data/single_model_final_test_Z1_metals.h5",
    )
    parser.add_argument(
        "--show-values", action="store_true",
        help="Print all heating and cooling array values, not only summaries.",
    )
    return parser.parse_args()


def decode(value):
    """Decode scalar HDF5 byte strings while leaving other values unchanged."""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def read_metal_pie_table(path):
    """Read the MetalPIE group and return its contents as Python objects."""
    with h5py.File(path, "r") as handle:
        if "MetalPIE" not in handle:
            raise ValueError("HDF5 file does not contain a MetalPIE group")
        table = handle["MetalPIE"]
        result = {
            "attributes": {key: decode(value) for key, value in table.attrs.items()},
            "axes": {name: np.asarray(table[f"axes/{name}"]) for name in AXES},
            "rates": {name: np.asarray(table[f"rates/{name}"]) for name in RATES},
            "cloudy": {},
        }
        for name in ("input_file", "abundance_file", "command", "version"):
            result["cloudy"][name] = decode(table[f"cloudy/{name}"][()])
    return result


def print_rate_summary(name, values, show_values):
    finite = values[np.isfinite(values)]
    print(f"{name}: shape={values.shape}, dtype={values.dtype}")
    if finite.size:
        print(f"  range=[{np.min(finite):.6e}, {np.max(finite):.6e}]")
    else:
        print("  no finite values")
    if show_values:
        print(values)


def main():
    args = parse_args()
    if not args.input.exists():
        raise SystemExit(f"File does not exist: {args.input}")
    try:
        table = read_metal_pie_table(args.input)
    except (OSError, KeyError, ValueError) as exc:
        raise SystemExit(f"Could not read {args.input}: {exc}") from exc

    print(f"File: {args.input}")
    print("\nMetalPIE attributes:")
    for name, value in table["attributes"].items():
        print(f"  {name} = {value}")

    print("\nCloudy setup:")
    print("--- input_file ---")
    print(table["cloudy"]["input_file"], end="")
    print("--- abundance_file ---")
    print(table["cloudy"]["abundance_file"], end="")
    print("--- command ---")
    print(table["cloudy"]["command"])
    print("--- version ---")
    print(table["cloudy"]["version"])

    print("\nAxes:")
    for name in AXES:
        print(f"  {name}: {table['axes'][name]}")

    print("\nRates:")
    for name in RATES:
        print_rate_summary(name, table["rates"][name], args.show_values)


if __name__ == "__main__":
    main()
