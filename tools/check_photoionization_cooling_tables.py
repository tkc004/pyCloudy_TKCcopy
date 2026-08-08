#!/usr/bin/env python3
"""Check the H/He and metal photoionization cooling HDF5 tables.

Example::

    python tools/check_photoionization_cooling_tables.py \
        --data-dir data --stem photoionization_cooling_30kK

The checker exits with status 1 for structural or numerical errors. Negative
metal contributions are reported as warnings because subtracting the H/He
baseline from two separately ionized Cloudy models can legitimately produce
small negative values in some cells.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


COORDINATES = (
    "temperature_K",
    "log10_temperature_K",
    "hydrogen_density_cm-3",
    "log10_hydrogen_density_cm-3",
    "metallicity_Zsun",
    "log10_ionization_parameter",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--stem", default="photoionization_cooling_30kK",
                        help="Base filename without _HHe.h5 or _metals.h5.")
    parser.add_argument("--hhe", type=Path, default=None)
    parser.add_argument("--metals", type=Path, default=None)
    return parser.parse_args()


def check_table(path, expected_component):
    errors = []
    warnings = []
    with h5py.File(path, "r") as handle:
        if "cooling_erg_cm-3_s" not in handle:
            errors.append("missing cooling_erg_cm-3_s dataset")
            return errors, warnings
        cooling = np.asarray(handle["cooling_erg_cm-3_s"])
        if cooling.ndim != 4:
            errors.append(f"cooling dataset has {cooling.ndim} dimensions; expected 4")
        if not np.all(np.isfinite(cooling)):
            errors.append("cooling dataset contains NaN or infinite values")
        if expected_component == "hydrogen+helium" and np.any(cooling < 0):
            errors.append("H/He cooling contains negative values")
        if expected_component == "metals" and np.any(cooling < 0):
            warnings.append("metal contribution contains negative values")

        for name in COORDINATES:
            if name not in handle:
                errors.append(f"missing coordinate dataset {name}")
        if all(name in handle for name in COORDINATES):
            if not np.allclose(handle["log10_temperature_K"],
                               np.log10(handle["temperature_K"])):
                errors.append("temperature and log-temperature axes disagree")
            if not np.allclose(handle["log10_hydrogen_density_cm-3"],
                               np.log10(handle["hydrogen_density_cm-3"])):
                errors.append("density and log-density axes disagree")
            if np.any(np.asarray(handle["temperature_K"]) <= 0):
                errors.append("temperature axis is not positive")
            if np.any(np.asarray(handle["hydrogen_density_cm-3"]) <= 0):
                errors.append("hydrogen-density axis is not positive")
            expected_shape = (
                len(handle["metallicity_Zsun"]),
                len(handle["temperature_K"]),
                len(handle["hydrogen_density_cm-3"]),
                len(handle["log10_ionization_parameter"]),
            )
            if cooling.shape != expected_shape:
                errors.append(f"cooling shape {cooling.shape} != {expected_shape}")

        component = handle.attrs.get("component", "")
        if component != expected_component:
            errors.append(f"component metadata is {component!r}, expected {expected_component!r}")
        if handle.attrs.get("cooling_units") != "erg cm^-3 s^-1":
            errors.append("unexpected cooling units metadata")

        print(f"{path}: shape={cooling.shape}, min={cooling.min():.3e}, max={cooling.max():.3e}")
    return errors, warnings


def main():
    args = parse_args()
    hhe_path = args.hhe or args.data_dir / f"{args.stem}_HHe.h5"
    metals_path = args.metals or args.data_dir / f"{args.stem}_metals.h5"
    all_errors = []
    all_warnings = []
    for path, component in ((hhe_path, "hydrogen+helium"), (metals_path, "metals")):
        if not path.exists():
            all_errors.append(f"missing file: {path}")
            continue
        errors, warnings = check_table(path, component)
        all_errors.extend(f"{path}: {error}" for error in errors)
        all_warnings.extend(f"{path}: {warning}" for warning in warnings)

    for warning in all_warnings:
        print(f"WARNING: {warning}")
    for error in all_errors:
        print(f"ERROR: {error}")
    if all_errors:
        return 1
    print("Cooling tables passed structural and numerical checks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
