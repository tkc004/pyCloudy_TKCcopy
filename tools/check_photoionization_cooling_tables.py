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
                        help="Base filename before _Z<label>_HHe.h5 or _Z<label>_metals.h5.")
    parser.add_argument("--metallicity", type=float, default=1.0,
                        help="Metallicity label used when --hhe/--metals are omitted.")
    parser.add_argument("--hhe", type=Path, default=None)
    parser.add_argument("--metals", type=Path, default=None)
    parser.add_argument("--total", type=Path, default=None,
                        help="Optional HM12 total H/He+metal MetalPIE table.")
    return parser.parse_args()


def check_table(path, expected_component):
    errors = []
    warnings = []
    with h5py.File(path, "r") as handle:
        if "cooling_erg_cm-3_s" not in handle:
            errors.append("missing cooling_erg_cm-3_s dataset")
            return errors, warnings
        cooling = np.asarray(handle["cooling_erg_cm-3_s"])
        if "heating_erg_cm-3_s" not in handle:
            errors.append("missing heating_erg_cm-3_s dataset")
            return errors, warnings
        heating = np.asarray(handle["heating_erg_cm-3_s"])
        if cooling.ndim != 4:
            errors.append(f"cooling dataset has {cooling.ndim} dimensions; expected 4")
        if heating.ndim != 4:
            errors.append(f"heating dataset has {heating.ndim} dimensions; expected 4")
        if heating.shape != cooling.shape:
            errors.append(f"heating shape {heating.shape} != cooling shape {cooling.shape}")
        if np.any(np.isinf(cooling)):
            errors.append("cooling dataset contains infinite values")
        if np.any(np.isnan(cooling)):
            warnings.append("cooling dataset contains excluded cells stored as NaN")
        if np.any(np.isinf(heating)):
            errors.append("heating dataset contains infinite values")
        if np.any(np.isnan(heating)):
            warnings.append("heating dataset contains excluded cells stored as NaN")
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
            if heating.shape != expected_shape:
                errors.append(f"heating shape {heating.shape} != {expected_shape}")

        component = handle.attrs.get("component", "")
        if component != expected_component:
            errors.append(f"component metadata is {component!r}, expected {expected_component!r}")
        if handle.attrs.get("cooling_units") != "erg cm^-3 s^-1":
            errors.append("unexpected cooling units metadata")
        if handle.attrs.get("heating_units") != "erg cm^-3 s^-1":
            errors.append("unexpected heating units metadata")

        finite = cooling[np.isfinite(cooling)]
        finite_heating = heating[np.isfinite(heating)]
        print(
            f"{path}: shape={cooling.shape}, finite={finite.size}/{cooling.size}, "
            f"cooling=[{np.min(finite) if finite.size else np.nan:.3e}, "
            f"{np.max(finite) if finite.size else np.nan:.3e}], "
            f"heating=[{np.min(finite_heating) if finite_heating.size else np.nan:.3e}, "
            f"{np.max(finite_heating) if finite_heating.size else np.nan:.3e}]"
        )
    return errors, warnings


def check_metal_pie_table(path, expected_type="photoionization_equilibrium_metals"):
    """Check the MetalPIE grouped schema."""
    errors = []
    warnings = []
    with h5py.File(path, "r") as handle:
        if "MetalPIE" not in handle:
            return check_table(path, "metals")
        table = handle["MetalPIE"]
        hm12 = "redshift" in table.get("axes", {})
        required_attrs = {
            "table_type": expected_type,
            "cooling_units": "erg cm^-3 s^-1",
            "heating_units": "erg cm^-3 s^-1",
            "abundance_reference": "solar",
            "spectrum_type": (
                "Haardt-Madau 2012 UV background" if hm12 else "blackbody"
            ),
            "axis_order": (
                "temperature,density,redshift,metallicity" if hm12
                else "temperature,density,ionization_parameter,metallicity"
            ),
        }
        for name, expected in required_attrs.items():
            if table.attrs.get(name) != expected:
                errors.append(f"MetalPIE attribute {name!r} is not {expected!r}")
        axes = table.get("axes")
        rates = table.get("rates")
        if axes is None or rates is None:
            errors.append("MetalPIE must contain axes and rates groups")
            return errors, warnings
        axis_names = (
            "log10_temperature_K", "log10_hydrogen_density_cm-3",
            "redshift" if hm12 else "log10_ionization_parameter",
            "metallicity_Zsun",
        )
        if expected_type == "photoionization_equilibrium_total":
            rate_names = ("photoheating_erg_cm3_s", "cooling_erg_cm3_s")
        else:
            rate_names = ("metal_photoheating_erg_cm3_s", "metal_cooling_erg_cm3_s")
        if any(name not in axes for name in axis_names):
            errors.extend(f"missing MetalPIE/axes/{name}" for name in axis_names if name not in axes)
            return errors, warnings
        if any(name not in rates for name in rate_names):
            errors.extend(f"missing MetalPIE/rates/{name}" for name in rate_names if name not in rates)
            return errors, warnings
        expected_shape = tuple(len(axes[name]) for name in axis_names)
        for name in rate_names:
            values = np.asarray(rates[name])
            if values.shape != expected_shape:
                errors.append(f"{name} shape {values.shape} != {expected_shape}")
            if np.any(~np.isfinite(values)):
                errors.append(f"{name} contains non-finite values")
            if np.any(values < 0):
                errors.append(f"{name} contains negative values")
        print(
            f"{path}: MetalPIE shape={expected_shape}, "
            f"heating range=[{np.min(rates[rate_names[0]]):.3e}, "
            f"{np.max(rates[rate_names[0]]):.3e}], "
            f"cooling range=[{np.min(rates[rate_names[1]]):.3e}, "
            f"{np.max(rates[rate_names[1]]):.3e}]"
        )
    return errors, warnings


def main():
    args = parse_args()
    z_label = f"{args.metallicity:g}".replace("-", "m").replace(".", "p")
    hhe_path = args.hhe or args.data_dir / f"{args.stem}_Z{z_label}_HHe.h5"
    metals_path = args.metals or args.data_dir / f"{args.stem}_Z{z_label}_metals.h5"
    all_errors = []
    all_warnings = []
    paths = [
        (hhe_path, "hydrogen+helium", "photoionization_equilibrium_metals"),
        (metals_path, "metals", "photoionization_equilibrium_metals"),
    ]
    if args.total is not None:
        paths.append((args.total, "total", "photoionization_equilibrium_total"))
    for path, component, expected_type in paths:
        if not path.exists():
            all_errors.append(f"missing file: {path}")
            continue
        if component == "metals":
            errors, warnings = check_metal_pie_table(path, expected_type)
        elif component == "total":
            errors, warnings = check_metal_pie_table(path, expected_type)
        else:
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
