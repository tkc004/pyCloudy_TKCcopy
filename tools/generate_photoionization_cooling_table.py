#!/usr/bin/env python3
"""Generate a photoionization-equilibrium cooling table with Cloudy.

The table is evaluated on a grid of metallicity, gas temperature, hydrogen
density, and ionization parameter.  The temperature is held fixed with
Cloudy's ``constant temperature`` command, so the result is a cooling
function sampled at the requested temperatures while Cloudy solves the
photoionization state.

The two HDF5 datasets have shape
``(metallicity, temperature, hydrogen_density, ionization_parameter)`` and
contain volumetric cooling rates in erg cm^-3 s^-1.  The H/He file contains a
Cloudy run with metals suppressed.  The metal file contains full-metal
cooling minus that H/He baseline.  Each value is the arithmetic mean over the
Cloudy zones.

Example::

    python tools/generate_photoionization_cooling_table.py \
        --output photoionization_cooling_30kK.h5 \
        --teff 30000 \
        --logT-min 3.7 --logT-max 4.5 --nT 33 \
        --lognH-min 0 --lognH-max 4 --nnH 17 \
        --logU-min -4 --logU-max -1 --nU 13 \
        --metallicities 0.1 0.3 1.0 2.0 \
        --workers 4

Requirements: ``numpy``, ``h5py``, and a working Cloudy executable.
"""

from __future__ import annotations

import argparse
import itertools
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import h5py
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pyCloudy as pc  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a photoionization cooling table in HDF5 format."
    )
    parser.add_argument("--output", type=Path, default=Path("photoionization_cooling.h5"),
                        help="Base HDF5 filename; split files are written under --data-dir.")
    parser.add_argument("--data-dir", type=Path, default=Path("data"),
                        help="Directory for generated HDF5 tables, relative to the run directory.")
    parser.add_argument("--teff", type=float, default=30000.0,
                        help="Blackbody effective temperature in K.")
    parser.add_argument("--logT-min", type=float, default=3.0)
    parser.add_argument("--logT-max", type=float, default=5.0)
    parser.add_argument("--nT", type=int, default=101)
    parser.add_argument("--lognH-min", type=float, default=-2.0,
                        help="Minimum log10 hydrogen density in cm^-3.")
    parser.add_argument("--lognH-max", type=float, default=6.0,
                        help="Maximum log10 hydrogen density in cm^-3.")
    parser.add_argument("--nnH", type=int, default=81)
    parser.add_argument("--logU-min", type=float, default=-4.0)
    parser.add_argument("--logU-max", type=float, default=-1.0)
    parser.add_argument("--nU", type=int, default=31)
    parser.add_argument("--metallicities", type=float, nargs="+",
                        default=[0.1, 0.3, 1.0, 2.0],
                        help="Metallicity values in Z/Zsun.")
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of concurrent Cloudy runs.")
    parser.add_argument("--work-dir", type=Path, default=Path("cooling_models"),
                        help="Directory where individual Cloudy models are kept.")
    parser.add_argument("--cloudy-exe", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite the output HDF5 file if it exists.")
    parser.add_argument("--resume", action="store_true",
                        help="Reuse completed Cloudy models in --work-dir.")
    return parser.parse_args()


def default_cloudy_executable() -> Path:
    """Find the newest bundled Cloudy executable or one available on PATH."""
    bundled = sorted((REPO_ROOT / "Cloudy_exe" / "Cloudy").glob("*/source/cloudy.exe"))
    if bundled:
        return bundled[-1]
    found = shutil.which("cloudy.exe") or shutil.which("cloudy")
    if found:
        return Path(found)
    raise RuntimeError("Cloudy executable not found; pass --cloudy-exe.")


def make_cloudy_input(model_path, temperature, nH, metallicity, log_u, teff):
    """Write one isolated Cloudy input file."""
    model = pc.CloudyInput(str(model_path))
    model.set_BB(
        Teff=teff,
        lumi_unit="ionization parameter",
        lumi_value=log_u,
    )
    model.set_cste_density(np.log10(nH))
    model.set_cste_temperature(temperature)
    # Cloudy's ``metals`` argument is logarithmic.  The CLI and HDF5 coordinate
    # use the clearer physical quantity Z/Zsun, hence the conversion here.
    model.set_abund(
        predef="hii region",
        metals=np.log10(metallicity),
        nograins=True,
    )
    model.set_stop("neutral column density 20")
    model.set_comment(
        f"T={temperature:g} K, nH={nH:g} cm-3, Z/Zsun={metallicity:g}, logU={log_u:g}"
    )
    model.print_input(to_file=True, verbose=False)


def cooling_file_is_usable(path):
    """Return whether a saved cooling file contains at least one Ctot value."""
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        values = np.genfromtxt(path, comments="#", usecols=(3,), dtype=float)
    except (OSError, ValueError):
        return False
    return bool(np.any(np.isfinite(np.atleast_1d(values))))


def run_cloudy_rate(model_name, temperature, nH, metallicity, log_u, args):
    """Run or resume one isolated Cloudy model and return mean Ctot."""
    model_path = args.work_dir / model_name
    output_path = model_path.with_suffix(".out")
    cooling_path = model_path.with_suffix(".cool")
    reusable = args.resume and output_path.exists() and cooling_file_is_usable(cooling_path)
    if not reusable:
        make_cloudy_input(model_path, temperature, nH, metallicity, log_u, args.teff)

    if not reusable:
        completed = subprocess.run(
            [str(args.cloudy_exe), "-p", model_name],
            cwd=model_path.parent,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            log_path = model_path.with_suffix(".cloudy.log")
            log_path.write_text(completed.stdout or "", encoding="utf-8")
            raise RuntimeError(
                f"Cloudy failed for {model_name} with exit code {completed.returncode}; "
                f"see {log_path}"
            )

    # The fourth numeric column in ``save last cooling`` is Ctot.  Read this
    # file directly instead of loading every Cloudy extension through
    # CloudyModel; this also tolerates runs that do not produce auxiliary
    # radius/heat files under heavy parallel load.
    cooling = np.genfromtxt(
        cooling_path,
        comments="#",
        usecols=(3,),
        dtype=float,
    )
    cooling = np.atleast_1d(np.asarray(cooling, dtype=float))
    cooling = cooling[np.isfinite(cooling)]
    if cooling.size == 0:
        raise RuntimeError(f"No finite cooling values were produced for {model_name}")
    return float(np.mean(cooling))


def run_one(task, args):
    """Run total and H/He models and return their separated rates."""
    index, temperature, nH, metallicity, log_u = task
    total = run_cloudy_rate(
        f"total_model_{index:07d}", temperature, nH, metallicity, log_u, args
    )
    hhe = run_cloudy_rate(
        # Cloudy c22 asserts for an exactly zero-metal initialization.  A
        # 1e-10 solar metal floor is numerically safe and negligible for the
        # H/He baseline.
        f"hhe_model_{index:07d}", temperature, nH, 1.0e-10, log_u, args
    )
    return index, hhe, total - hhe


def compute_grid(temperatures, densities, metallicities, log_us, args):
    """Compute the grid, returning shape (nZ, nT, nnH, nU)."""
    shape = (len(metallicities), len(temperatures), len(densities), len(log_us))
    cooling_hhe = np.full(shape, np.nan, dtype=float)
    cooling_metals = np.full(shape, np.nan, dtype=float)
    tasks = []
    index = 0
    for temperature, nH, metallicity, log_u in itertools.product(
        temperatures, densities, metallicities, log_us
    ):
        tasks.append((index, temperature, nH, metallicity, log_u))
        index += 1

    coordinates = {
        (temperature, nH, metallicity, log_u): (z, t, d, u)
        for z, metallicity in enumerate(metallicities)
        for t, temperature in enumerate(temperatures)
        for d, nH in enumerate(densities)
        for u, log_u in enumerate(log_us)
    }

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(run_one, task, args) for task in tasks]
        for completed, future in enumerate(as_completed(futures), start=1):
            index, hhe_rate, metal_rate = future.result()
            _, temperature, nH, metallicity, log_u = tasks[index]
            z, t, d, u = coordinates[(temperature, nH, metallicity, log_u)]
            cooling_hhe[z, t, d, u] = hhe_rate
            cooling_metals[z, t, d, u] = metal_rate
            print(f"  completed {completed}/{len(tasks)}", flush=True)
    return cooling_hhe, cooling_metals


def write_hdf5(filename, temperatures, densities, metallicities, log_us,
               cooling, args, component):
    """Write coordinates, cooling values, and reproducibility metadata."""
    with h5py.File(filename, "w") as handle:
        handle.create_dataset("temperature_K", data=temperatures)
        handle.create_dataset("log10_temperature_K", data=np.log10(temperatures))
        handle.create_dataset("hydrogen_density_cm-3", data=densities)
        handle.create_dataset("log10_hydrogen_density_cm-3", data=np.log10(densities))
        handle.create_dataset("metallicity_Zsun", data=metallicities)
        handle.create_dataset("log10_ionization_parameter", data=log_us)
        handle.create_dataset(
            "cooling_erg_cm-3_s",
            data=cooling,
            compression="gzip",
            compression_opts=4,
        )
        handle.attrs["description"] = (
            f"Cloudy photoionization {component} cooling table at fixed gas temperature."
        )
        handle.attrs["spectrum"] = "blackbody"
        handle.attrs["stellar_Teff_K"] = args.teff
        handle.attrs["cooling_units"] = "erg cm^-3 s^-1"
        handle.attrs["temperature_units"] = "K"
        handle.attrs["hydrogen_density_units"] = "cm^-3"
        handle.attrs["metallicity_units"] = "Z/Zsun"
        handle.attrs["ionization_parameter_definition"] = "U = Phi(H)/(nH*c)"
        handle.attrs["cooling_definition"] = "arithmetic mean of local Cloudy zone cooling rates"
        handle.attrs["component"] = component
        handle.attrs["axis_order"] = (
            "cooling_erg_cm-3_s[metallicity, temperature, hydrogen_density, logU]"
        )


def main():
    args = parse_args()
    if args.workers < 1:
        raise SystemExit("ERROR: --workers must be at least 1")
    if args.nT < 1 or args.nnH < 1 or args.nU < 1:
        raise SystemExit("ERROR: nT, nnH, and nU must be positive")
    if args.teff <= 0 or any(value <= 0 for value in args.metallicities):
        raise SystemExit("ERROR: Teff and metallicities must be positive")
    args.cloudy_exe = args.cloudy_exe or default_cloudy_executable()
    if not args.cloudy_exe.exists():
        raise SystemExit(f"ERROR: Cloudy executable does not exist: {args.cloudy_exe}")
    args.work_dir.mkdir(parents=True, exist_ok=True)

    temperatures = np.logspace(args.logT_min, args.logT_max, args.nT)
    densities = np.logspace(args.lognH_min, args.lognH_max, args.nnH)
    log_us = np.linspace(args.logU_min, args.logU_max, args.nU)
    metallicities = np.asarray(args.metallicities, dtype=float)

    print("Photoionization cooling-table generation")
    print("==========================================")
    print(f"Cloudy executable = {args.cloudy_exe}")
    print(f"Output file       = {args.output}")
    print(f"Blackbody Teff    = {args.teff:g} K")
    print(f"Grid shape        = ({len(metallicities)}, {len(temperatures)}, {len(densities)}, {len(log_us)})")
    print(f"Workers           = {args.workers}")

    cooling_hhe, cooling_metals = compute_grid(
        temperatures, densities, metallicities, log_us, args
    )
    args.data_dir.mkdir(parents=True, exist_ok=True)
    stem = args.data_dir / args.output.stem
    hhe_output = stem.with_name(stem.name + "_HHe.h5")
    metals_output = stem.with_name(stem.name + "_metals.h5")
    if not args.overwrite and (hhe_output.exists() or metals_output.exists()):
        raise SystemExit(
            f"ERROR: split output exists; use --overwrite: {hhe_output}, {metals_output}"
        )
    write_hdf5(
        hhe_output, temperatures, densities, metallicities, log_us,
        cooling_hhe, args, "hydrogen+helium",
    )
    write_hdf5(
        metals_output, temperatures, densities, metallicities, log_us,
        cooling_metals, args, "metals",
    )
    print(f"Wrote: {hhe_output}")
    print(f"Wrote: {metals_output}")
    print("Main dataset: cooling_erg_cm-3_s[metallicity, temperature, hydrogen_density, logU]")


if __name__ == "__main__":
    main()
