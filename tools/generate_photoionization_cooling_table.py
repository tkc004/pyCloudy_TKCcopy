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
import csv
import datetime
import errno
import itertools
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

import h5py
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pyCloudy as pc  # noqa: E402


HHE_ELEMENT_COMMANDS = tuple(
    f"element {element} abundance -30"
    for element in (
        "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen",
        "Fluorine", "Neon", "Sodium", "Magnesium", "Aluminum", "Silicon",
        "Phosphorus", "Sulphur", "Chlorine", "Argon", "Potassium", "Calcium",
        "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron",
        "Cobalt", "Nickel", "Copper", "Zinc",
    )
)
HHE_OFF_COMMANDS = tuple(command.replace(" abundance -30", " off")
                          for command in HHE_ELEMENT_COMMANDS)


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
    parser.add_argument("--iterations", type=int, default=3,
                        help="Cloudy ionization iterations per model.")
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
    parser.add_argument("--keep-cloudy-files", action="store_true",
                        help="Keep each Cloudy model directory and its output files.")
    parser.add_argument("--cloudy-exe", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite the output HDF5 file if it exists.")
    parser.add_argument("--resume", action="store_true",
                        help="Reuse completed Cloudy models in --work-dir.")
    parser.add_argument("--retry-rejected", action="store_true",
                        help="Retry cells previously rejected for ionization failures.")
    return parser.parse_args()


def executable_matches_host(path):
    """Return whether a native executable appears compatible with this host."""
    try:
        magic = path.read_bytes()[:4]
    except OSError:
        return False
    system = platform.system()
    if system == "Linux":
        return magic == b"\x7fELF" or magic.startswith(b"#!")
    if system == "Darwin":
        return magic in (b"\xfe\xed\xfa\xce", b"\xfe\xed\xfa\xcf",
                         b"\xce\xfa\xed\xfe", b"\xcf\xfa\xed\xfe")
    if system == "Windows":
        return magic[:2] == b"MZ"
    return True


def default_cloudy_executable() -> Path:
    """Find a host-compatible bundled Cloudy or one available on PATH."""
    bundled = sorted((REPO_ROOT / "Cloudy_exe" / "Cloudy").glob("*/source/cloudy.exe"))
    found = shutil.which("cloudy.exe") or shutil.which("cloudy")
    if found:
        return Path(found)
    compatible = [path for path in bundled if executable_matches_host(path)]
    if compatible:
        return compatible[-1]
    raise RuntimeError("Cloudy executable not found; pass --cloudy-exe.")


def make_cloudy_input(model_path, temperature, nH, metallicity, log_u, teff,
                      iterations, hhe=False):
    """Write one isolated Cloudy input file."""
    model = pc.CloudyInput(str(model_path))
    model.set_BB(
        Teff=teff,
        lumi_unit="ionization parameter",
        lumi_value=log_u,
    )
    model.set_cste_density(np.log10(nH))
    model.set_cste_temperature(temperature, others="linear")
    model.set_abund(predef="hii region", nograins=True)
    if hhe:
        model.set_other(HHE_ELEMENT_COMMANDS + HHE_OFF_COMMANDS + (
            "metals 1e-30",
        ))
    else:
        # Cloudy's ``metals`` argument is logarithmic.  The CLI and HDF5
        # coordinate use the clearer physical quantity Z/Zsun.
        model.set_abund(metals=np.log10(metallicity))
    # Molecular chemistry is not part of this optically thin ionized-gas
    # table.  Apply this to both the total-metal and H/He input recipes.
    model.set_other(("no molecules",))
    model.set_stop("zone 1")
    model.set_iterate(to_convergence=True)
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


class CloudyCellRejected(RuntimeError):
    """A Cloudy cell produced an ionization-convergence failure."""


def ionization_failure_count(output_path):
    """Read Cloudy's ionization-failure count from the .out summary."""
    if not output_path.exists():
        return None
    text = output_path.read_text(encoding="utf-8", errors="replace")
    matches = re.findall(r"Failures:\s*\d+\s+thermal,\s*\d+\s+pressure,\s*(\d+)\s+ionization", text)
    return int(matches[-1]) if matches else None


def cleanup_model_artifacts(model_name, model_dir, work_dir, keep=False):
    """Remove Cloudy files after their cooling rate has been extracted."""
    if keep:
        return
    if model_dir != work_dir and model_dir.is_dir():
        last_error = None
        for attempt in range(5):
            try:
                shutil.rmtree(model_dir)
                return
            except OSError as exc:
                last_error = exc
                time.sleep(0.5 * (attempt + 1))
        # The cooling value is already safely in memory.  Do not discard a
        # completed grid point merely because NFS delayed directory updates.
        print(
            f"WARNING: could not fully remove {model_dir} after 5 attempts: "
            f"{last_error}",
            flush=True,
        )
        return
    for path in work_dir.glob(f"{model_name}.*"):
        if path.is_file():
            try:
                path.unlink()
            except OSError as exc:
                print(f"WARNING: could not remove {path}: {exc}", flush=True)


def run_cloudy_rate(model_name, temperature, nH, metallicity, log_u, args,
                    hhe=False):
    """Run or resume one isolated Cloudy model and return mean Ctot."""
    # Keep every new Cloudy invocation in a private directory.  Unique
    # prefixes alone are insufficient because Cloudy can also create shared
    # temporary files in its current working directory.  Retain compatibility
    # with completed models from older runs, which stored files directly in
    # --work-dir.
    legacy_path = args.work_dir / model_name
    legacy_out = legacy_path.with_suffix(".out")
    legacy_cool = legacy_path.with_suffix(".cool")
    legacy_reusable = (
        args.resume and legacy_out.exists() and cooling_file_is_usable(legacy_cool)
    )
    if legacy_reusable:
        model_dir = args.work_dir
        model_path = legacy_path
    else:
        model_dir = args.work_dir / model_name
        model_dir.mkdir(parents=True, exist_ok=True)
        model_path = model_dir / model_name
    output_path = model_path.with_suffix(".out")
    cooling_path = model_path.with_suffix(".cool")
    reusable = args.resume and output_path.exists() and cooling_file_is_usable(cooling_path)
    if not reusable:
        if model_dir != args.work_dir and model_dir.exists():
            shutil.rmtree(model_dir)
            model_dir.mkdir(parents=True, exist_ok=True)
        make_cloudy_input(
            model_path, temperature, nH, metallicity, log_u,
            args.teff, args.iterations, hhe,
        )

    if not reusable:
        try:
            completed = subprocess.run(
                [str(args.cloudy_exe), "-p", model_name],
                cwd=model_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                env={**os.environ, "OMP_NUM_THREADS": "1"},
            )
        except OSError as exc:
            if exc.errno == errno.ENOEXEC:
                raise RuntimeError(
                    f"Cloudy executable is not native to this host: {args.cloudy_exe}. "
                    "Compile Cloudy on the cluster or pass a compatible Linux path "
                    "with --cloudy-exe."
                ) from exc
            raise
        failures = ionization_failure_count(output_path)
        if failures and failures > 0:
            cleanup_model_artifacts(
                model_name, model_dir, args.work_dir, args.keep_cloudy_files
            )
            raise CloudyCellRejected(
                f"{model_name}: Cloudy reported {failures} ionization failures"
            )
        if completed.returncode != 0:
            log_path = model_path.with_suffix(".cloudy.log")
            log_path.write_text(completed.stdout or "", encoding="utf-8")
            raise RuntimeError(
                f"Cloudy failed for {model_name} with exit code {completed.returncode}; "
                f"see {log_path}"
            )

    failures = ionization_failure_count(output_path)
    if failures and failures > 0:
        cleanup_model_artifacts(
            model_name, model_dir, args.work_dir, args.keep_cloudy_files
        )
        raise CloudyCellRejected(
            f"{model_name}: Cloudy reported {failures} ionization failures"
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
    rate = float(np.mean(cooling))
    cleanup_model_artifacts(
        model_name, model_dir, args.work_dir, args.keep_cloudy_files
    )
    return rate


def run_one(task, args):
    """Run total and H/He models and return their separated rates."""
    index, temperature, nH, metallicity, log_u = task
    total = run_cloudy_rate(
        f"total_model_{index:07d}", temperature, nH, metallicity, log_u, args
    )
    hhe = run_cloudy_rate(
        f"hhe_model_{index:07d}", temperature, nH, metallicity, log_u, args,
        hhe=True,
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

    completed_rows = {}
    rejected_indices = set()
    checkpoint = args.checkpoint
    if args.resume and checkpoint.exists():
        with checkpoint.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                row_index = int(row["index"])
                if row.get("status", "complete") == "rejected":
                    rejected_indices.add(row_index)
                else:
                    completed_rows[row_index] = row
        for index, row in completed_rows.items():
            if 0 <= index < len(tasks):
                _, temperature, nH, metallicity, log_u = tasks[index]
                if np.isclose(float(row["temperature_K"]), temperature) and np.isclose(
                    float(row["nH_cm-3"]), nH
                ) and np.isclose(float(row["metallicity_Zsun"]), metallicity) and np.isclose(
                    float(row["log10_ionization_parameter"]), log_u
                ):
                    z, t, d, u = coordinates_for_task(
                        index, tasks, metallicities, temperatures, densities, log_us
                    )
                    cooling_hhe[z, t, d, u] = float(row["hhe_cooling_erg_cm-3_s"])
                    cooling_metals[z, t, d, u] = float(row["metal_cooling_erg_cm-3_s"])

    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_exists = checkpoint.exists()
    checkpoint_handle = checkpoint.open("a", newline="", encoding="utf-8")
    checkpoint_writer = csv.DictWriter(
        checkpoint_handle,
        fieldnames=[
            "index", "temperature_K", "nH_cm-3", "metallicity_Zsun",
            "log10_ionization_parameter", "hhe_cooling_erg_cm-3_s",
            "metal_cooling_erg_cm-3_s", "status", "error",
        ],
    )
    if not checkpoint_exists:
        checkpoint_writer.writeheader()

    pending_tasks = [
        task for task in tasks
        if task[0] not in completed_rows
        and not (args.resume and task[0] in rejected_indices and not args.retry_rejected)
    ]

    coordinates = {
        (temperature, nH, metallicity, log_u): (z, t, d, u)
        for z, metallicity in enumerate(metallicities)
        for t, temperature in enumerate(temperatures)
        for d, nH in enumerate(densities)
        for u, log_u in enumerate(log_us)
    }

    completed = 0
    task_iterator = iter(pending_tasks)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        in_flight = {}
        for _ in range(args.workers):
            try:
                task = next(task_iterator)
            except StopIteration:
                break
            in_flight[executor.submit(run_one, task, args)] = task

        # Keep only ``workers`` tasks submitted at once.  This prevents a huge
        # grid from eagerly creating all model inputs/futures before threads
        # have capacity to run them.
        while in_flight:
            finished, _ = wait(in_flight, return_when=FIRST_COMPLETED)
            for future in finished:
                task = in_flight.pop(future)
                _, temperature, nH, metallicity, log_u = task
                try:
                    index, hhe_rate, metal_rate = future.result()
                except CloudyCellRejected as exc:
                    index = task[0]
                    checkpoint_writer.writerow({
                        "index": index,
                        "temperature_K": temperature,
                        "nH_cm-3": nH,
                        "metallicity_Zsun": metallicity,
                        "log10_ionization_parameter": log_u,
                        "hhe_cooling_erg_cm-3_s": "",
                        "metal_cooling_erg_cm-3_s": "",
                        "status": "rejected",
                        "error": str(exc),
                    })
                    checkpoint_handle.flush()
                    completed += 1
                    print(f"  rejected {completed}/{len(pending_tasks)}: {exc}", flush=True)
                    try:
                        next_task = next(task_iterator)
                    except StopIteration:
                        continue
                    in_flight[executor.submit(run_one, next_task, args)] = next_task
                    continue
                z, t, d, u = coordinates[(temperature, nH, metallicity, log_u)]
                cooling_hhe[z, t, d, u] = hhe_rate
                cooling_metals[z, t, d, u] = metal_rate
                checkpoint_writer.writerow({
                    "index": index,
                    "temperature_K": temperature,
                    "nH_cm-3": nH,
                    "metallicity_Zsun": metallicity,
                    "log10_ionization_parameter": log_u,
                    "hhe_cooling_erg_cm-3_s": hhe_rate,
                    "metal_cooling_erg_cm-3_s": metal_rate,
                    "status": "complete",
                    "error": "",
                })
                checkpoint_handle.flush()
                completed += 1
                print(f"  completed {completed}/{len(pending_tasks)} pending models", flush=True)
                try:
                    next_task = next(task_iterator)
                except StopIteration:
                    continue
                in_flight[executor.submit(run_one, next_task, args)] = next_task
    checkpoint_handle.close()
    return cooling_hhe, cooling_metals


def coordinates_for_task(index, tasks, metallicities, temperatures, densities, log_us):
    """Return array indices for a task without storing a coordinate dictionary."""
    _, temperature, nH, metallicity, log_u = tasks[index]
    z = int(np.flatnonzero(np.isclose(metallicities, metallicity))[0])
    t = int(np.flatnonzero(np.isclose(temperatures, temperature))[0])
    d = int(np.flatnonzero(np.isclose(densities, nH))[0])
    u = int(np.flatnonzero(np.isclose(log_us, log_u))[0])
    return z, t, d, u


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
        handle.attrs["pycloudy_version"] = pc.__version__
        handle.attrs["cloudy_version"] = cloudy_version_from_path(args.cloudy_exe)
        handle.attrs["cloudy_executable"] = str(args.cloudy_exe)
        handle.attrs["spectrum_description"] = "blackbody"
        handle.attrs["spectrum_Teff_K"] = args.teff
        handle.attrs["hhe_abundance_recipe"] = (
            "all metals abundance -30, all metals off, metals 1e-30, no molecules"
        )
        handle.attrs["python_version"] = platform.python_version()
        handle.attrs["platform"] = platform.platform()
        handle.attrs["run_timestamp_utc"] = datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat()
        handle.attrs["command_line"] = shlex.join(sys.argv)
        handle.attrs["grid_logT_min"] = args.logT_min
        handle.attrs["grid_logT_max"] = args.logT_max
        handle.attrs["grid_nT"] = args.nT
        handle.attrs["grid_lognH_min"] = args.lognH_min
        handle.attrs["grid_lognH_max"] = args.lognH_max
        handle.attrs["grid_nnH"] = args.nnH
        handle.attrs["grid_logU_min"] = args.logU_min
        handle.attrs["grid_logU_max"] = args.logU_max
        handle.attrs["grid_nU"] = args.nU
        handle.attrs["workers"] = args.workers
        handle.attrs["cloudy_iterations"] = args.iterations
        handle.attrs["geometry"] = "plane-parallel slab"
        handle.attrs["stop_criterion"] = "zone 1"
        handle.attrs["iterate_mode"] = "to convergence"
        handle.attrs["temperature_command"] = "constant temperature T K linear"
        handle.attrs["axis_order"] = (
            "cooling_erg_cm-3_s[metallicity, temperature, hydrogen_density, logU]"
        )


def metallicity_label(value):
    """Create a filesystem-safe metallicity label."""
    label = f"{value:g}".replace("-", "m").replace(".", "p")
    return label


def cloudy_version_from_path(path):
    """Extract a version such as c22.02 from a Cloudy executable path."""
    match = re.search(r"(?:^|[/_])c(\d+\.\d+)(?:[/_.]|$)", str(path), re.IGNORECASE)
    return match.group(1) if match else "unknown"


def main():
    args = parse_args()
    if args.workers < 1:
        raise SystemExit("ERROR: --workers must be at least 1")
    if args.iterations < 1:
        raise SystemExit("ERROR: --iterations must be at least 1")
    if args.nT < 1 or args.nnH < 1 or args.nU < 1:
        raise SystemExit("ERROR: nT, nnH, and nU must be positive")
    if args.teff <= 0 or any(value <= 0 for value in args.metallicities):
        raise SystemExit("ERROR: Teff and metallicities must be positive")
    args.cloudy_exe = args.cloudy_exe or default_cloudy_executable()
    if not args.cloudy_exe.exists():
        raise SystemExit(f"ERROR: Cloudy executable does not exist: {args.cloudy_exe}")
    if not executable_matches_host(args.cloudy_exe):
        raise SystemExit(
            f"ERROR: Cloudy executable is not compatible with {platform.system()}: "
            f"{args.cloudy_exe}. Compile Cloudy on this host or pass --cloudy-exe."
        )
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.data_dir.mkdir(parents=True, exist_ok=True)
    args.checkpoint = args.data_dir / f"{args.output.stem}.checkpoint.csv"
    if not args.resume and args.checkpoint.exists():
        args.checkpoint.unlink()

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
    stem = args.data_dir / args.output.stem
    for z_index, metallicity in enumerate(metallicities):
        z_label = metallicity_label(metallicity)
        hhe_output = stem.with_name(stem.name + f"_Z{z_label}_HHe.h5")
        metals_output = stem.with_name(stem.name + f"_Z{z_label}_metals.h5")
        if not args.overwrite and (hhe_output.exists() or metals_output.exists()):
            raise SystemExit(
                f"ERROR: output exists; use --overwrite: {hhe_output}, {metals_output}"
            )
        # Keep a one-element metallicity axis in every file, making the files
        # directly compatible with the existing checker and lookup code.
        write_hdf5(
            hhe_output, temperatures, densities, metallicities[z_index:z_index + 1],
            log_us, cooling_hhe[z_index:z_index + 1], args, "hydrogen+helium",
        )
        write_hdf5(
            metals_output, temperatures, densities, metallicities[z_index:z_index + 1],
            log_us, cooling_metals[z_index:z_index + 1], args, "metals",
        )
        print(f"Wrote: {hhe_output}")
        print(f"Wrote: {metals_output}")
    print("Main dataset: cooling_erg_cm-3_s[metallicity, temperature, hydrogen_density, logU]")


if __name__ == "__main__":
    main()
