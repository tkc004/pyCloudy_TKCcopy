#!/usr/bin/env python3
"""Generate MetalPIE heating/cooling tables using the Cloudy HM12 background.

The table axes are temperature, hydrogen density, HM12 redshift, and
metallicity.  One ``*_metals.h5`` file is written per metallicity.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import shlex
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

import h5py
import numpy as np

import generate_photoionization_cooling_table as base


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("metal_pie_hm12.h5"))
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--logT-min", type=float, default=1.0)
    parser.add_argument("--logT-max", type=float, default=8.0)
    parser.add_argument("--nT", type=int, default=101)
    parser.add_argument("--lognH-min", type=float, default=-6.0)
    parser.add_argument("--lognH-max", type=float, default=2.0)
    parser.add_argument("--nnH", type=int, default=81)
    parser.add_argument("--redshift-min", type=float, default=0.0)
    parser.add_argument("--redshift-max", type=float, default=15.93)
    parser.add_argument("--n-redshift", type=int, default=32)
    parser.add_argument("--metallicities", type=float, nargs="+",
                        default=[0.0, 0.01, 0.1, 1.0, 2.0])
    parser.add_argument("--teff", type=float, default=0.0,
                        help="Unused compatibility option; HM12 supplies the spectrum.")
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--work-dir", type=Path, default=Path("cooling_models_hm12"))
    parser.add_argument("--keep-cloudy-files", action="store_true")
    parser.add_argument("--cloudy-exe", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--retry-rejected", action="store_true")
    return parser.parse_args()


def make_hm12_input(model_path, temperature, nH, metallicity, redshift,
                    teff, iterations, hhe=False):
    """Write a Cloudy input using ``table HM12 redshift`` and return its text."""
    model = base.pc.CloudyInput(str(model_path))
    model.set_other((f"table HM12 redshift {redshift:g}",))
    model.set_cste_density(np.log10(nH))
    model.set_cste_temperature(temperature, others="linear")
    model.set_abund(predef="hii region", nograins=True)
    if hhe:
        model.set_other(base.HHE_ELEMENT_COMMANDS + base.HHE_OFF_COMMANDS + (
            "metals 1e-30",
        ))
    else:
        if metallicity <= 0:
            raise ValueError("zero metallicity does not require a total-metal model")
        model.set_other((f"metals {np.log10(metallicity):g}",))
    model.set_other(("no molecules",))
    model.set_stop("zone 1")
    model.set_iterate(to_convergence=True)
    model.set_comment(
        f"T={temperature:g} K, nH={nH:g} cm-3, Z/Zsun={metallicity:g}, z={redshift:g}"
    )
    model.print_input(to_file=True, verbose=False)
    return model_path.with_suffix(".in").read_text(encoding="utf-8")


def configure_base(args):
    """Adapt the shared Cloudy runner to the HM12 input recipe."""
    base.make_cloudy_input = make_hm12_input
    args.teff = 0.0


def compute_grid(temperatures, densities, redshifts, metallicities, args):
    shape = (len(metallicities), len(temperatures), len(densities), len(redshifts))
    cooling = np.full(shape, np.nan, dtype=float)
    heating = np.full(shape, np.nan, dtype=float)
    total_cooling = np.full(shape, np.nan, dtype=float)
    total_heating = np.full(shape, np.nan, dtype=float)
    tasks = [
        (index, temperature, nH, metallicity, redshift)
        for index, (temperature, nH, metallicity, redshift) in enumerate(
            itertools.product(temperatures, densities, metallicities, redshifts)
        )
    ]
    coordinates = {
        (temperature, nH, metallicity, redshift): (z, t, d, r)
        for z, metallicity in enumerate(metallicities)
        for t, temperature in enumerate(temperatures)
        for d, nH in enumerate(densities)
        for r, redshift in enumerate(redshifts)
    }
    completed_rows = {}
    rejected_indices = set()
    checkpoint = args.checkpoint
    if args.resume and checkpoint.exists():
        with checkpoint.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                index = int(row["index"])
                if row.get("status", "complete") == "rejected":
                    rejected_indices.add(index)
                elif row.get("metal_heating_erg_cm-3_s", "") and row.get(
                    "metal_cooling_erg_cm-3_s", ""
                ) and row.get("total_heating_erg_cm-3_s", "") and row.get(
                    "total_cooling_erg_cm-3_s", ""
                ):
                    completed_rows[index] = row
        for index, row in completed_rows.items():
            if 0 <= index < len(tasks):
                _, temperature, nH, metallicity, redshift = tasks[index]
                if all((
                    np.isclose(float(row["temperature_K"]), temperature),
                    np.isclose(float(row["nH_cm-3"]), nH),
                    np.isclose(float(row["metallicity_Zsun"]), metallicity),
                    np.isclose(float(row["redshift"]), redshift),
                )):
                    z, t, d, r = coordinates[(temperature, nH, metallicity, redshift)]
                    cooling[z, t, d, r] = float(row["metal_cooling_erg_cm-3_s"])
                    heating[z, t, d, r] = float(row["metal_heating_erg_cm-3_s"])
                    total_cooling[z, t, d, r] = float(row["total_cooling_erg_cm-3_s"])
                    total_heating[z, t, d, r] = float(row["total_heating_erg_cm-3_s"])

    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    exists = checkpoint.exists()
    checkpoint_handle = checkpoint.open("a", newline="", encoding="utf-8")
    writer = csv.DictWriter(
        checkpoint_handle,
        fieldnames=[
            "index", "temperature_K", "nH_cm-3", "redshift", "metallicity_Zsun",
            "metal_heating_erg_cm-3_s", "metal_cooling_erg_cm-3_s", "status", "error",
            "total_heating_erg_cm-3_s", "total_cooling_erg_cm-3_s",
        ],
    )
    if not exists:
        writer.writeheader()

    pending = [
        task for task in tasks
        if task[0] not in completed_rows
        and not (args.resume and task[0] in rejected_indices and not args.retry_rejected)
    ]
    first_input_file = None
    completed = 0
    task_iterator = iter(pending)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        in_flight = {}
        for _ in range(args.workers):
            try:
                task = next(task_iterator)
            except StopIteration:
                break
            in_flight[executor.submit(base.run_one, task, args)] = task
        while in_flight:
            finished, _ = wait(in_flight, return_when=FIRST_COMPLETED)
            for future in finished:
                task = in_flight.pop(future)
                index, temperature, nH, metallicity, redshift = task
                try:
                    (
                        _, hhe_rate, metal_rate, hhe_heating_rate,
                        metal_heating_rate, input_file,
                    ) = future.result()
                    if first_input_file is None and input_file:
                        first_input_file = input_file
                    z, t, d, r = coordinates[(temperature, nH, metallicity, redshift)]
                    cooling[z, t, d, r] = max(float(metal_rate), 0.0)
                    heating[z, t, d, r] = max(float(metal_heating_rate), 0.0)
                    total_cooling[z, t, d, r] = max(
                        float(hhe_rate) + float(metal_rate), 0.0
                    )
                    total_heating[z, t, d, r] = max(
                        float(hhe_heating_rate) + float(metal_heating_rate), 0.0
                    )
                    writer.writerow({
                        "index": index, "temperature_K": temperature, "nH_cm-3": nH,
                        "redshift": redshift, "metallicity_Zsun": metallicity,
                        "metal_heating_erg_cm-3_s": heating[z, t, d, r],
                        "metal_cooling_erg_cm-3_s": cooling[z, t, d, r],
                        "total_heating_erg_cm-3_s": total_heating[z, t, d, r],
                        "total_cooling_erg_cm-3_s": total_cooling[z, t, d, r],
                        "status": "complete", "error": "",
                    })
                except base.CloudyCellRejected as exc:
                    writer.writerow({
                        "index": index, "temperature_K": temperature, "nH_cm-3": nH,
                        "redshift": redshift, "metallicity_Zsun": metallicity,
                        "metal_heating_erg_cm-3_s": "", "metal_cooling_erg_cm-3_s": "",
                        "total_heating_erg_cm-3_s": "", "total_cooling_erg_cm-3_s": "",
                        "status": "rejected", "error": str(exc),
                    })
                checkpoint_handle.flush()
                completed += 1
                print(f"  completed {completed}/{len(pending)} pending models", flush=True)
                try:
                    next_task = next(task_iterator)
                except StopIteration:
                    continue
                in_flight[executor.submit(base.run_one, next_task, args)] = next_task
    checkpoint_handle.close()
    return cooling, heating, total_cooling, total_heating, first_input_file


def write_hm12_table(path, temperatures, densities, redshifts, metallicities,
                     cooling, heating, args, input_file, component):
    cooling = np.maximum(np.transpose(cooling, (1, 2, 3, 0)), 0.0)
    heating = np.maximum(np.transpose(heating, (1, 2, 3, 0)), 0.0)
    with h5py.File(path, "w") as handle:
        table = handle.create_group("MetalPIE")
        table.attrs["table_type"] = (
            "photoionization_equilibrium_metals"
            if component == "metals" else "photoionization_equilibrium_total"
        )
        table.attrs["schema_version"] = 1
        table.attrs["cooling_units"] = "erg cm^-3 s^-1"
        table.attrs["heating_units"] = "erg cm^-3 s^-1"
        table.attrs["abundance_reference"] = "solar"
        table.attrs["cloudy_version"] = base.cloudy_version_from_path(args.cloudy_exe)
        table.attrs["spectrum_type"] = "Haardt-Madau 2012 UV background"
        table.attrs["radiation_background"] = "table HM12 redshift"
        table.attrs["axis_order"] = "temperature,density,redshift,metallicity"
        table.attrs["pycloudy_version"] = base.pc.__version__
        table.attrs["cloudy_executable"] = str(args.cloudy_exe)
        table.attrs["command_line"] = shlex.join(__import__("sys").argv)
        table.attrs["geometry"] = "plane-parallel slab"
        table.attrs["stop_criterion"] = "zone 1"
        table.attrs["iterate_mode"] = "to convergence"
        table.attrs["temperature_command"] = "constant temperature T K linear"
        table.attrs["molecular_chemistry"] = "disabled (no molecules)"
        table.attrs["component"] = component
        if component == "metals":
            table.attrs["rate_definition"] = (
                "total HM12 Cloudy rate minus H/He-only HM12 Cloudy rate; residuals clipped to zero"
            )
            heating_name = "metal_photoheating_erg_cm3_s"
            cooling_name = "metal_cooling_erg_cm3_s"
        else:
            table.attrs["rate_definition"] = "H/He plus metal volumetric Cloudy rates"
            heating_name = "photoheating_erg_cm3_s"
            cooling_name = "cooling_erg_cm3_s"
        axes = table.create_group("axes")
        axes.create_dataset("log10_temperature_K", data=np.log10(temperatures))
        axes.create_dataset("log10_hydrogen_density_cm-3", data=np.log10(densities))
        axes.create_dataset("redshift", data=redshifts)
        axes.create_dataset("metallicity_Zsun", data=metallicities)
        rates = table.create_group("rates")
        rates.create_dataset(heating_name, data=heating,
                             compression="gzip", compression_opts=4)
        rates.create_dataset(cooling_name, data=cooling,
                             compression="gzip", compression_opts=4)
        cloudy = table.create_group("cloudy")
        cloudy.create_dataset("input_file", data=np.bytes_(input_file or ""))
        cloudy.create_dataset(
            "abundance_file",
            data=np.bytes_("abundances hii region no grains\nmetals = log10(Z/Zsun)\nno molecules\n"),
        )
        cloudy.create_dataset("command", data=np.bytes_(f"{args.cloudy_exe} -p <model>"))
        cloudy.create_dataset("version", data=np.bytes_(base.cloudy_version_from_path(args.cloudy_exe)))


def main():
    args = parse_args()
    if args.workers < 1 or args.iterations < 1:
        raise SystemExit("--workers and --iterations must be positive")
    if args.nT < 1 or args.nnH < 1 or args.n_redshift < 1:
        raise SystemExit("nT, nnH, and n-redshift must be positive")
    if args.redshift_min < 0 or args.redshift_max > 15.93 or args.redshift_min > args.redshift_max:
        raise SystemExit("HM12 redshift must be within 0 <= z <= 15.93")
    if any(value < 0 for value in args.metallicities):
        raise SystemExit("metallicities must be non-negative")
    args.cloudy_exe = args.cloudy_exe or base.default_cloudy_executable()
    if not args.cloudy_exe.exists() or not base.executable_matches_host(args.cloudy_exe):
        raise SystemExit(f"Cloudy executable is missing or incompatible: {args.cloudy_exe}")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.data_dir.mkdir(parents=True, exist_ok=True)
    args.checkpoint = args.data_dir / f"{args.output.stem}.checkpoint.csv"
    if not args.resume and args.checkpoint.exists():
        args.checkpoint.unlink()
    temperatures = np.logspace(args.logT_min, args.logT_max, args.nT)
    densities = np.logspace(args.lognH_min, args.lognH_max, args.nnH)
    redshifts = np.linspace(args.redshift_min, args.redshift_max, args.n_redshift)
    metallicities = np.asarray(args.metallicities, dtype=float)
    configure_base(args)
    print("HM12 photoionization cooling-table generation")
    print("===============================================")
    print(f"Cloudy executable = {args.cloudy_exe}")
    print(f"Grid shape = ({len(metallicities)}, {len(temperatures)}, {len(densities)}, {len(redshifts)})")
    cooling, heating, total_cooling, total_heating, input_file = compute_grid(
        temperatures, densities, redshifts, metallicities, args
    )
    stem = args.data_dir / args.output.stem
    for z_index, metallicity in enumerate(metallicities):
        label = base.metallicity_label(metallicity)
        metal_output = stem.with_name(stem.name + f"_Z{label}_metals.h5")
        total_output = stem.with_name(stem.name + f"_Z{label}_total.h5")
        if not args.overwrite and (metal_output.exists() or total_output.exists()):
            raise SystemExit(
                f"ERROR: output exists; use --overwrite: {metal_output}, {total_output}"
            )
        write_hm12_table(
            metal_output, temperatures, densities, redshifts,
            metallicities[z_index:z_index + 1],
            cooling[z_index:z_index + 1], heating[z_index:z_index + 1],
            args, input_file, "metals",
        )
        write_hm12_table(
            total_output, temperatures, densities, redshifts,
            metallicities[z_index:z_index + 1],
            total_cooling[z_index:z_index + 1], total_heating[z_index:z_index + 1],
            args, input_file, "hydrogen+helium+metals",
        )
        print(f"Wrote: {metal_output}")
        print(f"Wrote: {total_output}")


if __name__ == "__main__":
    main()
