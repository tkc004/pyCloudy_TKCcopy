#!/usr/bin/env python
# coding: utf-8

"""Quickstart example converted from the README.

The Cloudy executable is resolved relative to the repository so the script can
run without editing paths. Set CLOUDY_EXE in the environment to override it.
"""

from pathlib import Path

import pyCloudy as pc


def main():
    script_dir = Path(__file__).resolve().parent
    emis_tab = [
        "H  1  4861.33A",
        "O  3  5006.84A",
    ]

    cloudy_exe = None
    for base_dir in (script_dir, *script_dir.parents):
        candidate = base_dir / "Cloudy_exe" / "Cloudy" / "c22.02" / "source" / "cloudy.exe"
        if candidate.exists():
            cloudy_exe = candidate
            break
    if cloudy_exe is None:
        raise FileNotFoundError("Could not find Cloudy_exe/Cloudy/c22.02/source/cloudy.exe")

    pc.config.cloudy_exe = str(cloudy_exe)

    model_dir = script_dir / "temp_models"
    model_dir.mkdir(exist_ok=True)
    output_file = script_dir / "output.txt"

    c_input = pc.CloudyInput(str(model_dir / "M17_quickstart"))
    c_input.set_BB(Teff=40000, lumi_unit="q(H)", lumi_value=47)
    c_input.set_cste_density(2.0)
    c_input.set_radius(r_in=17.3)
    c_input.set_abund(predef="ism")
    c_input.set_sphere(True)
    c_input.set_iterate()
    c_input.set_emis_tab(emis_tab)
    c_input.set_distance(dist=1.0, unit="kpc", linear=True)
    c_input.print_input(to_file=True, verbose=False)
    c_input.run_cloudy()

    model = pc.CloudyModel(str(model_dir / "M17_quickstart"))
    lines = [
        f"zones: {model.n_zones}",
        f"te[:5]: {model.te[:5]}",
        f"Hb_SB: {model.get_Hb_SB()}",
    ]
    output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
