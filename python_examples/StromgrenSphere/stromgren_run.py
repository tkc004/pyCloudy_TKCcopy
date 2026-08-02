from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pyCloudy as pc

from _example_utils import find_cloudy_exe
from stromgren_config import ELESTATE, load_common_config, load_mode_config, mode_file_specs, resolve_abund, write_roman


def save_profile_tables(mod, mode_name, data_dir):
    spec = mode_file_specs(mode_name)
    x = mod.radius / 3.086e21
    for ele in spec["elements"]:
        for state in range(ELESTATE[ele]):
            romanstate = write_roman(state + 1)
            y = np.log10(mod.get_ionic(ele, state))
            xy = np.vstack((x, y)).T
            filename = spec["x_template"].format(ele=ele, roman=romanstate, mode=mode_name)
            np.savetxt(Path(data_dir) / filename, xy, delimiter=",")
    y = np.log10(mod.te)
    xy = np.vstack((x, y)).T
    np.savetxt(Path(data_dir) / spec["t_template"].format(mode=mode_name), xy, delimiter=",")


def run_model(mode_name, script_dir):
    script_dir = Path(script_dir)
    data_dir = script_dir / "data"
    figure_dir = script_dir / "figures"
    model_dir = script_dir / "temp_models"
    options_dir = script_dir / "options"
    data_dir.mkdir(exist_ok=True)
    figure_dir.mkdir(exist_ok=True)
    model_dir.mkdir(exist_ok=True)

    pc.log_.level = 3

    common_config = load_common_config(options_dir)
    mode_config = load_mode_config(options_dir, mode_name)

    model_name = f"model_{mode_name}"
    full_model_name = str(model_dir / model_name)
    dens = float(common_config["dens"])
    Teff = float(common_config["Teff"])
    qH = float(common_config["qH"])
    r_min = float(common_config["r_min"])
    dist = float(common_config["dist"])
    emis_tab = common_config["emis_tab"]
    options = tuple(mode_config["options"])
    abund = resolve_abund(mode_config)

    c_input = pc.CloudyInput(full_model_name)
    c_input.set_BB(Teff=Teff, lumi_unit="q(H)", lumi_value=qH)
    c_input.set_cste_density(dens)
    c_input.set_radius(r_in=np.log10(r_min))
    c_input.set_abund(ab_dict=abund, nograins=True)
    c_input.set_other(options)
    c_input.set_iterate(5)
    c_input.set_sphere()
    c_input.set_emis_tab(emis_tab)
    c_input.set_distance(dist=dist, unit="kpc", linear=True)
    c_input.print_input(to_file=True, verbose=False)

    pc.log_.message(f"Running {model_name}", calling="test1")
    pc.config.cloudy_exe = str(find_cloudy_exe(script_dir))
    pc.log_.timer("Starting Cloudy", quiet=True, calling="test1")
    c_input.run_cloudy()
    pc.log_.timer("Cloudy ended after seconds:", calling="test1")

    mod = pc.CloudyModel(full_model_name)

    dir(mod)
    mod.n_ions
    mod.print_stats()
    mod.print_lines()
    mod.get_ab_ion_vol_ne("H", 1)
    mod.get_T0_ion_vol_ne("O", 2)
    mod.log_U_mean
    mod.log_U_mean_ne
    print("T0 = {0:7.1f}K, t2 = {1:6.4f}".format(mod.T0, mod.t2))
    print(
        "Hbeta Equivalent width = {0:6.1f}, Hbeta Surface Brightness = {1:4.2e}".format(
            mod.get_Hb_EW(), mod.get_Hb_SB()
        )
    )

    for line in mod.emis_labels:
        print("{0} {1:10.3e} {2:7.2f}".format(line, mod.get_emis_vol(line), mod.get_emis_vol(line) / mod.get_emis_vol("H__1_486133A") * 100.0))

    save_profile_tables(mod, mode_name, data_dir)
    context = {
        "script_dir": script_dir,
        "data_dir": data_dir,
        "figure_dir": figure_dir,
        "model_dir": model_dir,
        "options_dir": options_dir,
        "full_model_name": full_model_name,
        "model_name": model_name,
    }
    return mod, abund, context
