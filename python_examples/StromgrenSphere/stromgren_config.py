from collections import OrderedDict
from pathlib import Path

import numpy as np

from _example_utils import load_yaml


AMU = {
    "He": 4.0,
    "C": 12.011,
    "N": 14.0067,
    "O": 15.9994,
    "Ne": 20.179,
    "S": 32.06,
    "Ar": 39.948,
    "Fe": 55.847,
    "Cl": 35.453,
    "Mg": 24.305,
}

ELESTATE = {"H": 2, "He": 3, "O": 4, "N": 4}
DEFAULT_MODE = "HHe_fixedT"


def write_roman(num):
    roman = OrderedDict()
    roman[1000] = "M"
    roman[900] = "CM"
    roman[500] = "D"
    roman[400] = "CD"
    roman[100] = "C"
    roman[90] = "XC"
    roman[50] = "L"
    roman[40] = "XL"
    roman[10] = "X"
    roman[9] = "IX"
    roman[5] = "V"
    roman[4] = "IV"
    roman[1] = "I"

    def roman_num(value):
        for r in roman.keys():
            x, y = divmod(value, r)
            yield roman[r] * x
            value -= r * x
            if value <= 0:
                break

    return "".join(a for a in roman_num(num))


def load_common_config(options_dir):
    return load_yaml(Path(options_dir) / "common.yml")


def load_mode_config(options_dir, mode_name):
    config_path = Path(options_dir) / f"{mode_name}.yml"
    if not config_path.exists():
        raise FileNotFoundError(f"Missing StromgrenSphere mode config: {config_path}")
    return load_yaml(config_path)


def resolve_abund(config):
    if "abund" in config and config["abund"] is not None:
        return config["abund"]
    if "massfrac" in config and config["massfrac"] is not None:
        massfrac = config["massfrac"]
        if "H" not in massfrac:
            raise KeyError("massfrac config must include hydrogen")
        abund = {}
        for elem, atomic_weight in AMU.items():
            if elem not in massfrac:
                continue
            temp_abund = massfrac[elem] / massfrac["H"] / atomic_weight
            abund[elem] = np.log10(temp_abund) if temp_abund > 1e-12 else -12.0
        return abund
    raise KeyError("Mode config must define either 'abund' or 'massfrac'")


def mode_file_specs(mode_name):
    specs = {
        "H": {
            "elements": ["H"],
            "x_template": "x{ele}{roman}_Cloudy_Stromgren_real_{mode}.txt",
            "t_template": "T_Cloudy_Stromgren_real_{mode}.txt",
        },
        "H_fixedT": {
            "elements": ["H"],
            "x_template": "x{ele}{roman}_Cloudy_Stromgren_real_{mode}_T4.txt",
            "t_template": "T_Cloudy_Stromgren_{mode}_T4.txt",
        },
        "HHe": {
            "elements": ["H", "He"],
            "x_template": "x{ele}{roman}_Cloudy_Stromgren_real_{mode}.txt",
            "t_template": "T_Cloudy_Stromgren_real_{mode}.txt",
        },
        "HHe_fixedT": {
            "elements": ["H", "He"],
            "x_template": "x{ele}{roman}_Cloudy_Stromgren_real_{mode}_T4.txt",
            "t_template": "T_Cloudy_Stromgren_real_{mode}_T4.txt",
        },
        "HHe_fixedT_caseB": {
            "elements": ["H", "He"],
            "x_template": "x{ele}{roman}_Cloudy_Stromgren_real_{mode}_T4_caseB.txt",
            "t_template": "T_Cloudy_Stromgren_real_{mode}_T4_caseB.txt",
        },
        "HHeZ": {
            "elements": ["H", "He", "O", "N"],
            "x_template": "x{ele}{roman}_Cloudy_Stromgren_real_{mode}.txt",
            "t_template": "T_Cloudy_Stromgren_real_{mode}.txt",
        },
    }
    if mode_name not in specs:
        raise ValueError(f"Unsupported mode: {mode_name}")
    return specs[mode_name]
