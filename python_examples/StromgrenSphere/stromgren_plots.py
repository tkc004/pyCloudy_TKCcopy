from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _example_utils import save_fig
from stromgren_config import write_roman


def save_temperature_profile(mod, figure_dir):
    fig = plt.figure(figsize=(10, 10))
    plt.plot(mod.radius / 3.086e18, np.log10(mod.te), label="Te")
    plt.ylim([3.0, 5.0])
    plt.xlim([0, 10])
    plt.legend(loc=3)
    save_fig(fig, Path(figure_dir) / "temperature_profile.png")


def save_te_ne_scatter(mod, figure_dir):
    fig = plt.figure(figsize=(10, 10))
    plt.scatter(mod.te / 1e3, mod.ne / 1e4, c=mod.depth / np.max(mod.depth), edgecolors="none")
    plt.colorbar()
    plt.xlabel("Te [kK]")
    plt.ylabel(r"Ne [$10^4$ cm$^{-3}$]")
    save_fig(fig, Path(figure_dir) / "te_ne_scatter.png")


def save_continuum(mod, figure_dir):
    fig = plt.figure(figsize=(10, 10))
    plt.loglog(mod.get_cont_x(unit="Ang"), mod.get_cont_y(cont="incid", unit="Jy"), label="Incident")
    plt.loglog(mod.get_cont_x(unit="Ang"), mod.get_cont_y(cont="diffout", unit="Jy"), label="Diff Out")
    plt.loglog(mod.get_cont_x(unit="Ang"), mod.get_cont_y(cont="ntrans", unit="Jy"), label="Net Trans")
    plt.xlim((100, 100000))
    plt.ylim((1e-9, 1e1))
    plt.xlabel("Angstrom")
    plt.ylabel("Jy")
    plt.legend(loc=4)
    save_fig(fig, Path(figure_dir) / "continuum.png")


def save_ionic_fractions(mod, figure_dir, abund):
    fig = plt.figure(figsize=(10, 10))
    plt.plot(mod.radius / 3.086e18, mod.get_ionic("H", 1), label="H+")
    plt.plot(mod.radius / 3.086e18, mod.get_ionic("H", 0), label="H0")
    plt.plot(mod.radius / 3.086e18, mod.get_ionic("He", 2) * 10 ** abund["He"], label="He++")
    plt.plot(mod.radius / 3.086e18, mod.get_ionic("He", 1) * 10 ** abund["He"], label="He+")
    plt.plot(mod.radius / 3.086e18, mod.get_ionic("He", 0) * 10 ** abund["He"], label="He0")
    plt.yscale("log")
    plt.ylim([1e-4, 1])
    plt.xlim([0, 15.0])
    plt.legend(loc=3)
    save_fig(fig, Path(figure_dir) / "ionic_fractions.png")


def save_nitrogen_ion_fractions(mod, figure_dir):
    element = "N"
    totionic = sum(mod.get_ionic(element, i) for i in range(4))
    fig = plt.figure(figsize=(10, 10))
    for i in range(4):
        plt.plot(mod.radius / 3.086e21, mod.get_ionic(element, i) / totionic, label=element + write_roman(i + 1))
    plt.legend()
    save_fig(fig, Path(figure_dir) / "nitrogen_ion_fractions.png")


def save_all_plots(mod, mode_name, figure_dir, abund):
    save_temperature_profile(mod, figure_dir)
    save_te_ne_scatter(mod, figure_dir)
    save_continuum(mod, figure_dir)
    save_ionic_fractions(mod, figure_dir, abund)
    if mode_name == "HHeZ":
        save_nitrogen_ion_fractions(mod, figure_dir)
