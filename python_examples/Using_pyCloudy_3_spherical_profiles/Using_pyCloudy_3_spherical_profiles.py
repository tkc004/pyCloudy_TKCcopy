#!/usr/bin/env python
# coding: utf-8
"""Spherical pyCloudy model driven by an external radial profile file.

The CSV file uses radius in pc, velocity in km/s, and hydrogen density in atom
cm^-3. Cloudy receives the density profile through ``dlaw table radius``; C3D
receives the velocity profile through its user velocity function.
"""

import os
import sys
import warnings
from pathlib import Path

import numpy as np

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parents[0]))
from _example_utils import find_cloudy_exe, save_fig

temp_model_dir = script_dir / "temp_models"
temp_model_dir.mkdir(exist_ok=True)
fig_dir = script_dir / "figures"
fig_dir.mkdir(exist_ok=True)
mpl_config_dir = temp_model_dir / ".mplconfig"
cache_dir = temp_model_dir / ".cache"
mpl_config_dir.mkdir(exist_ok=True)
cache_dir.mkdir(exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(mpl_config_dir)
os.environ["XDG_CACHE_HOME"] = str(cache_dir)

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import pyCloudy as pc


PROFILE_FILE = script_dir / "radial_profiles.csv"

MODEL_NAME = "M3D_spherical_profiles"
DIM = 101
PROJ_AXIS = 0
EMIS_TAB = [
    "H  1  4861.32A", "H  1  6562.80A", "Ca B  5875.64A",
    "N  2  6583.45A", "O  1  6300.30A", "O  2  3726.03A",
    "O  2  3728.81A", "O  3  5006.84A", "O  3  4363.21A",
    "O 3R  4363.00A", "O 3C  4363.00A", "S  2  6716.44A",
    "S  2  6730.82A", "Cl 3  5517.71A", "Cl 3  5537.87A",
    "O  1  63.1679m", "O  1  145.495m", "C  2  157.636m",
]


def validate_profiles(radius_pc, velocity_kms, density_cm3):
    radius_pc = np.asarray(radius_pc, dtype=float)
    velocity_kms = np.asarray(velocity_kms, dtype=float)
    density_cm3 = np.asarray(density_cm3, dtype=float)
    if not (radius_pc.ndim == velocity_kms.ndim == density_cm3.ndim == 1):
        raise ValueError("Radial profiles must be one-dimensional arrays")
    if not (len(radius_pc) == len(velocity_kms) == len(density_cm3)):
        raise ValueError("Radial profiles must have the same number of points")
    if len(radius_pc) < 2:
        raise ValueError("At least two radial profile points are required")
    if np.any(radius_pc <= 0) or np.any(density_cm3 <= 0):
        raise ValueError("Radii and hydrogen densities must be positive")
    if np.any(np.diff(radius_pc) <= 0):
        raise ValueError("Radius values must be strictly increasing")
    return radius_pc, velocity_kms, density_cm3


def read_profiles(path):
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    required_columns = ("RADIUS_PC", "VELOCITY_KMS", "DENSITY_CM3")
    if data.dtype.names is None or any(column not in data.dtype.names for column in required_columns):
        raise ValueError(
            "Profile file must have header: RADIUS_PC,VELOCITY_KMS,DENSITY_CM3"
        )
    return validate_profiles(
        data["RADIUS_PC"],
        data["VELOCITY_KMS"],
        data["DENSITY_CM3"],
    )


def plot_input_profiles(
    radius_pc,
    velocity_kms,
    density_cm3,
    title="Input spherical radial profiles",
    filename="input_radial_profiles.png",
    velocity_label="Radial velocity",
):
    figure, density_axis = plt.subplots(figsize=(10, 6))
    velocity_axis = density_axis.twinx()
    density_line = density_axis.plot(
        radius_pc,
        density_cm3,
        "o-",
        color="tab:blue",
        label="Hydrogen density",
    )
    velocity_line = velocity_axis.plot(
        radius_pc,
        velocity_kms,
        "s--",
        color="tab:orange",
        label=velocity_label,
    )
    density_axis.set_xlabel("Radius [pc]")
    density_axis.set_ylabel(r"Hydrogen density [cm$^{-3}$]", color="tab:blue")
    velocity_axis.set_ylabel("Radial velocity [km/s]", color="tab:orange")
    density_axis.set_ylim(bottom=0.0)
    velocity_axis.set_ylim(bottom=0.0)
    density_axis.tick_params(axis="y", labelcolor="tab:blue")
    velocity_axis.tick_params(axis="y", labelcolor="tab:orange")
    lines = density_line + velocity_line
    density_axis.legend(lines, [line.get_label() for line in lines], loc="best")
    density_axis.set_title(title)
    figure.tight_layout()
    save_fig(figure, fig_dir / filename)


def polynomial_velocity_profile(radius_pc, params):
    """Evaluate C3D's polynomial radial velocity law in km/s."""
    radius_pc = np.asarray(radius_pc, dtype=float)
    normalized_radius = radius_pc / radius_pc[-1]
    velocity_kms = sum(
        parameter * normalized_radius**power
        for power, parameter in enumerate(params)
    )
    # C3D explicitly sets the velocity vector to zero at the origin.
    velocity_kms[radius_pc == 0.0] = 0.0
    return velocity_kms


def density_table_commands(radius_pc, density_cm3):
    """Return Cloudy commands for a log-radius/log-density table."""
    # Extend the table slightly so Cloudy never evaluates exactly outside it.
    table_radius = np.r_[radius_pc[0] * 1.0e-6, radius_pc, radius_pc[-1] * 1.0e6]
    table_density = np.r_[density_cm3[0], density_cm3, density_cm3[-1]]
    pairs = [
        (np.log10(radius * pc.CST.PC), np.log10(density))
        for radius, density in zip(table_radius, table_density)
    ]
    commands = ["dlaw table radius"]
    commands.extend("continue {0:.8f} {1:.8f}".format(*pair) for pair in pairs)
    commands.append("end of dlaw")
    return commands


def build_model(model_path, radius_pc, density_cm3):
    model = pc.CloudyInput(str(model_path))
    model.set_BB(80000.0, "q(H)", 47.3)
    model.set_grains()
    model.set_radius(r_in=np.log10(radius_pc[0] * pc.CST.PC))
    model.set_stop("radius {0:.8f}".format(np.log10(radius_pc[-1] * pc.CST.PC)))
    table_commands = density_table_commands(radius_pc, density_cm3)
    model.set_dlaw("table radius")
    # Cloudy requires table continuation rows immediately after the dlaw line.
    model._filling_factor = None
    model.set_other(table_commands[1:])
    model.set_emis_tab(EMIS_TAB)
    model.set_sphere()
    model.print_input(to_file=True, verbose=False)


def set_user_velocity(m3d, radius_pc, velocity_kms):
    def velocity_function(params):
        profile_radius, profile_velocity, cub_coord = params
        radius = cub_coord.r / pc.CST.PC
        speed = np.interp(
            radius,
            profile_radius,
            profile_velocity,
            left=0.0,
            right=profile_velocity[-1],
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            direction = np.divide(
                1.0,
                cub_coord.r,
                out=np.zeros_like(cub_coord.r),
                where=cub_coord.r != 0.0,
            )
        return speed * cub_coord.x * direction, speed * cub_coord.y * direction, speed * cub_coord.z * direction

    m3d.set_velocity(
        velocity_law="user",
        params=[radius_pc, velocity_kms, m3d.cub_coord],
        user_function=velocity_function,
    )
    m3d.config_profile(size_spectrum=41, vel_max=max(25.0, float(np.max(np.abs(velocity_kms))) * 1.2), v_turb=0.01)


def plot_profiles(m3d, x_pos, y_pos, title):
    axis = plt.gca()
    plt.plot(m3d.vel_tab, m3d.get_profile("H__1_486132A", axis="x")[:, x_pos, y_pos] * 5, label=r"H$\beta$")
    plt.plot(m3d.vel_tab, m3d.get_profile("N__2_658345A", axis="x")[:, x_pos, y_pos] * 5, label=r"[NII]$\lambda$6584")
    plt.plot(m3d.vel_tab, m3d.get_profile("O__3_500684A", axis="x")[:, x_pos, y_pos], label=r"[OIII]$\lambda$5007")
    axis.set_title(title)
    axis.set_xlabel("Velocity [km/s]")
    axis.set_ylabel(r"Scaled emissivity [erg s$^{-1}$ cm$^{-3}$]")
    plt.legend()


def safe_divide(num, den):
    return np.divide(num, den, out=np.zeros_like(num), where=den != 0.0)


def show_image(data, title, colorbar_label, **kwargs):
    image = plt.imshow(data, **kwargs)
    plt.title(title)
    colorbar = plt.colorbar(image)
    colorbar.set_label(colorbar_label)
    return image


def add_rgb_colorbars(figure):
    figure.subplots_adjust(bottom=0.22)
    channels = (
        ("[NII] 6584 (red channel)", "black", "red"),
        ("[OIII] 5007 (green channel)", "black", "lime"),
        (r"H$\beta$ 4861 (blue channel)", "black", "blue"),
    )
    for index, (label, low_color, high_color) in enumerate(channels):
        bar_y = 0.14 - index * 0.055
        colorbar_axis = figure.add_axes((0.30, bar_y, 0.60, 0.025))
        colormap = LinearSegmentedColormap.from_list(label, [low_color, high_color])
        mappable = ScalarMappable(norm=Normalize(0, 255), cmap=colormap)
        mappable.set_array([])
        colorbar = figure.colorbar(mappable, cax=colorbar_axis, orientation="horizontal")
        figure.text(0.28, bar_y + 0.0125, label, ha="right", va="center", fontsize=10)


def show_rgb_with_colorbars(image, extent):
    figure, axis = plt.subplots(figsize=(10, 12))
    axis.imshow(image, extent=extent, origin="lower", interpolation="nearest", vmin=0, vmax=255)
    axis.set_title("RGB emission image")
    axis.set_xlabel("Projected x [pc]")
    axis.set_ylabel("Projected z [pc]")
    add_rgb_colorbars(figure)
    return figure


def add_profiles_to_rgb(m3d, image_axis, ref, nx=20, ny=20):
    """Overlay profile panels using the RGB axes position as the coordinate frame."""
    profiles = m3d.get_profile(ref, axis="x")
    size_x, size_y = profiles.shape[1:]
    sx = int(size_x / nx)
    sy = int(size_y / ny)
    dx = int((size_x - sx * nx) / 2)
    dy = int((size_y - sy * ny) / 2)
    figure = image_axis.figure
    figure.canvas.draw()
    image_bbox = image_axis.get_position()
    for ix in range(nx):
        for iy in range(ny):
            profile = profiles[
                :, dx + ix * sx:dx + (ix + 1) * sx,
                dy + iy * sy:dy + (iy + 1) * sy,
            ].sum(axis=1).sum(axis=1)
            if profile.sum() <= 0.0:
                continue
            profile /= np.max(profile)
            panel_position = (
                image_bbox.x0 + (dy + iy * sy) * image_bbox.width / size_y,
                image_bbox.y0 + (dx + ix * sx) * image_bbox.height / size_x,
                sy * image_bbox.width / size_y,
                sx * image_bbox.height / size_x,
            )
            profile_axis = figure.add_axes(panel_position)
            profile_axis.plot(profile, color="yellow")
            profile_axis.set_ylim(0.0, 1.05)
            profile_axis.set_xticks([])
            profile_axis.set_yticks([])
            profile_axis.patch.set_alpha(0.0)


def other_plots(m3d, proj_axis, n_cut):
    plt.subplot(331)
    show_image(
        m3d.get_emis("H__1_486132A").sum(axis=proj_axis) * m3d.cub_coord.cell_size,
        "Hb",
        r"Emissivity [erg s$^{-1}$ cm$^{-2}$]",
    )

    plt.subplot(332)
    show_image(
        m3d.get_emis("N__2_658345A").sum(axis=proj_axis) * m3d.cub_coord.cell_size,
        "[NII]",
        r"Emissivity [erg s$^{-1}$ cm$^{-2}$]",
    )

    plt.subplot(333)
    show_image(
        m3d.get_emis("O__3_500684A").sum(axis=proj_axis) * m3d.cub_coord.cell_size,
        "[OIII]",
        r"Emissivity [erg s$^{-1}$ cm$^{-2}$]",
    )

    hb = m3d.get_emis("H__1_486132A").sum(axis=proj_axis)
    nii = m3d.get_emis("N__2_658345A").sum(axis=proj_axis)
    oiii = m3d.get_emis("O__3_500684A").sum(axis=proj_axis)
    plt.subplot(334)
    show_image(safe_divide(nii, hb), "[NII]/Hb", "Ratio [dimensionless]")
    plt.subplot(335)
    show_image(safe_divide(oiii, hb), "[OIII]/Hb", "Ratio [dimensionless]")

    plt.subplot(336)
    show_image(m3d.get_ionic("O", 1)[n_cut, :, :], "O+ cut", "Ionic fraction [dimensionless]")

    oplus = m3d.get_ionic("O", 1)
    nplus = m3d.get_ionic("N", 1)
    ratio = safe_divide(nplus, oplus)
    plt.subplot(337)
    plt.scatter(oplus.ravel(), ratio.ravel(), c=np.abs(m3d.cub_coord.theta.ravel()), edgecolors="none")
    plt.title("Colored by |Theta|")
    plt.xlabel("O+ / O")
    plt.ylabel("N+/O+ / N/O")
    plt.colorbar(label=r"$|\Theta|$ [deg]")

    plt.subplot(338)
    plt.scatter(oplus.ravel(), ratio.ravel(), c=m3d.relative_depth.ravel(), vmin=0, vmax=1, edgecolors="none")
    plt.title("Colored by position in the nebula")
    plt.xlabel("O+ / O")
    plt.ylabel("N+/O+ / N/O")
    plt.colorbar(label="Relative depth [dimensionless]")

    n2 = m3d.get_ionic("N", 2)
    plt.subplot(339)
    weighted_map = safe_divide((ratio * n2).sum(axis=proj_axis), n2.sum(axis=proj_axis))
    show_image(weighted_map, "N+/O+ / N/O weighted by NII", "Ratio [dimensionless]")
    plt.contour(weighted_map, levels=[1.0])


radius_pc, velocity_kms, density_cm3 = read_profiles(PROFILE_FILE)
plot_input_profiles(radius_pc, velocity_kms, density_cm3)
poly_params = [20.0, 60.0]
poly_velocity_kms = polynomial_velocity_profile(radius_pc, poly_params)
plot_input_profiles(
    radius_pc,
    poly_velocity_kms,
    density_cm3,
    title="Input density with C3D polynomial velocity (params=[20, 60])",
    filename="input_radial_profiles_wpolyv.png",
    velocity_label="C3D polynomial velocity",
)
cloudy_exe = find_cloudy_exe(script_dir)
pc.config.cloudy_exe = str(cloudy_exe)
model_path = temp_model_dir / MODEL_NAME
build_model(model_path, radius_pc, density_cm3)
pc.print_make_file(dir_=str(temp_model_dir) + "/")
pc.run_cloudy(dir_=str(temp_model_dir) + "/", n_proc=6, model_name=MODEL_NAME, use_make=True)
c_output = pc.CloudyModel(str(model_path))
m3d = pc.C3D(c_output, dims=DIM, center=True, n_dim=1)

m3d.set_velocity(params=poly_params)
m3d.config_profile(size_spectrum=51, vel_max=50.0, v_turb=0.01)
n_cut = (DIM - 1) // 2

plt.figure(figsize=(10, 10))
plot_profiles(
    m3d,
    n_cut,
    n_cut,
    "Line profiles: default polynomial velocity law (params=[20, 60])",
)
save_fig(plt.gcf(), fig_dir / "profile_default.png")

set_user_velocity(m3d, radius_pc, velocity_kms)

plt.figure(figsize=(10, 10))
plot_profiles(
    m3d,
    n_cut,
    n_cut,
    "Line profiles: user velocity profile from radial_profiles.csv",
)
save_fig(plt.gcf(), fig_dir / "profile_user_velocity.png")

plt.figure(figsize=(15, 15))
other_plots(m3d, PROJ_AXIS, n_cut)
save_fig(plt.gcf(), fig_dir / "derived_maps.png")

image = m3d.get_RGB(list_emis=["N__2_658345A", "O__3_500684A", "H__1_486132A"])
rgb_extent = [
    m3d.cub_coord.x_vec[0] / pc.CST.PC,
    m3d.cub_coord.x_vec[-1] / pc.CST.PC,
    m3d.cub_coord.z_vec[0] / pc.CST.PC,
    m3d.cub_coord.z_vec[-1] / pc.CST.PC,
]
save_fig(show_rgb_with_colorbars(image, rgb_extent), fig_dir / "rgb_compact.png")

rgb_with_profiles_figure, rgb_axis = plt.subplots(figsize=(15, 15))
rgb_image = rgb_axis.imshow(
    image,
    extent=rgb_extent,
    origin="lower",
    interpolation="nearest",
    vmin=0,
    vmax=255,
)
rgb_axis.set_title("RGB emission image with profiles")
rgb_axis.set_xlabel("Projected x [pc]")
rgb_axis.set_ylabel("Projected z [pc]")
add_rgb_colorbars(rgb_with_profiles_figure)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    add_profiles_to_rgb(m3d, rgb_axis, ref=3, nx=20, ny=20)
save_fig(rgb_with_profiles_figure, fig_dir / "rgb_with_profiles.png")

f, ax = plt.subplots()
n2map = m3d.get_emis("N__2_658345A").sum(axis=PROJ_AXIS)
hbmap = m3d.get_emis("H__1_486132A").sum(axis=PROJ_AXIS)
o3map = m3d.get_emis("O__3_500684A").sum(axis=PROJ_AXIS)
mask = np.logical_and.reduce([mapl > 0.01 * mapl.max() for mapl in (hbmap, o3map, n2map)])
ax.scatter(np.log10(safe_divide(n2map, hbmap)[mask]), np.log10(safe_divide(o3map, hbmap)[mask]))
ax.set_xlabel("log10([NII]/Hb)")
ax.set_ylabel("log10([OIII]/Hb)")
save_fig(f, fig_dir / "diagnostic_scatter.png")
