#!/usr/bin/env python
# coding: utf-8
"""Spherical pyCloudy model driven by an external RHD wind profile file.

The CSV file uses radius in pc, velocity in km/s, hydrogen density in atom
cm^-3, and temperature in K. Cloudy receives radial density and temperature
profiles through ``dlaw table radius`` and ``tlaw table radius``; C3D receives
the velocity profile through its user velocity function.
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np

script_dir = Path(__file__).resolve().parent
# Import the repository checkout rather than an older installed pyCloudy.
sys.path.insert(0, str(script_dir.parents[1]))
sys.path.insert(1, str(script_dir.parents[0]))
from _example_utils import find_cloudy_exe, save_fig

base_temp_model_dir = script_dir / "temp_models"
base_fig_dir = script_dir / "figures"
temp_model_dir = base_temp_model_dir
fig_dir = base_fig_dir
base_temp_model_dir.mkdir(exist_ok=True)
base_fig_dir.mkdir(exist_ok=True)
base_mpl_config_dir = base_temp_model_dir / ".mplconfig"
base_cache_dir = base_temp_model_dir / ".cache"
base_mpl_config_dir.mkdir(exist_ok=True)
base_cache_dir.mkdir(exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(base_mpl_config_dir)
os.environ["XDG_CACHE_HOME"] = str(base_cache_dir)

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import pyCloudy as pc


PROFILE_DIR = script_dir / "radial_profiles"

MODEL_NAME_PREFIX = "M3D_spherical_profiles_rhd_wind_temperature"
DIM = 151
PROJ_AXIS = 0
PROFILE_DISPLAY_LIMITS = (-40.0, 40.0)
XRAY_FILE_SUFFIX = ".xray.cont"
XRAY_MIN_KEV = 0.1
XRAY_MAX_KEV = 10.0
XRAY_ZONE_INTERVAL = 10
EMIS_TAB = [
    "H  1  4861.32A", "H  1  6562.80A", "Ca B  5875.64A",
    "N  2  6583.45A", "O  1  6300.30A", "O  2  3726.03A",
    "O  2  3728.81A", "O  3  5006.84A", "O  3  4363.21A",
    "O 3R  4363.00A", "O 3C  4363.00A", "S  2  6716.44A",
    "S  2  6730.82A", "Cl 3  5517.71A", "Cl 3  5537.87A",
    "O  1  63.1679m", "O  1  145.495m", "C  2  157.636m",
]


def configure_output_dirs(profile_file):
    """Create isolated model, cache, and figure directories for one CSV."""
    global temp_model_dir, fig_dir
    profile_stem = profile_file.stem
    temp_model_dir = base_temp_model_dir / profile_stem
    fig_dir = base_fig_dir / profile_stem
    mpl_config_dir = temp_model_dir / ".mplconfig"
    cache_dir = temp_model_dir / ".cache"
    for directory in (temp_model_dir, fig_dir, mpl_config_dir, cache_dir):
        directory.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_config_dir)
    os.environ["XDG_CACHE_HOME"] = str(cache_dir)


def validate_profiles(radius_pc, velocity_kms, density_cm3, temperature_k):
    radius_pc = np.asarray(radius_pc, dtype=float)
    velocity_kms = np.asarray(velocity_kms, dtype=float)
    density_cm3 = np.asarray(density_cm3, dtype=float)
    temperature_k = np.asarray(temperature_k, dtype=float)
    if not (radius_pc.ndim == velocity_kms.ndim == density_cm3.ndim == temperature_k.ndim == 1):
        raise ValueError("Radial profiles must be one-dimensional arrays")
    if not (len(radius_pc) == len(velocity_kms) == len(density_cm3) == len(temperature_k)):
        raise ValueError("Radial profiles must have the same number of points")
    if len(radius_pc) < 2:
        raise ValueError("At least two radial profile points are required")
    if np.any(radius_pc <= 0) or np.any(density_cm3 <= 0) or np.any(temperature_k <= 0):
        raise ValueError("Radii, hydrogen densities, and temperatures must be positive")
    if np.any(np.diff(radius_pc) <= 0):
        raise ValueError("Radius values must be strictly increasing")
    return radius_pc, velocity_kms, density_cm3, temperature_k


def read_profiles(path):
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    required_columns = ("RADIUS_PC", "VELOCITY_KMS", "DENSITY_CM3", "TEMP_K")
    if data.dtype.names is None or any(column not in data.dtype.names for column in required_columns):
        raise ValueError(
            "Profile file must have header: RADIUS_PC,VELOCITY_KMS,DENSITY_CM3,TEMP_K"
        )
    return validate_profiles(
        data["RADIUS_PC"],
        data["VELOCITY_KMS"],
        data["DENSITY_CM3"],
        data["TEMP_K"],
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
        np.where(velocity_kms > 0.0, velocity_kms, np.nan),
        "s--",
        color="tab:orange",
        label=velocity_label,
    )
    density_axis.set_xlabel("Radius [pc]")
    density_axis.set_ylabel(r"Hydrogen density [cm$^{-3}$]", color="tab:blue")
    velocity_axis.set_ylabel("Radial velocity [km/s]", color="tab:orange")
    density_axis.set_ylim(bottom=0.0)
    velocity_axis.set_yscale("log")
    velocity_axis.set_ylim(bottom=0.5)
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


def temperature_table_commands(radius_pc, temperature_k):
    """Return Cloudy commands using log10 radius in cm and log10 temperature."""
    table_radius = np.r_[radius_pc[0] * 1.0e-6, radius_pc, radius_pc[-1] * 1.0e6]
    # Cloudy's tlaw table stores log10 temperature, while TEMP_K is in K.
    table_temperature = np.log10(
        np.r_[temperature_k[0], temperature_k, temperature_k[-1]]
    )
    pairs = [
        (np.log10(radius * pc.CST.PC), temperature)
        for radius, temperature in zip(table_radius, table_temperature)
    ]
    commands = ["tlaw table radius"]
    commands.extend("continue {0:.8f} {1:.8f}".format(*pair) for pair in pairs)
    commands.append("end of tlaw")
    return commands


def build_model(model_path, radius_pc, density_cm3, temperature_k, include_xray=False):
    model = pc.CloudyInput(str(model_path))
    model.set_BB(80000.0, "q(H)", 49.0)
    model.set_grains()
    model.set_radius(r_in=np.log10(radius_pc[0] * pc.CST.PC))
    model.set_stop("radius {0:.8f}".format(np.log10(radius_pc[-1] * pc.CST.PC)))
    table_commands = density_table_commands(radius_pc, density_cm3)
    model.set_dlaw("table radius")
    # Cloudy requires table continuation rows immediately after the dlaw line.
    model._filling_factor = None
    model.set_other(table_commands[1:] + temperature_table_commands(radius_pc, temperature_k))
    if include_xray:
        model.set_other(
            'save continuum "{0}" every {1}'.format(
                XRAY_FILE_SUFFIX,
                XRAY_ZONE_INTERVAL,
            )
        )
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
    # Resolve the velocity window shown in the user-profile figures rather
    # than spreading a small number of bins over the full wind speed range.
    m3d.config_profile(size_spectrum=1001, vel_max=100.0, v_turb=0.01)


def plot_profiles(m3d, x_pos, y_pos, title, velocity_limits=None):
    axis = plt.gca()
    plt.plot(m3d.vel_tab, m3d.get_profile("H__1_486132A", axis="x")[:, x_pos, y_pos] * 5, label=r"H$\beta$")
    plt.plot(m3d.vel_tab, m3d.get_profile("N__2_658345A", axis="x")[:, x_pos, y_pos] * 5, label=r"[NII]$\lambda$6584")
    plt.plot(m3d.vel_tab, m3d.get_profile("O__3_500684A", axis="x")[:, x_pos, y_pos], label=r"[OIII]$\lambda$5007")
    axis.set_title(title)
    axis.set_xlabel("Velocity [km/s]")
    axis.set_ylabel(r"Scaled emissivity [erg s$^{-1}$ cm$^{-3}$]")
    if velocity_limits is not None:
        axis.set_xlim(*velocity_limits)
    plt.legend()


def plot_xray_luminosity_profile(c_output, continuum_path):
    """Plot the radial 0.1--10 keV luminosity from Cloudy's per-zone continuum."""
    continuum = np.loadtxt(continuum_path, comments="#", usecols=(0, 6))
    n_zones = len(c_output.radius)
    n_points = len(c_output.get_cont_x())
    if continuum.shape[0] % n_points != 0:
        raise ValueError("Per-zone X-ray continuum output has an incomplete spectrum")
    n_saved_zones = continuum.shape[0] // n_points
    continuum = continuum.reshape(n_saved_zones, n_points, 2)
    zone_indices = np.arange(0, n_zones, XRAY_ZONE_INTERVAL)[:n_saved_zones]
    if len(zone_indices) != n_saved_zones:
        raise ValueError("Cloudy X-ray zone sampling exceeds the model zones")
    energy_kev = continuum[0, :, 0] * pc.CST.RYD_EV / 1000.0
    xray_mask = (energy_kev >= XRAY_MIN_KEV) & (energy_kev <= XRAY_MAX_KEV)
    if not np.any(xray_mask):
        raise ValueError("Cloudy continuum does not contain the requested X-ray band")
    luminosity = np.trapz(
        continuum[:, xray_mask, 1],
        x=np.log(energy_kev[xray_mask]),
        axis=1,
    )
    positive = luminosity > 0.0
    figure, axis = plt.subplots(figsize=(10, 6))
    axis.plot(
        c_output.radius[zone_indices][positive] / pc.CST.PC,
        luminosity[positive],
        color="tab:purple",
    )
    axis.set_title("Cloudy radial X-ray luminosity (0.1--10 keV)")
    axis.set_xlabel("Radius [pc]")
    axis.set_ylabel(r"X-ray luminosity [erg s$^{-1}$]")
    axis.set_yscale("log")
    axis.grid(True, which="both", alpha=0.25)
    figure.tight_layout()
    save_fig(figure, fig_dir / "xray_luminosity_profile.png")


def plot_cloudy_radial_profiles(c_output, continuum_path=None):
    """Save Cloudy electron-temperature and line-emissivity profiles."""
    radius_pc = c_output.radius / pc.CST.PC

    temperature_figure, temperature_axis = plt.subplots(figsize=(10, 6))
    temperature_axis.plot(radius_pc, c_output.te, color="tab:red")
    temperature_axis.set_title("Cloudy radial electron temperature")
    temperature_axis.set_xlabel("Radius [pc]")
    temperature_axis.set_ylabel("Electron temperature [K]")
    temperature_axis.set_yscale("log")
    temperature_axis.grid(True, which="both", alpha=0.25)
    temperature_figure.tight_layout()
    save_fig(temperature_figure, fig_dir / "temperature_profile.png")

    line_figure, line_axis = plt.subplots(figsize=(10, 6))
    lines = (
        ("H__1_486132A", r"H$\beta$ 4861", "tab:blue"),
        ("N__2_658345A", "[NII] 6584", "tab:orange"),
        ("O__3_500684A", "[OIII] 5007", "tab:green"),
    )
    for ref, label, color in lines:
        emissivity = np.asarray(c_output.get_emis(ref))
        positive = emissivity > 0.0
        line_axis.plot(radius_pc[positive], emissivity[positive], label=label, color=color)
    line_axis.set_title("Cloudy radial line emissivities")
    line_axis.set_xlabel("Radius [pc]")
    line_axis.set_ylabel(r"Emissivity [erg s$^{-1}$ cm$^{-3}$]")
    line_axis.set_yscale("log")
    line_axis.grid(True, which="both", alpha=0.25)
    line_axis.legend()
    line_figure.tight_layout()
    save_fig(line_figure, fig_dir / "emissivity_radial_profiles.png")
    if continuum_path is not None:
        plot_xray_luminosity_profile(c_output, continuum_path)


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


def make_rgb_image(red, green, blue):
    """Normalize three projected emissivity maps into an RGB image."""
    channels = []
    for channel in (red, green, blue):
        maximum = np.nanmax(channel)
        channels.append(
            np.zeros_like(channel, dtype=np.uint8)
            if maximum <= 0.0
            else np.clip(channel / maximum * 255.0, 0.0, 255.0).astype(np.uint8)
        )
    return np.stack(channels, axis=-1)


def project_spherical_emissivity(
    c_output, m3d, ref, proj_axis, pixel_samples=4, ray_samples=2049
):
    """Project a 1D spherical emissivity profile onto the C3D image plane.

    The C3D grid coordinates are pixel locations, while Cloudy reports values
    at many unequal radial zone centers.  Supersampling each image pixel and
    integrating along the ray prevents narrow radial zones from becoming
    conspicuous numerical rings in the projected image.
    """
    coordinates = [m3d.cub_coord.x_vec, m3d.cub_coord.y_vec, m3d.cub_coord.z_vec]
    projected_axes = [axis for axis in range(3) if axis != proj_axis]
    first, second = (coordinates[axis] for axis in projected_axes)
    first_grid, second_grid = np.meshgrid(first, second, indexing="ij")
    first_step = np.median(np.diff(first)) if len(first) > 1 else 0.0
    second_step = np.median(np.diff(second)) if len(second) > 1 else 0.0
    subpixel_offsets = np.linspace(-0.5, 0.5, pixel_samples, endpoint=False)
    subpixel_offsets += 0.5 / pixel_samples

    impact_parameter = np.concatenate(
        [
            np.sqrt(
                (first_grid + first_step * offset_first) ** 2
                + (second_grid + second_step * offset_second) ** 2
            ).ravel()
            for offset_first in subpixel_offsets
            for offset_second in subpixel_offsets
        ]
    )
    radius = np.asarray(c_output.radius, dtype=float)
    emissivity = np.asarray(c_output.get_emis(ref), dtype=float)
    # Cloudy values are zone-center values. Extend them to the origin and
    # outer edge so the interpolation does not create a hollow center.
    radius = np.r_[0.0, radius, radius[-1]]
    emissivity = np.r_[emissivity[0], emissivity, emissivity[-1]]
    ray_coordinate = np.linspace(-radius[-1], radius[-1], ray_samples)
    ray_step = ray_coordinate[1] - ray_coordinate[0]
    projected = np.empty(impact_parameter.size)
    chunk_size = 512
    for start in range(0, impact_parameter.size, chunk_size):
        stop = start + chunk_size
        ray_radius = np.sqrt(
            impact_parameter[start:stop, None] ** 2 + ray_coordinate[None, :] ** 2
        )
        ray_emissivity = np.interp(
            ray_radius, radius, emissivity, left=0.0, right=0.0
        )
        projected[start:stop] = np.trapz(ray_emissivity, dx=ray_step, axis=1)

    n_subpixels = pixel_samples**2
    image = projected.reshape(n_subpixels, *first_grid.shape).mean(axis=0)
    return image


def add_profiles_to_rgb(m3d, image_axis, ref, image_proj_axis=1, nx=20, ny=20):
    """Overlay profile panels using the RGB axes position as the coordinate frame."""
    profile_axis = ("x", "y", "z")[image_proj_axis]
    profiles = m3d.get_profile(ref, axis=profile_axis)
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
            profile_axis.plot(m3d.vel_tab, profile, color="yellow")
            profile_axis.set_ylim(0.0, 1.05)
            profile_axis.set_xlim(*PROFILE_DISPLAY_LIMITS)
            profile_axis.set_xticks([])
            profile_axis.set_yticks([])
            profile_axis.patch.set_alpha(0.0)


def other_plots(m3d, proj_axis, n_cut):
    c_output = m3d.m[0]
    hbmap = project_spherical_emissivity(c_output, m3d, "H__1_486132A", proj_axis)
    niimap = project_spherical_emissivity(c_output, m3d, "N__2_658345A", proj_axis)
    oiiimap = project_spherical_emissivity(c_output, m3d, "O__3_500684A", proj_axis)

    plt.subplot(331)
    show_image(
        hbmap,
        "Hb",
        r"Emissivity [erg s$^{-1}$ cm$^{-2}$]",
    )

    plt.subplot(332)
    show_image(
        niimap,
        "[NII]",
        r"Emissivity [erg s$^{-1}$ cm$^{-2}$]",
    )

    plt.subplot(333)
    show_image(
        oiiimap,
        "[OIII]",
        r"Emissivity [erg s$^{-1}$ cm$^{-2}$]",
    )

    hb = hbmap
    nii = niimap
    oiii = oiiimap
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


def run_profile(profile_file, include_xray=False):
    """Run the complete workflow for one external profile CSV."""
    configure_output_dirs(profile_file)
    profile_stem = profile_file.stem
    model_name = f"{MODEL_NAME_PREFIX}_{profile_stem}"
    radius_pc, velocity_kms, density_cm3, temperature_k = read_profiles(profile_file)
    plot_input_profiles(
        radius_pc,
        velocity_kms,
        density_cm3,
        title=f"Input spherical radial profiles: {profile_file.name}",
    )
    poly_params = [20.0, 60.0]
    poly_velocity_kms = polynomial_velocity_profile(radius_pc, poly_params)
    plot_input_profiles(
        radius_pc,
        poly_velocity_kms,
        density_cm3,
        title=f"Input density with C3D polynomial velocity: {profile_file.name}",
        filename="input_radial_profiles_wpolyv.png",
        velocity_label="C3D polynomial velocity",
    )
    cloudy_exe = find_cloudy_exe(script_dir)
    pc.config.cloudy_exe = str(cloudy_exe)
    model_path = temp_model_dir / model_name
    build_model(model_path, radius_pc, density_cm3, temperature_k, include_xray=include_xray)
    pc.print_make_file(dir_=str(temp_model_dir) + "/")
    pc.run_cloudy(dir_=str(temp_model_dir) + "/", n_proc=6, model_name=model_name, use_make=True)
    c_output = pc.CloudyModel(str(model_path))
    xray_path = model_path.with_name(model_path.name + XRAY_FILE_SUFFIX) if include_xray else None
    plot_cloudy_radial_profiles(c_output, xray_path)
    m3d = pc.C3D(c_output, dims=DIM, center=True, n_dim=1)

    m3d.set_velocity(params=poly_params)
    m3d.config_profile(size_spectrum=51, vel_max=50.0, v_turb=0.01)
    n_cut = (DIM - 1) // 2

    plt.figure(figsize=(10, 10))
    plot_profiles(
        m3d,
        n_cut,
        n_cut,
        f"Line profiles: default polynomial velocity law ({profile_stem})",
    )
    save_fig(plt.gcf(), fig_dir / "profile_default.png")

    set_user_velocity(m3d, radius_pc, velocity_kms)

    plt.figure(figsize=(10, 10))
    plot_profiles(
        m3d,
        n_cut,
        n_cut,
        f"Line profiles: user velocity profile from {profile_file.name}",
        velocity_limits=PROFILE_DISPLAY_LIMITS,
    )
    save_fig(plt.gcf(), fig_dir / "profile_user_velocity.png")

    plt.figure(figsize=(15, 15))
    other_plots(m3d, PROJ_AXIS, n_cut)
    save_fig(plt.gcf(), fig_dir / "derived_maps.png")

    rgb_maps = [
        project_spherical_emissivity(c_output, m3d, ref, proj_axis=1)
        for ref in ("N__2_658345A", "O__3_500684A", "H__1_486132A")
    ]
    image = make_rgb_image(*rgb_maps)
    rgb_extent = [
        m3d.cub_coord.x_vec[0] / pc.CST.PC,
        m3d.cub_coord.x_vec[-1] / pc.CST.PC,
        m3d.cub_coord.z_vec[0] / pc.CST.PC,
        m3d.cub_coord.z_vec[-1] / pc.CST.PC,
    ]
    save_fig(show_rgb_with_colorbars(image, rgb_extent), fig_dir / "rgb_compact.png")

    rgb_with_profiles_figure, rgb_axis = plt.subplots(figsize=(15, 15))
    rgb_axis.imshow(
        image,
        extent=rgb_extent,
        origin="lower",
        interpolation="nearest",
        vmin=0,
        vmax=255,
    )
    rgb_axis.set_title(f"RGB emission image with [NII] 6584 line profiles: {profile_stem}")
    rgb_axis.set_xlabel("Projected x [pc]")
    rgb_axis.set_ylabel("Projected z [pc]")
    add_rgb_colorbars(rgb_with_profiles_figure)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        add_profiles_to_rgb(m3d, rgb_axis, ref=3, image_proj_axis=1, nx=20, ny=20)
    save_fig(rgb_with_profiles_figure, fig_dir / "rgb_with_profiles.png")

    f, ax = plt.subplots()
    n2map = m3d.get_emis("N__2_658345A").sum(axis=PROJ_AXIS)
    hbmap = m3d.get_emis("H__1_486132A").sum(axis=PROJ_AXIS)
    o3map = m3d.get_emis("O__3_500684A").sum(axis=PROJ_AXIS)
    mask = np.logical_and.reduce(
        [line_map > 0.01 * line_map.max() for line_map in (hbmap, o3map, n2map)]
    )
    ax.scatter(
        np.log10(safe_divide(n2map, hbmap)[mask]),
        np.log10(safe_divide(o3map, hbmap)[mask]),
    )
    ax.set_xlabel("log10([NII]/Hb)")
    ax.set_ylabel("log10([OIII]/Hb)")
    save_fig(f, fig_dir / "diagnostic_scatter.png")
    plt.close("all")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--profile",
    type=Path,
    help="run one profile CSV instead of every CSV in radial_profiles/",
)
parser.add_argument(
    "--xray",
    action="store_true",
    help="save and plot the radial 0.1--10 keV X-ray luminosity for each profile",
)
args = parser.parse_args()
profile_files = (
    [args.profile]
    if args.profile is not None
    else sorted(
        profile_file
        for profile_file in PROFILE_DIR.glob("*.csv")
        if profile_file.stem != "radial_profile_0Myr"
    )
)
if not profile_files:
    raise FileNotFoundError(f"No profile CSV files found in {PROFILE_DIR}")
for profile_file in profile_files:
    profile_file = profile_file.expanduser().resolve()
    if not profile_file.is_file():
        raise FileNotFoundError(profile_file)
    run_profile(profile_file, include_xray=args.xray)
