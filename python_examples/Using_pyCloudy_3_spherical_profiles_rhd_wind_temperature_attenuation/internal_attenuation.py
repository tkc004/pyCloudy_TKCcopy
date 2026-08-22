"""Line-of-sight dust attenuation for C3D emission diagnostics.

Definitions used here:

* ``m3d.nH`` is the dust-bearing hydrogen-equivalent density in cm^-3; the
  dust-to-gas ratio is included in ``A_V/N_H``.
* ``CCM 89`` supplies the wavelength dependence, with ``R_V=3.1``.
* The observer is on the positive side of the selected Cartesian axis.
* Cell path lengths are the coordinate spacing along that axis, in cm.
* Empty cells (``nH <= 0``) and regions beyond the cube have zero opacity.
  Emission-cell optical depth uses half of the local cell column, while the
  external boundary has no additional column.
"""

import numpy as np

import pyCloudy as pc
from pyCloudy.utils.physics import atomic_mass
from pyCloudy.utils import misc


DEFAULT_AV_PER_NH = 5.3e-22  # mag cm^2 H atom^-1


def sigma_lambda_per_h(wavelength_angstrom, law="CCM 89", r_v=3.1,
                       av_per_nh=DEFAULT_AV_PER_NH):
    """Return dust extinction cross-section per H atom in cm^2."""
    red_corr = pc.RedCorr(E_BV=1.0, R_V=r_v, law=law)
    a_lambda_over_av = np.asarray(red_corr.X(wavelength_angstrom), dtype=float) / r_v
    return av_per_nh / 1.086 * a_lambda_over_av


def _axis_number(axis):
    if axis in ("x", 0):
        return 0
    if axis in ("y", 1):
        return 1
    if axis in ("z", 2):
        return 2
    raise ValueError("axis must be one of 'x', 'y', 'z', 0, 1, or 2")


def transmission_cube(m3d, wavelength_angstrom, axis="x", law="CCM 89",
                      r_v=3.1, av_per_nh=DEFAULT_AV_PER_NH,
                      observer_side="positive"):
    """Compute cell transmission toward the selected observer boundary.

    ``observer_side`` is currently ``"positive"`` (the high-coordinate end
    of the selected axis) or ``"negative"``. This explicit choice prevents
    the attenuation direction from being hidden in an array operation.
    """
    axis = _axis_number(axis)
    if observer_side not in ("positive", "negative"):
        raise ValueError("observer_side must be 'positive' or 'negative'")
    coordinates = (m3d.cub_coord.x_vec, m3d.cub_coord.y_vec, m3d.cub_coord.z_vec)
    coordinate = np.asarray(coordinates[axis], dtype=float)
    ds = float(np.median(np.diff(coordinate))) if len(coordinate) > 1 else 0.0
    hydrogen_density = np.nan_to_num(np.asarray(m3d.nH, dtype=float), nan=0.0)
    hydrogen_density = np.maximum(hydrogen_density, 0.0)
    opacity = hydrogen_density * sigma_lambda_per_h(
        wavelength_angstrom, law, r_v, av_per_nh
    ) * ds
    if observer_side == "positive":
        tau_from_observer = np.flip(
            np.cumsum(np.flip(opacity, axis=axis), axis=axis), axis=axis
        )
    else:
        tau_from_observer = np.cumsum(opacity, axis=axis)
    # Emission is represented at the cell centre, so exclude half of its own
    # optical depth from the foreground column.
    tau_after_cell = np.maximum(tau_from_observer - 0.5 * opacity, 0.0)
    return np.exp(-tau_after_cell)


def attenuated_projection(m3d, ref, wavelength_angstrom, axis="x", **kwargs):
    """Integrate a line emissivity cube with internal dust transmission."""
    axis = _axis_number(axis)
    emissivity = np.asarray(m3d.get_emis(ref), dtype=float)
    return (emissivity * transmission_cube(
        m3d, wavelength_angstrom, axis=axis, **kwargs
    )).sum(axis=axis)


def attenuated_profile(m3d, ref, wavelength_angstrom, axis="x", **kwargs):
    """Calculate a velocity profile with cell-by-cell internal attenuation."""
    axis = _axis_number(axis)
    line_ref = m3d.m[0]._l_emis(ref)
    elem = misc.get_elem_ion(line_ref)[0].capitalize()
    emissivity = np.asarray(m3d.get_emis(line_ref), dtype=float)
    transmission = transmission_cube(m3d, wavelength_angstrom, axis=axis, **kwargs)
    emissivity *= transmission
    coeff1 = np.sqrt(2 * pc.CST.BOLTZMANN * 1e4 / pc.CST.HMASS) / 1e5
    zeta_0 = np.sqrt(
        m3d.v_turb**2 + coeff1**2 * m3d.te / 1e4 / atomic_mass(elem)
    )
    velocity = (m3d.cub_coord.vel_x, m3d.cub_coord.vel_y, m3d.cub_coord.vel_z)[axis]
    result_shape = list(emissivity.shape)
    del result_shape[axis]
    result = np.zeros((m3d.size_spectrum, *result_shape)).squeeze()
    for index, velocity_bin in enumerate(m3d.vel_tab):
        delta_v = velocity + velocity_bin
        result[index] = (
            emissivity * m3d.profile_function(x=delta_v, zeta_0=zeta_0)
        ).sum(axis=axis)
    return result
