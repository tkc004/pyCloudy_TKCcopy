# C3D Reference

This page documents `pyCloudy.C3D`, the helper that builds a pseudo-3D model
from a list of 1D `CloudyModel` objects.

The implementation lives in `pyCloudy/c3d/model_3d.py`.

## What It Takes As Input

`C3D` is constructed from a list of already-loaded `CloudyModel` objects.
Those models should represent different viewing angles or sectors of the same
physical object.

Typical constructor arguments:

- `list_of_models`: list of `CloudyModel` objects
- `dims`: 1D, 2D, or 3D grid size; it may be a single integer or an array-like shape
- `center`: whether the coordinate origin is centered in the cube
- `angles`: rotation angles in degrees, used to orient the cube
- `n_dim`: interpolation mode, usually `2` for angular interpolation or `3` for angular plus azimuthal interpolation
- `file_coeffs`: optional file containing previously saved interpolation coefficients
- `interp_method`: interpolation method for angular matching
- `plan_sym`: use mirror symmetry across the equatorial plane when interpolating
- `r_max`: override the 3D cube size inferred from the models
- `r_interp_method`: radial interpolation method, such as `numpy` or `scipy`

In practice, `C3D` is built from models that already have their viewing angle
set, then it interpolates between those models onto a 3D Cartesian grid.

## What It Produces

Once created, `C3D` exposes 3D arrays and integrated diagnostics derived from
the input 1D models.

Common outputs include:

- `nH`
- `ne`
- `te`
- `ff`
- `log_U`
- `log_U_mean`
- `emis_labels`
- emissivity cubes returned by `get_emis(ref)`
- ionic fraction cubes returned by `get_ionic(elem, ion)`
- integrated masses and line luminosities
- cached line profiles, ionic cubes, and emissivity cubes

## Lifecycle

1. Read or build a list of `CloudyModel` objects.
2. Create `C3D(models, ...)`.
3. Inspect the 3D fields or generate maps and profiles.
4. Plot or export the derived 3D quantities.

## Constructor

### `__init__(list_of_models, dims=51, center=True, angles=None, n_dim=2, file_coeffs=None, interp_method=None, plan_sym=False, r_max=None, r_interp_method='numpy')`

Create a 3D object from the input models.

Notes:

- `n_dim=2` is the common case for interpolating between models in angle.
- `n_dim=3` enables 3D angular interpolation.
- `r_interp_method` controls how values are interpolated radially inside each model.

## Main Properties

When these methods say `ref`, they mean a line reference: either a line label
such as `H__1_486133A` or a zero-based line index such as `0`.

Examples:

- `c3d.get_emis("H__1_486133A")`
- `c3d.get_emis_vol("O__3_500684A")`
- `c3d.get_emis(0)`

When these methods say `elem`, they mean an element symbol recognized by the
input models, such as `H`, `He`, `C`, `N`, `O`, `Ne`, `S`, `Ar`, `Cl`, `Fe`,
or `Si`.

When they say `ion`, they mean the zero-based ion stage for that element:

- `0` means neutral, for example `H0` or `O0`
- `1` means singly ionized, for example `H+` or `O+`
- `2` means doubly ionized, for example `He++` or `O++`

### Physical fields

- `nH` - hydrogen density cube
- `ne` - electron density cube
- `te` - electron temperature cube
- `ff` - filling-factor cube
- `log_U` - ionization parameter cube
- `log_U_mean` - volume-weighted mean ionization parameter

### Geometry and coordinates

- `cub_coord` - internal `CubCoord` geometry object
- `angles` - current rotation angles
- `r_max` - cube size used for the coordinate grid
- `x_unit` - output unit for map axes

### Model content

- `emis_labels` - available emissivity labels
- `get_emis_list(available=False)` - emissivity cubes that have been computed, or all available labels if `available=True`
- `get_ionic_list()` - ionic cubes that have been computed
- `del_emis(ref)` - remove one cached emissivity cube
- `del_ionic(elem, ion)` - remove one cached ionic cube
- `print_all_emis_vol(norm=None)` - print all line luminosities, optionally normalized to a reference line

## Main Methods

### Emissivities

- `get_emis(ref)` - interpolate a line emissivity cube for a line label or index
- `get_emis_vol(ref, at_earth=False)` - volume-integrated line luminosity, optionally normalized to Earth distance
- `get_T0_emis(ref)` - emissivity-weighted temperature for the selected line
- `get_t2_emis(ref)` - emissivity-weighted temperature fluctuation for the selected line
- `get_vel_emis(ref)` - emissivity-weighted velocity for a line, when a velocity field is configured

### Ions and abundances

- `get_ionic(elem, ion)` - interpolate an ionic fraction cube for one element and ion stage
- `get_ab_ion_vol(elem=None, ion=None)` - volume-weighted ionic abundance
- `get_ab_ion_vol_ne(elem, ion)` - ionic abundance weighted by hydrogen density
- `get_T0_ion_vol(elem, ion)` - temperature weighted by ionic abundance
- `get_T0_ion_vol_ne(elem, ion)` - density-weighted temperature diagnostic for one ion
- `get_t2_ion_vol(elem, ion)` - temperature fluctuation for one ion
- `get_t2_ion_vol_ne(elem, ion)` - density-weighted temperature fluctuation for one ion
- `get_vel_ionic(elem, ion)` - velocity cube associated with one ionic fraction
- `get_vel_emis(ref)` - emissivity-weighted velocity cube for one line

### Profiles and maps

- `set_velocity(*args, **kwargs)` - pass velocity-law settings through to `CubCoord`
- `config_profile(size_spectrum=21, vel_max=20., v_turb=5., profile_function='gaussian')` - configure line-profile generation
- `get_profile(ref, axis='x')` - extract a line profile along the chosen axis
- `get_profile_list()` - list the configured profiles
- `del_profile(ref=None, axis='x')` - clear one cached profile or reset all of them
- `plot_profiles(Nx=10, Ny=10, ref=None, axis='x', normalized=True, i_fig=None, ...)` - plot profiles in a grid
- `get_RGB(list_emis=[0, 1, 2], axes=1)` - build RGB composite cubes or maps from three emissivity channels

### Saving and setup

- `save_coeffs(file_coeffs)` - save interpolation coefficients for reuse
- changing `angles` reconfigures the interpolation mesh and recomputes the internal interpolation coefficients

### Notes on cached results

`C3D` caches expensive results such as emissivity cubes, ionic cubes, and line
profiles. If you change the geometry, velocity field, or profile settings, the
relevant caches are cleared automatically. If you want to free a specific
cached result, use `del_emis(...)`, `del_ionic(...)`, or `del_profile(...)`.

## Practical Example

```python
import pyCloudy as pc

models = [
    pc.CloudyModel("models/angle_1"),
    pc.CloudyModel("models/angle_2"),
]
c3d = pc.C3D(models, dims=51, n_dim=2)
print(c3d.log_U_mean)
print(c3d.get_emis_vol("H__1_486133A"))
```

```python
fig = c3d.plot_profiles(ref="H__1_486133A")
```

## Notes

- `C3D` is an interpolation and analysis layer, not a Cloudy runner.
- It depends on `CloudyModel` objects as inputs.
- The exact results depend on how similar the input models are and on the interpolation method chosen.
- If you want a quick visual check of a 3D result, `get_RGB(...)` can combine
  three emissivity channels into a false-color composite.
