# pyCloudy_TKCcopy API Reference

This page summarizes the main public objects exposed by `pyCloudy_TKCcopy`.
It is intentionally compact and points to the most commonly used classes and
helpers. For implementation details, see the source files in `pyCloudy/`.

## Top-Level Imports

The package root exports the most common entry points:

- `pyCloudy.config` - runtime configuration object
- `pyCloudy.log_` - shared logger instance
- `pyCloudy.CloudyInput` - create Cloudy input files
- `pyCloudy.CloudyModel` - read and analyze Cloudy outputs
- `pyCloudy.load_models` - load several Cloudy models
- `pyCloudy.print_make_file` - write a `make` helper file for Cloudy runs
- `pyCloudy.run_cloudy` - launch Cloudy from Python
- `pyCloudy.C3D` - build 3D models from 1D runs
- `pyCloudy.CubCoord` - geometry helper for 3D models
- `pyCloudy.MdB` - 3MdB database interface
- `pyCloudy.CST` - physical constants and unit helpers
- `pyCloudy.RedCorr` - reddening correction helper
- `pyCloudy.sextract`, `pyCloudy.save`, `pyCloudy.restore` - common utility functions

## `pyCloudy.c1d`

### `CloudyInput`

Create Cloudy input files and optionally run Cloudy.

Common methods:

- `set_BB(...)` - define a blackbody source
- `set_star(...)` - define a custom source
- `set_cste_density(...)` - define a constant density model
- `set_dlaw(...)` - define a density law
- `set_radius(...)` - set inner and outer radii
- `set_abund(...)` - set abundances
- `set_other(...)` - add raw Cloudy commands
- `set_iterate(...)` - configure iterations
- `set_sphere(...)` - select spherical or open geometry
- `set_grains(...)` - configure dust grains
- `set_stop(...)` - define stopping criteria
- `set_distance(...)` - set distance and normalization
- `set_emis_tab(...)` - define an emission line table
- `print_input(...)` - write the input file
- `run_cloudy(...)` - execute Cloudy

### `CloudyModel`

Load a completed Cloudy model and access derived data.

Frequently used properties and methods:

- `zones`, `n_zones`
- `depth`, `radius`, `thickness`
- `ne`, `nH`, `te`, `ff`
- `cool`, `heat`
- `abunds`
- `get_ionic(elem, ion)`
- `get_line(ref)`
- `get_emis(ref)`
- `get_emis_vol(ref, at_earth=False)`
- `get_T0_ion_vol(elem=None, ion=None)`
- `get_ab_ion_vol(elem=None, ion=None)`
- `get_t2_ion_vol(elem=None, ion=None)`
- `get_T0_emis(ref)`
- `get_t2_emis(ref)`
- `get_cont_x(...)`
- `get_cont_y(...)`
- `get_integ_spec(...)`
- `get_interp_cont(...)`
- `plot_spectrum(...)`
- `print_lines(...)`
- `print_stats()`

### Module Helpers

- `load_models(...)` - load one or more models by name or from a model list
- `print_make_file(...)` - generate a makefile for running Cloudy
- `run_cloudy(...)` - run a set of models from a directory

## `pyCloudy.c3d`

### `CubCoord`

Coordinate and geometry helper for 3D constructions.

Useful capabilities:

- build a 3D coordinate grid
- apply rotation angles
- define a velocity law
- expose `x`, `y`, `z`, `r`, `theta`, `phi`

### `C3D`

Construct a 3D object from a list of 1D Cloudy models.

Common methods and properties:

- `nH`, `ne`, `te`, `ff`
- `log_U`, `log_U_mean`
- `get_emis(ref)`
- `get_emis_list(...)`
- `get_emis_vol(ref, at_earth=False)`
- `get_ionic(elem, ion)`
- `get_ionic_list()`
- `get_T0_emis(ref)`
- `get_t2_emis(ref)`
- `get_T0_ion_vol(elem, ion)`
- `get_t2_ion_vol(elem, ion)`
- `get_ab_ion_vol(elem=None, ion=None)`
- `get_vel_ionic(elem, ion)`
- `get_vel_emis(ref)`
- `config_profile(...)`
- `get_profile(ref, axis='x')`
- `plot_profiles(...)`
- `get_RGB(...)`

## `pyCloudy.db`

### `MdB`

Database access wrapper used for 3MdB-style workflows.

Typical operations:

- connect to the database
- select a database or temporary database
- inspect tables
- execute SQL statements
- read query results into Python structures

### `MdB_subproc`

Fallback database helper that uses a subprocess-based access path.

## `pyCloudy.utils`

### `physics`

- `CST` - physical constants
- `flux_convert(...)`
- `planck(...)`
- `atomic_mass(...)`
- `get_metallicity(...)`
- `get_abund_nicholls(...)`
- `get_abunds_Ni17_G24(...)`

### `misc`

General-purpose utilities used throughout the project.

Notable helpers include:

- `sextract(...)`
- `save(...)`
- `restore(...)`
- `convert_label(...)`
- `int_to_roman(...)`
- `roman_to_int(...)`
- `get_elem_ion(...)`
- `pyneb2cloudy(...)`
- `cloudy2pyneb(...)`
- `make_mask(...)`
- `quiet_divide(...)`
- `quiet_log10(...)`
- `fill_from_file(...)`
- `write_cols(...)`
- `read_atm_ascii(...)`
- `convert_c13_c17(...)`
- `convert_c17_c13(...)`
- `correc_He1(...)`
- `mytrapz(...)`

### `red_corr`

- `RedCorr` - reddening correction and law management

### `astro`

- `conv_arc(...)` - convert projected distance, distance, and angular size

## Typical Usage

```python
import pyCloudy as pc

pc.config.cloudy_exe = "/usr/local/Cloudy/c23.01/source/cloudy.exe"

c_input = pc.CloudyInput("models/example")
c_input.set_BB(Teff=40000, lumi_unit="q(H)", lumi_value=47)
c_input.set_cste_density(2.0)
c_input.set_radius(r_in=17.3)
c_input.set_abund(predef="ism")
c_input.print_input(to_file=True)

model = pc.CloudyModel("models/example")
print(model.get_Hb_SB())
```

## Source Pointers

If you want the full implementation details, start here:

- `pyCloudy/c1d/cloudy_model.py`
- `pyCloudy/c3d/model_3d.py`
- `pyCloudy/db/MdB.py`
- `pyCloudy/utils/misc.py`
- `pyCloudy/utils/physics.py`
- `pyCloudy/utils/red_corr.py`
- `pyCloudy/utils/astro.py`

