# pyCloudy_TKCcopy API Reference

This page summarizes the main public objects exposed by `pyCloudy_TKCcopy`.
It is intentionally compact and points to the most commonly used classes and
helpers. For implementation details, see the source files in `pyCloudy/`.

For detailed model-building and model-reading references, see:

- [CloudyInput.md](CloudyInput.md)
- [CloudyModel.md](CloudyModel.md)
- [C3D.md](C3D.md)

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

Common entry points include `set_BB(...)`, `set_cste_density(...)`,
`set_abund(...)`, `print_input(...)`, and `run_cloudy(...)`. For the full
method-by-method reference and parameter notes, see the dedicated
[CloudyInput.md](CloudyInput.md) page.

### `CloudyModel`

Load a completed Cloudy model and access the output files and diagnostics.

If you need multiple models, use `load_models(model_name=None, mod_list=None,
n_sample=None, verbose=False, **kwargs)`. It searches for matching `.out`
files when `model_name` is given, or reads the paths in `mod_list` when you
already have a list of models. The optional `n_sample` argument randomly
selects a subset, `verbose=True` prints progress, and extra keyword arguments
are forwarded to `CloudyModel`.

See the dedicated [CloudyModel.md](CloudyModel.md) page for the full
output-reading and analysis reference.

### Module Helpers

- `load_models(model_name=None, mod_list=None, n_sample=None, verbose=False, **kwargs)` - load one or more models by base name or from an explicit file list
- `print_make_file(dir_=None)` - generate a Makefile that points to the configured Cloudy executable
- `run_cloudy(dir_=None, n_proc=1, use_make=True, model_name=None, precom="", cloudy_version=None)` - run Cloudy from Python, either through `make` or by calling the executable directly

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

See the dedicated [C3D.md](C3D.md) page for the full input/output and method
reference.

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
