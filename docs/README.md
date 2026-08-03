# pyCloudy_TKCcopy Documentation

This directory is the documentation hub for `pyCloudy_TKCcopy`, the fork of
pyCloudy maintained in this repository. The project already ships with a number
of notebooks, PDFs, and example scripts under `pyCloudy/docs/`, and this file
gives the top-level structure so new users can find the right starting point
quickly.

`pyCloudy_TKCcopy` is a Python library for working with the input and output
files of Cloudy, the photoionization and plasma simulation code by Gary Ferland
and collaborators.

Originally version (pyCloudy) is developed by Christophe Morisset
(* Ch. Morisset (Instituto de Astronomia, Universidad Nacional Autonoma de Mexico))

This version (pyCloudy_TKCcopy) is a fork of the original pyCloudy, maintained
by Tsang Keung Chan, with some modifications.

## What pyCloudy Does

pyCloudy_TKCcopy is a companion library for the Cloudy photoionization code. It
helps you:

- build Cloudy input files from Python
- run Cloudy from scripts or notebooks
- read and analyze Cloudy output
- combine multiple 1D models into 3D structures
- query 3MdB-style databases
- apply astronomy and reddening utilities used in the workflow

## Package Map

### 1D models

Use `pyCloudy.c1d` when you are preparing or reading a single Cloudy run.

- `CloudyInput` writes Cloudy input files and can launch Cloudy
- `CloudyModel` loads output files and exposes derived quantities
- `load_models` helps load several models at once

Useful file types handled by the 1D layer include the standard Cloudy outputs
such as `.rad`, `.cont`, `.phy`, `.ovr`, `.heat`, `.cool`, `.opd`, `.emis`,
and element/ionic abundance files.

### 3D models

Use `pyCloudy.c3d` when you want to build a 3D nebula from a set of 1D runs.

- `C3D` interpolates physical quantities onto a 3D grid
- `CubCoord` manages the coordinate system, rotation, and geometry

The 3D layer can return integrated quantities, emission maps, profiles, and
volume-weighted diagnostics.

### Database support

Use `pyCloudy.db` if you work with the 3MdB database or similar SQL-backed
archives.

- `MdB` provides direct database access
- `MdB_subproc` provides a subprocess-based fallback

### Utilities

Use `pyCloudy.utils` for common helpers.

- `physics` for constants, conversions, and abundance helpers
- `misc` for file handling, geometry, labels, and convenience functions
- `red_corr` for reddening laws and correction factors
- `astro` for astronomy-specific conversions

## Suggested Reading Order

If you are new to this fork, this is the best path:

1. Read the top-level `README.rst`
2. Open `pyCloudy/docs/Using_pyCloudy_1.ipynb`
3. Continue through `Using_pyCloudy_2.ipynb`, `Using_pyCloudy_3.ipynb`, and
   `Using_pyCloudy_4.ipynb`
4. Review `Using_pyCloudy_with_PyNeb.ipynb` if you use PyNeb
5. Review `Using_pyCloudy_MdB.ipynb` if you use the database tools

For a spherical 3D workflow with user-supplied radius, velocity, and density
profiles, see [`Using_pyCloudy_3_spherical_profiles.md`](Using_pyCloudy_3_spherical_profiles.md).

## Examples

The original example notebooks are stored in [`pyCloudy/docs/`](../pyCloudy/docs/).
Where available, each example below links to its Python script, notebook, and
exported PDF.

- [`quickstart.py`](../python_examples/quickstart/quickstart.py): Minimal introduction to building, running, and reading a 1D Cloudy model.
- [`Using_pyCloudy_1.py`](../python_examples/Using_pyCloudy_1/Using_pyCloudy_1.py) | [`notebook`](../pyCloudy/docs/Using_pyCloudy_1.ipynb) | [`PDF`](../pyCloudy/docs/Using_pyCloudy_1.pdf): Basic 1D workflow for creating, running, and analyzing a spherical model.
- [`Using_pyCloudy_2.py`](../python_examples/Using_pyCloudy_2/Using_pyCloudy_2.py) | [`notebook`](../pyCloudy/docs/Using_pyCloudy_2.ipynb) | [`PDF`](../pyCloudy/docs/Using_pyCloudy_2.pdf): Grid workflow for running multiple 1D models with varied physical parameters.
- [`Using_pyCloudy_3.py`](../python_examples/Using_pyCloudy_3/Using_pyCloudy_3.py) | [`notebook`](../pyCloudy/docs/Using_pyCloudy_3.ipynb) | [`PDF`](../pyCloudy/docs/Using_pyCloudy_3.pdf): Pseudo-3D workflow built from angular 1D models, including line profiles, maps, and RGB images.
- [`Using_pyCloudy_4.py`](../python_examples/Using_pyCloudy_4/Using_pyCloudy_4.py) | [`notebook`](../pyCloudy/docs/Using_pyCloudy_4.ipynb) | [`PDF`](../pyCloudy/docs/Using_pyCloudy_4.pdf): C3D workflow for extracting slit or aperture observations from a spherical model.
- [`Using_pyCloudy_3_spherical_profiles.py`](../python_examples/Using_pyCloudy_3_spherical_profiles/Using_pyCloudy_3_spherical_profiles.py): Spherical pseudo-3D workflow driven by external radial velocity and density profiles.
- [`Using_pyCloudy_3_spherical_profiles_rhd.py`](../python_examples/Using_pyCloudy_3_spherical_profiles_rhd/Using_pyCloudy_3_spherical_profiles_rhd.py): RHD-profile variant using `radial_profile_rhd.csv` and an ionizing photon rate of `10^49 s^-1`; see [`Using_pyCloudy_3_spherical_profiles_rhd.md`](Using_pyCloudy_3_spherical_profiles_rhd.md).
- [`Using_pyCloudy_with_PyNeb.py`](../python_examples/Using_pyCloudy_with_PyNeb/Using_pyCloudy_with_PyNeb.py) | [`notebook`](../pyCloudy/docs/Using_pyCloudy_with_PyNeb.ipynb) | [`PDF`](../pyCloudy/docs/Using_pyCloudy_with_PyNeb.pdf): PyNeb workflow for comparing atomic data sets and deriving temperature and density diagnostics.
- [`Using_pyCloudy_MdB.py`](../python_examples/Using_pyCloudy_MdB/Using_pyCloudy_MdB.py) | [`notebook`](../pyCloudy/docs/Using_pyCloudy_MdB.ipynb) | [`PDF`](../pyCloudy/docs/Using_pyCloudy_MdB.pdf): Database workflow for querying 3MdB models and plotting diagnostic diagrams.
- [`shocks.py`](../python_examples/shocks/shocks.py) | [`notebook`](../pyCloudy/docs/shocks.ipynb): Shock and elliptical-nebula example producing RGB and density diagnostic maps.

## Quick Examples

### Build a Cloudy input file

```python
import pyCloudy as pc

pc.config.cloudy_exe = "/usr/local/Cloudy/c23.01/source/cloudy.exe"

c_input = pc.CloudyInput("models/example")
c_input.set_BB(Teff=40000, lumi_unit="q(H)", lumi_value=47)
c_input.set_cste_density(2.0)
c_input.set_radius(r_in=17.3)
c_input.set_abund(predef="ism")
c_input.set_iterate()
c_input.set_sphere(True)
c_input.print_input(to_file=True)
```

### Read an existing model

```python
import pyCloudy as pc

model = pc.CloudyModel("models/example")
print(model.get_Hb_SB())
print(model.get_emis("O  3  5006.84A"))
```

### Build a 3D object

```python
import pyCloudy as pc

models = [pc.CloudyModel("models/example")]
c3d = pc.C3D(models, dims=51)
print(c3d.get_emis_list())
```

## API Reference

See [`API.md`](API.md) for a compact reference to the main public classes and
helpers exported by the package.

See [`CloudyInput.md`](CloudyInput.md) for the full `CloudyInput` method reference.

See [`CloudyModel.md`](CloudyModel.md) for the main model-reading and output-analysis reference.

See [`C3D.md`](C3D.md) for the 3D interpolation and map-building reference.

## Working With Cloudy

pyCloudy_TKCcopy expects a local Cloudy executable. You can set it in either of these
ways:

- assign `pc.config.cloudy_exe` in Python
- export the `CLOUDY_EXE` environment variable before importing pyCloudy

If a model run fails, check that the executable path is correct and that the
Cloudy build you are using matches the model files you are generating.

## Data and Examples

The repository includes sample model outputs in `tests/models/` and example
scripts in `pyCloudy/3MdB_17/`.

The documentation notebooks in `pyCloudy/docs/` are the best place to see the
package in action. They show:

- basic model creation
- line and continuum analysis
- 3D visualization workflows
- PyNeb integration
- 3MdB workflows

## Testing

Run the tests with:

```bash
pytest
```

Some tests may require a local Cloudy installation. If you are working on a
machine with a different Cloudy path, update the executable path before running
the tests.

## Building The Documentation

The local documentation entry point lives in `docs/`.

## Notes

- pyCloudy_TKCcopy is a companion library, not a replacement for Cloudy.
- Many features are useful only when optional dependencies are installed.
- The repository includes older examples and notebooks for reference; they are
  still helpful for understanding the package structure even if you write new
  code in a different style.
