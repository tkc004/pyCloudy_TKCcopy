# CloudyModel Reference

This page documents `pyCloudy.CloudyModel`, the class used to read Cloudy
output files and analyze a finished model.

The implementation lives in `pyCloudy/c1d/cloudy_model.py`.

## What It Reads

`CloudyModel` reads the output files produced by a Cloudy run and exposes them
through Python properties and helper methods.

If you want to load many models at once, use `pyCloudy.load_models(...)`.
That helper searches for matching output files and returns a list of
`CloudyModel` objects.

Common files include:

- `.rad` for the radius grid
- `.phy` for physical quantities
- `.ovr` for overview data
- `.emis` for emissivities
- `.cont` for continuum data
- `.ele_*` for ionic abundances by element
- `.heat`, `.cool`, `.opd`, and grain-related files when present

## Lifecycle

1. Run Cloudy for a model.
2. Create `CloudyModel("models/example")`.
3. Access the physical arrays, integrated values, and diagnostics.
4. Print summaries or generate plots from the loaded model.

## Constructor

### `__init__(model_name, ...)`

Load a Cloudy model from the files written by the run.

- `model_name`: base path of the model, without the output extension
- Optional flags control which extensions are read

Typical options include:

- `read_all_ext`
- `read_emis`
- `read_grains`
- `read_cont`
- `list_elem`
- `distance`
- `line_is_log`
- `emis_is_log`

## Core Data

Frequently used properties include:

- `radius`
- `depth`
- `thickness`
- `ne`
- `nH`
- `te`
- `ff`
- `abunds`
- `cool`
- `heat`
- `T0`
- `t2`
- `H_mass`
- `Hbeta`

## Useful Methods

When these methods say `ref`, they mean a line reference: either a line label
such as `H__1_486133A` or a zero-based line index such as `0`.

Examples:

- `model.get_emis("H__1_486133A")`
- `model.get_line("O__3_500684A")`
- `model.get_emis(0)`

When these methods say `elem`, they mean an element symbol recognized by the
model, such as `H`, `He`, `C`, `N`, `O`, `Ne`, `S`, `Ar`, `Cl`, `Fe`, or `Si`.

When they say `ion`, they mean the zero-based ion stage for that element:

- `0` means neutral, for example `H0` or `O0`
- `1` means singly ionized, for example `H+` or `O+`
- `2` means doubly ionized, for example `He++` or `O++`
- In practice, `elem` and `ion` are the pair you pass to ion-sensitive methods
  such as `get_ionic(...)`, `get_ab_ion_vol(...)`, and `get_T0_ion_vol(...)`.

The continuum helpers use these accepted values. The `unit` choice controls
what kind of x-axis or y-axis quantity you get back, while `cont` selects
which continuum component is being sampled.

- `unit` for `get_cont_x(...)`: `Ryd`, `eV`, `Ang`, `mu`, `cm-1`, `Hz`, `kHz`, `MHz`, `GHz`
- `cont` for `get_cont_y(...)`: `incid`, `trans`, `diffout`, `ntrans`, `reflec`, `total`
- `unit` for `get_cont_y(...)`: `esc`, `ec3`, `es`, `esA`, `esHz`, `esAc`, `ec3A`, `esHzc`, `Jy`, `Q`, `Wcmu`, `WmHz`, `WmA`, `phs`, `phsmu`, `phsc`, `phsmuc`
- `dist_norm` for `get_cont_y(...)` and `get_interp_cont(...)`: `at_earth`, `r_out`, or a numeric distance in cm

In short:

- use `cont='incid'` for the incident source spectrum
- use `cont='trans'` for transmitted radiation
- use `cont='diffout'` for diffuse outward radiation
- use `cont='ntrans'` for net transmitted radiation
- use `cont='reflec'` for reflected radiation
- use `cont='total'` for the total continuum

And for `dist_norm`:

- `at_earth` applies the model distance, when one was set with `set_distance(...)`
- `r_out` normalizes at the model outer radius
- a numeric value lets you supply a custom distance in centimeters

### Ionization and abundances

- `get_ionic(elem, ion)` - ionic fraction for an element and ionization stage
- `get_ab_ion_vol(elem, ion)` - volume-weighted ionic fraction
- `get_ab_ion_vol_ne(elem, ion)` - ionic fraction weighted by hydrogen density
- `get_T0_ion_vol(elem, ion)` - volume-weighted electron temperature for an ion
- `get_T0_ion_vol_ne(elem, ion)` - density-weighted version of `get_T0_ion_vol(...)`
- `get_t2_ion_vol_ne(elem, ion)` - density-weighted temperature fluctuation

### Lines and emissivities

- `get_line(ref)` - line intensity
- `get_emis(ref)` - emissivity for a line
- `get_emis_vol(ref, at_earth=False)` - volume-integrated emissivity
- `get_emis_rad(ref)` - radial integral of the emissivity
- `get_T0_emis(ref)` - emissivity-weighted electron temperature
- `get_T0_emis_rad(ref)` - radial emissivity-weighted electron temperature
- `get_ne_emis(ref)` - emissivity-weighted electron density
- `get_t2_emis(ref)` - emissivity-weighted temperature fluctuation
- `print_lines(...)` - print line intensities
- `print_stats()` - print a compact summary of the model

### PyNeb and custom emissivities

These helpers let you update or extend the emissivity table after the model
has already been read.

- `emis_from_pyneb(emis_labels=None, atoms=None)` - replace selected
  emissivities with values computed from PyNeb atoms
- `add_emis_from_pyneb(new_label, pyneb_atom, label=None, wave=None)` - add a
  new emissivity line computed from a PyNeb `Atom` or `RecAtom`
- `copy_line(new_label, old_label)` - duplicate an existing emissivity under a
  new label

`emis_from_pyneb(...)` is useful when you want to recompute line emissivities
from a PyNeb atomic dataset but keep the Cloudy ionization structure and
physical conditions. `add_emis_from_pyneb(...)` is useful when Cloudy did not
write a line you want to track, but PyNeb can compute it from an atom object
you already have.

Example:

```python
import pyneb as pn

o2 = pn.Atom("O", 2)
model.add_emis_from_pyneb("O__2R_4639", o2, label="4638.86")
```

### Continuum

- `get_cont_x(unit=...)` - continuum x-axis in wavelength, frequency, energy, or wavenumber
- `get_cont_y(cont=..., unit=..., dist_norm=...)` - continuum flux or intensity in the requested units
- `get_integ_spec(...)` - integrated spectrum over a wavelength range
- `get_interp_cont(...)` - interpolated continuum value
- `plot_spectrum(...)` - plot the spectrum

### Thermal diagnostics

- `get_T0_emis(ref)` - emissivity-weighted temperature
- `get_t2_emis(ref)` - emissivity-weighted temperature fluctuation
- `get_Hb_SB()` - Hbeta surface brightness
- `get_Hb_EW()` - Hbeta equivalent width

## Practical Example

```python
import pyCloudy as pc

model = pc.CloudyModel("models/example")
print(model.print_stats())
print(model.get_Hb_EW())
print(model.get_emis_vol("H__1_486133A"))
```

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
model.plot_spectrum(ax=ax, cont="ntrans", xunit="Ang", yunit="Jy")
fig.savefig("spectrum.png")
```

## Notes

- The class is a reader and analysis layer, not a runner.
- The exact output files available depend on the Cloudy run options.
- If you need a specific file, first make sure your input requested it.
