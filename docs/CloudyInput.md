# CloudyInput Reference

This page documents `pyCloudy.CloudyInput`, the helper used to assemble Cloudy input files.
The implementation lives in `pyCloudy/c1d/cloudy_model.py`.

## What It Does

`CloudyInput` collects the model ingredients for a single Cloudy run:

- source spectrum
- density and geometry
- abundances and grains
- stopping criteria and emissivity tables
- extra Cloudy commands

Use it to build an `.in` file, then optionally launch Cloudy from Python.

## Lifecycle

1. Create the object with a model name.
2. Configure the source, density, radius, abundances, and other commands.
3. Write the input with `print_input(...)`.
4. Run Cloudy with `run_cloudy(...)` if desired.

## Methods

### `__init__(model_name=None)`

Create a new input builder.

- `model_name`: base name used for the generated input and output files.

### `set_save_str(save='save')`

Choose whether generated save commands use `save` or `punch`.

- Invalid values are replaced with `save`.

### `set_radius(r_in=None, r_out=None)`

Set the inner and optional outer radius in log(cm).

- `r_in`: inner radius
- `r_out`: optional outer radius

### `set_BB(Teff=None, lumi_unit=None, lumi_value=None)`

Define a blackbody source.

- `Teff`: effective temperature in K
- `lumi_unit`: Cloudy luminosity unit such as `q(H)` or `logU`
- `lumi_value`: numeric value for the chosen unit

### `set_star(SED=None, SED_params=None, lumi_unit=None, lumi_value=None)`

Define a custom source spectrum.

- Use this for table-based or non-blackbody SEDs.
- `SED_params` may be a scalar, string, list, or tuple.

### `set_cste_density(dens=None, ff=None, others=None)`

Set a constant density model.

- `dens`: density in log(cm-3)
- `ff`: optional filling factor
- `others`: additional text appended to the Cloudy density command

### `set_dlaw(dlaw_params, ff=None)`

Set a user-defined density law.

- `dlaw_params`: one value or a sequence of values
- `ff`: optional filling factor

### `set_fudge(fudge_params=None)`

Add Cloudy `fudge factors`.

- Useful for custom correction terms or special setups.

### `set_sphere(sphere=True)`

Enable or disable spherical geometry.

- `True` adds `sphere`
- `False` removes it

### `set_iterate(n_iter=None, to_convergence=False)`

Configure iteration behavior.

- `n_iter=None` prints a plain `iterate`
- `n_iter=0` removes the iterate directive
- `n_iter>0` sets an explicit iteration count
- `to_convergence=True` requests iteration to convergence

### `set_grains(grains=None)`

Set grain commands.

- `None` clears the grain list
- A scalar appends one grain command
- A list or tuple appends several grain commands

### `set_stop(stop_criter=None)`

Add stopping criteria.

- `None` clears the list
- A scalar or sequence appends one or more stopping commands

### `read_emis_file(emis_file, N_char=14, emergent=False)`

Load emissivity labels from a text file.

- `N_char` controls how many characters are read from each line
- `emergent=True` switches to emergent emissivities

### `set_emis_tab(emis_tab_str=None, emergent=False, absolute=False)`

Set the emissivity label list directly in Python.

- `emis_tab_str` is a sequence of line labels
- `emergent=True` requests emergent emissivities
- `absolute=True` requests absolute emissivities

### `import_file(file_=None)`

Append a text file directly into the generated input.

- `None` clears the imported block

### `set_line_file(line_file=None, absolute=False, emergent=False)`

Attach a line-list file for Cloudy output.

- `absolute` and `emergent` control the output style

### `set_theta_phi(theta=None, phi=None)`

Set orientation angles for 3D-style inputs.

- `None` for both values clears the angles
- Either angle may be provided independently

### `set_abund(predef=None, elem=None, value=None, nograins=True, ab_dict=None, metals=None, metalsgrains=None, metalsdeplete=None)`

Define abundances.

- `predef`: Cloudy predefined abundance set
- `elem` and `value`: set one abundance directly
- `ab_dict`: dictionary of element abundances
- `nograins`: disable grain-related abundance handling
- `metals`, `metalsgrains`, and `metalsdeplete`: optional metal-scaling controls

### `set_other(other_str=None)`

Add raw Cloudy commands.

- Accepts one string or a sequence of strings
- `None` clears the command list

### `set_comment(comment=None)`

Add human-readable comments to the generated input.

- Accepts one string or a sequence of strings
- `None` clears the list

### `set_C3D_comment(comment=None)`

Add `#C3D` comments for 3D workflows.

- Accepts one string or a sequence of strings
- `None` clears the list

### `set_distance(dist=None, unit='kpc', linear=True)`

Set the distance used for model normalization.

- `unit` may be `kpc`, `Mpc`, `parsecs`, or `cm`
- `linear=True` writes the linear-distance form

### `set_heat_cooling(cextra=None, hextra=None)`

Store custom heating and cooling additions.

- `cextra` affects cooling
- `hextra` affects heating

### `print_input(to_file=True, verbose=False)`

Write the Cloudy input file.

- `to_file=True` writes `<model_name>.in`
- `verbose=True` echoes the generated input to the console

### `run_cloudy(dir_=None, n_proc=1, use_make=True, model_name=None, precom='')`

Run Cloudy for the configured model.

- `dir_`: directory holding the model files
- `n_proc`: number of processes
- `use_make=True`: use the make-based runner when possible
- `model_name`: optional model override
- `precom`: optional prefix text before the Cloudy command

## Practical Example

```python
import pyCloudy as pc

c_input = pc.CloudyInput("models/example")
c_input.set_BB(Teff=40000, lumi_unit="q(H)", lumi_value=47)
c_input.set_cste_density(2.0)
c_input.set_radius(r_in=17.3)
c_input.set_abund(predef="ism")
c_input.set_grains(["silicate", "graphite"])
c_input.set_stop(["temperature 4000 K", "ionization parameter -2"])
c_input.print_input(to_file=True)
```

## Notes

- This page is a practical summary, not the source of truth.
- If a method changes in `pyCloudy/c1d/cloudy_model.py`, the implementation wins.
