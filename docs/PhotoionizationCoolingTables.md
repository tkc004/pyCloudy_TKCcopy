# Photoionization equilibrium cooling tables

`tools/generate_photoionization_cooling_table.py` generates optically thin
photoionization-equilibrium cooling tables by running Cloudy on a grid of
temperature, hydrogen density, ionization parameter, and metallicity.

The script uses a blackbody spectrum by default. The current command-line
interface exposes its effective temperature with `--teff`; for example,
`--teff 30000` uses a 30,000 K blackbody.

## Requirements

- Python with `numpy`, `h5py`, and the local `pyCloudy` package available.
- A Cloudy executable compiled for the host machine. Pass its path with
  `--cloudy-exe` when it is not found automatically.
- Enough disk space for the HDF5 tables and, if requested, retained Cloudy
  output files.

Run the command from the project directory so that the default `data/` and
`cooling_models/` directories are created there.

## Example: one-cell test

Use a one-cell run to check the Cloudy executable and input recipe before
starting a large grid:

```bash
python tools/generate_photoionization_cooling_table.py \
  --output single_model_test.h5 \
  --teff 30000 \
  --logT-min 4 --logT-max 4 --nT 1 \
  --lognH-min 0 --lognH-max 0 --nnH 1 \
  --logU-min -2 --logU-max -2 --nU 1 \
  --metallicities 1.0 \
  --workers 1 \
  --overwrite
```

This writes two files under `data/`:

```text
data/single_model_test_Z1_HHe.h5
data/single_model_test_Z1_metals.h5
```

The metallicity axis is retained as a one-element axis in each file. Heating
and cooling rates are stored in `erg cm^-3 s^-1` in the datasets
`heating_erg_cm-3_s` and `cooling_erg_cm-3_s`, respectively, with axis order
`[metallicity, temperature, hydrogen_density, logU]`.

## Full-grid example

```bash
python tools/generate_photoionization_cooling_table.py \
  --output photoionization_cooling_30kK.h5 \
  --teff 30000 \
  --logT-min 3.7 --logT-max 4.5 --nT 64 \
  --lognH-min 0 --lognH-max 4 --nnH 32 \
  --logU-min -4 --logU-max -1 --nU 32 \
  --metallicities 0.1 0.3 1.0 2.0 \
  --workers 4 \
  --cloudy-exe /path/to/cloudy.exe
```

Each metallicity produces an H/He file and a metal-cooling file. The metal
rate is calculated as the total rate minus the H/He rate. The grid is
thread-parallelized, while each Cloudy process is restricted to one OpenMP
thread through `OMP_NUM_THREADS=1`.

## Cloudy model setup

Each cell is a constant-temperature, plane-parallel slab. The generated input
uses:

- a blackbody spectrum and the requested ionization parameter;
- constant hydrogen density and `constant temperature T K linear`;
- `stop zone 1` and `iterate to convergence`;
- H II-region abundances with no grains;
- `no molecules` for both the full-metal and H/He calculations.

The H/He calculation sets elements from lithium through zinc to abundance
`-30`, turns those elements off, and adds `metals 1e-30`. This suppresses
metal cooling while preserving the separate hydrogen and helium contribution.

## Resume and retained Cloudy files

The script writes a checkpoint CSV in `data/` while running. It stores H/He
and metal heating and cooling columns for each completed cell. Resume a
stopped run with the same grid and `--resume`:

```bash
python tools/generate_photoionization_cooling_table.py \
  --output photoionization_cooling_30kK.h5 \
  ... \
  --resume
```

By default, a cell's `.in`, `.out`, `.cool`, `.phy`, `.heat`, `.cont`, and
other Cloudy files are deleted after its cooling rate has been read. Use
`--keep-cloudy-files` to retain the model directories for inspection:

```bash
python tools/generate_photoionization_cooling_table.py \
  --output single_model_keep_test.h5 \
  --teff 30000 \
  --logT-min 4 --logT-max 4 --nT 1 \
  --lognH-min 0 --lognH-max 0 --nnH 1 \
  --logU-min -2 --logU-max -2 --nU 1 \
  --metallicities 1.0 \
  --workers 1 \
  --overwrite \
  --keep-cloudy-files
```

Retained files are placed in `cooling_models/<model-name>/`. Ionization
failure cells are excluded from the HDF5 table as `NaN` and recorded in the
checkpoint. `--retry-rejected` retries such cells on a subsequent resumed
run.

## Checking generated tables

Use the companion checker after a run:

```bash
python tools/check_photoionization_cooling_tables.py \
  --data-dir data \
  --stem single_model_test \
  --metallicity 1.0
```

The checker verifies the HDF5 structure, coordinate axes, finite values, and
basic cooling-rate sanity. It accepts `NaN` cells because those represent
Cloudy cells explicitly excluded after ionization-convergence failures.

## Run metadata

Each HDF5 file stores reproducibility metadata as HDF5 attributes, including
the pyCloudy version, Cloudy version and executable path, spectrum and
effective temperature, command line, grid definition, geometry, stop
criterion, iteration mode, and H/He abundance recipe.
