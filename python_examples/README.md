# Python Examples

This directory contains plain Python exports of the example notebooks stored in
`pyCloudy/docs/`.

## Source Map

Each example lives in its own directory named after the script.

- `Using_pyCloudy_1/Using_pyCloudy_1.py` from `pyCloudy/docs/Using_pyCloudy_1.ipynb`
- `Using_pyCloudy_2/Using_pyCloudy_2.py` from `pyCloudy/docs/Using_pyCloudy_2.ipynb`
- `Using_pyCloudy_3/Using_pyCloudy_3.py` from `pyCloudy/docs/Using_pyCloudy_3.ipynb`
- `Using_pyCloudy_4/Using_pyCloudy_4.py` from `pyCloudy/docs/Using_pyCloudy_4.ipynb`
- `Using_pyCloudy_MdB/Using_pyCloudy_MdB.py` from `pyCloudy/docs/Using_pyCloudy_MdB.ipynb`
- `Using_pyCloudy_with_PyNeb/Using_pyCloudy_with_PyNeb.py` from `pyCloudy/docs/Using_pyCloudy_with_PyNeb.ipynb`
- `shocks/shocks.py` from `pyCloudy/docs/shocks.ipynb`
- `quickstart/quickstart.py` from the README quickstart example

## How To Run

Run any example with:

```bash
python python_examples/Using_pyCloudy_1/Using_pyCloudy_1.py
```

For the quickstart example:

```bash
python python_examples/quickstart/quickstart.py
```

Many of the examples expect:

- a working Cloudy installation
- `pc.config.cloudy_exe` set to the correct Cloudy binary
- optional packages such as `matplotlib`, `PyNeb`, or `pandas`

## Notes

- The converted scripts are meant as readable command-line examples.
- The original notebooks remain in `pyCloudy/docs/` if you prefer the
  interactive version.
- Some scripts write files to temporary or model directories defined inside the
  example code, so you may want to review the paths before running them.
