pyCloudy
========

pyCloudy is a Python library for working with the input and output files of
Cloudy, the photoionization and plasma simulation code by Gary Ferland and
collaborators.

The package focuses on three practical tasks:

* building Cloudy input files from Python
* reading and analyzing Cloudy output files
* combining multiple 1D runs into 3D nebula models

It also includes helpers for the 3MdB database, reddening corrections, common
astronomy utilities, and small numerical tools used throughout the workflow.

Project home
------------

* Website: https://sites.google.com/site/pycloudy/
* Repository: https://github.com/Morisset/pyCloudy
* Manual and notebooks: see the bundled documentation in ``pyCloudy/docs/``
* Documentation landing page: ``docs/README.md``

Requirements
------------

pyCloudy depends on:

* Python 3.8 or newer
* numpy
* scipy

Optional features become available when extra packages are installed:

* matplotlib for plotting and triangulation support
* astropy for FITS support
* Pillow for RGB image helpers
* pandas for table handling
* PyMySQL or MySQLdb for database access
* pyneb for atomic and emission line utilities

You also need a working Cloudy installation if you want to run models from
pyCloudy.

Installation
------------

Install from PyPI:

.. code-block:: bash

   pip install pyCloudy

Install with optional features:

.. code-block:: bash

   pip install "pyCloudy[full]"

Or install the development version from GitHub:

.. code-block:: bash

   pip install -U git+https://github.com/Morisset/pyCloudy.git

If you already have the source tree, you can also install it in editable mode:

.. code-block:: bash

   pip install -e .

Configuring Cloudy
------------------

pyCloudy looks for the Cloudy executable in ``pc.config.cloudy_exe``. You can
set it in Python:

.. code-block:: python

   import pyCloudy as pc

   pc.config.cloudy_exe = "/usr/local/Cloudy/c23.01/source/cloudy.exe"

Or define the ``CLOUDY_EXE`` environment variable before importing pyCloudy.

The package also contains a small set of preset paths in ``pc.config.cloudy_dict``
for common Cloudy versions.

Quick start
-----------

The two main entry points are ``CloudyInput`` for building a model and
``CloudyModel`` for reading the results.

Create and run a model:

.. code-block:: python

   import pyCloudy as pc

   pc.config.cloudy_exe = "/usr/local/Cloudy/c23.01/source/cloudy.exe"

   c_input = pc.CloudyInput("models/M17")
   c_input.set_BB(Teff=40000, lumi_unit="q(H)", lumi_value=47)
   c_input.set_cste_density(2.0)
   c_input.set_radius(r_in=17.3)
   c_input.set_abund(predef="ism")
   c_input.set_sphere(True)
   c_input.set_iterate()
   c_input.set_distance(dist=1.0, unit="kpc", linear=True)
   c_input.print_input(to_file=True)
   c_input.run_cloudy()

Read a finished model:

.. code-block:: python

   import pyCloudy as pc

   model = pc.CloudyModel("models/M17")
   print(model.n_zones)
   print(model.te[:5])
   print(model.get_Hb_SB())

If you have several 1D models sampled along different directions or physical
conditions, you can assemble them into a 3D representation:

.. code-block:: python

   import pyCloudy as pc

   c3d = pc.C3D([pc.CloudyModel("models/M17")], dims=51)
   print(c3d.get_emis_list())

Main modules
------------

``pyCloudy.c1d``
  Cloudy input generation and 1D model analysis.

  * ``CloudyInput`` writes Cloudy input files and can run Cloudy.
  * ``CloudyModel`` reads the standard Cloudy output files and exposes
    derived quantities such as temperatures, densities, ionic fractions,
    emissivities, line intensities, continua, and equivalent widths.
  * ``load_models`` loads multiple models from a list or pattern.

``pyCloudy.c3d``
  3D model construction from a set of 1D Cloudy runs.

  * ``C3D`` builds interpolated 3D fields and can return volume-integrated
    physical quantities, emission maps, ion distributions, and line profiles.
  * ``CubCoord`` manages the geometry, coordinate grid, and velocity laws.

``pyCloudy.db``
  3MdB database helpers.

  * ``MdB`` provides a database access layer for querying the 3MdB tables.
  * ``MdB_subproc`` offers a subprocess-based fallback in environments where
    direct connectors are not available.

``pyCloudy.utils``
  Shared utilities.

  * ``physics`` contains physical constants and common conversions.
  * ``misc`` contains general helpers for I/O, geometry, labels, and data
    transformations.
  * ``red_corr`` implements reddening correction laws.
  * ``astro`` contains astronomy-specific conversions.

Typical workflow
----------------

1. Prepare a ``CloudyInput`` instance.
2. Set the ionizing source, density law, abundances, geometry, and stopping
   criteria.
3. Write the input file and run Cloudy.
4. Load the results with ``CloudyModel``.
5. Inspect line emissivities, ionic fractions, temperatures, and continua.
6. Optionally combine multiple 1D models into a 3D cube with ``C3D``.

Working with emission lines
---------------------------

The package supports both Cloudy output lines and user-defined emission lines.
Common helpers include:

* ``get_line``
* ``get_emis``
* ``get_emis_vol``
* ``print_lines``
* ``emis_from_pyneb``
* ``add_emis_from_pyneb``

If you are using PyNeb, install the optional dependency and then import the
line list or atomic data you need before adding lines to the model.

Database support
----------------

The ``pyCloudy.db`` package was designed for the 3MdB database workflow. The
main class is ``MdB``. It can connect to a database, execute SQL queries, and
return results in a Python-friendly format.

The bundled ``pyCloudy/3MdB_17`` directory contains example scripts and notes
for working with Cloudy version 17 datasets.

Documentation and examples
--------------------------

The repository bundles a more extensive set of notebooks and example PDFs in
``pyCloudy/docs/``. These cover:

* basic 1D usage
* PyNeb integration
* 3D model construction
* 3MdB database workflows
* shock-related examples

If you are looking for a starting point, open the notebooks in that directory
in this order:

1. ``Using_pyCloudy_1.ipynb``
2. ``Using_pyCloudy_2.ipynb``
3. ``Using_pyCloudy_3.ipynb``
4. ``Using_pyCloudy_4.ipynb``
5. ``Using_pyCloudy_with_PyNeb.ipynb``
6. ``Using_pyCloudy_MdB.ipynb``

Testing
-------

The test suite uses example model outputs stored in ``tests/models``.
Some tests expect a local Cloudy installation and may require you to update the
Cloudy executable path before running them.

Run the tests with:

.. code-block:: bash

   pytest

Notes and limitations
---------------------

* pyCloudy does not replace Cloudy. It is a companion library that helps you
  build, run, and analyze Cloudy models.
* Some features depend on optional scientific Python packages.
* Database and model-run helpers assume you have local access to the required
  infrastructure, such as a Cloudy binary or a MySQL-compatible server.

Warranty
--------

pyCloudy is provided as is, without warranty. Use it for science, but verify
the results with your own checks and domain knowledge.

Acknowledgements
----------------

This project is partly supported by grants DGAPA/PAPIIT-107215 and
CONACyT-CB2015-254132.
