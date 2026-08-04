RHD stellar-wind radial-profile example with temperature
=========================================================

This example is based on
``Using_pyCloudy_3_spherical_profiles_rhd_wind`` and reads the local
``radial_profile_rhd_wind.csv`` file. It uses the ``RADIUS_PC``,
``VELOCITY_KMS``, ``DENSITY_CM3``, and ``TEMP_K`` columns. The density and
temperature profiles are passed to Cloudy as radial ``dlaw table radius`` and
``tlaw table radius`` laws, respectively. ``TEMP_K`` values are converted to
Cloudy's logarithmic temperature-table values.

No Cloudy ``wind`` command is added. The velocity profile is passed to C3D as
a user velocity field, as in the wind example.

Run it from the repository root with:

.. code-block:: bash

   python python_examples/Using_pyCloudy_3_spherical_profiles_rhd_wind_temperature/Using_pyCloudy_3_spherical_profiles_rhd_wind_temperature.py

The X-ray continuum calculation is disabled by default. Enable it explicitly
with ``--xray`` to save the per-zone continuum and write the X-ray figure.

Figures are written under ``figures/`` in this example directory.
