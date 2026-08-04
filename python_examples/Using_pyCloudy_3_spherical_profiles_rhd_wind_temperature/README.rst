RHD stellar-wind radial-profile example with temperature
=========================================================

This example is based on
``Using_pyCloudy_3_spherical_profiles_rhd_wind`` and reads CSV files from the
``radial_profiles/`` directory. Each file uses the ``RADIUS_PC``,
``VELOCITY_KMS``, ``DENSITY_CM3``, and ``TEMP_K`` columns. The density and
temperature profiles are passed to Cloudy as radial ``dlaw table radius`` and
``tlaw table radius`` laws, respectively. ``TEMP_K`` values are converted to
Cloudy's logarithmic temperature-table values. The ionizing photon rate is
``q(H)=10^49 s^-1``.

No Cloudy ``wind`` command is added. The velocity profile is passed to C3D as
a user velocity field, as in the wind example.

Run it from the repository root with:

.. code-block:: bash

python python_examples/Using_pyCloudy_3_spherical_profiles_rhd_wind_temperature/Using_pyCloudy_3_spherical_profiles_rhd_wind_temperature.py

By default, all CSV files in ``radial_profiles/`` are run except
``radial_profile_0Myr.csv``. Run one file explicitly with:

.. code-block:: bash

   python python_examples/Using_pyCloudy_3_spherical_profiles_rhd_wind_temperature/Using_pyCloudy_3_spherical_profiles_rhd_wind_temperature.py --profile python_examples/Using_pyCloudy_3_spherical_profiles_rhd_wind_temperature/radial_profiles/radial_profile_1Myr.csv

The X-ray continuum calculation is disabled by default. Enable it explicitly
with ``--xray`` to save the sampled per-zone continuum and write an X-ray
figure for each profile.

Each CSV has isolated output directories under ``figures/<csv-stem>/`` and
``temp_models/<csv-stem>/``. The figures include input profiles, Cloudy
temperature and emissivity radial profiles, line profiles, derived maps, RGB
images, and diagnostic scatter plots.


Parameters of the simulations:
  radiative_transfer_source_photon_rate:
    value: 1.0e49
    unit: 1/s
  wind_mass_loss_rate:
    value: 1.0e-6
    unit: Msun/yr
  wind_velocity:
    value: 1000.0
    unit: km/s
  wind_temperature:
    value: 100.0
    unit: K
  hydrogen_number_density:
    value: 100.0
    unit: 1/cm**3
  initial_temperature:
    value: 100.0
    unit: K
