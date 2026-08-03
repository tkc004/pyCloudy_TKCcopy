RHD stellar-wind radial-profile example
=======================================

This example is identical to
``Using_pyCloudy_3_spherical_profiles_rhd`` but reads
``radial_profile_rhd_wind.csv`` from that example's directory. The supplied
profile represents a stellar wind with a mass-loss rate of
``10^-6 Msun/yr`` and a wind speed of ``1000 km/s``.

Only the ``RADIUS_PC``, ``VELOCITY_KMS``, and ``DENSITY_CM3`` columns are used.
The additional ``TEMP_K`` column in the CSV is intentionally ignored for now.
No Cloudy wind command is added; the radial velocity and density fields come
from the external profile file.

Run it from the repository root with:

.. code-block:: bash

   python python_examples/Using_pyCloudy_3_spherical_profiles_rhd_wind/Using_pyCloudy_3_spherical_profiles_rhd_wind.py

Figures are written under ``figures/`` in this example directory. In addition
to the standard C3D diagnostics, the script writes Cloudy radial electron
temperature and line-emissivity profiles.
