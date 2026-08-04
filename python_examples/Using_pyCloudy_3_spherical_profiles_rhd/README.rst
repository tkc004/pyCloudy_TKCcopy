RHD spherical radial-profile example
====================================

This example is based on ``Using_pyCloudy_3_spherical_profiles`` but reads
``radial_profile_rhd.csv`` and uses an ionizing photon rate of
``q(H) = 10^49 s^-1``.

The CSV radial profile comes from the 1D RadHydropy
``DynamicStromgrenSpherePhotoheating20pc1D`` simulation:
https://github.com/tkc004/RadHydropy/tree/main/example/DynamicStromgrenSpherePhotoheating20pc1D

Run it from the repository root with:

.. code-block:: bash

   python python_examples/Using_pyCloudy_3_spherical_profiles_rhd/Using_pyCloudy_3_spherical_profiles_rhd.py

The CSV file must contain the columns ``RADIUS_PC``, ``VELOCITY_KMS``, and
``DENSITY_CM3``. The script writes Cloudy model files under ``temp_models/``
and diagnostic figures under ``figures/`` in this example directory.

For the complete spherical-profile workflow description, see
``docs/Using_pyCloudy_3_spherical_profiles.md``.
