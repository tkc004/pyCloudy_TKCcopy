# RHD Spherical Radial Profiles Example

This example is based on
[`Using_pyCloudy_3_spherical_profiles`](../python_examples/Using_pyCloudy_3_spherical_profiles/Using_pyCloudy_3_spherical_profiles.py),
but reads the RHD profile file
[`radial_profile_rhd.csv`](../python_examples/Using_pyCloudy_3_spherical_profiles_rhd/radial_profile_rhd.csv).

It uses the same spherical density and velocity workflow, diagnostic maps, RGB
figures, and line profiles as the original example. The ionizing photon rate
is set to:

```text
q(H) = 10^49 s^-1
```

## Running

From the repository root:

```bash
python python_examples/Using_pyCloudy_3_spherical_profiles_rhd/Using_pyCloudy_3_spherical_profiles_rhd.py
```

The script requires a working Cloudy executable and writes model files under
`python_examples/Using_pyCloudy_3_spherical_profiles_rhd/temp_models/` and
figures under its `figures/` directory.

## Profile File

The CSV header is:

```text
RADIUS_PC,VELOCITY_KMS,DENSITY_CM3
```

Radius is in pc, velocity is in km/s, and hydrogen density is in cm^-3. The
profile is interpolated onto the C3D grid for the user-velocity line profile.

The line-profile figure uses the velocity profile supplied by the RHD CSV,
interpolated onto the C3D grid.
