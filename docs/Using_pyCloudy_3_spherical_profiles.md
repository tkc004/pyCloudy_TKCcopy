# Spherical Radial Profiles Example

This example extends [`Using_pyCloudy_3`](../python_examples/Using_pyCloudy_3/Using_pyCloudy_3.py)
to a spherically symmetric nebula controlled by external radial profiles.

The example reads:

- radius in pc
- radial velocity in km/s
- hydrogen density in atom cm^-3

from [`radial_profiles.csv`](../python_examples/Using_pyCloudy_3_spherical_profiles/radial_profiles.csv).

## Profile File

The CSV header must contain these exact column names:

```text
RADIUS_PC,VELOCITY_KMS,DENSITY_CM3
```

Example:

```text
RADIUS_PC,VELOCITY_KMS,DENSITY_CM3
0.05,0.0,1200.0
0.10,8.0,1100.0
0.20,18.0,900.0
0.40,30.0,650.0
0.70,42.0,450.0
1.00,50.0,300.0
```

Radius values must be positive and strictly increasing. Density values must be
positive. Velocity values may be positive or negative.

## Running

From the repository root, run:

```bash
python python_examples/Using_pyCloudy_3_spherical_profiles/Using_pyCloudy_3_spherical_profiles.py
```

The example requires a working Cloudy executable and the plotting dependencies
used by pyCloudy. It writes Cloudy model files under `temp_models/` and figures
under `figures/` in the example directory.

## How Profiles Are Used

The density profile is passed to Cloudy with a radial `dlaw table`. Cloudy
interpolates the supplied log-radius/log-density pairs when solving the
spherical model.

The velocity profile is interpolated onto the C3D radius grid and converted to
Cartesian velocity components. The velocity is zero at the origin and is
held at the final supplied value beyond the last profile point.

## Outputs

The example produces the input profile plot plus the same diagnostic outputs as
the original 3D example:

- [`input_radial_profiles.png`](../python_examples/Using_pyCloudy_3_spherical_profiles/figures/input_radial_profiles.png)
- [`profile_default.png`](../python_examples/Using_pyCloudy_3_spherical_profiles/figures/profile_default.png)
- [`profile_user_velocity.png`](../python_examples/Using_pyCloudy_3_spherical_profiles/figures/profile_user_velocity.png)
- [`derived_maps.png`](../python_examples/Using_pyCloudy_3_spherical_profiles/figures/derived_maps.png)
- [`rgb_compact.png`](../python_examples/Using_pyCloudy_3_spherical_profiles/figures/rgb_compact.png)
- [`rgb_with_profiles.png`](../python_examples/Using_pyCloudy_3_spherical_profiles/figures/rgb_with_profiles.png)
- [`diagnostic_scatter.png`](../python_examples/Using_pyCloudy_3_spherical_profiles/figures/diagnostic_scatter.png)

The RGB figures use `[NII] 6584` for red, `[OIII] 5007` for green, and H-beta
4861 for blue. The RGB-with-profiles figure overlays spatially aligned line
profiles on the projected image; its axes are in pc.
