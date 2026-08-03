from pathlib import Path

import numpy as np
import pyCloudy as pc
import pytest


def test_set_tlaw_table_writes_temperature_pairs(tmp_path):
    model_path = Path(tmp_path) / "temperature_table"
    cloudy_input = pc.CloudyInput(str(model_path))
    cloudy_input.set_tlaw([(17.0, 10000.0), (18.0, 8000.0)])
    cloudy_input.print_input(to_file=True)

    lines = model_path.with_suffix(".in").read_text().splitlines()
    start = lines.index("tlaw table radius")
    assert lines[start + 1:start + 4] == [
        "continue 17.0 10000.0",
        "continue 18.0 8000.0",
        "end of tlaw",
    ]


def test_set_tlaw_direct_command(tmp_path):
    model_path = Path(tmp_path) / "temperature_law"
    cloudy_input = pc.CloudyInput(str(model_path))
    cloudy_input.set_tlaw("DB96")
    cloudy_input.print_input(to_file=True)

    assert "tlaw DB96" in model_path.with_suffix(".in").read_text().splitlines()


def test_set_cste_temperature_writes_constant_temperature(tmp_path):
    model_path = Path(tmp_path) / "constant_temperature"
    cloudy_input = pc.CloudyInput(str(model_path))
    cloudy_input.set_cste_temperature(10000.0)
    cloudy_input.print_input(to_file=True)

    assert "constant temperature 10000 K" in model_path.with_suffix(".in").read_text().splitlines()


def test_cloudy_constant_temperature_output(tmp_path):
    cloudy_exe = (
        Path(__file__).resolve().parents[1]
        / "Cloudy_exe" / "Cloudy" / "c22.02" / "source" / "cloudy.exe"
    )
    if not cloudy_exe.exists():
        pytest.skip("Bundled Cloudy executable is not available")

    model_path = Path(tmp_path) / "constant_temperature_model"
    cloudy_input = pc.CloudyInput(str(model_path))
    cloudy_input.set_BB(40000.0, "q(H)", 47.0)
    cloudy_input.set_radius(r_in=17.0)
    cloudy_input.set_cste_density(2.0)
    cloudy_input.set_cste_temperature(10000.0)
    cloudy_input.set_stop("radius 18")
    cloudy_input.print_input(to_file=True)
    pc.config.cloudy_exe = str(cloudy_exe)
    cloudy_input.run_cloudy(use_make=False)

    physical_conditions = np.loadtxt(
        model_path.with_suffix(".phy"),
        comments="#",
        usecols=(0, 1, 2, 3, 6),
    )
    radii, temperatures, hydrogen_densities, electron_densities, filling_factors = physical_conditions.T
    assert temperatures.size > 0
    assert np.allclose(temperatures, 10000.0, rtol=0.0, atol=1.0)
    assert np.all(np.diff(radii) > 0.0)
    assert np.allclose(hydrogen_densities, 100.0, rtol=0.0, atol=0.1)
    assert np.all(np.isfinite(electron_densities))
    assert np.all(electron_densities > 0.0)
    assert np.allclose(filling_factors, 1.0, rtol=0.0, atol=1.0e-6)
