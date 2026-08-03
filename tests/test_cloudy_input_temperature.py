from pathlib import Path

import pyCloudy as pc


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
