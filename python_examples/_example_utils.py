from pathlib import Path

import yaml


def load_yaml(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def find_cloudy_exe(start_dir):
    target = Path("Cloudy_exe/Cloudy/c22.02/source/cloudy.exe")
    start_dir = Path(start_dir)
    for current in [start_dir, *start_dir.parents]:
        candidate = current / target
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find Cloudy executable under {start_dir}")


def save_fig(fig, path, dpi=150):
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
