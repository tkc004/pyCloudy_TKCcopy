#!/usr/bin/env python
# coding: utf-8

import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir))
sys.path.insert(0, str(script_dir.parent))

from stromgren_plots import save_all_plots
from stromgren_run import run_model
from stromgren_config import DEFAULT_MODE


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_MODE
    mod, abund, context = run_model(mode, script_dir)
    save_all_plots(mod, mode, context["figure_dir"], abund)


if __name__ == "__main__":
    main()
