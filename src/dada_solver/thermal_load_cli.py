"""Command-line reporting for external cooling-cell loads."""

from __future__ import annotations

import argparse
from pathlib import Path

from dada_solver.thermal_load_configuration import load_cooling_cell_scenario
from dada_solver.thermal_load_reporting import format_cooling_cell_requirements


def main(arguments: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Calculate an ideal external cooling-cell thermal load."
    )
    parser.add_argument("configuration", type=Path, help="Cooling-load TOML file")
    options = parser.parse_args(arguments)
    try:
        scenario = load_cooling_cell_scenario(options.configuration)
        print(format_cooling_cell_requirements(scenario))
        return 0
    except (OSError, ValueError) as error:
        print(f"thermal_load_error = {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
