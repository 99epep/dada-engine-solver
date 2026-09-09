"""Command-line entry point for coupled cycle/exchanger sizing."""

from __future__ import annotations

import argparse
from pathlib import Path

from dada_solver.exchangers.coupled_configuration import (
    load_coupled_exchanger_problem,
)
from dada_solver.exchangers.coupled_reporting import format_coupled_exchanger_report
from dada_solver.exchangers.coupling import CoupledExchangerSizer, CoupledSizingStatus


def main(arguments: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Iterate a DADA periodic cycle with geometric exchanger sizing."
    )
    parser.add_argument("configuration", type=Path, help="Coupled sizing TOML file")
    options = parser.parse_args(arguments)
    try:
        loaded = load_coupled_exchanger_problem(options.configuration)
        result = CoupledExchangerSizer(
            loaded.cold_search,
            loaded.hot_search,
            loaded.settings,
            progress_callback=lambda message: print(message, flush=True),
        ).solve(loaded.simulation)
        print(format_coupled_exchanger_report(result))
        return 0 if result.status is CoupledSizingStatus.CONVERGED else 2
    except (OSError, ValueError) as error:
        print(f"coupled_exchanger_error = {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
