"""Command-line entry point for geometric exchanger screening."""

from __future__ import annotations

import argparse
from pathlib import Path

from dada_solver.exchangers.configuration import (
    load_exchanger_screening_configuration,
)
from dada_solver.exchangers.reporting import (
    format_exchanger_optimization_report,
    format_exchanger_screening_report,
)
from dada_solver.exchangers.sizing import screen_parallel_channel_designs
from dada_solver.exchangers.optimization import (
    ContinuousGeometryBounds,
    ExchangerObjective,
    optimize_parallel_channel_geometry,
)


def main(arguments: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Screen DADA parallel-channel heat-exchanger geometries."
    )
    parser.add_argument("configuration", type=Path, help="Exchanger TOML file")
    parser.add_argument(
        "--optimize",
        choices=[item.value for item in ExchangerObjective],
        help="Refine feasible grid seeds for the selected objective",
    )
    options = parser.parse_args(arguments)
    try:
        config = load_exchanger_screening_configuration(options.configuration)
        common = dict(
            channel_counts=config.channel_counts,
            operating_point=config.operating_point,
            gas=config.gas,
            transport=config.transport,
            thermal_resistances=config.thermal_resistances,
            requirements=config.requirements,
            surface_roughness=config.surface_roughness,
            minor_loss_coefficient=config.minor_loss_coefficient,
        )
        if options.optimize is not None:
            result = optimize_parallel_channel_geometry(
                **common,
                bounds=ContinuousGeometryBounds(
                    min(config.channel_widths),
                    max(config.channel_widths),
                    min(config.channel_heights),
                    max(config.channel_heights),
                    min(config.channel_lengths),
                    max(config.channel_lengths),
                ),
                seed_widths=config.channel_widths,
                seed_heights=config.channel_heights,
                seed_lengths=config.channel_lengths,
                objective=ExchangerObjective(options.optimize),
            )
            print(format_exchanger_optimization_report(result))
            return 0 if result.success else 2
        results = screen_parallel_channel_designs(
            **common,
            channel_widths=config.channel_widths,
            channel_heights=config.channel_heights,
            channel_lengths=config.channel_lengths,
        )
        print(format_exchanger_screening_report(results))
        return 0 if any(item.feasible for item in results) else 2
    except (OSError, ValueError) as error:
        print(f"exchanger_error = {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
