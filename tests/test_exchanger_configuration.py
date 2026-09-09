from pathlib import Path

from dada_solver.exchangers.configuration import (
    load_exchanger_screening_configuration,
)
from dada_solver.exchangers.coupled_configuration import (
    load_coupled_exchanger_problem,
)


def test_load_exchanger_screening_example() -> None:
    path = Path(__file__).parents[1] / "examples" / "exchanger_screening_example.toml"

    loaded = load_exchanger_screening_configuration(path)

    assert loaded.example_data
    assert loaded.channel_counts == (20, 40, 80)
    assert loaded.requirements.maximum_mach_number == 0.3
    assert loaded.requirements.maximum_pressure_drop_fraction == 0.05


def test_load_coupled_exchanger_example() -> None:
    path = Path(__file__).parents[1] / "examples" / "exchanger_coupled_example.toml"

    loaded = load_coupled_exchanger_problem(path)

    assert loaded.example_data
    assert loaded.cold_search.channel_counts == (20, 40, 80)
    assert loaded.settings.maximum_iterations == 6
