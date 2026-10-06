from pathlib import Path

from dada_solver.sizing.configuration import load_sizing_problem
from dada_solver.sizing.design import DesignParameter
from dada_solver.sizing.objectives import MinimizeTotalSweptVolume


def test_sizing_configuration_resolves_relative_base(tmp_path):
    base = tmp_path / "machine.toml"
    base.write_bytes((Path(__file__).parent / "data" / "sizing_machine.toml").read_bytes())
    path = tmp_path / "sizing.toml"
    path.write_text('''[problem]
base_configuration = "machine.toml"
[[variables]]
parameter = "small_swept_volume"
lower_bound = 0.0001
upper_bound = 0.001
initial_value = 0.0002
[objective]
type = "minimize_total_swept_volume"
[[constraints]]
type = "maximum_pressure"
limit = 200000.0
[optimizer]
objective_scale = 0.001
unavailable_objective_penalty = 1000.0
unavailable_constraint_margin = -1.0
maximum_iterations = 5
function_tolerance = 0.000001
[optimizer.constraint_scales]
maximum_pressure = 100000.0
''')

    loaded = load_sizing_problem(path)

    assert loaded.base_configuration_path == base
    assert isinstance(loaded.problem.objective, MinimizeTotalSweptVolume)
    assert tuple(v.parameter for v in loaded.problem.variables) == (
        DesignParameter.SMALL_SWEPT_VOLUME,
    )
    assert loaded.problem.variables[0].initial_value == 0.0002
    assert tuple(c.name for c in loaded.problem.constraints) == ("maximum_pressure",)
    assert loaded.optimization_settings.constraint_scales == {"maximum_pressure": 100000.0}
    assert loaded.cooling_cell_scenario is None
