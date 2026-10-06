"""Tests exercise maintained code without loading application scripts or data."""
import ast
from pathlib import Path


def test_python_tests_do_not_depend_on_application_directory():
    directory = 'exam' + 'ples'
    script_modules = {
        'experimental_periodic_interpolation', 'diagnose_periodic_map',
        'numba_wall_prototype', 'compare_slider_crank_motion',
        'search_compact_motor_motion', 'compact_coupler_geometry',
        'search_small_rocker_e', 'prepare_pedal_cell_valve_regions',
        'benchmark_solver_acceleration', 'refine_motor_four_stage_hx9d_variable_gas',
        'optimize_motor_hybrid_compact_260k', 'optimize_sixbar_pairs_thermo5d_3952_hlat25',
        'search_small_sixbar_stage2f', 'search_large_sixbar_stageL1',
    }
    violations = []
    for path in Path(__file__).parent.rglob('*.py'):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                # A bare forbidden-path token in a documentary guard is not a dependency.
                parts = node.value.replace('\\', '/').split('/')
                if node.value == directory or (directory in parts and parts[-1]):
                    violations.append((path.name, node.lineno, node.value))
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                modules = ([alias.name for alias in node.names] if isinstance(node, ast.Import)
                           else [node.module or ''])
                for module in modules:
                    if module.split('.')[0] in script_modules | {directory}:
                        violations.append((path.name, node.lineno, module))
    assert not violations, violations
