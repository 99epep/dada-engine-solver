"""Run both motion families through the unchanged periodic thermodynamic evaluator.

This validates architecture, not optimality or a useful-power hardware design.
Run from the repository root with PYTHONPATH=src.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
from dada_solver.configuration import load_simulation_configuration
from dada_solver.sizing.evaluator import evaluate_configuration

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--maximum-cycles', type=int, default=50)
parser.add_argument('--output', type=Path, default=Path('outputs/pluggable_kinematics_smoke.json'))
args = parser.parse_args()
results = {}
for family, filename in [('four_bar', 'motor_demonstrator_original_325c.toml'),
                         ('free', 'motor_free_kinematics.toml')]:
    config = load_simulation_configuration(Path('examples')/filename)
    config = replace(config, numerical=replace(config.numerical, maximum_cycles=args.maximum_cycles))
    result = evaluate_configuration(config)
    results[family] = dict(status=result.status.value, message=result.periodic.message,
        cycles=len(result.periodic.history),
        indicated_power_w=None if result.performance is None else result.performance.gas_power,
        configuration=str(Path('examples')/filename),
        scope='constant_UA_orifice_architecture_smoke; not_hardware_optimization')
    print(f'{family}: {results[family]}', flush=True)
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(results, indent=2))
