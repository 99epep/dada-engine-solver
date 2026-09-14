"""Evaluate independent Stage 2F/L1 six-bars with unchanged Stage K2 physics.

Run: PYTHONPATH=src python3 examples/evaluate_motor_champion_sixbar.py
Lengths remain dimensionless; this example does not infer physical crank radii.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import time

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.evaluator import json_values
from dada_solver.six_bar import IndependentSixBarVolumeKinematics, load_six_bar_mechanism
from compare_motor_motion_laws_stage7A5 import (
    A5_CAMPAIGN, ROOT, _candidate_design_and_mass, _evaluate, _same_inventory_design,
)
from evaluate_motor_champion_four_bar_k2 import (
    K2_REPORT, _motion_fit, _rescaled_initial_state,
)
from optimize_motor_exchanger_asymmetry_stageK1 import _hardware_metrics, _scaled_exchanger

SMALL = ROOT / 'outputs/small_sixbar_stage2f_r4_freeH_tightaxis.json'
LARGE = ROOT / 'outputs/large_sixbar_stageL1.json'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--hydraulic-loss-multiplier', type=float, default=1.0,
                        help='Scale tube, header and valve loss coefficients at fixed flow.')
    args = parser.parse_args()
    factor = args.hydraulic_loss_multiplier
    if not math.isfinite(factor) or factor <= 0:
        parser.error('Hydraulic loss multiplier must be finite and positive.')
    if args.output is None:
        suffix = '' if factor == 1 else f'_hydraulic_x{factor:g}'
        args.output = ROOT/f'outputs/motor_champion_sixbar_k2{suffix}.json'
    definition = CampaignDefinition(A5_CAMPAIGN)
    base, mass, saved, _, candidate_id = _candidate_design_and_mass(definition)
    k2 = json.loads(K2_REPORT.read_text())['best_feasible']
    limits = base.configuration.machine_volumes
    kinematics = IndependentSixBarVolumeKinematics(
        load_six_bar_mechanism(SMALL, restart=0, second_branch=1),
        load_six_bar_mechanism(LARGE), limits.small_cylinder, limits.large_cylinder,
    )
    inventory = _same_inventory_design(base, mass, kinematics)
    design = replace(inventory,
        heat_in=_scaled_exchanger(inventory.heat_in, float(k2['k_i'])),
        heat_out=_scaled_exchanger(inventory.heat_out, float(k2['k_o'])))
    if factor != 1:
        def scaled_losses(exchanger):
            # dp = a*m_dot + b*m_dot**2. Valve b is proportional to 1/CdA**2.
            return replace(exchanger, inputs=replace(exchanger.inputs,
                core_loss_multiplier=exchanger.inputs.core_loss_multiplier*factor,
                header_loss_coefficient=exchanger.inputs.header_loss_coefficient*factor),
                outlet_valve_cda_m2=exchanger.outlet_valve_cda_m2/math.sqrt(factor))
        design = replace(design, heat_in=scaled_losses(design.heat_in),
                         heat_out=scaled_losses(design.heat_out))
    initial = _rescaled_initial_state(saved, base, design)
    report = {
        'experiment': 'Independent six-bar S/L mechanisms with Stage K2 thermodynamics',
        'angle_convention': 'Study angle; original stored phases; build_model reverses motor direction once',
        'sources': {name: {'path': str(path.relative_to(ROOT)),
                          'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                          'contents': json.loads(path.read_text())}
                    for name, path in [('small', SMALL), ('large', LARGE)]},
        'selection': {'small': {'restart': 0, 'second_branch': 1}, 'large': 'global best'},
        'kinematics': asdict(kinematics),
        'configuration': asdict(design.configuration),
        'hydraulic_sensitivity': {
            'loss_multiplier': factor,
            'scope': 'Internal tube, header and valve losses; external-air losses excluded',
            'heat_in_inputs': asdict(design.heat_in.inputs),
            'heat_out_inputs': asdict(design.heat_out.inputs),
            'heat_in_valve_cda_m2': design.heat_in.outlet_valve_cda_m2,
            'heat_out_valve_cda_m2': design.heat_out.outlet_valve_cda_m2,
        },
        'wall_numerical_settings': asdict(definition.wall_numerical_settings),
        'thermodynamic_basis': {
            'stage_a5_candidate_id': candidate_id,
            'total_working_gas_mass_kg': mass,
            'stage_k2_k_i': k2['k_i'], 'stage_k2_k_o': k2['k_o'],
            'reference_efficiency': k2['indicated_thermal_efficiency'],
            'reference_indicated_power_w': k2['indicated_power_w'],
            'heat_in': _hardware_metrics(design.heat_in),
            'heat_out': _hardware_metrics(design.heat_out),
        },
        'initial_guess_state': initial.tolist(),
        'motion_fit_against_stage_f6_target': _motion_fit(kinematics),
        'efficiency_boundary': 'Indicated gas work / external-source heat input',
        'useful_mechanical_power_w': None,
        'mechanical_losses': 'unknown',
        'diagnostic_limitations': [
            'WallDiagnosticCycle does not reconstruct valve events; its generic '
            'non_nominal event-sequence label is not a measured valve chronology.',
            'C and H are legacy internal labels: C is motor heat-in Hi; H is heat-out Ho.',
        ],
    }
    started = time.monotonic()
    result, last_state = _evaluate('motor_champion_sixbar_k2', design, definition, initial_state=initial)
    report.update(result=result, elapsed_seconds=time.monotonic()-started,
                  last_complete_state=last_state.tolist() if last_state is not None else None)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(json_values(report), indent=2)+'\n')
    print(json.dumps(json_values({'result': result, 'elapsed_seconds': report['elapsed_seconds']}), indent=2))


if __name__ == '__main__':
    main()
