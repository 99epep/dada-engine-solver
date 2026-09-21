"""Evaluate the frozen 300 K motor champion with all passive-valve placements."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import time
import tomllib

import numpy as np

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.wall_backend import WallBackendSettings
from compare_motor_motion_laws_stage7A5 import A5_CAMPAIGN, ROOT, _evaluate
from optimize_motor_temperature_point_fixed import (
    COLD_K, MAX_FLOW, MAX_P, MAX_T, base_geometry, build_design, feasibility,
)
from refine_motor_four_stage_hx9d_variable_gas import _load_basis

SOURCE = ROOT/'outputs/motor_temperature_map/dT_300K/report.json'
SETTINGS = ROOT/'examples/motor_valve_placement_screen_300k.toml'
OUTPUT = ROOT/'outputs/motor_valve_placement_screen_300K.json'
HOT_K = 598.15
POWER_FLOOR_W = 40.0
PLACEMENTS = (
    ('DD', 'downstream', 'downstream'),
    ('UD', 'upstream', 'downstream'),
    ('DU', 'downstream', 'upstream'),
    ('UU', 'upstream', 'upstream'),
)


def extrema(result, group):
    values = result.get('diagnostics', {}).get(group, {}).values()
    return max((float(item['maximum']) for item in values), default=None)


def write_report(rows, source_sha):
    converged = [r for r in rows if r['converged']]
    feasible = [r for r in rows if r['feasible_under_300K_constraints']]
    payload = dict(
        study='frozen 300 K motor champion passive check-valve placement screen',
        source_report=str(SOURCE.relative_to(ROOT)), source_sha256=source_sha,
        cold_source_temperature_k=COLD_K, hot_source_temperature_k=HOT_K,
        nominal_flow='S -> H_i -> L -> H_o -> S',
        constraints=dict(minimum_motor_power_w=POWER_FLOOR_W,
            maximum_pressure_pa=MAX_P, maximum_temperature_k=MAX_T,
            maximum_absolute_mass_flow_kg_s=MAX_FLOW,
            valid_thermodynamic_exchanger_model=True),
        results=rows,
        best_converged_efficiency=(max(converged,key=lambda r:r['indicated_thermal_efficiency'])['topology_code'] if converged else None),
        best_feasible_efficiency=(max(feasible,key=lambda r:r['indicated_thermal_efficiency'])['topology_code'] if feasible else None),
        complete=len(rows)==4,
    )
    OUTPUT.write_text(json.dumps(payload,indent=2)+'\n')
    return payload


def main():
    backend_data=tomllib.loads(SETTINGS.read_text())['wall_backend']
    if backend_data.get('name')!='numba':
        raise ValueError('The valve-placement screen requires the Numba backend.')
    source=json.loads(SOURCE.read_text()); champion=source['champion']
    initial=np.asarray(champion['last_complete_state'],dtype=float)
    if initial.shape!=(10,): raise ValueError('The 300 K champion lacks a ten-state endpoint.')
    campaign=CampaignDefinition(A5_CAMPAIGN)
    _,base=_load_basis(); geometry=base_geometry(base)
    parameters=champion['parameters']
    backend=WallBackendSettings('numba')
    definition=SimpleNamespace(wall_numerical_settings=campaign.wall_numerical_settings,
        wall_backend=backend,exact_kinematics_cache=True,shared_replay=True)
    rows=[];source_sha=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    print('topology | converged | feasible | eta % | power W | cycles | elapsed s | backend')
    for code,heat_in,heat_out in PLACEMENTS:
        design=build_design(base,geometry,parameters,HOT_K)
        design=replace(design,configuration=replace(design.configuration,
            heat_in_valve_placement=heat_in,heat_out_valve_placement=heat_out))
        captures=[];started=time.perf_counter()
        result,state=_evaluate('valve_placement_'+code,design,definition,
            initial_state=initial.copy(),periodic_observer=captures.append)
        elapsed=time.perf_counter()-started
        periodic=captures[0]; stats=periodic.backend_statistics
        actual=stats.get('actual_backend')
        if actual!='numba' or stats.get('fallback_calls',0):
            raise RuntimeError(f'{code} did not run entirely on Numba: {stats}')
        converged=result.get('status')=='converged'
        feasible,reasons=feasibility(result,POWER_FLOOR_W) if converged else (False,[])
        eta=result.get('indicated_thermal_efficiency')
        row=dict(topology_code=code,heat_in_valve_placement=heat_in,
            heat_out_valve_placement=heat_out,status=result.get('status'),
            converged=converged,periodic_cycle_count=len(periodic.history),
            elapsed_seconds=elapsed,indicated_thermal_efficiency=eta,
            indicated_power_w=result.get('indicated_power_w'),
            heat_input_w=result.get('heat_input_w'),heat_out_w=result.get('heat_out_w'),
            fraction_of_carnot=(eta/(1-COLD_K/HOT_K) if eta is not None else None),
            maximum_pressure_pa=result.get('maximum_pressure_pa') or extrema(result,'pressure_extrema'),
            maximum_temperature_k=result.get('maximum_temperature_k') or extrema(result,'temperature_extrema'),
            maximum_absolute_mass_flow_kg_s=result.get('maximum_absolute_mass_flow_kg_s'),
            validity_verdict=result.get('validity',{}).get('verdict'),
            failed_validity_criteria=result.get('validity',{}).get('failed_criteria',[]),
            maximum_tube_reynolds=result.get('maximum_tube_reynolds'),
            maximum_tube_mach=result.get('maximum_tube_mach_number'),
            final_periodic_state=state.tolist() if state is not None else None,
            requested_rhs_backend=stats.get('requested_backend'),actual_rhs_backend=actual,
            fallback_count=stats.get('fallback_calls',0),fallback_reason=stats.get('fallback_reason'),
            feasible_under_300K_constraints=feasible,
            failed_300K_constraints=reasons,
            backend_statistics=stats)
        if code == 'DD':
            fields=('indicated_thermal_efficiency','indicated_power_w','heat_input_w',
                    'heat_out_w','maximum_pressure_pa','maximum_temperature_k',
                    'maximum_absolute_mass_flow_kg_s')
            reference=champion['result']
            comparisons={name:dict(reference=reference[name],replayed=row[name],
                difference=row[name]-reference[name]) for name in fields}
            row['historical_DD_comparison']=comparisons
            row['historical_DD_regression_passed']=all(np.isclose(
                item['replayed'],item['reference'],rtol=5e-5,atol=1e-9)
                for item in comparisons.values())
        rows.append(row);write_report(rows,source_sha)
        if code == 'DD' and not row['historical_DD_regression_passed']:
            raise RuntimeError('DD replay did not reproduce the historical 300 K result.')
        print(f"{code:>2} | {str(converged):>9} | {str(feasible):>8} | "
              f"{100*eta if eta is not None else float('nan'):7.3f} | "
              f"{row['indicated_power_w'] if row['indicated_power_w'] is not None else float('nan'):8.3f} | "
              f"{len(periodic.history):6d} | {elapsed:9.3f} | {actual}",flush=True)
    report=write_report(rows,source_sha)
    print('best converged efficiency:',report['best_converged_efficiency'])
    print('best feasible efficiency:',report['best_feasible_efficiency'])


if __name__=='__main__': main()
