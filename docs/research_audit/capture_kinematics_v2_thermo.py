"""Freeze original machine inputs/warm starts and stored metrics without solving.

Run from the checkout with PYTHONPATH=src:examples. Existing motion references
were captured separately before extraction. These portable thermal fixtures
remove historical examples and large histories from the regression runtime.
"""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'tests/fixtures/research_v2/thermodynamic.json'


def main():
    if DEST.exists(): raise RuntimeError('Refusing to replace a reviewed reference.')
    from dada_solver.campaign.definition import CampaignDefinition
    from compare_motor_motion_laws_stage7A5 import A5_CAMPAIGN,_candidate_design_and_mass,_same_inventory_design
    from optimize_motor_exchanger_asymmetry_stageK1 import _scaled_exchanger
    from evaluate_motor_champion_four_bar_k2 import _rescaled_initial_state
    from compact_coupler_geometry import load_compact_coupler
    from evaluate_slider_crank_k2 import SliderMotion,SliderCrankKinematics
    from refine_motor_fourier_c2_260k import _load_basis,base_geometry,build_design
    from optimize_motor_temperature_point import uniform_metadata,warm_start
    from optimize_motor_free_spline_260k_v3 import _fixed_inventory_design,_make_kinematics
    sources={}
    def read(relative):
        p=ROOT/relative;sources[relative]=hashlib.sha256(p.read_bytes()).hexdigest()
        return json.loads(p.read_text())
    def predecessor(folder,cid):
        p=ROOT/'outputs'/folder/'history.jsonl'
        sources[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
        for line in p.open():
            r=json.loads(line)
            if r['candidate_id']==cid:return r
        raise ValueError('Missing historical predecessor')
    cases={}
    def save(family,design,definition,initial,stored,backend):
        wrapper=design.build()
        metrics={k:stored[k] for k in ('indicated_power_w','indicated_thermal_efficiency','heat_input_w','maximum_absolute_mass_flow_kg_s')}
        for metric,field in [('maximum_pressure_pa','pressure_extrema'),('maximum_temperature_k','temperature_extrema')]:
            metrics[metric]=max(v['maximum'] for v in stored['diagnostics'][field].values())
        data=dict(schema_version=2,configuration=asdict(design.configuration),heat_in=asdict(design.heat_in),heat_out=asdict(design.heat_out),geometry={},
            wall_settings=asdict(definition.wall_numerical_settings),warm_start=dict(values=np.asarray(initial).tolist(),source_candidate_id='historical_reference_initial_state',wall_capacities_j_k=[wrapper.heat_in.wall_capacity_j_k,wrapper.heat_out.wall_capacity_j_k]),
            provenance=dict(purpose='Historical thermal parity with original warm start',sources=dict(sources)))
        cases[family]=dict(basis=data,metrics=metrics,cycles=stored['convergence']['cycles_completed'],backend=backend)
    cd=CampaignDefinition(A5_CAMPAIGN)
    base,mass,saved,*_=_candidate_design_and_mass(cd)
    k2=read('outputs/motor_exchanger_asymmetry_stageK2.json')['best_feasible']
    hardware=replace(base,heat_in=_scaled_exchanger(base.heat_in,k2['k_i']),heat_out=_scaled_exchanger(base.heat_out,k2['k_o']))
    raw=read('outputs/slider_crank_target_comparison.json')
    kin=SliderCrankKinematics(*(SliderMotion(**{k:r[k] for k in SliderMotion.__dataclass_fields__}) for r in [raw['sides'][side]['inverted_offset'] for side in ('small','large')]),base.configuration.machine_volumes.small_cylinder,base.configuration.machine_volumes.large_cylinder)
    design=_same_inventory_design(hardware,mass,kin)
    save('slider_crank',design,cd,_rescaled_initial_state(saved,base,design),read('outputs/slider_crank_k2.json')['results']['inverted_offset'],'python')
    read('outputs/compact_motor_coupler_search.json')
    kin=load_compact_coupler(ROOT/'outputs/compact_motor_coupler_search.json',base.configuration.machine_volumes.small_cylinder,base.configuration.machine_volumes.large_cylinder)
    design=_same_inventory_design(hardware,mass,kin)
    save('four_bar',design,cd,_rescaled_initial_state(saved,base,design),read('outputs/compact_motor_coupler_k2.json')['result'],'python')
    cd,base=_load_basis();g=base_geometry(base)
    fb=read('outputs/motor_fourier_c2_8h_refine/report.json')['best_feasible']
    design=build_design(base,g,fb['thermo'],fb['small_coefficients'],fb['large_coefficients'],fb['harmonics'])
    uniform,hw,_=uniform_metadata(design)
    warm=predecessor('motor_fourier_c2_8h_refine',fb['warm_start_source'])
    initial=warm_start(warm['last_complete_state'],warm['hardware'],uniform,hw)
    # The pressure-charged historical case becomes the exact resolved inventory.
    design=_fixed_inventory_design(design,float(uniform[:8:2].sum()))
    save('fourier_c2',design,cd,initial,fb['result'],'numba')
    sb=read('outputs/motor_spline_from_linear_260k/report.json')['best_feasible']
    design=_fixed_inventory_design(design,float(fb['result']['total_mass_kg']))
    design=replace(design,kinematics=_make_kinematics(design.configuration.machine_volumes,sb['small_controls'],sb['large_controls'],sb['small_phase_deg'],sb['large_phase_deg']))
    assert sb['warm_start_mode']=='nearest_periodic_state' and not sb['safe_retry_used']
    initial=predecessor('motor_spline_from_linear_260k',sb['warm_start_source'])['last_complete_state']
    save('free_spline',design,cd,initial,sb['result'],'numba')
    DEST.write_text(json.dumps(cases,indent=2)+'\n')
    print('Captured four portable historical thermal cases without integration.')

if __name__=='__main__':main()
