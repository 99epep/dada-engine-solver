"""Freeze pre-extraction trajectories, source identities and historical inputs.

Run once from the checkout with PYTHONPATH=src:examples. No optimization or
integration is performed, and existing references are never overwritten.
"""
from dataclasses import asdict, fields
import hashlib
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'tests/fixtures/research_v2'


def main():
    if (DEST/'reference.json').exists():
        raise RuntimeError('Reference already exists; capture into a separate reviewed location.')
    from dada_solver.research.sixbar import load_basis
    from dada_solver.six_bar import SixBarCylinderMechanism, IndependentSixBarVolumeKinematics
    from dada_solver.kinematics import HarmonicVolumeKinematics, IdealPiecewiseLinearVolumeKinematics
    from evaluate_slider_crank_k2 import SliderMotion, SliderCrankKinematics
    from compact_coupler_geometry import load_compact_coupler
    from optimize_motor_fourier_c2_260k import FourierVolumeKinematics
    from optimize_motor_hybrid_c2_15p_260k import StructuredKinematics15
    from optimize_motor_free_spline_260k_v3 import _make_kinematics
    sources = {}
    def read(path):
        path = ROOT/path
        sources[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        return json.loads(path.read_text())
    bpath = ROOT/'src/dada_solver/research/data/rank01_basis.json'
    basis = load_basis(bpath, hashlib.sha256(bpath.read_bytes()).hexdigest())
    limits = basis.configuration.machine_volumes
    s,l = limits.small_cylinder, limits.large_cylinder
    slider = read('outputs/slider_crank_target_comparison.json')
    fourbar_path = 'outputs/compact_motor_coupler_search.json'
    read(fourbar_path)
    fourbar = load_compact_coupler(ROOT/fourbar_path,s,l)
    fourier = read('outputs/motor_fourier_c2_8h_refine/report.json')['best_feasible']
    structured = read('outputs/motor_hybrid_c2_15p_260k/candidate_3952_motion_target.json')
    spline = read('outputs/motor_spline_from_linear_260k/report.json')['best_feasible']
    inputs = dict(limits=asdict(limits), slider={side:{k:v for k,v in slider['sides'][side]['inverted_offset'].items()
        if k in SliderMotion.__dataclass_fields__} for side in ('small','large')},
        fourbar={f.name:asdict(getattr(fourbar,f.name)) if f.name.endswith('_assembly') else getattr(fourbar,f.name)
            for f in fields(fourbar) if f.init and 'volume_limits' not in f.name},
        fourier={k:fourier[k] for k in ('harmonics','small_coefficients','large_coefficients')},
        structured=structured, spline={k:spline[k] for k in ('small_controls','large_controls','small_phase_deg','large_phase_deg')},
        sixbar_families={})
    motions = dict(harmonic=HarmonicVolumeKinematics(s,l,math.radians(249)),
        slider_crank=SliderCrankKinematics(*(SliderMotion(**inputs['slider'][side]) for side in ('small','large')),s,l),
        four_bar=fourbar,
        fourier_c2=FourierVolumeKinematics(s,l,**inputs['fourier']),
        structured_c2=StructuredKinematics15(limits,structured['parameters']),
        free_spline=_make_kinematics(limits,**dict(small=inputs['spline']['small_controls'],large=inputs['spline']['large_controls'],
            small_phase_deg=spline['small_phase_deg'],large_phase_deg=spline['large_phase_deg'])),
        ideal_piecewise=IdealPiecewiseLinearVolumeKinematics(s,l,.12,.14,.21))
    for rank in (1,4,12,50):
        raw=read(f'outputs/sixbar_thermo_coupled_3952_fine_hlat25/rank_{rank:02d}/best_pair.json')
        pair={}
        for side in ('small','large'):
            p={k.lower():v for k,v in raw[side]['parameters'].items()}
            p['primary_branch']=raw[side]['primary']['assembly_branch']
            p['second_branch']=raw[side]['second_branch']
            pair[side]={f.name:p[f.name] for f in fields(SixBarCylinderMechanism) if f.init}
        inputs['sixbar_families'][str(rank)]=dict(geometry=pair, historical_metrics={side:{k:v for k,v in raw[side].items()
            if isinstance(v,(float,int))} for side in ('small','large')})
        motions[f'six_bar_{rank}']=IndependentSixBarVolumeKinematics(SixBarCylinderMechanism(**pair['small']),SixBarCylinderMechanism(**pair['large']),s,l)
    theta=np.linspace(-2*math.pi,4*math.pi,4321)
    arrays={'theta':theta}
    for name,motion in motions.items():
        arrays[name]=np.array([motion.cylinder_volumes_and_derivatives(float(a)) for a in theta])
        if name in ('structured_c2','free_spline'):
            if name=='structured_c2':
                acc=motion.motion.normalized((-theta/(2*math.pi))%1)[2]/(2*math.pi)**2
                arrays[name+'_second']=acc.T*np.array([s.swept,l.swept])
            else:
                arrays[name+'_second']=np.array([motion.small_cylinder_volume_second_derivative(theta),motion.large_cylinder_volume_second_derivative(theta)]).T
    for path in ['examples/evaluate_slider_crank_k2.py','examples/compare_slider_crank_motion.py',
                 'examples/compact_coupler_geometry.py','examples/optimize_motor_fourier_c2_260k.py',
                 'examples/fit_motor_hybrid_c2_15p_to_source16.py','examples/optimize_motor_hybrid_c2_15p_260k.py',
                 'examples/optimize_motor_free_spline_260k_v3.py']:
        sources[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    DEST.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(DEST/'trajectories.npz',**arrays)
    (DEST/'reference.json').write_text(json.dumps(dict(schema_version=1,inputs=inputs,sources=sources,
        purpose='Pre-extraction historical dense-grid parity; no search-order guarantee'),indent=2)+'\n')
    print(f'Captured {len(motions)} motions and {len(theta)} angles per motion.')


if __name__=='__main__': main()
