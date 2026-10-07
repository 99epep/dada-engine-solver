"""Structured geometric refit, current spline regression and exact study handoff."""
from dataclasses import asdict
from pathlib import Path
import json
import math
import tomllib
from types import SimpleNamespace
import numpy as np
import pytest
from dada_solver.free_kinematics import FreeMotionDefinition
from dada_solver.hybrid_compact_kinematics import HybridCompactKinematics
from dada_solver.geometry import CylinderVolumeLimits
from dada_solver.research.motion_target import MotionTarget, _hybrid_events
from dada_solver.research.motion_refit import MotionRefitRequest, MotionRefitResult, RefitPolicy, structured_model, monotonicity, plot_refit
from dada_solver.research.refit_study import refit_study
from dada_solver.research.families import STRUCTURED_DEFAULTS
from dada_solver.research.presets import initialize_kinematics
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.study_io import dumps

POLICY=RefitPolicy(dense_samples=360,maximum_evaluations=160,curvature_scales=(1.,),kink_width_starts=(.25,.65))

@pytest.mark.parametrize('count', [4,8,15,24])
def test_shape_chart_inverse_and_positive_affine_invariance(count):
    rng = np.random.default_rng(251)
    for _ in range(12):
        controls = rng.normal(size=count)
        definition = FreeMotionDefinition(tuple(controls),1.,2.)
        z = definition.to_shape_coordinates()
        rebuilt = FreeMotionDefinition.from_shape_coordinates(z,1.,2.)
        assert len(z) == count-2
        np.testing.assert_allclose(rebuilt.control_values,definition.control_values,atol=2e-14,rtol=2e-14)
        transformed = FreeMotionDefinition(tuple(7.+3.*controls),1.,2.)
        np.testing.assert_allclose(transformed.to_shape_coordinates(),z,atol=3e-14,rtol=3e-14)
    pole = np.r_[np.full(count-1,1/math.sqrt((count-1)*count)), -(count-1)/math.sqrt((count-1)*count)]
    with pytest.raises(ValueError,match='pole'):
        FreeMotionDefinition(tuple(pole),1.,2.).to_shape_coordinates()



def target_of(model,events=None):
    return MotionTarget.from_kinematics(model,source={'test':'production study angle'},samples=2881,events=events)


@pytest.fixture(scope='module')
def structured_fit():
    p=dict(STRUCTURED_DEFAULTS,small_max_deg=147.,small_down_duration_deg=169.,large_down_duration_deg=196.,
           small_down_kink_u=.45,small_down_kink_q=.43,large_up_kink_u=.57,large_up_kink_q=.55,
           small_max_curvature=25.,small_min_curvature=21.,large_max_curvature=23.,large_min_curvature=19.)
    model=structured_model(p)
    assert monotonicity(model.motion)['valid']
    target=target_of(model)
    return model,target,MotionRefitRequest(target,policy=POLICY).execute()


def test_exact_structured_refit_and_angle_convention(structured_fit):
    model,target,result=structured_fit
    raw=result.data['scientific'];rebuilt=structured_model(raw['parameters'])
    assert raw['parameter_count']==15 and raw['destination_family']=='structured_c2_15p'
    for side in ('small','large'):
        d=raw['sides'][side]['diagnostics']
        assert d['position_rms']<2e-6
        assert d['extrema_count']==2 and raw['monotonicity']['valid']
        angles=np.linspace(0,2*math.pi,1501)
        np.testing.assert_allclose(getattr(rebuilt,side+'_cylinder_volume')(angles),getattr(model,side+'_cylinder_volume')(angles),atol=2e-5)
        np.testing.assert_allclose(getattr(rebuilt,side+'_cylinder_volume')(angles+2*math.pi),getattr(rebuilt,side+'_cylinder_volume')(angles),atol=1e-14)
    assert rebuilt.large_cylinder_volume(0.)==2.
    assert rebuilt.small_cylinder_volume(-math.radians(raw['parameters']['small_max_deg']))==pytest.approx(2.)
    assert rebuilt.large_cylinder_volume_derivative(.2)<0
    assert MotionRefitResult.from_data(result.data).data==result.data
    print('Structured RMS',raw['combined']['position_rms'])


@pytest.mark.parametrize('variant',[0,1,2])
def test_hybrid_fit_is_monotone_and_parameters_are_free(variant):
    model=HybridCompactKinematics(CylinderVolumeLimits(1.,2.),CylinderVolumeLimits(1.,2.),
        small_max_deg=140.+5*variant,small_down_duration_deg=175.-5*variant,large_down_duration_deg=190.+5*variant,
        large_down_rounding=.12+.02*variant,small_up_rounding=.14+.02*variant,
        small_down_kink_u=.47,small_down_kink_q=.48,large_up_kink_u=.53,large_up_kink_q=.52)
    target=target_of(model,{s:_hybrid_events(model,s) for s in ('small','large')})
    result=MotionRefitRequest(target,policy=POLICY).execute();raw=result.data['scientific']
    assert raw['monotonicity']['valid'] and raw['combined']['position_rms']<.02
    assert len(raw['parameters'])==15 and len(raw['starts'])==2
    assert abs(raw['parameters']['small_down_duration_deg']-raw['initial_parameters']['small_down_duration_deg'])>1e-3
    for side in ('small','large'):
        d=raw['sides'][side]['diagnostics']
        assert d['extrema_count']==2 and not d['target_acceleration_used']
        assert d['events'] and d['velocity_rms'] is not None
    print('Hybrid',variant,{s:raw['sides'][s]['diagnostics']['position_rms'] for s in ('small','large')})
    if variant==0:
        assert MotionRefitRequest(target,policy=POLICY).execute().payload_json==result.payload_json


def test_acceleration_is_not_used_and_velocity_can_be_absent(structured_fit):
    _,target,result=structured_fit
    raw=target.scientific
    for side in ('small','large'):raw['sides'][side]['second_derivative']=[1e6]*len(raw['angles_rad'])
    changed=MotionTarget.create(raw['angles_rad'],raw['sides'],source=raw['source'])
    assert MotionRefitRequest(changed,policy=POLICY).execute().data['scientific']['parameters']==result.data['scientific']['parameters']
    for side in ('small','large'):
        raw['sides'][side]['first_derivative']=None;raw['sides'][side]['second_derivative']=None
    changed=MotionTarget.create(raw['angles_rad'],raw['sides'],source=raw['source'])
    other=MotionRefitRequest(changed,policy=POLICY).execute().data['scientific']
    assert all(other['sides'][s]['diagnostics']['velocity_rms'] is None for s in ('small','large'))


def no_thermodynamics(monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('Refit must not integrate'))


def test_candidate_generates_exact_machine_with_fifteen_active_coordinates(tmp_path,monkeypatch):
    from tests.test_research_rescale import snapshot
    no_thermodynamics(monkeypatch)
    source=initialize_kinematics(tmp_path/'source.toml','hybrid_compact','hybrid_compact')
    raw=tomllib.loads(source.read_text());row=next(r for r in raw['parameters'] if r['name']=='operation.frequency_hz')
    v=row.pop('value');row.update(initial=v,lower=v/2,upper=v*2,kind='continuous',transform='linear');source.write_text(dumps(raw))
    study,definition,record=snapshot(source,tmp_path/'campaign',{'operation.frequency_hz':v*1.25})
    output,result=refit_study(tmp_path/'campaign',tmp_path/'next/study.toml',candidate=record['candidate_id'],policy=POLICY,radius=.05)
    new=load_study(output);assert len(new.space.parameters)==15
    assert {p.name for p in new.space.parameters}=={f'kinematics.{"small" if n.startswith("small_") else "large"}.{n}' for n in STRUCTURED_DEFAULTS}
    for side in ('small','large'):assert new.settings[side]=={'family':'structured_c2_15p'}
    original=definition.adapter.build(dict(study.fixed_parameters,**record['physical']))
    rebuilt=compile_study(new).adapter.build(dict(new.fixed_parameters,**{p.name:p.initial for p in new.space.parameters}))
    assert asdict(original.configuration)==asdict(rebuilt.configuration)
    assert asdict(original.heat_in)==asdict(rebuilt.heat_in)
    assert asdict(original.heat_out)==asdict(rebuilt.heat_out)
    for key in ('constraints','objective','numerical','mechanical_constraints'):assert new.data[key]==study.data[key]
    assert new.basis.data['provenance']['motion_refit']['radius_fraction']==.05
    from dada_solver.research.cli import main
    assert main(['validate',str(output)])==0


def test_cli_report_plot_and_removed_spline_options(tmp_path,monkeypatch,structured_fit):
    no_thermodynamics(monkeypatch)
    _,target,result=structured_fit
    monkeypatch.setattr(MotionRefitRequest,'execute',lambda self:result)
    source=tmp_path/'target.json';target.save(source)
    from dada_solver.research.cli import main
    assert main(['motion-refit',str(source),'--validate-only'])==0
    assert main(['motion-refit',str(source),'--report',str(tmp_path/'refit.json'),'--plot','--no-show'])==0
    for option in ('--count','--family','--shape-radius','--phase-radius-fraction'):
        with pytest.raises(SystemExit):main(['motion-refit',str(source),option,'15'])
    with pytest.raises(SystemExit):main(['motion-refit',str(source),'--output',str(tmp_path/'no-machine.toml')])
    import matplotlib.pyplot as plt
    fig=plot_refit(result);assert len(fig.axes)==6
    assert 'structured_c2_15p refit' in fig._suptitle.get_text()
    assert not any(c.get_label()=='15 spline nodes' for ax in fig.axes for c in ax.collections)
    fig.canvas.draw();plt.close(fig)


def test_reference_inventory_stays_fixed(tmp_path,monkeypatch):
    from tests.test_research_reference_charge import derived
    no_thermodynamics(monkeypatch)
    source=initialize_kinematics(tmp_path/'source.toml','hybrid_compact','hybrid_compact')
    raw=tomllib.loads(source.read_text());derived(raw);source.write_text(dumps(raw))
    original=load_study(source);mass=compile_study(original).adapter.build(original.fixed_parameters).configuration.charge.total_mass
    output,_=refit_study(source,tmp_path/'structured.toml',policy=POLICY)
    new=load_study(output)
    assert new.data['policies']['charge']=='explicit_inventory'
    assert new.fixed_parameters['charge.total_mass_kg']==mass


def test_nonmonotone_structured_is_rejected():
    assert not monotonicity(structured_model(dict(STRUCTURED_DEFAULTS,small_max_curvature=2000.)).motion)['valid']


def test_free_spline_family_remains_loadable_and_evaluable(tmp_path,monkeypatch):
    import dada_solver.campaign.evaluator as evaluator
    from dada_solver.integration import IntegrationInterrupted
    from dada_solver.research.schema import candidate_for_values
    source=initialize_kinematics(tmp_path/'spline.toml','free_spline','free_spline')
    study=load_study(source);definition=compile_study(study)
    definition.adapter.build(study.fixed_parameters).build()
    calls=[]
    def bounded_solver(*args,**kwargs):
        calls.append(args)
        raise IntegrationInterrupted('Bounded regression check')
    monkeypatch.setattr(evaluator,'solve_periodic_wall_machine',bounded_solver)
    result=evaluator.MachineEvaluator(definition).evaluate(candidate_for_values(definition,{}))
    assert len(calls)==1 and result['status']=='budget_exhausted'
