"""Model domains are distinct from study requirements and search bounds."""
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace as NS
import tomllib

import pytest

from dada_solver.configuration import ValidityThresholds
from dada_solver.campaign.evaluator import finalize_microtube_validity
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel, MicrotubeDomainError
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.research.cli import initialize
from dada_solver.research.presets import initialize_v2, initialize_v3
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.schema_v2 import V2Adapter
from dada_solver.research.study_io import dumps
from dada_solver.sizing.configuration import load_sizing_problem
from dada_solver.sizing.constraints import MaximumMachNumber, RequireValidThermodynamicModel
from dada_solver.validity import assess_cycle_validity, ValidityVerdict
from tests.test_periodic import create_static_cycle


@pytest.mark.parametrize('initializer',[initialize_v2,initialize_v3])
def test_generic_presets_do_not_inherit_motor_design_guards(tmp_path,initializer):
    study=load_study(initializer(tmp_path/'study.toml'))
    forbidden={'minimum_motor_power','maximum_pressure','maximum_temperature','maximum_absolute_mass_flow'}
    assert not forbidden.intersection(row['type'] for row in study.data['constraints'])
    assert study.data['screening']=={'samples':1440}
    assert set(study.basis.data['configuration']['validity'])=={'maximum_cp_variation','maximum_compressibility_deviation'}


def test_explicit_motor_constraints_are_available(tmp_path):
    path=initialize_v2(tmp_path/'study.toml');raw=tomllib.loads(path.read_text())
    raw['constraints'] += [dict(type='minimum_motor_power',required_power=25.,unit='W'),
        dict(type='maximum_pressure',limit=1.2e6,unit='Pa'),dict(type='maximum_temperature',limit=850.,unit='K'),
        dict(type='maximum_absolute_mass_flow',limit=.08,unit='kg/s'),dict(type='maximum_mach_number',limit=.2,unit='1')]
    path.write_text(dumps(raw));definition=compile_study(load_study(path))
    evaluation=NS(performance=NS(motor_power=24.),validity=NS(maximum_mach_number=.24),
        diagnostics=NS(pressure_extrema={'S':NS(maximum=1.3e6)},temperature_extrema={'S':NS(maximum=851.)},
            mass_flow_extrema={'port':NS(minimum=-.09,maximum=.01)}))
    for constraint in definition.constraints:
        if constraint.name=='valid_thermodynamic_model': continue
        value=constraint.evaluate(evaluation)
        assert value.available and not value.satisfied and value.margin<0


def test_optional_motor_volume_ceiling_may_exceed_66_litres(tmp_path):
    path=initialize_v2(tmp_path/'study.toml');raw=tomllib.loads(path.read_text())
    raw['screening']['maximum_large_enclosed_volume_m3']=.2
    path.write_text(dumps(raw));study=load_study(path)
    V2Adapter(study).build(study.fixed_parameters)
    raw['screening'].pop('maximum_large_enclosed_volume_m3');path.write_text(dumps(raw))
    assert load_study(path).study_id!=study.study_id


def test_legacy_validity_keys_are_ignored_not_serialized_and_pressure_is_observable(ideal_gas):
    thresholds=ValidityThresholds(maximum_pressure_equalization_error=.05,
        maximum_mach_number=.2, maximum_compressibility_deviation=.01,maximum_cp_variation=.01)
    assert set(asdict(thresholds))=={'maximum_cp_variation','maximum_compressibility_deviation'}
    model,_,cycle=create_static_cycle(ideal_gas)
    states=cycle.states.copy();states[7,-1]*=1.2
    report=assess_cycle_validity(replace(cycle,states=states),model,thresholds)
    assert report.maximum_pressure_equalization_error>.05
    assert not report.failed_criteria
    final=finalize_microtube_validity(report,500.,.24,requires_laminar=False)
    assert RequireValidThermodynamicModel().evaluate(NS(validity=final)).satisfied
    # Obsolete values cannot reject even if a historical file used tiny limits.
    limits=ValidityThresholds(1e-12,1e-12,.01,.01)
    assert assess_cycle_validity(replace(cycle,states=states),model,limits).failed_criteria==()


@pytest.mark.parametrize('mach,valid',[(.24,True),(.301,False)])
def test_microtube_own_mach_domain_and_independent_design_limit(ideal_gas,mach,valid):
    model=MicrotubeGasModel(maximum_mach=.3)
    bank=MicrotubeBank(1000,.8,.0003,.0001,.0006,.001)
    # Obtain the exact same local Mach conversion as the model, without changing
    # any pressure/transport/domain criterion to manufacture a passing result.
    unit=model.diagnose(bank,.001,2e5,1.99e5,300.)
    diagnostic=model.diagnose(bank,.001*mach/unit.mach,2e5,1.99e5,300.)
    assert diagnostic.mach==pytest.approx(mach)
    assert ('high_mach' not in diagnostic.issues)==valid
    if valid: model.require(diagnostic)
    else:
        with pytest.raises(MicrotubeDomainError,match='high_mach'): model.require(diagnostic)
    thermo,_,cycle=create_static_cycle(ideal_gas)
    generic=assess_cycle_validity(cycle,thermo,ValidityThresholds(.05,.2,.01,.01))
    final=finalize_microtube_validity(generic,diagnostic.reynolds,diagnostic.mach,
        requires_laminar=False,domain_failures=diagnostic.issues)
    assert (final.verdict is ValidityVerdict.VALID)==valid
    assert not MaximumMachNumber(.2).evaluate(NS(validity=final)).satisfied


@pytest.mark.parametrize('family',['four_bar','six_bar'])
def test_mechanical_quality_limits_belong_to_new_study_not_seed(tmp_path,family):
    path=initialize_v2(tmp_path/'study.toml',family,family);study=load_study(path)
    assert study.mechanical_constraints
    artifact_paths=list(tmp_path.glob('*.mechanism.json'))
    before={p:p.read_bytes() for p in artifact_paths}
    assert all(not artifact.scientific['constraints'] for artifact in study.artifacts.values())
    raw=study.data;raw['mechanical_constraints']=[];path.write_text(dumps(raw))
    relaxed=load_study(path)
    assert not relaxed.mechanical_constraints and relaxed.study_id!=study.study_id
    assert all(p.read_bytes()==before[p] for p in artifact_paths)
    assert all(relaxed.artifacts[s].content_hash==study.artifacts[s].content_hash for s in ('small','large'))


def test_historical_v1_preset_is_explicitly_preserved(tmp_path):
    study=load_study(initialize(tmp_path/'historical.toml'))
    assert study.data['study']['purpose']=='historical_physical_parity'
    assert next(c for c in study.data['constraints'] if c['type']=='maximum_absolute_mass_flow')['limit']==.08


def test_obsolete_sizing_constraints_and_scales_load_but_are_inactive(tmp_path):
    source=Path('examples/sizing_controlled_example.toml')
    raw=tomllib.loads(source.read_text())
    for key in ('base_configuration','cooling_load_configuration'):
        raw['problem'][key]=str((source.parent/raw['problem'][key]).resolve())
    path=tmp_path/'sizing.toml';path.write_text(dumps(raw))
    loaded=load_sizing_problem(path)
    assert all(c.name != "maximum_pressure_equalization_error" for c in loaded.problem.constraints)
    assert 'maximum_pressure_equalization_error' not in loaded.optimization_settings.constraint_scales


def test_historical_embedded_mechanical_constraints_still_apply(tmp_path):
    from dada_solver.research.artifacts import MechanismArtifact
    path=initialize_v2(tmp_path/'study.toml','six_bar','six_bar');study=load_study(path)
    raw=study.data
    for side in ('small','large'):
        old=study.artifacts[side]
        constraints=[{k:v for k,v in row.items() if k!='side'} for row in raw['mechanical_constraints'] if row['side']==side]
        artifact=MechanismArtifact.create('six_bar',old.scientific['geometry'],
            settings=old.scientific['settings'],constraints=constraints,provenance={'description':'historical embedded constraint fixture'})
        artifact_path=tmp_path/(side+'.historical.json')
        raw['kinematics'][side]['artifact']=artifact_path.name
        artifact.save(artifact_path)
        raw['kinematics'][side]['sha256']=artifact.content_hash
    raw['mechanical_constraints']=[];path.write_text(dumps(raw))
    historical=load_study(path)
    assert len(historical.mechanical_constraints)==20
    assert all(a.scientific['constraints'] for a in historical.artifacts.values())
    assert historical.study_id==load_study(path).study_id


def test_caloric_and_eos_approximation_checks_remain_active(ideal_gas,monkeypatch):
    from dada_solver.state import ThermodynamicState
    model,_,cycle=create_static_cycle(ideal_gas)
    pressures=[ThermodynamicState.from_array(values).pressures(ideal_gas,model.volumes(float(angle)))
               for angle,values in zip(cycle.angles,cycle.states.T)]
    import numpy as np
    replay=NS(require=lambda **kwargs:None,pressures=np.asarray(pressures).T)
    monkeypatch.setattr(ThermodynamicState,'fluid_states',lambda *args:[
        NS(compressibility_factor=1.1,cp=1000.),NS(compressibility_factor=1.,cp=1100.)])
    report=assess_cycle_validity(cycle,NS(gas=NS(),volumes=model.volumes),
        ValidityThresholds(maximum_compressibility_deviation=.01,maximum_cp_variation=.01),replay=replay)
    assert report.verdict is ValidityVerdict.INVALID
    assert set(report.failed_criteria)=={'compressibility_factor','heat_capacity_variation'}


def test_kinematic_metrics_survive_without_any_mechanical_limit(tmp_path,monkeypatch):
    from dada_solver.campaign.evaluator import MachineEvaluator,EvaluationControl,rejected
    from dada_solver.research.schema import candidate_for_values
    study=load_study(initialize_v2(tmp_path/'study.toml','slider_crank','slider_crank'))
    assert not study.mechanical_constraints
    definition=compile_study(study)
    def stop_before_integration(self,candidate,design,wrapper,derived,direction,control):
        return dict(rejected('integration_failure','test stops before integration'),derived=derived)
    monkeypatch.setattr(MachineEvaluator,'_wall',stop_before_integration)
    result=MachineEvaluator(definition)._evaluate_with_control(candidate_for_values(definition,{}),EvaluationControl())
    assert len(result['derived']['kinematic_metrics'])==2
    assert not result['derived'].get('mechanical_constraints')
