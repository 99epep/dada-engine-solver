"""Method 6.4: primary topology/cadence, never piston-position matching."""
from dataclasses import replace
import math
import numpy as np
import pytest
from dada_solver.six_bar import velocity_stationary_points
from dada_solver.research.motion_target import MotionTarget, MotionEvent, PeriodicTargetSide, PERIOD
from dada_solver.research.synthesis_six_bar import (
    PrimaryCadencePolicy, PrimaryProjection, primary_target_features,
    primary_cadence_score, primary_evidence, PRIMARY_POLICY, SixBarStageSearch)
from dada_solver.research.synthesis_search import SearchPolicy
from tests.test_research_six_bar_synthesis import geometry, plan, policy, primary_library, target_for


class HarmonicProjection:
    def __init__(self,third=0.,sign=1.,phase=0.,second_sine=0.):
        self.third,self.sign,self.phase,self.second_sine=third,sign,phase,second_sine
        self.axis=np.array([sign,0.])
        self.roots=velocity_stationary_points(lambda t:self.value(t)-1,lambda t:self.value(t,1),360)
        self.metrics=dict(minimum_primary_transmission_sine=.8,E_axis_variance_fraction=.9)
    def value(self,t,order=0):
        t=np.asarray(t)-self.phase
        if order==0: return 1.5+self.sign*.5*(np.cos(t)+self.third*np.cos(3*t)+self.second_sine*np.sin(2*t))/(1+self.third+abs(self.second_sine))
        return self.sign*.5*(-np.sin(t)-3*self.third*np.sin(3*t)+2*self.second_sine*np.cos(2*t))/(1+self.third+abs(self.second_sine))


def harmonic_target(phase=0.,velocity=True):
    p=HarmonicProjection(phase=phase);angles=np.linspace(0.,PERIOD,721)
    row=dict(position=(p.value(angles)-1).tolist(),first_derivative=p.value(angles,1).tolist() if velocity else None,
        second_derivative=None,events=[MotionEvent('maximum',phase%PERIOD),MotionEvent('minimum',(phase+math.pi)%PERIOD)])
    return MotionTarget.create(angles,dict(small=row,large=row),source={'identity':'harmonic-cadence-test'})


def test_position_and_internal_long_branch_variations_are_not_primary_objectives():
    target=harmonic_target();p=PrimaryCadencePolicy(long_symmetry_weight=.45)
    features=primary_target_features(PeriodicTargetSide(target,'large'),p)
    assert features['cadence']==dict(available=False,reason='no_resolved_two_cadence_transition')
    simple=primary_cadence_score(HarmonicProjection(),features,p)
    varying=primary_cadence_score(HarmonicProjection(third=.15),features,p)
    assert varying['score']==pytest.approx(simple['score'],abs=1e-12)
    assert varying['long_mirror_asymmetry_rms']<1e-12
    angles=np.linspace(0.,PERIOD,721)
    assert np.sqrt(np.mean((HarmonicProjection(third=.15).value(angles)-HarmonicProjection().value(angles))**2))>.01
    # This symmetric branch has multiple speed peaks, but still two real turns.
    assert len(HarmonicProjection(third=.15).roots)==2
    assert not any('position' in name or 'profile' in name for name in varying['terms'])


@pytest.mark.parametrize('phase',[0.,.4,PERIOD-.2,PERIOD+.4])
def test_cyclic_phase_and_signed_turnaround_correspondence(phase):
    p=PrimaryCadencePolicy();target=harmonic_target(phase)
    features=primary_target_features(PeriodicTargetSide(target,'small'),p)
    good=primary_cadence_score(HarmonicProjection(phase=phase),features,p)
    wrong=primary_cadence_score(HarmonicProjection(sign=-1.,phase=phase),features,p)
    assert good['turning_rms_rad']<1e-8
    assert good['endpoint_zero_rms']<1e-8
    assert wrong['score']>good['score']+1
    assert wrong['monotonicity_wrong_sign_rms']>0


def test_additional_turnarounds_are_penalized_or_explicitly_rejected():
    p=PrimaryCadencePolicy();features=primary_target_features(PeriodicTargetSide(harmonic_target(),'small'),p)
    projection=HarmonicProjection(third=.4)
    evidence=primary_cadence_score(projection,features,p)
    assert evidence['extra_crossings']==4 and evidence['zero_crossing_count']==6
    assert evidence['terms']['extra_crossing']==4*p.extra_crossing_weight
    assert evidence['monotonicity_wrong_sign_rms']>0
    with pytest.raises(ValueError,match='Additional primary'):
        primary_cadence_score(projection,features,replace(p,reject_extra_turnarounds=True))


def test_unavailable_velocity_and_missing_cadence_have_explicit_behavior():
    source=PeriodicTargetSide(harmonic_target(velocity=False),'small')
    assert primary_target_features(source,PrimaryCadencePolicy())['cadence']['reason']=='target_velocity_unavailable'
    with pytest.raises(ValueError,match='cadence unavailable'):
        primary_target_features(source,PrimaryCadencePolicy(missing_cadence='error'))
    with pytest.raises(ValueError,match='cadence unavailable'):
        primary_target_features(PeriodicTargetSide(harmonic_target(),'small'),PrimaryCadencePolicy(missing_cadence='error'))


def test_event_supported_short_cadence_ratio_and_displacement_sharing():
    from scipy.integrate import cumulative_trapezoid
    from scipy.interpolate import PchipInterpolator
    # Smooth skewed speed, with a marked transition on a wrapped short branch.
    phase=5.8;span=2.;u=np.linspace(0.,1.,4001)
    speed=u*(1-u)*np.exp(-8*u)
    progress=cumulative_trapezoid(speed,u,initial=0);norm=progress[-1];progress/=norm
    position=PchipInterpolator(u,progress)
    angles=np.linspace(0.,PERIOD,1441);offset=(angles-phase)%PERIOD
    short=offset<span;q=np.empty_like(angles);v=np.empty_like(angles)
    q[short]=1-position(offset[short]/span)
    v[short]=-((offset[short]/span)*(1-offset[short]/span)*np.exp(-8*offset[short]/span))/norm/span
    long=(offset[~short]-span)/(PERIOD-span)
    q[~short]=.5*(1-np.cos(math.pi*long));v[~short]=.5*math.pi*np.sin(math.pi*long)/(PERIOD-span)
    q[-1]=q[0];v[-1]=v[0]
    events=[MotionEvent('maximum',phase),MotionEvent('minimum',(phase+span)%PERIOD),MotionEvent('kink',(phase+.3*span)%PERIOD)]
    row=dict(position=q.tolist(),first_derivative=v.tolist(),second_derivative=None,events=events)
    target=MotionTarget.create(angles,dict(small=row,large=row),source={'identity':'skewed-short-cadence'})
    p=PrimaryCadencePolicy();features=primary_target_features(PeriodicTargetSide(target,'small'),p)
    cadence=features['cadence']
    assert cadence['available'] and cadence['source']=='event:kink'
    assert cadence['split_offset_rad']==pytest.approx(.6)
    assert cadence['first_is_fast'] and cadence['fast_slow_speed_ratio']>2
    assert .5<cadence['fast_displacement_fraction']<1
    without_events=MotionTarget.create(angles,dict(small=dict(row,events=events[:2]),large=row),source={'identity':'numerical-short-cadence'})
    detected=primary_target_features(PeriodicTargetSide(without_events,'small'),p)['cadence']
    assert detected['available'] and detected['source']=='two_level_velocity_estimate'
    assert detected['explained_variance']>=p.minimum_split_explained_variance


def test_long_symmetry_defaults_to_historical_weight_and_can_be_disabled():
    features=primary_target_features(PeriodicTargetSide(harmonic_target(),'small'),PrimaryCadencePolicy())
    projection=HarmonicProjection(second_sine=.1)
    disabled=primary_cadence_score(projection,features,PrimaryCadencePolicy(long_symmetry_weight=0.))
    enabled=primary_cadence_score(projection,features,PrimaryCadencePolicy(long_symmetry_weight=.45))
    assert disabled['terms']['long_mirror_asymmetry']==0
    assert enabled['terms']['long_mirror_asymmetry']>0
    assert PrimaryCadencePolicy().long_symmetry_weight==.45
    assert enabled['terms']['long_mirror_asymmetry']==pytest.approx(.45*enabled['long_mirror_asymmetry_rms'])
    projection.metrics=dict(minimum_primary_transmission_sine=.1,E_axis_variance_fraction=.3)
    preferences=primary_cadence_score(projection,features,PrimaryCadencePolicy())
    assert preferences['terms']['transmission_preference']>0 and preferences['terms']['directionality_preference']>0
    assert preferences['score']==pytest.approx(sum(preferences['terms'].values()))


def test_real_pca_sign_choice_is_not_selected_by_position_rms():
    artifact=primary_library(target_for()).members[0]['mechanisms']['large']
    from dada_solver.research.artifacts import MechanismArtifact
    a=MechanismArtifact.from_data(artifact);projection=PrimaryProjection(a.reconstruct(),sign=-1.)
    angles=np.linspace(0.,PERIOD,721)
    extrema={r['kind']:r['angle_rad'] for r in projection.roots}
    row=dict(position=(projection.value(angles)-1).tolist(),first_derivative=projection.value(angles,1).tolist(),
        second_derivative=None,events=[MotionEvent(k,v) for k,v in extrema.items()])
    target=MotionTarget.create(angles,dict(small=row,large=row),source={'identity':'negative-PCA'})
    evidence=primary_evidence(a,target,'large',360)
    np.testing.assert_allclose(evidence['primary_projection']['axis'],projection.axis)
    assert evidence['primary_cadence']['turning_rms_rad']<1e-8
    assert evidence['primary_cadence']['policy_version']==PRIMARY_POLICY
    assert 'Diagnostic only' in evidence['primary_projection']['meaning']


@pytest.mark.parametrize('bad',[{'monotonicity_weight':-1},{'long_symmetry_weight':True},
    {'missing_cadence':'invent'},{'turning_scale_rad':0.},{'minimum_split_explained_variance':2.}])
def test_policy_validation(bad):
    with pytest.raises(ValueError):SearchPolicy(primary_cadence=bad)


def test_stationary_pause_is_not_a_parasitic_reversal():
    p=PrimaryCadencePolicy(reject_extra_turnarounds=True)
    features=primary_target_features(PeriodicTargetSide(harmonic_target(),'small'),p)
    projection=HarmonicProjection()
    projection.roots=(*projection.roots,dict(kind='stationary',angle_rad=1.,position=1.7))
    evidence=primary_cadence_score(projection,features,p)
    assert evidence['extra_crossings']==0
    assert evidence['stationary_point_count']==3


def test_primary_catalogue_shows_cadence_components_and_diagnostic_position(tmp_path):
    from dada_solver.research.synthesis_catalogue import render_synthesis_catalogue
    target=target_for();library=primary_library(target)
    path=render_synthesis_catalogue(library,target,tmp_path/'catalogue.html')
    text=path.read_text()
    for label in ('Primary topology / cadence','Internal cadence score','Turnarounds (study rad)',
                  'Fast/slow speed ratio','Fast displacement fraction','monotonicity',
                  'endpoint_zero','transmission_preference','directionality_preference',
                  'position and velocity errors are diagnostics only',PRIMARY_POLICY):
        assert label in text


def test_primary_constraints_remain_hard_and_effective_policy_is_recorded():
    target=target_for();p=policy(primary_cadence={'long_symmetry_weight':.45})
    engine=SixBarStageSearch(plan(target,'primary_discovery'),p,('large',),None,None)
    categories={'primary_branch':1};x=engine.encode(geometry())
    row=engine.evaluate(x,'large',categories,{'island':0,'seed':42},retain=True)
    assert row['valid']
    library=engine.final_library();artifact=library.members[0]['mechanisms']['large']
    assert artifact['provenance']['primary_cadence_policy']['long_symmetry_weight']==.45
    assert artifact['provenance']['primary_score_samples']==p.samples
    assert artifact['provenance']['primary_root_samples']==p.root_samples
    constrained=SixBarStageSearch(plan(target,'primary_discovery',constraints=(dict(
        metric='minimum_primary_transmission_sine',relation='minimum',limit=.999,unit='1'),)),p,('large',),None,None)
    assert not constrained.evaluate(x,'large',categories,{},retain=False)['valid']


@pytest.mark.parametrize('sine,admitted',[(.29,False),(.30,True),(.31,True)])
def test_primary_transmission_is_hard_with_inclusive_roundoff_boundary(sine,admitted,monkeypatch):
    import dada_solver.research.synthesis_six_bar as six
    g=dict(geometry(),primary_coupler=2.,primary_rocker=2.,primary_phase=0.,
           primary_ground=math.sqrt(8+8*math.sqrt(1-sine*sine))-1)
    target=harmonic_target();engine=SixBarStageSearch(plan(target,'primary_discovery'),policy(),('large',),None)
    if not admitted:
        monkeypatch.setattr(six,'primary_cadence_score',lambda *a,**k:pytest.fail('A good cadence cannot rescue transmission.'))
    row=engine.evaluate(engine.encode(g),'large',{'primary_branch':1},{})
    assert row['valid']==admitted
    assert row['metrics']['minimum_primary_transmission_sine']==pytest.approx(sine,abs=1e-14)
    record=next(r for r in row['primary_constraints'] if r['name']=='minimum_primary_transmission_sine')
    assert record['satisfied']==admitted
    if admitted:
        assert row['terms']['transmission_preference']>0  # 0.30 floor is not the 0.35 preference.
    else:
        assert row['reason']=='primary_transmission_below_minimum'
        assert not engine.archive.rows


def test_cartesian_excursion_rejects_before_pca_or_cadence(monkeypatch):
    import dada_solver.research.synthesis_six_bar as six
    g=dict(geometry(),primary_ground=.5,primary_coupler=1.,primary_rocker=1.,
           primary_e_along=.5,primary_e_normal=math.sqrt(.75),primary_phase=0.,primary_branch=-1)
    p=policy(bounds={'primary_ground':(.1,10.),'primary_coupler':(.5,12.),'primary_rocker':(.5,10.)})
    engine=SixBarStageSearch(plan(harmonic_target(),'primary_discovery'),p,('large',),None)
    monkeypatch.setattr(six,'PrimaryProjection',lambda *a,**k:pytest.fail('PCA must follow primary hard screening.'))
    row=engine.evaluate(engine.encode(g),'large',{'primary_branch':-1},{})
    assert not row['valid'] and row['reason']=='primary_e_span_below_minimum'
    assert row['metrics']['minimum_primary_transmission_sine']>.3
    assert row['metrics']['E_span_over_crank']<.75
    assert engine.rejections['primary_e_span_below_minimum']==1
    assert not engine.archive.rows


def test_cartesian_span_is_not_pca_span_and_screen_is_sign_independent():
    from dada_solver.research.mechanical_screen import PrimaryMechanicalPolicy
    from dada_solver.research.synthesis_six_bar import primary_cycle_data,primary_rejection_reason
    class EllipticPrimary:
        def joint_state(self,t):
            t=np.asarray(t)
            return dict(joints={'E':((.4*np.cos(t)-.01*np.sin(t))/math.sqrt(2),
                                    (.4*np.cos(t)+.01*np.sin(t))/math.sqrt(2))},
                        E_derivative=((- .4*np.sin(t)-.01*np.cos(t))/math.sqrt(2),
                                      (- .4*np.sin(t)+.01*np.cos(t))/math.sqrt(2)),
                        primary_transmission_sine=np.full_like(t,.8))
    cycle=primary_cycle_data(EllipticPrimary(),1440)
    plus=PrimaryProjection(EllipticPrimary(),sign=1.,cycle=cycle)
    minus=PrimaryProjection(EllipticPrimary(),sign=-1.,cycle=cycle)
    assert plus.metrics['E_projection_span_over_crank']>.75
    assert plus.metrics['E_span_over_crank']<.75
    p=PrimaryMechanicalPolicy()
    assert p.constraints(plus.metrics)==p.constraints(minus.metrics)
    assert primary_rejection_reason(p.constraints(plus.metrics))=='primary_e_span_below_minimum'
    relaxed=PrimaryMechanicalPolicy(minimum_e_span_over_crank=False)
    assert all(c['satisfied'] for c in relaxed.constraints(plus.metrics))


def test_hard_metrics_are_computed_once_before_two_projection_signs(monkeypatch):
    import dada_solver.research.synthesis_six_bar as six
    original=six.primary_cycle_data;calls=[]
    def counted(primary,samples):
        calls.append(samples);return original(primary,samples)
    monkeypatch.setattr(six,'primary_cycle_data',counted)
    engine=SixBarStageSearch(plan(harmonic_target(),'primary_discovery'),policy(),('large',),None)
    row=engine.evaluate(engine.encode(geometry()),'large',{'primary_branch':1},{})
    assert row['valid'] and calls==[360]


def test_final_primary_screen_is_recomputed_at_1440_and_cannot_be_compensated(monkeypatch):
    import dada_solver.research.synthesis_six_bar as six
    engine=SixBarStageSearch(plan(harmonic_target(),'primary_discovery'),policy(),('large',),None)
    assert engine.evaluate(engine.encode(geometry()),'large',{'primary_branch':1},{})['valid']
    original=six.primary_cycle_data;calls=[]
    def fail_final(primary,samples):
        calls.append(samples);state,points,metrics=original(primary,samples)
        if samples==1440:metrics=dict(metrics,minimum_primary_transmission_sine=.29,E_span_over_crank=.74)
        return state,points,metrics
    monkeypatch.setattr(six,'primary_cycle_data',fail_final)
    with pytest.raises(ValueError,match='No admissible'):engine.final_library()
    assert calls==[1440]
    assert engine.rejections['primary_transmission_below_minimum']==1
    assert engine.rejections['primary_e_span_below_minimum']==1
    assert engine.rejections['final_primary_transmission_below_minimum']==1
    assert engine.rejections['final_primary_e_span_below_minimum']==1
    assert all(r['margin']<0 for r in engine.rejection_evidence[-1]['constraints'])


def test_primary_profile_matches_archived_thresholds_without_complete_constraints(tmp_path,capsys):
    import tomllib
    from pathlib import Path
    from dada_solver.research.mechanical_screen import PrimaryMechanicalPolicy
    from dada_solver.research.cli import main
    archived=tomllib.loads(Path('docs/repro/mechanism_synthesis_search/policy.toml').read_text())
    p=PrimaryMechanicalPolicy();profile=p.profile()
    from importlib.resources import files
    resource=tomllib.loads(files('dada_solver.research').joinpath('data/six_bar_primary_design.toml').read_text())
    assert profile['thresholds']==resource['thresholds']
    for name,value in profile['thresholds'].items():assert value==archived['hard_constraints'][name]
    assert profile['final_mechanical_samples']==1440
    assert PrimaryCadencePolicy().long_symmetry_weight==archived['score_weights']['long_branch_mirror_asymmetry']
    target=harmonic_target();source=tmp_path/'target.json';target.save(source)
    args=['mechanism','synthesize',str(source),'--family','six_bar','--stage','primary_discovery','--validate-only']
    assert main(args)==0
    import json
    assert json.loads(capsys.readouterr().out)['primary_mechanical']==profile
    config=tmp_path/'search.toml';config.write_text('[primary_mechanical]\nminimum_primary_transmission_sine = false\nminimum_e_span_over_crank = 0.5\n[primary_cadence]\nlong_symmetry_weight = 0.0\n')
    assert main([*args,'--config',str(config)])==0
    thresholds=json.loads(capsys.readouterr().out)['primary_mechanical']['thresholds']
    assert thresholds==dict(minimum_primary_transmission_sine=None,minimum_e_span_over_crank=.5)
    engine=SixBarStageSearch(plan(target,'primary_discovery'),policy(),('large',),None)
    engine.evaluate(engine.encode(geometry()),'large',{'primary_branch':1},{})
    library=engine.final_library();member=library.members[0]
    assert member['mechanisms']['large']['scientific']['constraints']==[]
    assert member['mechanisms']['large']['provenance']['primary_mechanical_screen']==profile
    mechanical=member['metadata']['evidence']['large']['mechanical']
    assert mechanical['samples']==1440 and all(c['satisfied'] for c in mechanical['primary_constraints'])
    from dada_solver.research.synthesis_catalogue import render_synthesis_catalogue
    html=render_synthesis_catalogue(library,target,tmp_path/'catalogue.html').read_text()
    assert 'E_span_over_crank' in html and 'six_bar_primary_design_v1' in html
    assert '1440' in html and 'Weighted long symmetry contribution' in html
