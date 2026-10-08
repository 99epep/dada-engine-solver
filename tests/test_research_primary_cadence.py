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


def test_long_symmetry_is_opt_in_and_transmission_preferences_are_separate():
    features=primary_target_features(PeriodicTargetSide(harmonic_target(),'small'),PrimaryCadencePolicy())
    projection=HarmonicProjection(second_sine=.1)
    disabled=primary_cadence_score(projection,features,PrimaryCadencePolicy())
    enabled=primary_cadence_score(projection,features,PrimaryCadencePolicy(long_symmetry_weight=.45))
    assert disabled['terms']['long_mirror_asymmetry']==0
    assert enabled['terms']['long_mirror_asymmetry']>0
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
