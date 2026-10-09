"""Explicit family ownership and construction; no search policy or thermodynamics."""
from dataclasses import dataclass, fields
import math
from types import SimpleNamespace
import numpy as np
from dada_solver.composed_kinematics import CylinderLaw, ComposedKinematics
from dada_solver.four_bar import (FourBarLoop, FourBarSliderAssembly, SliderConstraint,
    RockerOutputPoint, CouplerOutputPoint, SharedCrankFourBarVolumeKinematics)
from dada_solver.six_bar import SixBarCylinderMechanism, IndependentSixBarVolumeKinematics
from dada_solver.slider_crank import SliderMotion, SliderCrankKinematics
from dada_solver.free_kinematics import FreeMotionDefinition, FreeKinematicsConfiguration, FreeKinematics
from dada_solver.phased_free_kinematics import PhaseShiftedFreeKinematics
from dada_solver.fourier_kinematics import FourierVolumeKinematics
from dada_solver.structured_kinematics import StructuredKinematics15
from dada_solver.hybrid_compact_kinematics import HybridCompactKinematics
from dada_solver.kinematics import HarmonicVolumeKinematics, IdealPiecewiseLinearVolumeKinematics
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.independent_four_stage_kinematics import IndependentFourStageVolumeKinematics
from dada_solver.mechanism_diagnostics import six_bar_metrics, four_bar_metrics, zero_crossings

FAMILIES=('harmonic','slider_crank','four_bar','six_bar','free_spline','fourier_c2',
          'structured_c2_15p','ideal_piecewise','four_stage','independent_four_stage','hybrid_compact')
PHYSICAL_FAMILIES=('slider_crank','four_bar','six_bar')
SIXBAR_CONTINUOUS=tuple(f.name for f in fields(SixBarCylinderMechanism) if f.init and 'branch' not in f.name)
PRIMARY_COORDINATES=SIXBAR_CONTINUOUS[:6]
DOWNSTREAM_COORDINATES=SIXBAR_CONTINUOUS[6:]
STRUCTURED_DEFAULTS=dict(small_max_deg=180.,small_down_duration_deg=180.,large_down_duration_deg=180.,
    small_max_curvature=16.,small_min_curvature=16.,large_max_curvature=16.,large_min_curvature=16.,
    small_up_bp_mid_q=.5,large_down_bp_mid_q=.5,small_down_kink_u=.5,small_down_kink_q=.5,
    small_down_kink_width_rel=.5,large_up_kink_u=.5,large_up_kink_q=.5,large_up_kink_width_rel=.5)
HYBRID_COMPACT_DEFAULTS = {'small_max_deg': 159.9896807151838, 'small_down_duration_deg': 154.34898715570617, 'large_down_duration_deg': 204.75091539779868, 'large_down_rounding': 0.09031610971409856, 'small_up_rounding': 0.1257187427943574, 'small_down_kink_u': 0.5641084556417698, 'small_down_kink_q': 0.7123342432191856, 'large_up_kink_u': 0.6898642219110528, 'large_up_kink_q': 0.445918856800329}


@dataclass(frozen=True)
class ParameterSpec:
    unit: str = '1'
    kind: str = 'continuous'
    positive: bool = False
    choices: tuple = ()

    def validate(self, value):
        if self.kind=='choice':
            if value not in self.choices: raise ValueError('Unknown categorical value.')
            return
        if self.kind=='boolean':
            if type(value) is not bool: raise ValueError('Expected a boolean scientific input.')
            return
        if self.kind=='branch':
            if type(value) is not int or value not in (-1,1): raise ValueError('Assembly branch/direction requires integer -1 or +1.')
            return
        if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value):
            raise ValueError('Expected a finite number.')
        if self.kind=='integer' and type(value) is not int: raise ValueError('Expected an integer count.')
        if self.positive and value<=0: raise ValueError('Expected a positive value.')


def parameter_specs(settings, side):
    family=settings.get('family')
    if family not in FAMILIES: raise ValueError(f'Unknown kinematic family: {family!r}.')
    spec=lambda names,unit='1':{k:ParameterSpec(unit) for k in names}
    if family=='harmonic': return spec(['phase_rad'],'rad')
    if family=='slider_crank':
        return dict(rod_over_crank=ParameterSpec('crank_radius',positive=True),offset_over_crank=ParameterSpec('crank_radius'),
            phase_rad=ParameterSpec('rad'),volume_increases_with_coordinate=ParameterSpec(kind='boolean'))
    if family=='six_bar':
        result={k:ParameterSpec('rad' if k in ('primary_phase','slider_axis_angle') else
            '1' if k.startswith('h_') else 'crank_radius',positive=k in ('primary_ground','primary_coupler','primary_rocker','link_ef','link_gf','piston_rod')) for k in SIXBAR_CONTINUOUS}
        result.update(primary_branch=ParameterSpec(kind='branch'),second_branch=ParameterSpec(kind='branch'))
        return result
    if family=='four_bar':
        if settings.get('output') not in ('rocker','coupler'): raise ValueError('Four-bar output must be rocker or coupler.')
        result=spec(['coupler','rocker','ground_x','ground_y','output_along','output_normal','rod_length','slider_origin_x','slider_origin_y'],'crank_radius')
        for k in ('coupler','rocker','rod_length'): result[k]=ParameterSpec('crank_radius',positive=True)
        result.update(spec(['axis_angle','phase_rad'],'rad'))
        result.update(loop_branch=ParameterSpec(kind='branch'),slider_branch=ParameterSpec(kind='branch'),
            crank_direction=ParameterSpec(kind='branch'),volume_increases_with_coordinate=ParameterSpec(kind='boolean'))
        return result
    if family=='free_spline':
        count=settings.get('count'); representation=settings.get('representation')
        if type(count) is not int or count<4: raise ValueError('Spline count must be an integer >= 4.')
        if representation not in ('controls','shape_coordinates'): raise ValueError('Choose controls or shape_coordinates.')
        prefix,n=('control',count) if representation=='controls' else ('shape',count-2)
        result={f'{prefix}_{i}':ParameterSpec(kind='fixed_continuous' if prefix=='control' else 'continuous') for i in range(n)}
        # The chart removes control offset/amplitude, not a continuous shift of
        # the spline knots. Phase remains an independently owned coordinate.
        result['phase_rad']=ParameterSpec('rad')
        return result
    if family=='fourier_c2':
        h=settings.get('harmonics')
        if type(h) is not int or h<1: raise ValueError('Fourier harmonics must be a positive integer.')
        return spec([f'coefficient_{i}' for i in range(2*h)])
    if family=='structured_c2_15p':
        return {k:ParameterSpec('deg' if k.endswith('_deg') else '1') for k in STRUCTURED_DEFAULTS if k.startswith(side+'_')}
    if family=='hybrid_compact':
        return {k:ParameterSpec('deg' if k.endswith('_deg') else '1') for k in HYBRID_COMPACT_DEFAULTS if k.startswith(side+'_')}
    if family=='ideal_piecewise': return spec(['small_lambda_target','large_lambda_target','adiabatic_sector_fraction'])
    if family=='four_stage': return spec(['t1','t2','t3',*(['a_s','b_s'] if side=='small' else ['a_l','b_l'])])
    return spec(['t0_s','t1_s','t2_s','t3_s','a_s','b_s'] if side=='small' else ['t1_l','t2_l','t3_l','a_l','b_l'])


def validate_settings(settings, side):
    family=settings.get('family')
    extra={'four_bar':{'output','envelope_frame_angle_rad','envelope_frame'},'free_spline':{'representation','count'},'fourier_c2':{'harmonics'}}.get(family,set())
    allowed={'family','crank_radius_m','artifact','sha256'}|extra
    if set(settings)-allowed: raise ValueError(f'Unknown {side} family settings: {sorted(set(settings)-allowed)}')
    if ('artifact' in settings)!=('sha256' in settings): raise ValueError('An artifact requires its content hash.')
    if 'artifact' in settings and family not in PHYSICAL_FAMILIES: raise ValueError('Mechanism artifacts require a physical family.')
    if 'crank_radius_m' in settings:
        if family not in PHYSICAL_FAMILIES: raise ValueError('An abstract motion has no physical crank scale.')
        ParameterSpec('m',positive=True).validate(settings['crank_radius_m'])
    if 'envelope_frame_angle_rad' in settings: ParameterSpec('rad').validate(settings['envelope_frame_angle_rad'])
    if 'envelope_frame' in settings:
        if settings['envelope_frame'] != 'slider_axis' or 'envelope_frame_angle_rad' in settings:
            raise ValueError('envelope_frame must be slider_axis, without a fixed frame angle.')
    return parameter_specs(settings,side)


def build_side(settings, parameters, side, limits):
    """Select one cylinder of an existing production backend without changing its law."""
    if settings.get('component')=='primary':
        raise ValueError('An incomplete primary has no piston law; synthesize its downstream linkage first.')
    family=settings['family']; p=parameters; backend_side=side; stroke=None; geometry=None
    if family=='harmonic':
        # Reuse the production small-side phase convention for either cylinder.
        model=HarmonicVolumeKinematics(limits,limits,p['phase_rad']); backend_side='small'
    elif family=='slider_crank':
        geometry=SliderMotion(**p); model=SliderCrankKinematics(geometry,geometry,limits,limits)
        stroke=geometry.stroke_over_crank
    elif family=='six_bar':
        geometry=SixBarCylinderMechanism(**p); model=IndependentSixBarVolumeKinematics(geometry,geometry,limits,limits)
        stroke=geometry.stroke_over_crank
    elif family=='four_bar':
        output=(RockerOutputPoint if settings['output']=='rocker' else CouplerOutputPoint)(p['output_along'],p['output_normal'])
        geometry=FourBarSliderAssembly(FourBarLoop(p['coupler'],p['rocker'],p['ground_x'],p['ground_y'],p['loop_branch']),
            output,SliderConstraint(p['slider_origin_x'],p['slider_origin_y'],p['axis_angle'],p['rod_length'],p['slider_branch']))
        model=SharedCrankFourBarVolumeKinematics(1.,geometry,geometry,limits,limits,
            p['volume_increases_with_coordinate'],p['volume_increases_with_coordinate'],p['phase_rad'],p['crank_direction'])
        stroke=model._small.stroke
    elif family=='free_spline':
        if settings['representation']=='controls':
            definition=FreeMotionDefinition(tuple(p[f'control_{i}'] for i in range(settings['count'])),limits.minimum,limits.maximum)
        else:
            definition=FreeMotionDefinition.from_shape_coordinates([p[f'shape_{i}'] for i in range(settings['count']-2)],limits.minimum,limits.maximum)
        model=PhaseShiftedFreeKinematics(FreeKinematics(FreeKinematicsConfiguration(definition,definition)),p['phase_rad'],p['phase_rad'])
    elif family=='fourier_c2':
        coefficients=tuple(p[f'coefficient_{i}'] for i in range(2*settings['harmonics']))
        model=FourierVolumeKinematics(limits,limits,coefficients,coefficients,settings['harmonics'])
    elif family=='structured_c2_15p':
        model=StructuredKinematics15(SimpleNamespace(small_cylinder=limits,large_cylinder=limits),dict(STRUCTURED_DEFAULTS,**p))
        # Only the selected piston contributes; the other helper branch is unused.
    elif family=='hybrid_compact':
        model=HybridCompactKinematics(limits,limits,**dict(HYBRID_COMPACT_DEFAULTS,**p))
    elif family=='ideal_piecewise': model=IdealPiecewiseLinearVolumeKinematics(limits,limits,**p)
    elif family=='four_stage':
        model=FourStageVolumeKinematics(limits,limits,**dict(dict(t1=.25,t2=.5,t3=.75,a_l=.5,b_l=.5,a_s=.5,b_s=.5),**p))
    else:
        model=IndependentFourStageVolumeKinematics(limits,limits,**dict(dict(t0_s=0.,t1_s=.25,t2_s=.5,t3_s=.75,t1_l=.25,t2_l=.5,t3_l=.75,a_l=.5,b_l=.5,a_s=.5,b_s=.5),**p))
    physical=None if stroke is None or 'crank_radius_m' not in settings else stroke*settings['crank_radius_m']
    return CylinderLaw(model,backend_side,physical),geometry


def side_metrics(settings, law, geometry, samples):
    family=settings['family']
    angles=np.linspace(0,2*math.pi,samples,endpoint=False)
    if family in ('slider_crank','four_bar','six_bar','free_spline'):
        velocity=np.asarray(law.value(angles,1)); position=np.asarray(law.value(angles))
    else:
        velocity=np.array([law.value(float(t),1) for t in angles])
        position=np.array([law.value(float(t)) for t in angles])
    limits=law.limits
    if np.min(position)<limits.minimum-1e-9*limits.swept or np.max(position)>limits.maximum+1e-9*limits.swept:
        raise ValueError('Motion exceeds its declared cylinder volume limits; no clipping is performed.')
    result=dict(zero_crossing_count=zero_crossings(velocity),maximum_absolute_first_derivative=float(np.max(np.abs(velocity))))
    try:
        result['maximum_absolute_second_derivative']=float(np.max(abs(law.value(angles,2)))) if family in ('slider_crank','free_spline') else max(abs(float(law.value(float(t),2))) for t in angles)
    except (AttributeError,NotImplementedError): pass
    if family=='free_spline':
        d=law.model.diagnostics[0]
        result.update(maximum_absolute_first_derivative=d.maximum_absolute_first_derivative,
                      maximum_absolute_second_derivative=d.maximum_absolute_second_derivative)
    elif family=='six_bar': result.update(six_bar_metrics(geometry,samples))
    elif family=='four_bar':
        frame = -float(geometry.slider.axis_angle) if settings.get('envelope_frame')=='slider_axis' else settings.get('envelope_frame_angle_rad',0.)
        result.update(four_bar_metrics(law.model,law.side,samples,frame))
        result['envelope_frame_angle_rad'] = frame
    elif family=='slider_crank':
        l,e=geometry.rod_over_crank,geometry.offset_over_crank
        result.update(stroke_over_crank=geometry.stroke_over_crank,closure_margin=l-1-abs(e),
                      minimum_rod_axis_cosine=math.sqrt(l*l-(1+abs(e))**2)/l)
    return result


MECHANICAL_METRICS={
    'maximum_absolute_first_derivative':('m^3/rad','maximum'),
    'maximum_absolute_second_derivative':('m^3/rad^2','maximum'),
    'zero_crossing_count':('1','equal'),
    'stroke_over_crank':('crank_radius',None),
    'minimum_primary_transmission_sine':('1','minimum'),
    'minimum_secondary_transmission_sine':('1','minimum'),
    'minimum_rod_axis_cosine':('1','minimum'),
    'stroke_over_envelope':('1','minimum'),
    'EH_over_crank':('crank_radius','maximum'),
    'H_axis_lateral_rms_over_stroke':('1','maximum'),
    'H_axis_lateral_span_over_stroke':('1','maximum'),
    'crank_axis_to_EFH_clearance_over_crank':('crank_radius','minimum'),
    'closure_margin':('crank_radius','minimum')}


def available_metrics(family):
    common={'maximum_absolute_first_derivative','zero_crossing_count'}
    if family not in ('four_bar','six_bar','ideal_piecewise','four_stage','independent_four_stage','hybrid_compact'): common.add('maximum_absolute_second_derivative')
    if family in PHYSICAL_FAMILIES: common.update(('stroke_over_crank','minimum_rod_axis_cosine'))
    if family in ('four_bar','six_bar'): common.add('minimum_primary_transmission_sine')
    if family=='four_bar': common.add('stroke_over_envelope')
    if family=='six_bar': common.update(('minimum_secondary_transmission_sine','EH_over_crank','H_axis_lateral_rms_over_stroke','H_axis_lateral_span_over_stroke','crank_axis_to_EFH_clearance_over_crank'))
    if family=='slider_crank': common.add('closure_margin')
    return common
