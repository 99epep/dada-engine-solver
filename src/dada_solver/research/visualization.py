"""Scientific motion/frame samples for future renderers, using solver geometry.

This is a kinematic export only. It never integrates or replays thermodynamics.
"""
import math
import numpy as np
from .schema import compile_study


def sample_motion(study, active_values=None, *, samples=361):
    if type(samples) is not int or samples<3: raise ValueError('At least three angle samples are required.')
    definition=compile_study(study)
    values={p.name:p.initial for p in study.space.parameters}
    if active_values is not None:
        study.space.encode(active_values); values=dict(active_values)
    values=dict(definition.fixed_parameters,**values)
    design=definition.adapter.build(values)
    angles=np.linspace(0,2*math.pi,samples)
    result=dict(schema_version=1,study_id=study.study_id,angle_domain='study_angle_before_operation_transform',
        angle_unit='rad',derivative_variable='study_angle_radians',angle=angles.tolist(),sides={},
        coupling=study.data.get('kinematics',{}).get('coupling','independent'),thermodynamic_replay=False)
    for side in ('small','large'):
        model=design.kinematics; limits=getattr(model,side+'_volume_limits')
        volume=np.array([getattr(model,side+'_cylinder_volume')(float(t)) for t in angles])
        derivative=np.array([getattr(model,side+'_cylinder_volume_derivative')(float(t)) for t in angles])
        try:
            second=np.array([getattr(model,side+'_cylinder_volume_second_derivative')(float(t)) for t in angles])/limits.swept
            second=second.tolist()
        except (AttributeError,NotImplementedError): second=None
        out=dict(volume_m3=volume.tolist(),normalized_position=((volume-limits.minimum)/limits.swept).tolist(),
            normalized_velocity_per_rad=(derivative/limits.swept).tolist(),normalized_acceleration_per_rad2=second,
            physical_stroke_m=getattr(model,side+'_physical_stroke'),joints=None,links=None,
            coordinate_convention='stored local frame; no presentation offset or implied physical shaft layout')
        law=getattr(model,side)
        backend=getattr(law,'model',model)
        family=study.settings[side]['family'] if hasattr(study,'settings') else 'six_bar'
        geometry=getattr(backend,getattr(law,'side',side),None)
        if family=='slider_crank':
            states=[geometry.joint_state(float(t)) for t in angles]
            out.update(joints={k:[list(s[k]) for s in states] for k in states[0]},links=[['A','B'],['B','P']])
        elif family=='six_bar':
            states=[geometry.joint_state(float(t)) for t in angles]
            out.update(joints={k:[list(s['joints'][k]) for s in states] for k in states[0]['joints']},
                links=[['A','B'],['B','C'],['C','D'],['B','E'],['C','E'],['E','F'],['F','G'],['E','H'],['F','H'],['H','P']])
        elif family=='four_bar':
            assembly=getattr(backend,law.side+'_assembly'); states=[assembly.evaluate(*backend._crank(float(t))) for t in angles]
            slider=assembly.slider
            joints={'A':[[0.,0.]]*samples,'D':[[assembly.loop.rocker_pivot_x,assembly.loop.rocker_pivot_y]]*samples,
                'B':[list(s.crank_pin) for s in states],'C':[list(s.coupler_joint) for s in states],
                'H':[list(s.output_point) for s in states],
                'P':[[slider.axis_origin_x+s.coordinate*math.cos(slider.axis_angle),slider.axis_origin_y+s.coordinate*math.sin(slider.axis_angle)] for s in states]}
            out.update(joints=joints,links=[['A','B'],['B','C'],['C','D'],['H','P']]+([['C','H'],['D','H']] if study.settings[side]['output']=='rocker' else [['B','H'],['C','H']]))
        out['joint_length_unit']='crank_radius' if out['joints'] else None
        out['crank_radius_m']=study.settings[side].get('crank_radius_m') if hasattr(study,'settings') else None
        result['sides'][side]=out
    return result


def sample_report_volumes(study, active_values, *, samples=721):
    """Sample physical volumes in solver cycle angle, with operation applied once."""
    if type(samples) is not int or samples < 3:
        raise ValueError('At least three angle samples are required.')
    definition = compile_study(study)
    study.space.encode(active_values)
    design = definition.adapter.build(dict(definition.fixed_parameters, **active_values))
    built = design.build()
    model = getattr(built, 'model', built)
    angles = np.linspace(0., 2*math.pi, samples)
    volumes = [model.volumes(float(angle)) for angle in angles]
    return dict(schema_version=1, kind='volumes', angle_domain='solver_cycle_angle',
        operation_transform='applied once by production factory', angle_unit='deg',
        volume_unit='m^3', angle=np.degrees(angles).tolist(),
        small=[v.small_cylinder for v in volumes], large=[v.large_cylinder for v in volumes],
        thermodynamic_replay=False, method='current production kinematics; no ODE integration')
